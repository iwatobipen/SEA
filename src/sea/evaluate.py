"""Scaffold holdout target retrieval; unmeasured pairs are never labelled inactive."""
import hashlib
import time
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from .data import load_dataset
from .statistics import calibrate, significance
from .chemistry import raw_score


def evaluate(dataset_id, max_queries=100):
    started=time.perf_counter()
    data,structures=load_dataset(dataset_id)
    def is_test(s):
        scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s),includeChirality=False)
        # All acyclic molecules remain in one group, avoiding cross-split empty scaffolds.
        return int(hashlib.sha256(scaffold.encode()).hexdigest()[:8],16)%5==0
    train={s:m for s,m in structures.items() if not is_test(s)}
    holdout=sorted(s for s in structures if is_test(s))[:max_queries]
    if len(train)<10 or not holdout:
        return {"status":"insufficient_scaffold_groups","dataset_id":dataset_id,"train_molecules":len(train),"test_molecules":len(holdout)}
    reference=[{**t,"smiles":[s for s in t["smiles"] if s in train]} for t in data["targets"]]
    reference=[t for t in reference if len(t["smiles"])>=10]
    if not reference:
        return {"status":"insufficient_training_targets","dataset_id":dataset_id}
    pool=list(train.values())
    model=calibrate(pool,repeats=32)
    truth={s:{t["id"] for t in data["targets"] if s in t["smiles"]} & {t["id"] for t in reference} for s in holdout}
    metrics={method:{"recall_at_1":[],"recall_at_5":[],"recall_at_10":[],"mrr":[]} for method in ("sea","max_tanimoto")}
    null_cache={}
    for s in holdout:
        relevant=truth[s]
        if not relevant:
            continue
        ranking=[]
        for t in reference:
            right=[train[x] for x in t["smiles"]]
            score,_,_=raw_score([structures[s]],right,model["threshold"],0)
            sig=significance(score,1,len(right),pool,model,len(reference),null_cache=null_cache)
            maximum=max(DataStructs.BulkTanimotoSimilarity(structures[s].fp,[m.fp for m in right]))
            ranking.append((t["id"],sig["log_e_value"],score,maximum))
        for method,key in (("sea",lambda r:(r[1],-r[2],r[0])),("max_tanimoto",lambda r:(-r[3],r[0]))):
            ids=[r[0] for r in sorted(ranking,key=key)]
            ranks=[i+1 for i,t in enumerate(ids) if t in relevant]
            for k in (1,5,10):
                metrics[method][f"recall_at_{k}"].append(sum(r<=k for r in ranks)/len(relevant))
            metrics[method]["mrr"].append(1/min(ranks))
    return {"status":"completed","dataset_id":dataset_id,"split":"global Murcko scaffold hash modulo 5; 20% holdout",
            "max_queries":max_queries,"evaluated_queries":len(metrics["sea"]["mrr"]),"training_molecules":len(train),
            "training_targets":len(reference),"calibration":model,"elapsed_seconds":time.perf_counter()-started,
            "metrics":{m:{k:float(np.mean(v)) if v else None for k,v in values.items()} for m,values in metrics.items()},
            "interpretation":"Known-target retrieval only; unmeasured targets are not negatives. Empirical p-value resolution can cause ties; raw score breaks ties."}
