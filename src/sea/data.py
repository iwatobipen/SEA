import hashlib
import json
import os
from pathlib import Path
import sqlite3
from datetime import datetime, timezone
from rdkit import rdBase
from .chemistry import molecule, Molecule, FP_SPEC, STANDARDIZATION
from rdkit import DataStructs
from functools import lru_cache
from .statistics import calibrate


def home():
    path = Path(os.environ.get("SEA_HOME", "var"))
    path.mkdir(parents=True, exist_ok=True)
    return path


def digest_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def save_dataset(targets, provenance, sizes=(1, 10, 30, 100), repeats=128, seed=42, progress=lambda *_: None, empirical_trials=999):
    unique = sorted({s for t in targets for s in t["smiles"]})
    pool = [molecule(s) for s in unique]
    if len(pool) < 10:
        raise ValueError("背景分子は10以上必要です")
    model = calibrate(pool, sizes, repeats, seed, progress)
    model["empirical_trials"] = empirical_trials
    manifest = {"created_at": datetime.now(timezone.utc).isoformat(), "fingerprint": FP_SPEC,
                "standardization": STANDARDIZATION, "rdkit_version": rdBase.rdkitVersion,
                "target_count": len(targets), "molecule_count": len(unique), **provenance}
    payload = {"manifest": manifest, "model": model, "targets": sorted(targets, key=lambda t: t["id"]),
               "fingerprints": {m.smiles: DataStructs.BitVectToBinaryText(m.fp).hex() for m in pool}}
    # Timestamp does not affect content-addressed identity.
    stable = {**payload, "manifest": {k:v for k,v in manifest.items() if k != "created_at"}}
    dataset_id = hashlib.sha256(json.dumps(stable, sort_keys=True).encode()).hexdigest()[:20]
    model_id = hashlib.sha256(json.dumps(model, sort_keys=True).encode()).hexdigest()[:20]
    payload["manifest"].update({"id": dataset_id, "model_id": model_id})
    directory = home()/"datasets"
    directory.mkdir(exist_ok=True)
    dest = directory/f"{dataset_id}.json"
    if not dest.exists():
        temp = dest.with_suffix(f".{os.getpid()}.tmp")
        temp.write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False))
        temp.replace(dest)
    manifests = home()/"manifests"
    manifests.mkdir(exist_ok=True)
    sidecar = manifests/f"{dataset_id}.json"
    if not sidecar.exists():
        sidecar.write_text(json.dumps(payload["manifest"],ensure_ascii=False))
    return dataset_id


def load_dataset(dataset_id):
    return _load_dataset(str(home().resolve()), dataset_id)


@lru_cache(maxsize=2)
def _load_dataset(root, dataset_id):
    if len(dataset_id) != 20 or any(c not in "0123456789abcdef" for c in dataset_id):
        raise ValueError("不正なdataset ID")
    data = json.loads((Path(root)/"datasets"/f"{dataset_id}.json").read_text())
    if data["manifest"]["rdkit_version"] != rdBase.rdkitVersion:
        raise ValueError("RDKit版が参照DBと異なります。DBを再構築してください")
    if "fingerprints" in data:
        structures = {s: Molecule(s,DataStructs.CreateFromBinaryText(bytes.fromhex(fp))) for s,fp in data["fingerprints"].items()}
    else:
        structures = {s: molecule(s) for t in data["targets"] for s in t["smiles"]}
    return data, structures


def list_datasets():
    manifests = {p.stem:json.loads(p.read_text()) for p in (home()/"manifests").glob("*.json")}
    for p in (home()/"datasets").glob("*.json"):
        if p.stem not in manifests:
            manifests[p.stem] = json.loads(p.read_text())["manifest"]
    return [manifests[k] for k in sorted(manifests)]


def build_chembl(source, release, settings, repeats=128, empirical_trials=999):
    con = sqlite3.connect(f"file:{Path(source).resolve()}?mode=ro", uri=True)
    version = con.execute("SELECT * FROM version").fetchall()
    if release.casefold() not in str(version).casefold():
        raise ValueError(f"指定版 {release} がSQLiteのversion情報 {version} と一致しません")
    query = """
    SELECT td.chembl_id, td.pref_name, cs.canonical_smiles, a.standard_type,
           a.pchembl_value, a.activity_id, a.assay_id, a.doc_id,
           a.standard_relation, a.data_validity_comment, a.potential_duplicate,
           md.chembl_id
    FROM activities a JOIN assays ass ON a.assay_id=ass.assay_id
    JOIN target_dictionary td ON ass.tid=td.tid
    JOIN molecule_dictionary md ON a.molregno=md.molregno
    LEFT JOIN molecule_hierarchy mh ON a.molregno=mh.molregno
    JOIN compound_structures cs ON cs.molregno=COALESCE(mh.parent_molregno,a.molregno)
    WHERE td.tax_id=? AND td.target_type=?
      AND ass.confidence_score=? AND ass.assay_type=?
      AND a.standard_type IN ({types})
    ORDER BY a.activity_id
    """
    kinds = settings.get("activity_types", ["Ki", "Kd", "IC50"])
    if not kinds or not set(kinds) <= {"Ki", "Kd", "IC50", "EC50"}:
        raise ValueError("Unsupported activity_types")
    if settings.get("standard_relation", "=") != "=":
        raise ValueError("Only exact measurements are supported")
    query = query.format(types=",".join("?" for _ in kinds))
    parameters = (settings.get("organism_tax_id", 9606), settings.get("target_type", "SINGLE PROTEIN"),
                  settings.get("confidence_score", 9), settings.get("assay_type", "B"), *kinds)
    import statistics
    from collections import Counter, defaultdict
    counts = Counter()
    measurements = defaultdict(list)
    names, evidence = {}, defaultdict(list)
    cache = {}
    for tid, name, smiles, kind, value, aid, assay, doc, relation, validity, duplicate, mid in con.execute(query, parameters):
        counts["candidate_rows"] += 1
        if relation != "=" or value is None or validity not in (None, "Manually validated") or duplicate:
            counts["quality_excluded"] += 1
            continue
        try:
            if smiles not in cache:
                cache[smiles] = molecule(smiles).smiles
            canonical = cache[smiles]
        except Exception:
            counts["invalid_structure"] += 1
            continue
        measurements[(tid, canonical, kind)].append(float(value))
        names[tid] = name
        evidence[(tid, canonical)].append({"activity_id": aid, "assay_id": assay, "doc_id": doc, "molecule_id": mid, "type": kind, "pchembl": value})
    con.close()
    actives = defaultdict(set)
    for (tid, smiles, kind), values in measurements.items():
        if statistics.median(values) >= settings.get("min_pchembl", 6):
            actives[tid].add(smiles)
    targets = [{"id": tid, "name": names[tid], "smiles": sorted(s),
                "evidence": {m: evidence[(tid,m)] for m in sorted(s)}}
               for tid,s in sorted(actives.items()) if len(s) >= settings.get("min_ligands", 10)]
    if not targets:
        raise ValueError("条件を満たす標的がありません")
    provenance = {"source": "ChEMBL", "release": release, "source_sha256": digest_file(source),
                  "filters": settings, "counts": dict(counts), "sqlite_version": version,
                  "attribution": f"ChEMBL data is from https://www.ebi.ac.uk/chembl — version {release}",
                  "data_license": "CC-BY-SA-3.0"}
    return save_dataset(targets, provenance, repeats=repeats, empirical_trials=empirical_trials)
