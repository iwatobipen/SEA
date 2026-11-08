import argparse
import json
from pathlib import Path
from .data import home, list_datasets, save_dataset, build_chembl
from .jobs import create_key, run_workers


def demo():
    # Synthetic chemical groups, deliberately not biological target annotations.
    targets = [
        {"id":"DEMO-ALCOHOL","name":"Alcohol set (synthetic demo)","smiles":["CO","CCO","CCCO","CCCCO","CCCCCO","CCCCCCO","CC(O)C","CCC(O)C","CC(C)CO","OCCO"]},
        {"id":"DEMO-ACID","name":"Carboxylic acid set (synthetic demo)","smiles":["CC(=O)O","CCC(=O)O","CCCC(=O)O","CCCCC(=O)O","CCCCCC(=O)O","CC(C)C(=O)O","O=C(O)c1ccccc1","O=C(O)CCc1ccccc1","O=C(O)Cc1ccccc1","O=C(O)C=CC"]},
        {"id":"DEMO-AROMATIC","name":"Aromatic set (synthetic demo)","smiles":["c1ccccc1","Cc1ccccc1","CCc1ccccc1","Oc1ccccc1","Nc1ccccc1","COc1ccccc1","Clc1ccccc1","Fc1ccccc1","Brc1ccccc1","c1ccncc1"]}]
    from .chemistry import normalize
    for t in targets:
        t["smiles"] = [m.smiles for m in normalize(t["smiles"])[0]]
    return save_dataset(targets,{"source":"synthetic-demo","release":"1","data_license":"CC0-1.0","attribution":"Synthetic examples; no biological target claims"},sizes=(1,10,20),repeats=32)


def download(release, checksum):
    """Download only an explicit official release; verify before extracting SQLite."""
    import hashlib
    import tarfile
    import httpx
    if not release.isdigit() or len(checksum)!=64 or any(c not in "0123456789abcdefABCDEF" for c in checksum):
        raise ValueError("Specify release as an integer and sha256 as 64 hexadecimal characters.")
    folder = home()/"downloads"
    folder.mkdir(exist_ok=True)
    name = f"chembl_{release}_sqlite.tar.gz"
    dest = folder/name
    url = f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{release}/{name}"
    temporary = dest.with_suffix(".part")
    digest = hashlib.sha256()
    with httpx.stream("GET",url,timeout=120,follow_redirects=True) as response:
        response.raise_for_status()
        with temporary.open("wb") as f:
            for chunk in response.iter_bytes(1024*1024):
                digest.update(chunk)
                f.write(chunk)
    if digest.hexdigest()!=checksum.lower():
        temporary.unlink()
        raise ValueError("SHA256 mismatch")
    temporary.replace(dest)
    with tarfile.open(dest) as archive:
        members = [m for m in archive if m.isfile() and m.name.endswith(".db")]
        if not members:
            raise ValueError("SQLite .db file missing")
        for member in members:
            target = folder/Path(member.name).name
            with archive.extractfile(member) as src,target.open("wb") as out:
                import shutil
                shutil.copyfileobj(src,out)
            print(target)


def main():
    parser = argparse.ArgumentParser(prog="sea")
    subs = parser.add_subparsers(dest="command",required=True)
    subs.add_parser("demo")
    subs.add_parser("datasets")
    key = subs.add_parser("key")
    key.add_argument("name")
    worker = subs.add_parser("worker")
    worker.add_argument("--count",type=int,default=2)
    build = subs.add_parser("build")
    build.add_argument("--sqlite",required=True)
    build.add_argument("--release",required=True)
    build.add_argument("--config",default="config/chembl.json")
    build.add_argument("--repeats",type=int,default=128)
    build.add_argument("--trials",type=int,default=999)
    fetch = subs.add_parser("download")
    fetch.add_argument("--release",required=True)
    fetch.add_argument("--sha256",required=True)
    calibration = subs.add_parser("calibrate")
    calibration.add_argument("dataset_id")
    calibration.add_argument("--sizes",default="1,10,30,100,300,1000")
    calibration.add_argument("--repeats",type=int,default=128)
    calibration.add_argument("--trials",type=int,default=999)
    annotation = subs.add_parser("annotate")
    annotation.add_argument("dataset_id")
    annotation.add_argument("--sqlite",required=True)
    evaluation = subs.add_parser("evaluate")
    evaluation.add_argument("dataset_id")
    evaluation.add_argument("--output",default="var/evaluation.json")
    args = parser.parse_args()
    if getattr(args,"trials",999) < 1:
        parser.error("trials must be positive")
    if args.command=="demo":
        print(demo())
    elif args.command=="key":
        print(create_key(args.name))
    elif args.command=="datasets":
        print(json.dumps(list_datasets(),ensure_ascii=False,indent=2))
    elif args.command=="worker":
        if not 1<=args.count<=32:
            parser.error("count must be 1..32")
        run_workers(args.count)
    elif args.command=="build":
        settings=json.loads(Path(args.config).read_text())
        print(build_chembl(args.sqlite,args.release,settings,args.repeats,args.trials))
    elif args.command=="download":
        download(args.release,args.sha256)
    elif args.command=="calibrate":
        from .data import load_dataset
        data,_=load_dataset(args.dataset_id)
        provenance={k:v for k,v in data["manifest"].items() if k not in ("id","model_id","created_at","fingerprint","standardization","rdkit_version","target_count","molecule_count")}
        print(save_dataset(data["targets"],provenance,sizes=tuple(map(int,args.sizes.split(","))),repeats=args.repeats,empirical_trials=args.trials))
    elif args.command=="annotate":
        from .annotations import annotate_dataset
        print(json.dumps(annotate_dataset(args.dataset_id,args.sqlite),ensure_ascii=False,indent=2))
    elif args.command=="evaluate":
        from .evaluate import evaluate
        result = evaluate(args.dataset_id)
        path=Path(args.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
        print(path)
