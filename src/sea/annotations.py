"""ChEMBL target descriptions, protein family hierarchy, and species metadata."""
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import sqlite3

TARGET_FIELDS = ("id", "name", "family", "families", "protein_classes", "organism", "tax_id")


def family_for_path(path):
    """Collapse authoritative ChEMBL ancestors into readable target families."""
    lowered = [label.casefold() for label in path]
    if any("g protein-coupled receptor" in label for label in lowered):
        return "GPCR"
    for label, family in (("kinase", "Kinase"), ("ion channel", "Ion channel"),
                          ("nuclear receptor", "Nuclear receptor"), ("transporter", "Transporter")):
        if label in lowered:
            return family
    # Enzyme subclasses (e.g. Protease) are more informative than "Enzyme".
    if path and path[0] == "Enzyme" and len(path) > 1:
        return path[1]
    return path[0] if path else "Unclassified"


def read_target_annotations(con, target_ids):
    classes = {r[0]: (r[1], r[2], r[3]) for r in con.execute(
        "SELECT protein_class_id,parent_id,pref_name,class_level FROM protein_classification")}
    def hierarchy(class_id):
        labels, seen = [], set()
        while class_id in classes and class_id not in seen:
            seen.add(class_id)
            parent, name, level = classes[class_id]
            if level and name:
                labels.append(name)
            class_id = parent
        return tuple(reversed(labels))
    ids = sorted(set(target_ids))
    targets = {}
    for start in range(0, len(ids), 500):
        batch = ids[start:start+500]
        marks = ",".join("?" for _ in batch)
        for tid, name, organism, tax_id in con.execute(
                f"SELECT chembl_id,pref_name,organism,tax_id FROM target_dictionary WHERE chembl_id IN ({marks})",batch):
            targets[tid] = {"id":tid,"name":name,"organism":organism,"tax_id":tax_id,
                            "family":"Unclassified","families":[],"protein_classes":[]}
        paths = {tid:set() for tid in batch}
        for tid, class_id in con.execute(f"""
            SELECT td.chembl_id,cc.protein_class_id FROM target_dictionary td
            JOIN target_components tc ON tc.tid=td.tid
            JOIN component_class cc ON cc.component_id=tc.component_id
            WHERE td.chembl_id IN ({marks})""", batch):
            path = hierarchy(class_id)
            if path:
                paths[tid].add(path)
        for tid, path_set in paths.items():
            if tid in targets:
                families = sorted({family_for_path(path) for path in path_set})
                targets[tid].update({"family":" / ".join(families) or "Unclassified",
                                     "families":families,"protein_classes":[list(p) for p in sorted(path_set)]})
    return targets


def validate_release(con, release):
    versions = con.execute("SELECT * FROM version").fetchall()
    if not any(str(row[0]).casefold() == release.casefold() for row in versions):
        raise ValueError(f"Release {release} does not match the SQLite version table")
    return versions


def annotate_dataset(dataset_id, source):
    from .data import home, digest_file
    annotation_path(dataset_id)  # Validate the ID before accessing any dataset file.
    data = json.loads((home()/"datasets"/f"{dataset_id}.json").read_text())
    manifest = data["manifest"]
    if manifest.get("source") != "ChEMBL":
        raise ValueError("Only ChEMBL datasets can be annotated from a ChEMBL SQLite database")
    ids = [t["id"] for t in data["targets"]]
    del data  # Do not retain fingerprints and evidence while hashing the source.
    with sqlite3.connect(f"file:{Path(source).resolve()}?mode=ro", uri=True) as con:
        validate_release(con, manifest["release"])
        checksum = digest_file(source)
        if checksum != manifest["source_sha256"]:
            raise ValueError("SQLite SHA256 does not match the dataset source")
        targets = read_target_annotations(con, ids)
    missing = set(ids)-targets.keys()
    if missing:
        raise ValueError(f"Source is missing {len(missing)} dataset targets")
    document = {"schema_version":1,"dataset_id":dataset_id,"source":"ChEMBL",
                "release":manifest["release"],"source_sha256":checksum,"targets":targets}
    document["annotation_id"] = hashlib.sha256(json.dumps(document,sort_keys=True).encode()).hexdigest()[:20]
    document["created_at"] = datetime.now(timezone.utc).isoformat()
    dest = annotation_path(dataset_id)
    dest.parent.mkdir(exist_ok=True)
    temporary = dest.with_suffix(f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(document,ensure_ascii=False,allow_nan=False))
    temporary.replace(dest)
    return {k:v for k,v in document.items() if k != "targets"} | {"target_count":len(targets)}


def annotation_path(dataset_id):
    from .data import home
    if len(dataset_id)!=20 or any(c not in "0123456789abcdef" for c in dataset_id):
        raise ValueError("Invalid dataset ID")
    return home()/"annotations"/f"{dataset_id}.json"


@lru_cache(maxsize=4)
def _read_annotation_file(path, modified, size):
    return json.loads(Path(path).read_text())


def load_annotations(dataset_id):
    path = annotation_path(dataset_id)
    try:
        status = path.stat()
    except FileNotFoundError:
        return None
    return _read_annotation_file(str(path.resolve()), status.st_mtime_ns, status.st_size)


def target_info(target, annotations=None):
    metadata = {"family":None,"families":[],"protein_classes":[],"organism":None,"tax_id":None,
                **{k:target[k] for k in TARGET_FIELDS if k in target}}
    if annotations:
        metadata.update(annotations["targets"].get(target["id"],{}))
    return metadata


def enrich_results(data):
    annotations = load_annotations(data["dataset"]["id"])
    for result in data["results"]:
        if result.get("target"):
            result["target"] = target_info(result["target"], annotations)
    if annotations:
        data["target_annotations"] = {k:v for k,v in annotations.items() if k != "targets"}
    return data
