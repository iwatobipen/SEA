from dataclasses import dataclass
from rdkit import Chem, DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize

FP_SPEC = {"kind": "RDKit-path", "bits": 2048, "min_path": 1, "max_path": 7}
STANDARDIZATION = "cleanup-fragment-parent-uncharge-canonical-isomeric-v1"


class InputError(ValueError):
    pass


@dataclass
class Molecule:
    smiles: str
    fp: object


def molecule(smiles: str) -> Molecule:
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > 10000:
        raise InputError("空または長すぎるSMILESです")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise InputError("不正なSMILESです")
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.FragmentParent(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    Chem.SanitizeMol(mol)
    if not any(a.GetAtomicNum() == 6 for a in mol.GetAtoms()):
        raise InputError("有機分子を入力してください")
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    return Molecule(canonical, Chem.RDKFingerprint(mol, minPath=1, maxPath=7, fpSize=2048))


def normalize(smiles: list[str]):
    if not smiles:
        raise InputError("集合が空です")
    unique, rows, errors = {}, [], []
    for i, value in enumerate(smiles, 1):
        try:
            item = molecule(value)
            duplicate = item.smiles in unique
            unique.setdefault(item.smiles, item)
            rows.append({"row": i, "input": value, "smiles": item.smiles, "duplicate": duplicate})
        except Exception as exc:
            errors.append({"row": i, "input": value, "message": str(exc)})
    if errors:
        raise InputError(errors)
    return list(unique.values()), rows


def raw_score(left, right, threshold, evidence_limit=20, check=lambda: None):
    score, count, evidence = 0.0, 0, []
    right_fps = [m.fp for m in right]
    for a in left:
        check()
        sims = DataStructs.BulkTanimotoSimilarity(a.fp, right_fps)
        for b, sim in zip(right, sims):
            if sim >= threshold:
                score += sim
                count += 1
                if evidence_limit:
                    evidence.append({"query_smiles": a.smiles, "reference_smiles": b.smiles, "tanimoto": sim})
        if evidence_limit and len(evidence) > evidence_limit * 2:
            evidence = sorted(evidence, key=lambda x: (-x["tanimoto"], x["query_smiles"], x["reference_smiles"]))[:evidence_limit]
    evidence.sort(key=lambda x: (-x["tanimoto"], x["query_smiles"], x["reference_smiles"]))
    return score, count, evidence[:evidence_limit]
