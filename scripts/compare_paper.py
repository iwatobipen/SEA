"""Reproducible fingerprint diagnosis and replay of the 2007 SEA supplement.

Run in the optional Pixi paper-validation environment. No service data is changed.
The published raw scores are inputs, not predictions reproduced from structures.
"""
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import zipfile

import httpx
import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from scipy import stats
import scipy

from sea.chemistry import FP_SPEC, STANDARDIZATION, Molecule, molecule, raw_score
from sea.statistics import original_tail, significance


QUERY = "COc1ccc2[nH]cc(c2c1)C(CC)CN"
REFERENCE = "COc1ccc2[nH]cc(CCN)c2c1"
SUPPLEMENT_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fnbt1284/"
    "MediaObjects/41587_2007_BFnbt1284_MOESM65_ESM.zip"
)
SUPPLEMENT_SHA256 = "1fa3c74bb4a5340d002ee712da815331152d58164d78e4e341111036d4e3ce1e"
WORKBOOK_NAME = "Shoichet Supplementary Data.xls"


def bit_comparison(a, b):
    """Independent binary-set definition of Tanimoto alongside RDKit's result."""
    left, right = set(a.GetOnBits()), set(b.GetOnBits())
    intersection, union = len(left & right), len(left | right)
    return {"query_bits": len(left), "reference_bits": len(right),
            "intersection_bits": intersection, "union_bits": union,
            "tanimoto": DataStructs.TanimotoSimilarity(a, b),
            "independent_tanimoto": intersection / union if union else 1.0}


def fingerprint_diagnosis():
    a, b = molecule(QUERY), molecule(REFERENCE)
    mols = [Chem.MolFromSmiles(m.smiles) for m in (a, b)]
    records = [{"spec": FP_SPEC, **bit_comparison(a.fp, b.fp)}]
    for radius in (1, 2, 3):
        for bits in (1024, 2048, 4096):
            generator = rdFingerprintGenerator.GetMorganGenerator(
                radius=radius, fpSize=bits, includeChirality=False, countSimulation=False,
                useBondTypes=True, includeRingMembership=True,
                includeRedundantEnvironments=False,
            )
            fps = [generator.GetFingerprint(m) for m in mols]
            records.append({"spec": {"kind": "Morgan", "radius": radius,
                                      "bits": bits, "chirality": False,
                                      "count_simulation": False}, **bit_comparison(*fps)})
    return {"query": QUERY, "reference": REFERENCE, "standardization": STANDARDIZATION,
            "canonical_structures": [a.smiles, b.smiles],
            "standardization_changed_structure": [
                Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=True) != m.smiles
                for s, m in zip((QUERY, REFERENCE), (a, b))],
            "fingerprints": records}


def raw_score_check():
    """Check production aggregation against explicit bit sets at several cutoffs."""
    inputs = [QUERY, REFERENCE, "CCO", "CCCO", "c1ccccc1O", "c1ccccc1N"]
    molecules = [molecule(s) for s in inputs]
    circular = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    checks = []
    for kind in ("RDKit-path", "Morgan-radius-2"):
        items = molecules if kind == "RDKit-path" else [
            Molecule(m.smiles, circular.GetFingerprint(Chem.MolFromSmiles(m.smiles)))
            for m in molecules]
        left, right = items[:4], items[1:]
        for threshold in (.30, .57, .80, 1.0):
            values = [bit_comparison(a.fp, b.fp)["independent_tanimoto"]
                      for a in left for b in right]
            contributing = [v for v in values if v >= threshold]
            expected = math.fsum(contributing)
            actual, count, _ = raw_score(left, right, threshold, 0)
            checks.append({"kind": kind, "threshold": threshold,
                           "production_score": actual, "independent_score": expected,
                           "absolute_error": abs(actual - expected),
                           "contributing_pairs": count,
                           "count_matches": count == len(contributing)})
    return checks


def independent_tail_check():
    """High-precision reference using Decimal exp/ln, without the production tail."""
    values = []
    with localcontext() as ctx:
        ctx.prec = 100
        pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286")
        gamma = Decimal("0.577215664901532860606512090082402431042159335939923598805767234884867726778")
        for z in (0, 20, 28, 119, 578, 1000):
            log_t = -Decimal(z) * pi / Decimal(6).sqrt() - gamma
            t = log_t.exp()
            # 1-exp(-t) loses precision for tiny t even in finite Decimal arithmetic.
            # ln(p) = ln(t) + ln(1-t/2+t^2/6-...), which remains stable.
            if t < Decimal("1e-20"):
                log_p = log_t + (1 - t / 2 + t * t / 6).ln()
            else:
                log_p = (1 - (-t).exp()).ln()
            _, actual = original_tail(z)
            values.append({"z": z, "production_log_p": actual,
                           "decimal_log_p": str(log_p),
                           "absolute_log_error": abs(actual - float(log_p))})
    return values


def read_workbook(content):
    try:
        import xlrd
    except ImportError as exc:
        raise SystemExit("Use: pixi run -e paper-validation python scripts/compare_paper.py ...") from exc
    digest = hashlib.sha256(content).hexdigest()
    if digest != SUPPLEMENT_SHA256:
        raise ValueError(f"Supplement ZIP checksum mismatch: {digest}")
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        workbook = xlrd.open_workbook(file_contents=archive.read(WORKBOOK_NAME))
    sheet = workbook.sheet_by_name("Sheet2")
    classes = {int(sheet.cell_value(i, 1)): {
        "size": int(sheet.cell_value(i, 2)), "name": sheet.cell_value(i, 3)}
        for i in range(4, sheet.nrows)}
    sheet = workbook.sheet_by_name("Sheet3")
    records = []
    for i in range(4, sheet.nrows):
        q, t, e, z, raw, identities, maximum = sheet.row_values(i)[1:]
        records.append({"query": int(q), "reference": int(t), "e": float(e),
                        "z": z, "raw": raw, "identities": identities, "max_tc": maximum})
    if len(classes) != 246 or len(records) != 246**2:
        raise ValueError("Unexpected supplement dimensions")
    pairs = {(r["query"], r["reference"]) for r in records}
    if pairs != {(a, b) for a in classes for b in classes}:
        raise ValueError("Missing or duplicate class comparisons")
    return classes, records


def summary(values):
    return {"count": len(values), **{
        label: float(np.quantile(values, quantile)) if values else None
        for label, quantile in (("median", .5), ("p95", .95), ("p99", .99), ("max", 1))}}


def replay(classes, records):
    comparisons = len(classes)**2
    factor = math.sqrt(6) / math.pi
    # This special model checks the production Gumbel branch with the original
    # standardized distribution. It is not the active ChEMBL calibration model.
    model = {"threshold": .57, "mean_a": .000424, "std_b": .00449, "std_c": .665,
             "gumbel_loc": -float(np.euler_gamma)*factor, "gumbel_scale": factor,
             "evd_usable": True, "sizes": [1, max(c["size"] for c in classes.values())]}
    errors, z_errors, production_errors = [], [], []
    regions = defaultdict(list)
    examples, by_query, zeros = [], defaultdict(list), []
    for row in records:
        a, b = classes[row["query"]], classes[row["reference"]]
        n = a["size"] * b["size"]
        z = (row["raw"] - model["mean_a"]*n) / (model["std_b"]*n**model["std_c"])
        row = {**row, "replayed_z": z}
        by_query[row["query"]].append(row)
        z_errors.append(abs(z - row["z"]))
        _, log_p = original_tail(row["z"])
        log_e = log_p + math.log(comparisons)
        if row["e"] > 0:
            error = abs(math.expm1(log_e - math.log(row["e"])))
            errors.append(error)
            region = "z_below_20" if row["z"] < 20 else (
                "z_20_to_30" if row["z"] < 30 else "z_at_least_30")
            regions[region].append(error)
        else:
            zeros.append(log_e)
        if min(a["size"], b["size"]) >= 10 and row["raw"] > 0:
            production = significance(row["raw"], a["size"], b["size"], [], model, comparisons)
            expected_log_e = original_tail(z)[1] + math.log(comparisons)
            production_errors.append(abs(production["log_e_value"] - expected_log_e))
        if row["query"] == 78389 and row["reference"] == 78362:
            examples.append({**row, "query_name": a["name"], "reference_name": b["name"],
                             "e_from_published_z": math.exp(log_e),
                             "e_from_rounded_coefficients": math.exp(original_tail(z)[1] + math.log(comparisons))})
    agreements = {"top_1_agreement": [], "mean_top_5_set_overlap": [],
                  "mean_top_10_set_overlap": [], "mean_spearman": []}
    for query, rows in by_query.items():
        rows = [r for r in rows if r["reference"] != query]
        # Published E=0 values are censored. Z retains their relative ordering.
        original = sorted(rows, key=lambda r: (-r["z"], r["reference"]))
        computed = sorted(rows, key=lambda r: (-r["replayed_z"], r["reference"]))
        for k in (1, 5, 10):
            overlap = {r["reference"] for r in original[:k]} & {r["reference"] for r in computed[:k]}
            metric = "top_1_agreement" if k == 1 else f"mean_top_{k}_set_overlap"
            agreements[metric].append(len(overlap)/min(k, len(rows)))
        agreements["mean_spearman"].append(float(stats.spearmanr(
            [r["z"] for r in rows], [r["replayed_z"] for r in rows]).statistic))
    return {"classes": len(classes), "comparisons": len(records),
            "published_z_to_e_relative_error": summary(errors),
            "published_z_to_e_by_z_region": {k: summary(v) for k, v in regions.items()},
            "nonzero_e_within_1_percent": sum(e <= .01 for e in errors),
            "published_zero_e_count": len(zeros),
            "published_zero_e_stable_log_e_range": [min(zeros), max(zeros)] if zeros else None,
            "rounded_raw_to_z_absolute_error": summary(z_errors),
            "production_gumbel_vs_original_tail_log_e_error": summary(production_errors),
            "ranking_replay": {"query_classes": len(by_query), "self_matches_excluded": True,
                "key": "descending Z, reference class ID for ties; published raw scores supplied",
                **{k: float(np.mean(v)) for k, v in agreements.items()}},
            "dhfr_ts_example": examples,
            "interpretation": "Statistical replay from published scores; not end-to-end biological replication."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplement", type=Path, default=Path(".cache/paper-comparison/data.zip"))
    parser.add_argument("--download", action="store_true", help="Download the pinned official ZIP if absent")
    parser.add_argument("--output", type=Path, default=Path("var/paper-comparison.json"))
    args = parser.parse_args()
    if not args.supplement.exists():
        if not args.download:
            parser.error("Supplement ZIP missing; use --download or supply --supplement PATH")
        response = httpx.get(SUPPLEMENT_URL, follow_redirects=True, timeout=120)
        response.raise_for_status()
        content = response.content
        # Verify before writing an untrusted response to the cache.
        if hashlib.sha256(content).hexdigest() != SUPPLEMENT_SHA256:
            raise ValueError("Downloaded supplement ZIP checksum mismatch")
        args.supplement.parent.mkdir(parents=True, exist_ok=True)
        args.supplement.write_bytes(content)
    classes, records = read_workbook(args.supplement.read_bytes())
    result = {"schema": 1, "environment": {"python": platform.python_version(),
                "rdkit": rdBase.rdkitVersion, "numpy": np.__version__, "scipy": scipy.__version__},
              "source": {"paper": "https://doi.org/10.1038/nbt1284", "url": SUPPLEMENT_URL,
                         "sha256": SUPPLEMENT_SHA256, "workbook": WORKBOOK_NAME},
              "fingerprint_diagnosis": fingerprint_diagnosis(),
              "independent_raw_score_checks": raw_score_check(),
              "independent_decimal_tail_checks": independent_tail_check(),
              "paper_statistical_replay": replay(classes, records)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), **result["paper_statistical_replay"]}, indent=2))


if __name__ == "__main__":
    main()
