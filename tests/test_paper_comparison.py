"""Offline checks for the fingerprint diagnosis and statistical replay tool."""
import importlib.util
import math
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/compare_paper.py"
spec = importlib.util.spec_from_file_location("paper_comparison", SCRIPT)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_reported_pair_depends_on_fingerprint():
    result = comparison.fingerprint_diagnosis()
    assert result["standardization_changed_structure"] == [False, False]
    path = result["fingerprints"][0]
    assert (path["intersection_bits"], path["union_bits"]) == (506, 616)
    assert path["tanimoto"] == pytest.approx(506 / 616)
    circular = [r for r in result["fingerprints"]
                if r["spec"]["kind"] == "Morgan" and r["spec"]["radius"] == 2]
    assert len(circular) == 3
    for row in circular:
        assert (row["intersection_bits"], row["union_bits"]) == (24, 44)
        assert row["tanimoto"] == pytest.approx(.5454545454545454)
    for row in result["fingerprints"]:
        assert row["tanimoto"] == pytest.approx(row["independent_tanimoto"])


def test_independent_aggregation_and_extreme_tail():
    for row in comparison.raw_score_check():
        assert row["absolute_error"] < 1e-12
        assert row["count_matches"]
    for row in comparison.independent_tail_check():
        assert math.isfinite(row["production_log_p"])
        assert row["absolute_log_error"] < 2e-12


def test_replay_uses_published_scores_and_excludes_self_matches():
    # The DHFR/TS published row; no chemical structures are available in the XLS.
    classes = {78389: {"size": 216, "name": "DHFR"},
               78362: {"size": 253, "name": "TS"},
               1: {"size": 20, "name": "Control"}}
    records = []
    for q in classes:
        for t in classes:
            raw = 772.25 if {q, t} == {78389, 78362} else (1000 if q == t else 1.)
            n = classes[q]["size"] * classes[t]["size"]
            z = (raw - .000424*n) / (.00449*n**.665)
            if {q, t} == {78389, 78362}:
                z = 117.57
            e = math.exp(comparison.original_tail(z)[1] + math.log(9))
            records.append({"query": q, "reference": t, "raw": raw, "z": z, "e": e})
    result = comparison.replay(classes, records)
    example = result["dhfr_ts_example"][0]
    assert example["replayed_z"] == pytest.approx(117.97623538114273)
    assert result["ranking_replay"]["self_matches_excluded"]
    assert result["ranking_replay"]["top_1_agreement"] == 1
    assert result["production_gumbel_vs_original_tail_log_e_error"]["max"] < 2e-12
