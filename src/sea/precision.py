"""Fixed-budget statistical profiles and pointwise simulation uncertainty."""
import hashlib
import json
from scipy.stats import beta

PROFILES = {"balanced": 199999, "research": 999999}
VERSION = "fixed-budget-v1"
NULL_DEFINITION = "independent-uniform-sets-without-replacement-v1"
SAMPLER_VERSION = "batched-exact-v1"
STREAM_POLICY = "pcg64-seed-sequence-seed-sorted-sizes-v1"


def profile(data, name):
    value = {"name": name, "version": VERSION, "planned_trials": PROFILES[name],
             "dataset_id": data["manifest"]["id"], "base_model_id": data["manifest"]["model_id"],
             "null_definition": NULL_DEFINITION, "seed": data["model"]["seed"],
             "sampler_version": SAMPLER_VERSION, "stream_policy": STREAM_POLICY,
             "testing_family": "per_query", "confidence_level": .95,
             "interval_method": "clopper-pearson-fixed-budget-pointwise"}
    value["id"] = hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:20]
    return value


def details(result, trial_count, phase, effective_profile):
    method = result["method"]
    info = {"phase": phase, "profile_id": effective_profile["id"],
            "planned_trials": effective_profile["planned_trials"], "trials": 0,
            "exceedances": None, "tail_probability_interval": None,
            "e_tail_interval": None, "confidence_level": .95,
            "interval_method": None, "decision_threshold_e": .05}
    if method == "pending":
        info["status"] = "pending"
    elif method == "monte_carlo":
        k = round(result["p_value"] * (trial_count + 1) - 1)
        lower = 0. if k == 0 else float(beta.ppf(.025, k, trial_count - k + 1))
        upper = 1. if k == trial_count else float(beta.ppf(.975, k + 1, trial_count - k))
        interval = [lower * result["comparisons"], upper * result["comparisons"]]
        info.update(trials=trial_count, exceedances=k,
                    tail_probability_interval=[lower, upper], e_tail_interval=interval,
                    interval_method="clopper-pearson-fixed-budget-pointwise")
        info["status"] = ("resolution_limited" if k == 0 else
                          "uncertain" if interval[0] <= .05 <= interval[1] else "evaluated")
    else:
        info["status"] = "exact" if method == "exact_zero" else "analytic"
        info["phase"] = "final"
    return info


def pending(comparisons):
    return {"z_score": None, "p_value": None, "log_p_value": None,
            "e_value": None, "log_e_value": None, "comparisons": comparisons,
            "method": "pending", "p_resolution": None}
