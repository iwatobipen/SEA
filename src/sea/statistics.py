"""SEA null calibration. Original MDDR constants are deliberately not defaults."""
import math
import numpy as np
from scipy import stats
from rdkit import DataStructs
from .chemistry import raw_score


def original_tail(z):
    """Standardized Gumbel survival in the original supplementary methods."""
    log_t = -z * math.pi / math.sqrt(6) - float(np.euler_gamma)
    if log_t > 40:
        return 1.0, 0.0
    if log_t < -35:
        return math.exp(log_t), log_t
    p = -math.expm1(-math.exp(log_t))
    return p, math.log(p)


def calibrate(pool, sizes=(1, 10, 30, 100), repeats=128, seed=42, progress=lambda *_: None):
    """Fit mu=a*n, sigma=b*n**c; select threshold by held-out CDF error."""
    sizes = sorted(set(s for s in sizes if s <= len(pool)))
    if len(sizes) < 2 or repeats < 16:
        raise ValueError("Calibration requires two sizes and at least 16 repeats")
    rng = np.random.default_rng(seed)
    thresholds = np.round(np.arange(.30, .801, .01), 2)
    products, blocks = [], []
    pairs = [(a, b) for a in sizes for b in sizes if a <= b]
    for step, (a, b) in enumerate(pairs):
        scores = []
        for _ in range(repeats):
            left = rng.choice(len(pool), a, replace=False)
            right = rng.choice(len(pool), b, replace=False)
            fps = [pool[i].fp for i in right]
            values = np.sort(np.concatenate([DataStructs.BulkTanimotoSimilarity(pool[i].fp, fps) for i in left]))
            cumulative = np.concatenate(([0.0], np.cumsum(values)))
            scores.append(cumulative[-1] - cumulative[np.searchsorted(values, thresholds, side="left")])
        blocks.append(np.array(scores))
        products.append(a*b)
        progress(step + 1, len(pairs))
    products = np.array(products, dtype=float)
    split = repeats // 2
    best = None
    diagnostics = []
    for k, threshold in enumerate(thresholds):
        train = np.array([b[:split, k] for b in blocks])
        means, stds = train.mean(axis=1), train.std(axis=1, ddof=1)
        # Singleton comparisons have a discrete distribution; do not fit them as EVD.
        valid = (stds > 1e-9) & np.array([min(a,b) >= 10 for a,b in pairs])
        if valid.sum() < 3:
            continue
        mean_a = float(np.dot(products[valid], means[valid]) / np.dot(products[valid], products[valid]))
        exponent, intercept = np.polyfit(np.log(products[valid]), np.log(stds[valid]), 1)
        sigma = np.exp(intercept) * products**exponent
        ztrain = ((train[valid] - mean_a*products[valid, None]) / sigma[valid, None]).ravel()
        loc, scale = stats.gumbel_r.fit(ztrain)
        validation = np.array([b[split:, k] for b in blocks])
        zvalid = ((validation[valid] - mean_a*products[valid, None]) / sigma[valid, None]).ravel()
        ks = float(stats.kstest(zvalid, "gumbel_r", args=(loc, scale)).statistic)
        candidate = {"threshold": float(threshold), "mean_a": mean_a,
                     "std_b": float(np.exp(intercept)), "std_c": float(exponent),
                     "gumbel_loc": float(loc), "gumbel_scale": float(scale), "validation_ks": ks}
        diagnostics.append({"threshold": float(threshold), "validation_ks": ks})
        if best is None or ks < best["validation_ks"]:
            best = candidate
    if best is None:
        best = {"threshold": .57, "mean_a": 0., "std_b": 0., "std_c": .5,
                "gumbel_loc": 0., "gumbel_scale": 1., "validation_ks": None}
    best.update({"seed": seed, "sizes": sizes, "repeats": repeats,
                 "evd_usable": best["validation_ks"] is not None and best["validation_ks"] <= .10,
                 "diagnostics": diagnostics, "empirical_trials": 999})
    return best


def significance(score, left_size, right_size, pool, model, comparisons=1, check=lambda: None, null_cache=None):
    n = left_size*right_size
    z = None
    sigma = model["std_b"] * n**model["std_c"]
    if sigma > 0:
        z = float((score - model["mean_a"]*n)/sigma)
    eligible = (model["evd_usable"] and min(left_size, right_size) >= 10
                and min(model["sizes"]) <= min(left_size, right_size)
                and max(left_size, right_size) <= max(model["sizes"]))
    resolution = None
    if score == 0:
        # Every null raw score is nonnegative, so P(null >= 0) is exactly 1.
        # This avoids sampling for targets with no contributing ligand pairs.
        p, log_p, method = 1.0, 0.0, "exact_zero"
    elif eligible and z is not None:
        normalized = (z-model["gumbel_loc"])/model["gumbel_scale"]
        # log survival is asymptotically -normalized; avoid subtraction cancellation.
        log_p = -normalized if normalized > 35 else float(stats.gumbel_r.logsf(normalized))
        p = math.exp(log_p)
        method = "gumbel"
    else:
        if max(left_size, right_size) > len(pool):
            raise ValueError("比較集合が背景分子数を超えています。より大きい背景DBを構築してください")
        trials = model["empirical_trials"]
        rng = np.random.default_rng(np.random.SeedSequence([model["seed"], min(left_size,right_size), max(left_size,right_size)]))
        key = (min(left_size,right_size), max(left_size,right_size))
        cached = null_cache.get(key) if null_cache is not None else None
        if cached is None:
            values = []
            for _ in range(trials):
                check()
                a = [pool[i] for i in rng.choice(len(pool), key[0], replace=False)]
                b = [pool[i] for i in rng.choice(len(pool), key[1], replace=False)]
                null, _, _ = raw_score(a, b, model["threshold"], 0, check)
                values.append(null)
            cached = np.sort(values)
            if null_cache is not None:
                null_cache[key] = cached
        exceed = trials - int(np.searchsorted(cached, score, side="left"))
        p = (exceed+1)/(trials+1)
        log_p = math.log(p)
        resolution = 1/(trials+1)
        method = "monte_carlo"
    return {"z_score": z, "p_value": p, "log_p_value": log_p,
            "e_value": p*comparisons, "log_e_value": log_p+math.log(comparisons),
            "comparisons": comparisons, "method": method, "p_resolution": resolution}
