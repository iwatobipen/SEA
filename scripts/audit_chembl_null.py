"""Independent uniform-null audit of a dataset's deployed Gumbel approximation."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy import stats

from sea.background import NullCache
from sea.data import load_dataset


def audit(dataset_id, directory, trials=8192, seed=20261109):
    started = time.perf_counter()
    data, structures = load_dataset(dataset_id)
    model = data['model']
    if not model['evd_usable'] or model['std_b'] <= 0:
        raise ValueError('This dataset does not enable a fitted Gumbel tail. Use the empirical-null audit instead.')
    pool = list(structures.values())
    sizes = [s for s in model['sizes'] if s >= 10]
    pairs = [(a, b) for a in sizes for b in sizes if a <= b]
    levels = [.1, .05, .01, .001]
    identity = hashlib.sha256(json.dumps({'model': model, 'seed': seed, 'trials': trials}, sort_keys=True).encode()).hexdigest()[:16]
    cache = NullCache(Path(directory) / identity, dataset_id)
    report = {'dataset_id': dataset_id, 'model_id': data['manifest']['model_id'],
              'threshold': model['threshold'], 'fingerprint': data['manifest']['fingerprint'],
              'original_seed': model['seed'], 'audit_seed': seed, 'trials_per_pair': trials,
              'declared_levels': levels, 'family_tests': len(pairs) * len(levels),
              'family_alpha': .05, 'null': 'independent uniform molecule sets; no replacement within sets; overlap across sets permitted',
              'pairs': []}
    for a, b in pairs:
        null = cache.sample(pool, {**model, 'seed': seed, 'empirical_trials': trials}, (a, b), lambda: None)
        sigma = model['std_b'] * (a * b)**model['std_c']
        z = (null - model['mean_a'] * a * b) / sigma
        p = stats.gumbel_r.sf((z - model['gumbel_loc']) / model['gumbel_scale'])
        p[null == 0] = 1.
        tails = []
        for level in levels:
            count = int(np.sum(p <= level))
            lower = 0. if count == 0 else float(stats.beta.ppf(.025, count, trials-count+1))
            upper = 1. if count == trials else float(stats.beta.ppf(.975, count+1, trials-count))
            test = float(stats.binomtest(count, trials, level, alternative='greater').pvalue)
            tails.append({'nominal_p': level, 'rejections': count, 'rate': count / trials,
                          'rate_pointwise_95_interval': [lower, upper], 'excess_test_p': test,
                          'excess_bonferroni_p': min(1., test * len(pairs) * len(levels))})
        record = {'sizes': [a, b], 'zero_fraction': float(np.mean(null == 0)),
                  'sample_mean': float(null.mean()), 'model_mean': model['mean_a'] * a * b,
                  'sample_sd': float(null.std(ddof=1)), 'model_sd': sigma,
                  'continuous_model_ks': float(stats.kstest(z, 'gumbel_r', args=(model['gumbel_loc'], model['gumbel_scale'])).statistic),
                  'tail_checks': tails}
        report['pairs'].append(record)
        print(json.dumps({'sizes': [a, b], 'ks': record['continuous_model_ks'],
                          'rates': {str(t['nominal_p']): t['rate'] for t in tails}}), flush=True)
    report['anti_conservative_detected'] = any(t['excess_bonferroni_p'] < .05 for r in report['pairs'] for t in r['tail_checks'])
    report['elapsed_seconds'] = time.perf_counter() - started
    report['interpretation'] = 'Independent simulations under the implemented uniform-set null, not biological inactive compounds. One-sided binomial excess tests use Bonferroni over declared size/level checks. Finite simulations do not validate arbitrarily small tails or all intermediate set sizes.'
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--trials', type=int, default=8192)
    parser.add_argument('--seed', type=int, default=20261109)
    args = parser.parse_args()
    if args.trials < 1000:
        parser.error('Use at least 1,000 independent trials per size pair.')
    result = audit(args.dataset_id, args.cache, args.trials, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'anti_conservative_detected': result['anti_conservative_detected'], 'elapsed_seconds': result['elapsed_seconds']}))
