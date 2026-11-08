"""Check empirical SEA tails using independent reference and query simulations."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from rdkit import DataStructs
from scipy import stats

from sea.background import NullCache
from sea.chemistry import MORGAN_SPEC, FP_SPEC, Molecule
from sea.data import load_dataset
from sea.parameter_validation import predictive_excess


def run(dataset_id, directory, representation, threshold):
    started = time.perf_counter()
    directory = Path(directory)
    data, structures = load_dataset(dataset_id)
    if representation == 'Morgan':
        packed = np.load(directory / 'morgan.npy', mmap_mode='r')
        pool = [Molecule(s, DataStructs.CreateFromBinaryText(row.tobytes())) for s, row in zip(structures, packed)]
        spec = MORGAN_SPEC
    else:
        pool, spec = list(structures.values()), FP_SPEC
    model = {**data['model'], 'threshold': threshold}
    identity = hashlib.sha256(json.dumps({'dataset_id': dataset_id, 'fingerprint': spec, 'threshold': threshold}, sort_keys=True).encode()).hexdigest()[:20]
    reference = NullCache(directory / 'empirical-audit' / identity / 'reference', dataset_id)
    observed = NullCache(directory / 'empirical-audit' / identity / 'observed', dataset_id)
    pairs = [(10, 10), (10, 30), (10, 100), (30, 30), (30, 100), (100, 100)]
    levels, b, m = [.1, .05, .01, .001], 8191, 8192
    output = []
    for a, c in pairs:
        reference.index = observed.index or reference.index
        background = reference.sample(pool, {**model, 'seed': 20261111, 'empirical_trials': b}, (a, c), lambda: None)
        observed.index = reference.index
        scores = observed.sample(pool, {**model, 'seed': 20261109, 'empirical_trials': m}, (a, c), lambda: None)
        p = (b - np.searchsorted(background, scores, side='left') + 1) / (b + 1)
        checks = []
        for level in levels:
            count = int(np.sum(p <= level))
            excess = predictive_excess(count, m, b, level)
            checks.append({'nominal_p': level, 'rejections': count, 'rate': count / m,
                           'predictive_excess_p': excess, 'bonferroni_p': min(1., excess * len(pairs) * len(levels))})
        output.append({'sizes': [a, c], 'zero_fraction': float(np.mean(scores == 0)), 'tail_checks': checks})
        print(json.dumps({'sizes': [a, c], 'rates': {str(t['nominal_p']): t['rate'] for t in checks}}), flush=True)
    return {'dataset_id': dataset_id, 'fingerprint': spec, 'threshold': threshold,
            'reference_trials': b, 'observed_trials': m, 'reference_seed': 20261111, 'observed_seed': 20261109,
            'pairs': output, 'anti_conservative_detected': any(t['bonferroni_p'] < .05 for row in output for t in row['tail_checks']),
            'elapsed_seconds': time.perf_counter()-started,
            'test': 'One-sided beta-binomial predictive tail accounts for the shared finite reference CDF. Bonferroni over 24 declared size/level checks. Ties use inclusive >= and remain conservative.',
            'limits': 'Uniform-set null only, not measured inactive compounds. Failure to detect excess does not prove biological calibration. Fixed-budget corrected Monte Carlo rank p-values also rely on exchangeability under the specified null.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--fingerprint', choices=['Morgan', 'RDKit-path'], required=True)
    parser.add_argument('--threshold', type=float, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.dataset_id, args.cache, args.fingerprint, args.threshold)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'anti_conservative_detected': result['anti_conservative_detected'], 'elapsed_seconds': result['elapsed_seconds']}))
