"""Measure fixed-budget Monte Carlo precision without modifying application data.

Example: pixi run python scripts/study_precision.py --dataset-id DATASET_ID
"""
import argparse
import json
from pathlib import Path
import tempfile
import time

import numpy as np
from scipy.stats import beta

from sea.background import FingerprintIndex, NullCache
from sea.chemistry import normalize, raw_score
from sea.data import load_dataset


def study(dataset_id, query, target_ids, trials):
    started = time.perf_counter()
    data, structures = load_dataset(dataset_id)
    load_seconds = time.perf_counter() - started
    left, records = normalize([query])
    pool = list(structures.values())
    targets = {target['id']: target for target in data['targets']}
    selected = [targets[target_id] for target_id in target_ids]
    started = time.perf_counter()
    index = FingerprintIndex(pool)
    index_seconds = time.perf_counter() - started
    measurements = []
    Path('.cache').mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='sea-precision-', dir='.cache') as temporary:
        for target in selected:
            right = [structures[s] for s in target['smiles']]
            score, _, _ = raw_score(left, right, data['model']['threshold'])
            for count in trials:
                print(f"Measuring {target['id']}, {count} trials...", flush=True)
                cache = NullCache(Path(temporary) / str(count), dataset_id)
                cache.index = index
                started = time.perf_counter()
                null = cache.sample(pool, {**data['model'], 'empirical_trials': count},
                                    (1, len(right)), lambda: None)
                seconds = time.perf_counter() - started
                exceed = count - int(np.searchsorted(null, score, side='left'))
                p = (exceed + 1) / (count + 1)
                # Pointwise fixed-budget intervals for the limiting null tail.
                # These are not simultaneous or anytime-valid intervals.
                lower = 0. if exceed == 0 else float(beta.ppf(.025, exceed, count - exceed + 1))
                upper = 1. if exceed == count else float(beta.ppf(.975, exceed + 1, count - exceed))
                measurements.append({
                    'target_id': target['id'], 'target_name': target['name'],
                    'reference_size': len(right), 'raw_score': score,
                    'trials': count, 'exceedances': exceed, 'seconds': seconds,
                    'p_value': p, 'e_value': p * len(targets),
                    'tail_probability_pointwise_95_interval': [lower, upper],
                })
                print(f"E={p * len(targets):.6g}; {seconds:.3f} s", flush=True)
    sizes = {len(t['smiles']) for t in targets.values()}
    return {'dataset_id': dataset_id, 'model_id': data['manifest']['model_id'],
            'rdkit_version': data['manifest']['rdkit_version'], 'numpy_version': np.__version__,
            'query': records, 'target_count': len(targets), 'unique_reference_sizes': len(sizes),
            'threshold': data['model']['threshold'], 'seed': data['model']['seed'],
            'load_seconds': load_seconds, 'packed_index_seconds': index_seconds,
            'measurements': measurements,
            'scope': 'Sampling timings share one packed index; exclude loading, index construction, API and job queue. The seed/size streams are shared between budgets, so budgets are not independent replicates.',
            'interpretation': 'Fixed-budget Monte Carlo precision/cost study, not biological validation. Tail intervals are pointwise and cannot justify optional stopping or simultaneous coverage.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True)
    parser.add_argument('--smiles', default='COc1ccc2[nH]cc(c2c1)C(CC)CN')
    parser.add_argument('--target-id', action='append', help='Repeat to measure several targets.')
    parser.add_argument('--trials', type=int, nargs='+', default=[999, 9999, 99999])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if any(n < 1 for n in args.trials):
        parser.error('Trial counts must be positive.')
    report = study(args.dataset_id, args.smiles, args.target_id or ['CHEMBL5043', 'CHEMBL1946'], args.trials)
    text = json.dumps(report, indent=2, allow_nan=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + '\n', encoding='utf-8')
    print(text)
