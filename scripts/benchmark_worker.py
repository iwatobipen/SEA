"""Compare scalar and optimized scoring without creating jobs or lab cache files.

Run with: pixi run python scripts/benchmark_worker.py --dataset-id DATASET_ID
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import tempfile
import time

import numpy as np
from rdkit import rdBase

from sea.background import NullCache
from sea.chemistry import normalize, raw_score
from sea.data import load_dataset
from sea.statistics import significance


def benchmark(dataset_id, smiles):
    started = time.perf_counter()
    data, molecules = load_dataset(dataset_id)
    loading = time.perf_counter() - started
    query, records = normalize(smiles)
    pool = list(molecules.values())
    targets = [(t, [molecules[s] for s in t['smiles']]) for t in data['targets']]

    def run(label, cache):
        print(f"Measuring {label}...", flush=True)
        started = time.perf_counter()
        results = []
        for target, ligands in targets:
            score, count, pairs = raw_score(query, ligands, data['model']['threshold'])
            stat = significance(score, len(query), len(ligands), pool, data['model'], len(targets), null_cache=cache)
            results.append({'target_id': target['id'], 'raw_score': score,
                            'contributing_pairs': count, 'similar_ligands': pairs, **stat})
        elapsed = time.perf_counter() - started
        print(f"{label}: {elapsed:.3f} s", flush=True)
        return results, elapsed

    # The ordinary dict selects the unchanged scalar reference sampler.
    reference_cache = {}
    reference, reference_seconds = run('scalar reference, cold', reference_cache)
    Path('.cache').mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='sea-benchmark-', dir='.cache') as temporary:
        cache = NullCache(temporary, dataset_id)
        measurements = []
        for label in ('cold', 'memory', 'disk'):
            if label == 'disk':
                cache = NullCache(temporary, dataset_id)
            results, seconds = run(f'optimized, {label}', cache)
            if results != reference:
                raise AssertionError(f'Optimized {label} results differ from the scalar reference.')
            if set(cache) != set(reference_cache) or any(
                not np.array_equal(values, cache[key]) for key, values in reference_cache.items()
            ):
                raise AssertionError(f'Optimized {label} background samples differ from the scalar reference.')
            measurements.append({'cache': label, 'seconds': seconds,
                                 'speedup_over_scalar_cold': reference_seconds / seconds,
                                 'exact_results': True, 'exact_background_samples': True})
    return {'dataset_id': dataset_id, 'model_id': data['manifest']['model_id'],
            'rdkit_version': rdBase.rdkitVersion, 'numpy_version': np.__version__,
            'query': records, 'target_count': len(targets), 'molecule_count': len(pool),
            'threshold': data['model']['threshold'], 'empirical_trials': data['model']['empirical_trials'],
            'background_distributions': len(reference_cache),
            'dataset_load_seconds': loading, 'scalar_cold_seconds': reference_seconds,
            'methods': dict(Counter(r['method'] for r in reference)), 'measurements': measurements,
            'scope': 'Scoring and significance only; excludes queueing, API validation, loading and job database writes.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True, help='Exact dataset ID printed by sea datasets.')
    parser.add_argument('--smiles', action='append', help='Repeat to benchmark one ensemble. Defaults to the example query.')
    parser.add_argument('--output', type=Path, help='Optional JSON report path.')
    args = parser.parse_args()
    report = benchmark(args.dataset_id, args.smiles or ['COc1ccc2[nH]cc(c2c1)C(CC)CN'])
    text = json.dumps(report, indent=2, allow_nan=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + '\n', encoding='utf-8')
    print(text)
