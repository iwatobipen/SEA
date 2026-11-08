"""Select SEA parameters on scaffold validation, then evaluate once on holdout.

Run prepare_parameter_validation.py first. Test results cannot select a second
candidate: a failed confirmation retains the existing model.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from rdkit import DataStructs, rdBase

from sea.data import load_dataset
from sea.parameter_validation import (PROTOCOL, adoption_passes, choose_queries,
    null_prefixes, paired_bootstrap, reference_split, retrieval, summarize)


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def run(dataset_id, directory, output):
    started = time.perf_counter()
    directory = Path(directory)
    protocol_file = directory / 'protocol.json'
    if not protocol_file.exists() or json.loads(protocol_file.read_text()) != PROTOCOL:
        raise ValueError('Write and review the fixed protocol before evaluating candidates.')
    data, structures = load_dataset(dataset_id)
    if data['manifest']['fingerprint']['kind'] != 'RDKit-path':
        raise ValueError('This fixed protocol compares an RDKit-path baseline with Morgan; choose the original path reference.')
    smiles = list(structures)
    identity = json.loads((directory / 'complete.json').read_text())
    assert identity['dataset_id'] == dataset_id
    assert identity['rdkit_version'] == rdBase.rdkitVersion
    assert identity['smiles_order_sha256'] == hashlib.sha256('\n'.join(smiles).encode()).hexdigest()
    groups = np.load(directory / 'groups.npy', mmap_mode='r')
    train, references, membership = reference_split(data, smiles, groups)
    selection_queries = choose_queries(smiles, groups, membership, references, 1, PROTOCOL['queries_per_partition'], PROTOCOL['seed'])
    test_queries = choose_queries(smiles, groups, membership, references, 0, PROTOCOL['queries_per_partition'], PROTOCOL['seed'])
    assert set(q['scaffold'] for q in selection_queries).isdisjoint(q['scaffold'] for q in test_queries)
    assert not set(groups[train]) & {q['scaffold'].encode() for q in selection_queries + test_queries}
    fingerprint_sets = {'RDKit-path': [m.fp for m in structures.values()]}
    morgan = np.load(directory / 'morgan.npy', mmap_mode='r')
    fingerprint_sets['Morgan'] = [DataStructs.CreateFromBinaryText(row.tobytes()) for row in morgan]
    train_hash = hashlib.sha256(train.tobytes()).hexdigest()[:12]
    sizes = [len(t['indices']) for t in references]
    report = {'dataset_id': dataset_id, 'base_model_id': data['manifest']['model_id'], 'protocol': PROTOCOL,
              'protocol_sha256': hashlib.sha256(json.dumps(PROTOCOL, sort_keys=True).encode()).hexdigest(),
              'rdkit_version': rdBase.rdkitVersion, 'source_threshold': data['model']['threshold'],
              'source_fingerprint': data['manifest']['fingerprint'],
              'training_molecules': len(train), 'training_targets': len(references),
              'all_molecules': len(smiles), 'all_targets': len(data['targets']),
              'partition_counts': {str(i): int(np.sum(membership == i)) for i in range(5)},
              'selection_queries': len(selection_queries), 'test_queries': len(test_queries),
              'unique_selection_scaffolds': len({q['scaffold'] for q in selection_queries}),
              'unique_test_scaffolds': len({q['scaffold'] for q in test_queries}),
              'training_indices_sha256': hashlib.sha256(train.tobytes()).hexdigest(),
              'queries_sha256': {phase: hashlib.sha256(json.dumps(queries, sort_keys=True).encode()).hexdigest()
                                 for phase, queries in [('selection', selection_queries), ('test', test_queries)]},
              'candidates': [], 'status': 'selecting'}
    write_json(output, report)
    print(f"Split: {len(train):,} training molecules; {len(references)} targets; {len(selection_queries)}/{len(test_queries)} selection/test queries.", flush=True)
    validation_rows = {}
    for name, fingerprints in fingerprint_sets.items():
        training = [fingerprints[i] for i in train]
        null, actual_sizes = null_prefixes(training, sizes, PROTOCOL['thresholds'], PROTOCOL['initial_null_trials'],
            PROTOCOL['seed'] + 11, directory, f'{dataset_id}-{train_hash}-{name}')
        scores, maximum = retrieval(fingerprints, train, references, selection_queries,
                                   PROTOCOL['thresholds'], null, actual_sizes)
        validation_rows[name] = scores
        report.setdefault('max_tanimoto_selection', {})[name] = summarize(maximum)
        for threshold, values in scores.items():
            report['candidates'].append({'fingerprint': name, 'threshold': float(threshold),
                                          'selection_metrics': summarize(values)})
        write_json(output, report)
    baseline = next(c for c in report['candidates'] if c['fingerprint'] == 'RDKit-path' and c['threshold'] == data['model']['threshold'])
    candidates = sorted(report['candidates'], key=lambda c: (-c['selection_metrics']['mrr'], -c['selection_metrics']['recall_at_10'],
                       c != baseline, c['fingerprint'], c['threshold']))
    winner = candidates[0]
    # Freeze the only contender before any test-query scoring.
    frozen = {'baseline': baseline, 'selection_winner': winner, 'frozen_before_test': True}
    write_json(directory / 'frozen-selection.json', frozen)
    report.update(frozen)
    report['status'] = 'confirming'
    write_json(output, report)
    print(f"Frozen selection: {winner['fingerprint']}, threshold {winner['threshold']}; MRR {winner['selection_metrics']['mrr']:.6f}", flush=True)
    paired = {'baseline': baseline, 'candidate': winner}
    confirmation = {}
    for representation in {c['fingerprint'] for c in paired.values()}:
        fingerprints = fingerprint_sets[representation]
        thresholds = sorted({c['threshold'] for c in paired.values() if c['fingerprint'] == representation})
        null, actual_sizes = null_prefixes([fingerprints[i] for i in train], sizes, thresholds,
            PROTOCOL['confirmatory_null_trials'], PROTOCOL['seed'] + 29,
            directory, f'{dataset_id}-{train_hash}-{representation}')
        for phase, queries in [('selection', selection_queries), ('test', test_queries)]:
            scores, maximum = retrieval(fingerprints, train, references, queries, thresholds, null, actual_sizes)
            for label, choice in paired.items():
                if choice['fingerprint'] == representation:
                    confirmation.setdefault(label, {})[phase] = scores[str(choice['threshold'])]
            report.setdefault('max_tanimoto_confirmatory', {}).setdefault(phase, {})[representation] = summarize(maximum)
    report['confirmation'] = {label: {phase: summarize(rows) for phase, rows in phases.items()}
                              for label, phases in confirmation.items()}
    comparison = paired_bootstrap(confirmation['candidate']['test'], confirmation['baseline']['test'],
                                  PROTOCOL['seed'], PROTOCOL['bootstrap_replicates'])
    report['independent_test_gain'] = dict(zip(['mrr', 'recall_at_1', 'recall_at_5', 'recall_at_10'],
        [{'mean': value, 'lower_95': low, 'upper_95': high} for value, low, high in
         zip(comparison['mean_gain'], comparison['lower_95'], comparison['upper_95'])]))
    report['retrieval_adoption_passed'] = adoption_passes(comparison) and winner != baseline
    report['status'] = 'completed'
    report['elapsed_seconds'] = time.perf_counter() - started
    report['limits'] = ['The candidate was frozen using selection data only; test data cannot select another candidate.',
        'Pointwise paired-bootstrap intervals resample unique query scaffolds; shared targets may still induce dependence.',
        'Study null streams have the same marginal definition but differ from production size-keyed streams.',
        'Study budgets are 8191/32767, below the production Balanced budget; rankings can change with Monte Carlo resolution.',
        PROTOCOL['interpretation'], 'A retrieval improvement does not establish calibrated Gumbel tails; validate those independently before adoption.']
    # Keep per-scaffold paired scores for independent interval recomputation.
    report['test_query_results'] = [
        {'scaffold_sha256': q['scaffold'], 'anchor_target_id': references[q['anchor']]['id'],
         'known_targets': [references[i]['id'] for i in q['truth']],
         'baseline': confirmation['baseline']['test'][i].tolist(),
         'candidate': confirmation['candidate']['test'][i].tolist()}
        for i, q in enumerate(test_queries)]
    write_json(output, report)
    print(json.dumps({k: report[k] for k in ['status', 'selection_winner', 'confirmation', 'independent_test_gain', 'retrieval_adoption_passed', 'elapsed_seconds']}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.dataset_id, args.cache, args.output)
