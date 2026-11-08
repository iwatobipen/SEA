"""Confirm the frozen candidate at the actual Balanced Monte Carlo budget."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from rdkit import DataStructs

from sea.data import load_dataset
from sea.parameter_validation import (PROTOCOL, adoption_passes, choose_queries,
    null_prefixes, paired_bootstrap, reference_split, retrieval, summarize)


def run(dataset_id, directory, output):
    started = time.perf_counter()
    directory = Path(directory)
    protocol = json.loads((directory / 'production-confirmation-protocol.json').read_text())
    assert protocol['null_trials'] == 199999 and protocol['declared_before_holdout_results_review']
    frozen = json.loads((directory / 'frozen-selection.json').read_text())
    data, structures = load_dataset(dataset_id)
    smiles = list(structures)
    groups = np.load(directory / 'groups.npy', mmap_mode='r')
    train, references, membership = reference_split(data, smiles, groups)
    test = choose_queries(smiles, groups, membership, references, 0, PROTOCOL['queries_per_partition'], PROTOCOL['seed'])
    choices = {'baseline': frozen['baseline'], 'candidate': frozen['selection_winner']}
    rows, maximum = {}, {}
    train_hash = hashlib.sha256(train.tobytes()).hexdigest()[:12]
    for representation in sorted({c['fingerprint'] for c in choices.values()}):
        if representation == 'RDKit-path':
            fps = [m.fp for m in structures.values()]
        else:
            packed = np.load(directory / 'morgan.npy', mmap_mode='r')
            fps = [DataStructs.CreateFromBinaryText(row.tobytes()) for row in packed]
        thresholds = sorted({c['threshold'] for c in choices.values() if c['fingerprint'] == representation})
        def progress(message, **kwargs):
            if 'null trials' not in message or int(message.split(': ')[-1].split('/')[0]) % 16384 == 0:
                print(message, flush=True)
        null, sizes = null_prefixes([fps[i] for i in train], [len(t['indices']) for t in references],
            thresholds, protocol['null_trials'], PROTOCOL['seed'] + 29, directory,
            f'{dataset_id}-{train_hash}-{representation}', progress)
        scores, max_rows = retrieval(fps, train, references, test, thresholds, null, sizes)
        maximum[representation] = summarize(max_rows)
        for label, choice in choices.items():
            if choice['fingerprint'] == representation:
                rows[label] = scores[str(choice['threshold'])]
        del fps, null
    comparison = paired_bootstrap(rows['candidate'], rows['baseline'], PROTOCOL['seed'], PROTOCOL['bootstrap_replicates'])
    result = {'dataset_id': dataset_id, 'protocol': protocol, 'choices': choices,
              'test_queries': len(test), 'training_targets': len(references),
              'null_trials': protocol['null_trials'], 'metrics': {k: summarize(v) for k, v in rows.items()},
              'max_tanimoto': maximum, 'paired_test_comparison': comparison,
              'production_adoption_passed': adoption_passes(comparison),
              'elapsed_seconds': time.perf_counter()-started,
              'timing_scope': 'Elapsed time covers this invocation only; reused null distributions may have been generated earlier. This is not an end-to-end sampling performance benchmark.',
              'limits': 'Same fixed scaffold holdout as the earlier comparison, re-evaluated at the production budget; not another independent biological replicate. Candidates remain frozen.'}
    result['test_query_results'] = [{'scaffold_sha256': q['scaffold'], 'anchor_target_id': references[q['anchor']]['id'],
        'baseline': rows['baseline'][i].tolist(), 'candidate': rows['candidate'][i].tolist()} for i, q in enumerate(test)]
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: result[k] for k in ['metrics', 'paired_test_comparison', 'production_adoption_passed', 'elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.dataset_id, args.cache, args.output)
