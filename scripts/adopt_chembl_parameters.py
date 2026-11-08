"""Publish an independently confirmed parameter variant into the chosen SEA_HOME.

Use an isolated staging SEA_HOME first. Existing source references and results
are never overwritten. Reports, runtime fingerprint agreement and cohort
identity are checked before saving a new content-addressed reference.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import DataStructs

from sea.annotations import annotation_path
from sea.chemistry import MORGAN_SPEC, FP_SPEC, molecule
from sea.data import load_dataset, save_parameter_variant
from sea.parameter_validation import approved_decision


def adopt(dataset_id, directory, report_path):
    directory = Path(directory)
    paths = {'retrieval': directory / 'retrieval-results.json', 'production': directory / 'production-budget-results.json',
             'legacy_gumbel_audit': directory / 'null-audit-results.json', 'empirical_null_audit': directory / 'empirical-null-results.json'}
    reports = {name: json.loads(path.read_text()) for name, path in paths.items()}
    candidate = approved_decision(reports['retrieval'], reports['production'], reports['empirical_null_audit'])
    if reports['retrieval']['dataset_id'] != dataset_id:
        raise ValueError('The requested dataset differs from the reports.')
    data, structures = load_dataset(dataset_id)
    identity = json.loads((directory / 'complete.json').read_text())
    assert identity['dataset_id'] == dataset_id
    assert identity['smiles_order_sha256'] == hashlib.sha256('\n'.join(structures).encode()).hexdigest()
    report = {'source_dataset_id': dataset_id, 'approved': True,
              'decision': 'Adopt the frozen representation/threshold candidate after both holdout budgets; use fixed-budget empirical null rather than an uncertified Gumbel approximation.',
              'source_report_sha256': {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()},
              **reports}
    encoded = json.dumps(report, indent=2, allow_nan=False) + '\n'
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(encoded)
    report_hash = hashlib.sha256(encoded.encode()).hexdigest()
    spec = MORGAN_SPEC if candidate['fingerprint'] == 'Morgan' else FP_SPEC
    if spec['kind'] == 'Morgan':
        packed = np.load(directory / 'morgan.npy', mmap_mode='r')
        for i in np.random.default_rng(20261108).choice(len(structures), min(64, len(structures)), replace=False):
            s = list(structures)[i]
            assert DataStructs.BitVectToBinaryText(molecule(s, spec).fp) == packed[i].tobytes()
        fingerprints = {s: row.tobytes().hex() for s, row in zip(structures, packed)}
    else:
        fingerprints = dict(data['fingerprints'])
    # Legacy power-law slots remain zero-valued compatibility fields. Z is
    # estimated from the actual size-specific null, never those placeholders.
    model = {'threshold': candidate['threshold'], 'mean_a': 0., 'std_b': 0., 'std_c': .5,
             'gumbel_loc': 0., 'gumbel_scale': 1., 'validation_ks': None, 'evd_usable': False,
             'sizes': [1, 10, 30, 100], 'seed': 42, 'repeats': 0, 'empirical_trials': 999,
             'empirical_z_score': True, 'fingerprint': spec,
             'calibration_method': 'independent-fixed-budget-empirical-null-v1',
             'analytic_tail_policy': 'disabled; no extrapolated Gumbel tail is certified',
             'parameter_validation_sha256': report_hash}
    validation = {'source_dataset_id': dataset_id, 'approved': True, 'protocol': 'chembl-parameters-v1',
                  'report_sha256': report_hash, 'report': 'docs/chembl37-parameter-validation-results.json',
                  'selected_fingerprint': spec, 'selected_threshold': model['threshold'],
                  'calibration_method': model['calibration_method'], 'gumbel_enabled': False,
                  'z_score_method': 'sample mean and standard deviation of the actual-size null',
                  'test_queries': reports['production']['test_queries'], 'null_trials': 199999,
                  'baseline_metrics': reports['production']['metrics']['baseline'],
                  'selected_metrics': reports['production']['metrics']['candidate'],
                  'scope': 'Known-active target retrieval on a global scaffold holdout; not inactive discrimination or binding-probability calibration.'}
    new_id = save_parameter_variant(dataset_id, model, spec, fingerprints, validation)
    old_annotations = annotation_path(dataset_id)
    if old_annotations.exists():
        annotations = json.loads(old_annotations.read_text())
        annotations['dataset_id'] = new_id
        dest = annotation_path(new_id)
        temporary = dest.with_suffix('.tmp')
        temporary.write_text(json.dumps(annotations, ensure_ascii=False, allow_nan=False))
        temporary.replace(dest)
    print(json.dumps({'dataset_id': new_id, 'parent_dataset_id': dataset_id, 'fingerprint': spec,
                      'threshold': model['threshold'], 'gumbel_enabled': False, 'report_sha256': report_hash}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    adopt(args.dataset_id, args.cache, args.report)
