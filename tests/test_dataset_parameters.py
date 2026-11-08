"""Dataset-defined representation and empirical statistics reach every job path."""
import hashlib
import json
import numpy as np
import pytest
from fastapi.testclient import TestClient
from rdkit import DataStructs

from sea.api import app
from sea.chemistry import FP_SPEC, MORGAN_SPEC, molecule, normalize, raw_score
from sea.cli import demo
from sea.data import apply_parameter_recipe, load_dataset, save_parameter_variant
from sea.jobs import claim, create_key, process
from sea.precision import PROFILES
from sea.statistics import significance


def empirical_model(base):
    return {**base, 'threshold': .45, 'evd_usable': False, 'empirical_z_score': True,
            'calibration_method': 'independent-fixed-budget-empirical-null-v1'}


def test_fingerprint_choice_is_explicit_and_legacy_is_unchanged():
    query = 'COc1ccc2[nH]cc(c2c1)C(CC)CN'
    reference = 'COc1ccc2[nH]cc(CCN)c2c1'
    for spec, expected in [(FP_SPEC, .8214285714285714), (MORGAN_SPEC, .5454545454545454)]:
        a, b = molecule(query, spec), molecule(reference, spec)
        assert DataStructs.TanimotoSimilarity(a.fp, b.fp) == pytest.approx(expected)
    assert molecule(query).fp == molecule(query, FP_SPEC).fp
    assert molecule(query, MORGAN_SPEC).smiles == molecule(query).smiles
    assert len(normalize(['CCO', 'OCC'], MORGAN_SPEC)[0]) == 1


def test_empirical_z_uses_actual_null_sample_mean_and_variance():
    pool = [molecule(s, MORGAN_SPEC) for s in ['CCO', 'CCCO', 'CCCCO', 'CCN', 'c1ccccc1']]
    model = {'threshold': .45, 'std_b': 100, 'std_c': .5, 'mean_a': 100,
             'evd_usable': False, 'seed': 42, 'empirical_trials': 31, 'sizes': [1, 2], 'empirical_z_score': True}
    cache = {}
    result = significance(1., 1, 2, pool, model, 3, null_cache=cache)
    assert result['z_score'] == (1. - np.mean(cache[(1, 2)])) / np.std(cache[(1, 2)], ddof=1)
    assert result['method'] == 'monte_carlo'
    # A null with zero variance cannot justify an infinite standardized score.
    constant = {(1, 2): np.zeros(31)}
    assert significance(1., 1, 2, pool, model, null_cache=constant)['z_score'] is None


@pytest.fixture
def variant(tmp_path, monkeypatch):
    monkeypatch.setenv('SEA_HOME', str(tmp_path))
    source = demo()
    data, pool = load_dataset(source)
    fingerprints = {s: DataStructs.BitVectToBinaryText(molecule(s, MORGAN_SPEC).fp).hex() for s in pool}
    validation = {'approved': True, 'source_dataset_id': source, 'report_sha256': 'test-fixture'}
    model = empirical_model(data['model'])
    selected = save_parameter_variant(source, model, MORGAN_SPEC, fingerprints, validation)
    return source, selected, model, fingerprints, validation


def test_variant_is_content_addressed_and_preserves_source(variant):
    source, selected, model, fingerprints, validation = variant
    old, old_pool = load_dataset(source)
    new, new_pool = load_dataset(selected)
    assert source != selected and new['targets'] == old['targets']
    assert set(old_pool) == set(new_pool)
    assert old['manifest']['fingerprint'] == FP_SPEC
    assert new['manifest']['fingerprint'] == MORGAN_SPEC
    assert new['manifest']['parent_dataset_id'] == source
    assert save_parameter_variant(source, model, MORGAN_SPEC, fingerprints, validation) == selected
    with pytest.raises(ValueError):
        save_parameter_variant(source, model, MORGAN_SPEC, fingerprints, {**validation, 'approved': False})
    with pytest.raises(ValueError):
        save_parameter_variant(source, model, MORGAN_SPEC, {}, validation)


def test_published_recipe_rebuilds_exact_variant_and_rejects_mismatches(variant, tmp_path):
    source, _, model, fingerprints, _ = variant
    candidate = {'fingerprint': 'Morgan', 'threshold': .45}
    comparison = {'lower_95': [.1, .1, .1, .1]}
    report = {
        'source_dataset_id': source,
        'retrieval': {'dataset_id': source, 'status': 'completed', 'retrieval_adoption_passed': True,
                      'selection_winner': candidate},
        'production': {'dataset_id': source, 'production_adoption_passed': True, 'null_trials': 199999,
                       'paired_test_comparison': comparison, 'test_queries': 512, 'choices': {'candidate': candidate}},
        'empirical_null_audit': {'dataset_id': source, 'anti_conservative_detected': False,
                                 'fingerprint': MORGAN_SPEC, 'threshold': .45},
    }
    report_file = tmp_path / 'report.json'
    report_file.write_text(json.dumps(report))
    digest = hashlib.sha256(report_file.read_bytes()).hexdigest()
    validation = {'approved': True, 'source_dataset_id': source, 'report_sha256': digest}
    model = {**model, 'parameter_validation_sha256': digest}
    expected = save_parameter_variant(source, model, MORGAN_SPEC, fingerprints, validation)
    recipe = {'source_dataset_id': source, 'expected_dataset_id': expected, 'report_file': 'report.json',
              'report_sha256': digest, 'model': model, 'fingerprint': MORGAN_SPEC, 'validation': validation,
              'fingerprints_sha256': hashlib.sha256(b''.join(bytes.fromhex(v) for v in fingerprints.values())).hexdigest()}
    path = tmp_path / 'recipe.json'
    path.write_text(json.dumps(recipe))
    assert apply_parameter_recipe(source, path) == expected
    with pytest.raises(ValueError, match='different reference cohort'):
        apply_parameter_recipe(expected, path)
    original_files = set((tmp_path / 'datasets').iterdir())
    for changed, message in [({'fingerprints_sha256': 'bad'}, 'Fingerprint generation'),
                             ({'expected_dataset_id': '0' * 20}, 'content identity'),
                             ({'report_sha256': 'bad'}, 'checksum'),
                             ({'fingerprint': {**MORGAN_SPEC, 'radius': 3}}, 'Fingerprint generation'),
                             ({'model': {**model, 'evd_usable': True}}, 'empirical P/E')]:
        path.write_text(json.dumps({**recipe, **changed}))
        with pytest.raises(ValueError, match=message):
            apply_parameter_recipe(source, path)
        assert set((tmp_path / 'datasets').iterdir()) == original_files


def test_empirical_model_disables_extrapolated_gumbel_for_ensembles():
    pool = [molecule('C' * n + 'O', MORGAN_SPEC) for n in range(1, 21)]
    model = {'threshold': .45, 'std_b': 0., 'std_c': .5, 'mean_a': 0.,
             'evd_usable': False, 'seed': 42, 'empirical_trials': 31, 'sizes': [1, 10, 100],
             'empirical_z_score': True, 'gumbel_loc': 0., 'gumbel_scale': 1.}
    result = significance(10., 10, 10, pool, model)
    assert result['method'] == 'monte_carlo' and result['p_resolution'] == 1/32


def test_preload_prefers_approved_reference_at_the_same_release(monkeypatch):
    import sea.data as data
    monkeypatch.setenv('SEA_PRELOAD_DATASETS', 'auto')
    loaded = []
    monkeypatch.setattr(data, 'load_dataset', loaded.append)
    monkeypatch.setattr(data, 'list_datasets', lambda: [
        {'id': 'unchecked', 'source': 'ChEMBL', 'release': 'ChEMBL_37', 'created_at': '2099-01-01'},
        {'id': 'approved', 'source': 'ChEMBL', 'release': 'ChEMBL_37', 'created_at': '2026-11-08', 'parameter_validation': {'approved': True}},
    ])
    assert data.preload_datasets() == ['approved'] and loaded == ['approved']


@pytest.mark.parametrize('profile', ['legacy', 'balanced'])
@pytest.mark.parametrize('kind', ['prediction', 'comparison'])
def test_workers_use_dataset_fingerprints_for_queries_and_ensembles(variant, monkeypatch, profile, kind):
    source, selected, model, _, _ = variant
    monkeypatch.setitem(PROFILES, 'balanced', 1999)
    client = TestClient(app)
    auth = {'X-API-Key': create_key('fingerprint-test')}
    payload = {'dataset_id': selected, 'statistics_profile': profile}
    if kind == 'prediction':
        payload['queries'] = [['CCCO']]
        endpoint = '/api/v1/predictions'
    else:
        payload.update(left=['CCCO'], right=['CCO', 'CCCCO'])
        endpoint = '/api/v1/comparisons'
    response = client.post(endpoint, headers=auth, json=payload)
    assert response.status_code == 202
    jid = response.json()['job_id']
    process(claim('fixture'))
    result = client.get(f'/api/v1/jobs/{jid}/result', headers=auth).json()
    data, pool = load_dataset(selected)
    for row in result['results']:
        left, _ = normalize(['CCCO'], MORGAN_SPEC)
        if kind == 'prediction':
            target = next(t for t in data['targets'] if t['id'] == row['target']['id'])
            right = [pool[s] for s in target['smiles']]
        else:
            right, _ = normalize(payload['right'], MORGAN_SPEC)
        expected = raw_score(left, right, .45)
        assert row['raw_score'] == expected[0] and row['similar_ligands'] == expected[2]
        assert row['method'] in ('exact_zero', 'monte_carlo')
    if profile == 'balanced':
        process(claim('refinement'))
        child = result['analysis']['refinement_job_id']
        assert client.get(f'/api/v1/jobs/{child}/result', headers=auth).json()['analysis']['complete']
