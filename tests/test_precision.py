"""Fixed budgets, exact checkpoint continuation and immutable linked jobs."""
import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from fastapi.testclient import TestClient
from scipy.stats import binom

from sea.api import app
from sea.chemistry import molecule
from sea.cli import demo
from sea.data import home, load_dataset
from sea.jobs import claim, connection, create_key, process
from sea.precision import PROFILES, details, profile
from sea.precision_cache import PrecisionCache
from sea.statistics import significance


@pytest.fixture
def sampling():
    pool = [molecule(s) for s in ['CCO', 'CCCO', 'CCCCO', 'CC(=O)O', 'c1ccccc1', 'CCN', 'CCCN', 'CCCCN']]
    model = {'threshold': .57, 'std_b': 0., 'std_c': .5, 'mean_a': 0.,
             'evd_usable': False, 'seed': 42, 'empirical_trials': 137, 'sizes': [1, 10]}
    return pool, model


@pytest.mark.parametrize('sizes', [(1, 1), (1, 8), (3, 7)])
def test_prefix_extension_and_smaller_budget_are_bit_exact(sampling, tmp_path, sizes):
    pool, model = sampling
    for budget in (67, 137, 31):
        cache = PrecisionCache(tmp_path, 'reference')
        effective = {**model, 'empirical_trials': budget}
        expected = {}
        assert significance(1., *sizes, pool, effective, 1517, null_cache=cache) == significance(
            1., *sizes, pool, effective, 1517, null_cache=expected)
        assert np.array_equal(cache[(*sizes, budget)], expected[sizes])
        if budget == 31:
            assert cache.index is None
    assert len(list(tmp_path.rglob('*.npz'))) == 1


def test_interrupted_block_preserves_rng_and_resumes_exactly(sampling, tmp_path, monkeypatch):
    from sea.background import FingerprintIndex
    pool, model = sampling
    original = FingerprintIndex.scores
    calls = []
    def interrupted(self, *args):
        calls.append(1)
        if len(calls) == 2:
            raise InterruptedError('Cancelled inside a block after RNG draws')
        return original(self, *args)
    monkeypatch.setattr(FingerprintIndex, 'scores', interrupted)
    with pytest.raises(InterruptedError):
        PrecisionCache(tmp_path, 'reference').sample(pool, model, (1, 8), lambda: None)
    with np.load(next(tmp_path.rglob('*.npz'))) as saved:
        assert len(saved['values']) == 64
    monkeypatch.setattr(FingerprintIndex, 'scores', original)
    resumed = PrecisionCache(tmp_path, 'reference').sample(pool, model, (1, 8), lambda: None)
    expected = {}
    significance(1., 1, 8, pool, model, 1517, null_cache=expected)
    assert np.array_equal(resumed, expected[(1, 8)])


def test_corruption_and_concurrent_sampling(sampling, tmp_path):
    pool, model = sampling
    def sample(_):
        cache = PrecisionCache(tmp_path, 'reference')
        return cache.sample(pool, model, (1, 8), lambda: None), cache.index
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(sample, range(2)))
    assert np.array_equal(results[0][0], results[1][0])
    assert sum(index is not None for _, index in results) == 1
    path = next(tmp_path.rglob('*.npz'))
    for contents in (b'truncated', path.read_bytes()[:100]):
        path.write_bytes(contents)
        regenerated, _ = sample(0)
        assert np.array_equal(regenerated, results[0][0])


def test_budget_and_seed_isolation_in_memory(sampling, tmp_path):
    pool, model = sampling
    cache = PrecisionCache(tmp_path, 'reference')
    a = cache.distribution(pool, {**model, 'empirical_trials': 31}, (1, 8), lambda: None)
    b = cache.distribution(pool, model, (1, 8), lambda: None)
    assert len(a) == 31 and len(b) == 137
    changed = PrecisionCache(tmp_path, 'reference').sample(pool, {**model, 'seed': 43}, (1, 8), lambda: None)
    assert len(list(tmp_path.rglob('*.npz'))) == 2
    assert not np.array_equal(b, changed)


def test_screening_reads_prefix_without_waiting_for_refinement_lock(sampling, tmp_path):
    import fcntl
    pool, model = sampling
    PrecisionCache(tmp_path, 'reference').sample(pool, model, (1, 8), lambda: None)
    path = next(tmp_path.rglob('*.lock'))
    with path.open('a+b') as locked:
        fcntl.flock(locked, fcntl.LOCK_EX)
        calls = []
        def check():
            calls.append(1)
            assert len(calls) <= 1, 'Screening waited for the refinement writer lock'
        reader = PrecisionCache(tmp_path, 'reference')
        values = reader.sample(pool, {**model, 'empirical_trials': 31}, (1, 8), check)
        assert len(values) == 31 and reader.index is None


def test_checksum_tampering_regenerates_samples(sampling, tmp_path):
    pool, model = sampling
    original = PrecisionCache(tmp_path, 'reference').sample(pool, model, (1, 8), lambda: None)
    path = next(tmp_path.rglob('*.npz'))
    with np.load(path) as archive:
        saved = dict(archive)
    saved['values'] = saved['values'] + 123
    np.savez(path, **saved)
    regenerated = PrecisionCache(tmp_path, 'reference').sample(pool, model, (1, 8), lambda: None)
    assert np.array_equal(original, regenerated)


@pytest.mark.parametrize('k,budget,n', [(0, 999, 1517), (5, 999, 3), (999, 999, 1), (1, 199999, 1517)])
def test_intervals_invert_binomial_tails(k, budget, n):
    stat = details({'method': 'monte_carlo', 'p_value': (k + 1) / (budget + 1), 'comparisons': n}, budget, 'final', {'id': 'p', 'planned_trials': budget})
    lo, hi = stat['tail_probability_interval']
    assert stat['exceedances'] == k and stat['trials'] == budget
    assert stat['e_tail_interval'] == [lo * n, hi * n]
    if k > 0:
        assert binom.sf(k - 1, budget, lo) == pytest.approx(.025)
    if k < budget:
        assert binom.cdf(k, budget, hi) == pytest.approx(.025)
    if k == 0:
        assert lo == 0 and stat['status'] == 'resolution_limited'


def test_production_profiles_have_reproducible_distinct_ids():
    assert PROFILES == {'balanced': 199999, 'research': 999999}
    data = {'manifest': {'id': 'dataset', 'model_id': 'model'}, 'model': {'seed': 42}}
    assert profile(data, 'balanced') == profile(data, 'balanced')
    assert profile(data, 'balanced')['id'] != profile(data, 'research')['id']
    modified = {**data, 'model': {'seed': 43}}
    assert profile(data, 'balanced')['id'] != profile(modified, 'balanced')['id']


@pytest.fixture
def phased(tmp_path, monkeypatch):
    monkeypatch.setenv('SEA_HOME', str(tmp_path))
    # Smaller budgets only for lifecycle tests. Production budgets are asserted
    # above and exercised by the browser smoke test without patching.
    monkeypatch.setitem(PROFILES, 'balanced', 1999)
    monkeypatch.setitem(PROFILES, 'research', 3999)
    dataset = demo()
    return TestClient(app), dataset, {'X-API-Key': create_key('alice')}, {'X-API-Key': create_key('bob')}


def submit(client, dataset, auth, **extra):
    return client.post('/api/v1/predictions', headers=auth, json={
        'dataset_id': dataset, 'queries': [['CCO'], ['CCN']], 'top_k': 1,
        'statistics_profile': 'balanced', **extra}).json()['job_id']


def test_evidence_screening_final_and_original_testing_scope(phased, monkeypatch):
    import sea.phased_jobs as pipeline
    client, dataset, a, b = phased
    jid = submit(client, dataset, a)
    original = pipeline.significance
    previews = []
    def inspect(*args, **kwargs):
        if not previews:
            response = client.get(f'/api/v1/jobs/{jid}/preview', headers=a)
            assert response.status_code == 200
            snapshot = response.json()
            assert snapshot['analysis']['phase'] == 'evidence'
            assert all(r['p_value'] is None and r['method'] == 'pending' for r in snapshot['results'])
            previews.append(snapshot)
        return original(*args, **kwargs)
    monkeypatch.setattr(pipeline, 'significance', inspect)
    process(claim('foreground', 'foreground'))
    source = client.get(f'/api/v1/jobs/{jid}/result', headers=a).json()
    assert source['analysis']['phase'] == 'screening' and not source['analysis']['complete']
    child = source['analysis']['refinement_job_id']
    assert claim('reserved', 'foreground') is None
    assert client.get(f'/api/v1/jobs/{child}', headers=b).status_code == 404
    assert client.get(f'/api/v1/jobs/{jid}/preview', headers=b).status_code == 404
    process(claim('background'))
    final = client.get(f'/api/v1/jobs/{child}/result', headers=a).json()
    assert final['analysis']['phase'] == 'final' and final['analysis']['complete']
    assert final['analysis']['total_comparisons'] == final['analysis']['evaluated_comparisons'] == 6
    assert len(final['results']) == 2  # top_k is presentation only.
    assert all(r['comparisons'] == 3 for r in final['results'])
    assert client.get(f'/api/v1/jobs/{jid}/result', headers=a).json() == source
    with connection() as con:
        all_rows = json.loads(con.execute('SELECT checkpoint FROM jobs WHERE id=?', (child,)).fetchone()[0])['results']
    data, structures = load_dataset(dataset)
    model = {**data['model'], 'empirical_trials': 1999}
    for row in all_rows:
        expected = significance(row['raw_score'], row['left_size'], row['right_size'], list(structures.values()), model, 3, null_cache={})
        assert all(row[key] == value for key, value in expected.items())
        assert row['statistics']['phase'] == 'final'
    exported = client.get(f'/api/v1/jobs/{child}/result?format=csv', headers=a)
    record = next(csv.DictReader(io.StringIO(exported.text)))
    assert record['analysis_complete'] == 'True' and record['planned_trials'] == '1999'
    assert record['source_job_id'] == jid and record['testing_family'] == 'per_query'


def test_idempotent_refinement_pause_resume_and_foreground_priority(phased):
    client, dataset, a, b = phased
    jid = submit(client, dataset, a)
    process(claim('screen'))
    child = client.get(f'/api/v1/jobs/{jid}', headers=a).json()['refinement_job_id']
    def refine(_):
        return client.post(f'/api/v1/jobs/{jid}/refine', headers=a, json={}).json()['job_id']
    with ThreadPoolExecutor(max_workers=2) as executor:
        assert list(executor.map(refine, range(2))) == [child, child]
    foreground = submit(client, dataset, a)
    job = claim('priority')
    assert job['id'] == foreground
    process(job)
    assert client.post(f'/api/v1/jobs/{child}/cancel', headers=b).status_code == 404
    assert client.post(f'/api/v1/jobs/{child}/cancel', headers=a).json()['status'] == 'cancelled'
    assert client.post(f'/api/v1/jobs/{child}/resume', headers=b).status_code == 404
    assert client.post(f'/api/v1/jobs/{child}/resume', headers=a).json()['status'] == 'queued'
    assert client.post(f'/api/v1/jobs/{child}/resume', headers=a).status_code == 409
    job = claim('refiner')
    assert job['id'] == child
    process(job)
    assert client.get(f'/api/v1/jobs/{child}/result', headers=a).json()['analysis']['complete']
    research = client.post(f'/api/v1/jobs/{jid}/refine', headers=a, json={'statistics_profile': 'research'}).json()['job_id']
    assert research != child
    assert client.get(f'/api/v1/jobs/{jid}', headers=a).json()['refinement_job_id'] == child
    assert client.get(f'/api/v1/jobs/{child}/result', headers=a).json()['analysis']['profile']['name'] == 'balanced'


def test_running_pause_retains_checkpoint_and_resumes(phased, monkeypatch):
    import sea.phased_jobs as pipeline
    client, dataset, a, b = phased
    jid = submit(client, dataset, a)
    process(claim('screen'))
    child = client.get(f'/api/v1/jobs/{jid}', headers=a).json()['refinement_job_id']
    original = pipeline.significance
    def pause(*args, **kwargs):
        # Cancellation during finalization must never publish complete=True.
        client.post(f'/api/v1/jobs/{child}/cancel', headers=a)
        return original(*args, **kwargs)
    monkeypatch.setattr(pipeline, 'significance', pause)
    process(claim('precision'))
    assert client.get(f'/api/v1/jobs/{child}', headers=a).json()['status'] == 'cancelled'
    preview = client.get(f'/api/v1/jobs/{child}/preview', headers=a).json()
    assert not preview['analysis']['complete']
    assert client.get(f'/api/v1/jobs/{child}/result', headers=a).status_code == 409
    monkeypatch.setattr(pipeline, 'significance', original)
    client.post(f'/api/v1/jobs/{child}/resume', headers=a)
    process(claim('resumed'))
    assert client.get(f'/api/v1/jobs/{child}/result', headers=a).json()['analysis']['complete']


def test_legacy_cannot_be_refined_and_invalid_profile_rejected(phased):
    client, dataset, a, b = phased
    jid = submit(client, dataset, a, statistics_profile='legacy')
    process(claim('legacy'))
    assert client.post(f'/api/v1/jobs/{jid}/refine', headers=a, json={}).status_code == 409
    assert client.post('/api/v1/predictions', headers=a, json={'dataset_id': dataset, 'queries': [['CCO']], 'statistics_profile': 'adaptive'}).status_code == 422
    assert client.post(f'/api/v1/jobs/{jid}/resume', headers=a).status_code == 409


def test_legacy_database_migration_preserves_completed_artifacts(tmp_path, monkeypatch):
    import sqlite3
    monkeypatch.setenv('SEA_HOME', str(tmp_path))
    with sqlite3.connect(tmp_path / 'service.sqlite') as con:
        con.execute('CREATE TABLE jobs (id TEXT PRIMARY KEY, owner TEXT, kind TEXT, payload TEXT, status TEXT, progress REAL DEFAULT 0, result TEXT, error TEXT, cancel INTEGER DEFAULT 0, created REAL, updated REAL, worker TEXT)')
        con.execute("INSERT INTO jobs(id,status,result) VALUES('old','completed','original bytes')")
    with connection() as con:
        row = con.execute("SELECT * FROM jobs WHERE id='old'").fetchone()
        assert row['profile'] == 'legacy' and row['result'] == 'original bytes'
        assert row['status'] == 'completed' and row['checkpoint'] is None


def test_preload_selection_is_bounded_and_can_be_disabled(monkeypatch):
    import sea.data as data
    manifests = [
        {'id': 'demo', 'source': 'synthetic-demo', 'created_at': '2026-12'},
        {'id': 'chembl', 'source': 'ChEMBL', 'release': 'ChEMBL_37', 'created_at': '2026-01'},
        {'id': 'older-release', 'source': 'ChEMBL', 'release': 'ChEMBL_36', 'created_at': '2026-02'},
    ]
    loaded = []
    monkeypatch.setattr(data, 'list_datasets', lambda: manifests)
    monkeypatch.setattr(data, 'load_dataset', loaded.append)
    monkeypatch.setenv('SEA_PRELOAD_DATASETS', 'auto')
    assert data.preload_datasets() == ['chembl'] and loaded == ['chembl']
    monkeypatch.setenv('SEA_PRELOAD_DATASETS', 'none')
    assert data.preload_datasets() == []
    monkeypatch.setenv('SEA_PRELOAD_DATASETS', 'one,two,three')
    with pytest.raises(ValueError):
        data.preload_datasets()
