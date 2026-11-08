"""Exact agreement with RDKit/scalar sampling, persistence and cancellation."""
import numpy as np
import pytest
from rdkit import DataStructs

from sea.background import FingerprintIndex, NullCache
from sea.chemistry import Molecule, raw_score
from sea.statistics import significance


@pytest.fixture
def pool():
    rng = np.random.default_rng(17)
    result = []
    for i, count in enumerate([0, 0, 1, 2, 3, 4, 20, 40, 50, 100, 150, 300, 600, 900, 1200, 1800, 2048] * 2):
        fp = DataStructs.ExplicitBitVect(2048)
        for position in rng.choice(2048, count, replace=False):
            fp.SetBit(int(position))
        result.append(Molecule(str(i), fp))
    return result


@pytest.mark.parametrize('threshold', [0., 1 / 3, .57, 1.])
def test_batched_scores_are_bit_exact_including_empty_fingerprints(pool, threshold):
    index = FingerprintIndex(pool)
    # Multiple trials, multiple left molecules, zero-bit and full fingerprints.
    left = np.array([[0, 1, 16], [2, 14, 4], [16, 13, 32]])
    right = np.array([[1, 0, 16, 20], [5, 14, 3, 2], [17, 32, 31, 33]])
    expected = [raw_score([pool[i] for i in a], [pool[i] for i in b], threshold, 0)[0]
                for a, b in zip(left, right)]
    assert np.array_equal(index.scores(left, right, threshold, lambda: None), expected)


def test_threshold_bound_retains_equality_and_excludes_next_float():
    pool = []
    for positions in ([0, 1], [0, 1, 2, 3]):
        fp = DataStructs.ExplicitBitVect(2048)
        fp.SetBitsFromList(positions)
        pool.append(Molecule(str(positions), fp))
    index = FingerprintIndex(pool)
    left, right = np.array([[0]]), np.array([[1]])
    assert index.scores(left, right, .5, lambda: None)[0] == .5
    assert index.scores(left, right, np.nextafter(.5, 1.), lambda: None)[0] == 0.


def test_cancellation_is_checked_between_intersection_chunks(pool):
    class Cancelled(Exception):
        pass
    checks = []
    def check():
        checks.append(1)
        if len(checks) == 2:
            raise Cancelled()
    index = FingerprintIndex(pool)
    full = np.full((1, 100), 16, dtype=np.int64)
    with pytest.raises(Cancelled):
        index.scores(full, full, .57, check)
    assert len(checks) == 2


@pytest.mark.parametrize('sizes', [(1, 1), (1, 8), (3, 7), (10, 10), (1, 34)])
def test_samples_and_statistics_match_scalar_reference_and_reload(pool, tmp_path, sizes):
    model = {'threshold': .57, 'std_b': 0., 'std_c': .5, 'mean_a': 0.,
             'evd_usable': False, 'seed': 42, 'empirical_trials': 131, 'sizes': [1, 10]}
    scalar = {}
    expected = significance(1., *sizes, pool, model, 1517, null_cache=scalar)
    cache = NullCache(tmp_path, 'dataset-one')
    actual = significance(1., *sizes, pool, model, 1517, null_cache=cache)
    assert actual == expected
    assert np.array_equal(cache[sizes], scalar[sizes])
    reloaded = NullCache(tmp_path, 'dataset-one')
    assert significance(1., *reversed(sizes), pool, model, 1517, null_cache=reloaded) == expected
    assert np.array_equal(reloaded[sizes], scalar[sizes])
    assert reloaded.index is None  # Disk hits need no packed fingerprint index.
    assert NullCache(tmp_path, 'dataset-two').directory != cache.directory
    assert NullCache(tmp_path / 'another-service', 'dataset-one').directory != cache.directory


@pytest.mark.parametrize('bad', [b'truncated', np.array([float('nan')] * 5), np.zeros(4), np.zeros(5, dtype=np.int32)])
def test_invalid_cache_is_regenerated(pool, tmp_path, bad):
    cache = NullCache(tmp_path, 'dataset')
    cache.directory.mkdir(parents=True)
    path = cache.directory / '1-2.npy'
    if isinstance(bad, bytes):
        path.write_bytes(bad)
    else:
        np.save(path, bad)
    model = {'seed': 42, 'threshold': .57, 'empirical_trials': 5}
    values = cache.sample(pool, model, (1, 2), lambda: None)
    assert values.shape == (5,) and np.isfinite(values).all()
    assert cache.index is not None
    assert np.array_equal(np.load(path, allow_pickle=False), values)


def test_cancelled_sampling_leaves_no_distribution(pool, tmp_path):
    class Cancelled(Exception):
        pass
    calls = []
    def check():
        calls.append(1)
        if len(calls) == 3:
            raise Cancelled()
    cache = NullCache(tmp_path, 'dataset')
    with pytest.raises(Cancelled):
        cache.sample(pool, {'seed': 42, 'threshold': .57, 'empirical_trials': 131}, (1, 8), check)
    assert not list(tmp_path.rglob('*.npy'))


def test_unwritable_disk_cache_still_returns_valid_scores(pool, tmp_path, monkeypatch):
    cache = NullCache(tmp_path, 'dataset')
    def blocked(*args, **kwargs):
        raise PermissionError('Read-only cache')
    monkeypatch.setattr(type(tmp_path), 'mkdir', blocked)
    values = cache.sample(pool, {'seed': 42, 'threshold': .57, 'empirical_trials': 5}, (1, 2), lambda: None)
    assert len(values) == 5 and np.isfinite(values).all()
