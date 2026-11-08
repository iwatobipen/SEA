"""Scientific validation safeguards, independent of reported ChEMBL outcomes."""
import hashlib

import numpy as np
import pytest
from rdkit import DataStructs

from sea.chemistry import molecule, raw_score
from sea.parameter_validation import (adoption_passes, buckets, choose_queries,
    null_prefixes, paired_bootstrap, raw_scores, reference_split, retrieval)


def test_global_split_removes_shared_scaffolds_from_every_target():
    groups = np.array([hashlib.sha256(str(i).encode()).hexdigest().encode() for i in range(100)])
    smiles = [f'molecule-{i}' for i in range(100)]
    data = {'targets': [{'id': 'A', 'smiles': smiles[:80]}, {'id': 'B', 'smiles': smiles[20:]}]}
    train, references, membership = reference_split(data, smiles, groups)
    training_groups = set(groups[train])
    for partition in (0, 1):
        queries = choose_queries(smiles, groups, membership, references, partition, 20, 42)
        assert len({q['scaffold'] for q in queries}) == len(queries)
        assert not training_groups & {q['scaffold'].encode() for q in queries}
        assert all(membership[q['index']] == partition and q['truth'] for q in queries)
    # Identical scaffold groups, including the empty-scaffold group, cannot
    # fall into different partitions even if multiple targets share ligands.
    assert buckets(np.array([groups[0], groups[0]])).tolist() == [membership[0]] * 2


@pytest.mark.parametrize('threshold', [.0, .3, .56, 1.])
def test_study_singleton_aggregation_matches_production_exactly(threshold):
    pool = [molecule(s) for s in ['CCO', 'CCCO', 'CCCCO', 'CCN', 'c1ccccc1']]
    refs = [{'indices': np.array([1, 2, 3])}, {'indices': np.array([2, 4, 0, 1])}]
    similarities = np.array(DataStructs.BulkTanimotoSimilarity(pool[0].fp, [m.fp for m in pool]))
    expected = [raw_score([pool[0]], [pool[i] for i in ref['indices']], threshold, 0)[0] for ref in refs]
    assert np.array_equal(raw_scores(similarities, refs, threshold), expected)


def test_prefix_null_matches_independent_uniform_draws_and_exact_scores(tmp_path):
    pool = [molecule(s) for s in ['CCO', 'CCCO', 'CCCCO', 'CCN', 'c1ccccc1']]
    thresholds, sizes, trials, seed = [.3, .56, 1.], [1, 3], 31, 123
    values, actual = null_prefixes([m.fp for m in pool], sizes, thresholds, trials, seed, tmp_path, 'fixture')
    rng = np.random.default_rng(seed)
    expected = np.zeros((len(thresholds), len(sizes), trials))
    for trial in range(trials):
        query = int(rng.integers(len(pool)))
        right = rng.choice(len(pool), max(sizes), replace=False)
        for i, threshold in enumerate(thresholds):
            for j, size in enumerate(sizes):
                expected[i, j, trial] = raw_score([pool[query]], [pool[k] for k in right[:size]], threshold, 0)[0]
    assert actual.tolist() == sizes
    assert np.array_equal(values, np.sort(expected, axis=2))
    reloaded, _ = null_prefixes([m.fp for m in pool], sizes, thresholds, trials, seed, tmp_path, 'fixture')
    assert np.array_equal(values, reloaded)


def test_retrieval_keeps_zero_scores_last_and_known_active_recall_only(tmp_path):
    pool = [molecule(s) for s in ['CCO', 'CCCO', 'c1ccccc1']]
    references = [{'id': 'alcohol', 'indices': np.array([0, 1])},
                  {'id': 'aromatic', 'indices': np.array([2])}]
    null, sizes = null_prefixes([m.fp for m in pool], [1, 2], [.56], 31, 12, tmp_path, 'fixture')
    scores, _ = retrieval([m.fp for m in pool], np.arange(3), references,
        [{'index': 0, 'truth': [0]}], [.56], null, sizes)
    assert scores['0.56'].tolist() == [[1., 1., 1., 1.]]


def test_bootstrap_is_paired_and_adoption_requires_gain_and_recall_guard():
    baseline = np.tile([.5, .2, .7, .8], (64, 1))
    equal = paired_bootstrap(baseline, baseline, repetitions=257)
    assert equal['lower_95'] == equal['upper_95'] == [0.] * 4
    assert not adoption_passes(equal)
    improved = paired_bootstrap(baseline + [.1, 0, 0, 0], baseline, repetitions=257)
    assert adoption_passes(improved)
    harmed = paired_bootstrap(baseline + [.1, 0, 0, -.1], baseline, repetitions=257)
    assert not adoption_passes(harmed)


def test_adoption_gate_rejects_mismatched_or_unconfirmed_parameters():
    from sea.parameter_validation import approved_decision
    candidate = {'fingerprint': 'Morgan', 'threshold': .45}
    retrieval = {'dataset_id': 'source', 'status': 'completed', 'retrieval_adoption_passed': True,
                 'selection_winner': candidate}
    production = {'dataset_id': 'source', 'production_adoption_passed': True,
                  'paired_test_comparison': {'lower_95': [.01, 0, 0, 0]},
                  'null_trials': 199999, 'test_queries': 512, 'choices': {'candidate': candidate}}
    audit = {'dataset_id': 'source', 'anti_conservative_detected': False,
             'fingerprint': {'kind': 'Morgan'}, 'threshold': .45}
    assert approved_decision(retrieval, production, audit) == candidate
    with pytest.raises(ValueError):
        approved_decision(retrieval, {**production, 'null_trials': 999}, audit)
    with pytest.raises(ValueError):
        approved_decision(retrieval, production, {**audit, 'anti_conservative_detected': True})
    with pytest.raises(ValueError):
        approved_decision(retrieval, production, {**audit, 'threshold': .4})
    with pytest.raises(ValueError):
        approved_decision({**retrieval, 'retrieval_adoption_passed': False}, production, audit)


def test_predictive_null_check_includes_shared_reference_uncertainty():
    from sea.parameter_validation import predictive_excess
    # B=1, three independent continuous query scores: the sole reference's
    # rank is uniform among four values, so the rejection count is uniform.
    assert predictive_excess(3, 3, 1, .5) == pytest.approx(.25)
    assert predictive_excess(2, 3, 1, .5) == pytest.approx(.5)
    assert predictive_excess(1, 3, 1, .5) == pytest.approx(.75)
    assert predictive_excess(0, 3, 1, .5) == pytest.approx(1.)
    assert predictive_excess(0, 3, 1, .01) == 1.  # Below the MC resolution.
