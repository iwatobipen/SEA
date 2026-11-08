"""Independent scaffold validation of SEA parameters; no inactive-label inference."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import DataStructs

THRESHOLDS = [.20, .25, .30, .35, .40, .45, .50, .55, .56, .57, .60, .65, .70, .75, .80]
PROTOCOL = {
    'version': 'chembl-parameters-v1', 'seed': 20261108,
    'split': 'global achiral Murcko scaffold SHA256 first 8 hex digits modulo 5',
    'training_buckets': [2, 3, 4], 'selection_bucket': 1, 'test_bucket': 0,
    'thresholds': THRESHOLDS, 'representations': ['RDKit-path', 'Morgan'],
    'selection_metric': 'mean reciprocal rank of first known target',
    'queries_per_partition': 512, 'one_query_per_scaffold': True,
    'initial_null_trials': 8191, 'confirmatory_null_trials': 32767,
    'adoption': 'selection winner only; independent test paired scaffold-bootstrap 95% lower MRR gain > 0 and Recall@10 lower gain >= -0.02; otherwise retain baseline',
    'bootstrap_replicates': 10000,
    'interpretation': 'Known-active target retrieval only; unmeasured targets are not inactive. No precision, false-positive, FDR or binding-probability claims.',
}


def buckets(groups):
    return np.array([int(value[:8], 16) % 5 for value in groups], dtype=np.uint8)


def reference_split(data, smiles, groups):
    membership = buckets(groups)
    train_indices = np.flatnonzero(membership >= 2)
    train_lookup = {smiles[i]: j for j, i in enumerate(train_indices)}
    references = []
    for target in data['targets']:
        indices = [train_lookup[s] for s in target['smiles'] if s in train_lookup]
        if len(indices) >= 10:
            references.append({**target, 'indices': np.array(indices, dtype=np.int32)})
    return train_indices, references, membership


def choose_queries(smiles, groups, membership, references, bucket, limit, seed):
    """One global scaffold per query, sampling shuffled target anchors in rounds."""
    index = {s: i for i, s in enumerate(smiles)}
    truth = defaultdict(set)
    candidates = []
    for tid, target in enumerate(references):
        values = []
        for s in target['smiles']:
            i = index[s]
            if membership[i] == bucket:
                values.append(i)
                truth[i].add(tid)
        candidates.append(values)
    rng = np.random.default_rng(np.random.SeedSequence([seed, bucket]))
    candidates = [rng.permutation(values).tolist() for values in candidates]
    used, selected = set(), []
    while len(selected) < limit:
        before = len(selected)
        for tid in rng.permutation(len(references)):
            values = candidates[tid]
            while values and groups[values[-1]] in used:
                values.pop()
            if values:
                i = values.pop()
                used.add(groups[i])
                selected.append({'index': i, 'anchor': int(tid), 'truth': sorted(truth[i]),
                                 'scaffold': groups[i].decode()})
                if len(selected) == limit:
                    break
        if len(selected) == before:
            break
    return selected


def raw_scores(similarities, references, threshold):
    """Bincount adds weights sequentially, matching singleton raw_score order."""
    indices = np.concatenate([t['indices'] for t in references])
    labels = np.repeat(np.arange(len(references)), [len(t['indices']) for t in references])
    values = similarities[indices]
    return np.bincount(labels, weights=np.where(values >= threshold, values, 0.), minlength=len(references))


def null_prefixes(fingerprints, sizes, thresholds, trials, seed, directory, name, progress=print):
    """Correct independent uniform sets, shared prefixes across target sizes.

    These study streams differ from production's size-keyed streams. Each
    size/threshold marginal has the same null. Correlated candidates make paired
    comparisons possible; they are never counted as independent replicates.
    """
    sizes = np.array(sorted(set(sizes)), dtype=np.int64)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    identity = {'name': name, 'sizes': sizes.tolist(), 'thresholds': list(thresholds),
                'trials': trials, 'seed': seed, 'pool_size': len(fingerprints), 'stream': 'independent-query-uniform-right-prefix-v1'}
    token = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]
    path = directory / f'null-{token}.npy'
    if path.exists():
        return np.load(path, mmap_mode='r'), sizes
    output = np.lib.format.open_memmap(path.with_suffix('.tmp.npy'), mode='w+', dtype=np.float64,
                                       shape=(len(thresholds), len(sizes), trials))
    rng = np.random.default_rng(seed)
    thresholds = np.asarray(thresholds)
    for trial in range(trials):
        query = int(rng.integers(len(fingerprints)))
        right = rng.choice(len(fingerprints), int(sizes[-1]), replace=False)
        similarities = np.array(DataStructs.BulkTanimotoSimilarity(fingerprints[query], [fingerprints[i] for i in right]))
        values = np.where(similarities[None, :] >= thresholds[:, None], similarities, 0.)
        output[:, :, trial] = np.add.accumulate(values, axis=1)[:, sizes - 1]
        if (trial + 1) % 2048 == 0:
            output.flush()
            progress(f'{name}: {trial+1}/{trials} null trials', flush=True)
    # Sorting occurs after stream generation; generation identity includes B.
    output.sort(axis=2)
    output.flush()
    del output
    path.with_suffix('.tmp.npy').replace(path)
    path.with_suffix('.json').write_text(json.dumps(identity, indent=2) + '\n')
    return np.load(path, mmap_mode='r'), sizes


def retrieval(fingerprints, train_indices, references, queries, thresholds, null, sizes, progress=print):
    training = [fingerprints[i] for i in train_indices]
    size_lookup = {int(n): i for i, n in enumerate(sizes)}
    target_sizes = np.array([len(t['indices']) for t in references])
    size_positions = np.array([size_lookup[int(n)] for n in target_sizes])
    edge_indices = np.concatenate([t['indices'] for t in references])
    edge_labels = np.repeat(np.arange(len(references)), target_sizes)
    offsets = np.cumsum(np.r_[0, target_sizes[:-1]])
    output = {str(t): [] for t in thresholds}
    baseline = []
    ids = np.array([t['id'] for t in references])

    def metrics(order, truth):
        positions = np.flatnonzero(np.isin(order, truth)) + 1
        return [1 / int(positions.min()), *[float(np.mean(positions <= k)) for k in (1, 5, 10)]]

    for qn, query in enumerate(queries):
        similarities = np.array(DataStructs.BulkTanimotoSimilarity(fingerprints[query['index']], training))
        edge_values = similarities[edge_indices]
        maximum = np.maximum.reduceat(edge_values, offsets)
        baseline.append(metrics(np.lexsort((ids, -maximum)), query['truth']))
        for ti, threshold in enumerate(thresholds):
            raw = np.bincount(edge_labels, weights=np.where(edge_values >= threshold, edge_values, 0.), minlength=len(references))
            exceed = np.array([len(null[ti, pos]) - np.searchsorted(null[ti, pos], value, side='left')
                               if value else null.shape[2] for pos, value in zip(size_positions, raw)])
            # N is the same for every target in a query. Sorting integer k is
            # exactly equivalent to sorting (k+1)/(B+1)*N without roundoff.
            ranking = np.lexsort((ids, -raw, exceed))
            output[str(threshold)].append(metrics(ranking, query['truth']))
        if (qn + 1) % 64 == 0:
            progress(f'Retrieval: {qn+1}/{len(queries)} queries', flush=True)
    return {t: np.array(rows) for t, rows in output.items()}, np.array(baseline)


def paired_bootstrap(candidate, baseline, seed=20261108, repetitions=10000):
    """Queries have distinct scaffolds; resample paired scaffold units."""
    difference = np.asarray(candidate) - np.asarray(baseline)
    rng = np.random.default_rng(seed)
    means = []
    for _ in range((repetitions + 127) // 128):
        indices = rng.integers(len(difference), size=(min(128, repetitions-len(means)), len(difference)))
        means.extend(difference[indices].mean(axis=1))
    interval = np.quantile(means, [.025, .975], axis=0)
    return {'mean_gain': difference.mean(axis=0).tolist(), 'lower_95': interval[0].tolist(), 'upper_95': interval[1].tolist()}


def adoption_passes(comparison):
    return comparison['lower_95'][0] > 0 and comparison['lower_95'][3] >= -.02


def summarize(rows):
    return dict(zip(['mrr', 'recall_at_1', 'recall_at_5', 'recall_at_10'], np.mean(rows, axis=0).tolist()))


def approved_decision(retrieval, production, null_audit):
    """The frozen contender must pass both budgets and the independent null audit."""
    if retrieval['status'] != 'completed' or not retrieval['retrieval_adoption_passed']:
        raise ValueError('The fixed selection/holdout study did not approve a candidate.')
    if not production['production_adoption_passed'] or not adoption_passes(production['paired_test_comparison']):
        raise ValueError('Production-budget holdout confirmation failed.')
    if production['null_trials'] != 199999 or production['test_queries'] < 256:
        raise ValueError('Production confirmation requires the declared budget and cohort.')
    candidate = retrieval['selection_winner']
    if candidate != production['choices']['candidate']:
        raise ValueError('The production candidate differs from the frozen selection.')
    if null_audit['anti_conservative_detected']:
        raise ValueError('Independent empirical-null checks detected excess rejection.')
    if candidate['fingerprint'] != null_audit['fingerprint']['kind'] or candidate['threshold'] != null_audit['threshold']:
        raise ValueError('Null audit does not match the selected parameters.')
    if len({r['dataset_id'] for r in (retrieval, production, null_audit)}) != 1:
        raise ValueError('Validation reports refer to different reference cohorts.')
    return candidate



def predictive_excess(count, queries, references, level):
    from scipy import stats
    # With a continuous null, the acceptance threshold's tail probability is
    # Beta(r, B+1-r). Ties make inclusive Monte Carlo p-values conservative.
    r = int(np.floor(level * (references + 1)))
    return 1. if r == 0 else float(stats.betabinom.sf(count - 1, queries, r, references + 1 - r))
