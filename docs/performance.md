# Worker performance

The worker accelerates the Monte Carlo fallback while keeping its statistical
definition unchanged. It still uses the dataset's fingerprint, threshold, random
seed, trial count, background pool order, and sampling without replacement.
Gumbel evaluation, normalization, evidence, target ranking, and API result fields
are unchanged.

## What changed

- Fingerprints are packed into a NumPy array on the first uncached Monte Carlo
  calculation. Batched bit intersections replace repeated Python molecule-list
  construction and individual RDKit bulk calls. The on-bit-count ratio rejects
  pairs that cannot reach the threshold. Threshold equality and RDKit's
  empty-fingerprint behavior are retained.
- The same alternating random draws and left-major sequential double additions
  retain exact agreement with the scalar reference sampler. Trial counts are
  not reduced, and no approximate nearest-neighbor search is involved.
- Sorted null scores are reused in worker memory and saved under
  `$SEA_HOME/cache/null-v1/DATASET_ID/SMALL_SIZE-LARGE_SIZE.npy`. Dataset IDs
  include calibration settings and the reference pool. Different data directories
  and datasets have separate caches. Files contain background scores only,
  without user queries, API keys, or job results.
- Cache writes use temporary files and atomic replacement. Another worker or a
  restarted worker can load a completed distribution. Concurrent initial misses
  may both compute the same deterministic distribution. Invalid or incomplete
  arrays are regenerated; a cache-write failure does not fail an analysis.
- Target ligand lists are constructed once per batch and shared across query
  iterations. Progress writes are limited to approximately five per second,
  with a final update at completion. Cancellation checks remain active during
  sampling and between fingerprint-intersection chunks.

The previous worker already reused null distributions in its own memory. The
new disk cache extends reuse across processes and restarts; cached timings
should not be interpreted as the speedup over the previous warm worker.

## Recorded ChEMBL 37 measurement

On November 8, 2026, one run on the development Linux x86_64 machine searched
1,517 targets against 441,860 unique reference molecules using the query
`COc1ccc2[nH]cc(c2c1)C(CC)CN`, threshold 0.56, and 999 Monte Carlo trials:

| Scoring path | Seconds |
| --- | ---: |
| Scalar reference, empty cache | 41.928 |
| Optimized, empty cache including packed-index construction | 11.094 |
| Optimized, worker-memory cache | 0.126 |
| Optimized, disk cache with a new cache instance | 0.132 |

The cold calculation was **3.78 times faster**. Dataset loading took a further
7.375 seconds and is excluded from the table. Warm measurements reuse all
required background distributions and the already loaded dataset. The disk
measurement exercises reloading distributions into a new cache instance; it
does not include starting a new Python process or loading the dataset again.

All 1,517 target score/statistic/evidence records and all 46,953 sorted null
scores (47 distributions × 999 trials) matched the scalar reference exactly.
The query produced 1,470 exact-zero evaluations and 47 Monte Carlo evaluations.
These are performance/equivalence checks, not a biological accuracy benchmark.

## Reproduce the benchmark

Select an exact dataset ID from `pixi run sea datasets`, then run:

```bash
pixi run python scripts/benchmark_worker.py \
  --dataset-id YOUR_20_CHARACTER_DATASET_ID \
  --smiles 'COc1ccc2[nH]cc(c2c1)C(CC)CN' \
  --output var/worker-performance.json
```

Replace `YOUR_20_CHARACTER_DATASET_ID` with the chosen ID. Repeat `--smiles`
to benchmark an ensemble. The script loads the local dataset and compares the
unchanged scalar sampler with optimized cold, memory-cache, and disk-cache
paths. It asserts exact equality of every target result and every sorted null
sample, then reports timings and software versions. It uses a temporary cache,
creates no jobs or API keys, and removes its cache afterward. No server or
worker service is needed, and existing application services are left running.

The benchmark measures scoring and significance across all reference targets.
It excludes queue waiting, API validation, dataset loading, query normalization,
and job-database writes. Dataset-loading time is recorded separately. The
machine-readable [recorded results](worker-performance-results.json) identify
the query, immutable dataset/model IDs, threshold, trial count, and versions.

## Operational considerations

The first request handled by a worker still loads the reference dataset.
Packed fingerprints are allocated only when a requested null distribution is
missing on disk. For 441,860 molecules with 2,048-bit fingerprints, the packed
array alone uses approximately 108 MiB per worker, in addition to the existing
dataset and temporary batch arrays. A worker with only disk-cache hits does not
need this array. The existing per-process dataset-cache limit is two datasets.

More workers can process more independent jobs concurrently and also use more
memory. Worker count does not divide one job among those processes. Inputs
with different ensemble sizes may require new distributions; additional
Monte Carlo trials also increase computation time. These timings do not measure
biological prediction accuracy or guarantee latency for other queries.

After a code update, restart the API and workers using the documented
`scripts/stop.sh` and `scripts/start.sh` commands. Derived cache files can be
regenerated from the immutable dataset; retain the dataset, API-key database,
and completed results independently of this cache.
