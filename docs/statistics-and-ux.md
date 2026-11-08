# Statistical precision and interactive search

Status: design recommendation supported by a limited executed cost study.
This document proposes future behavior; the current search defaults and
statistical calculations have not been changed by this study.

## Recommendation

Return molecular evidence and a clearly provisional screening result promptly,
then refine statistics in a separate background task. Reuse versioned null
distributions across queries and workers. Start with a predetermined final
Monte Carlo budget; introduce adaptive stopping only after its statistical
guarantees have been implemented and tested.

This separates time to first useful result from time to a reproducible final
analysis. Numerical precision, calibration under the statistical background,
and biological prediction accuracy are three separate validation requirements.
Increasing the trial count addresses only the first.

## What the current installation can establish

The ChEMBL 37 reference contains 1,517 targets and 441,860 unique molecules.
There are 600 distinct target ligand-set sizes. Single-molecule queries use
Monte Carlo evaluation even though this dataset has an eligible fitted Gumbel
model for some larger ensembles. A single-molecule distribution must not be
replaced by the ensemble Gumbel fit without separate validation.

The existing [performance benchmark](performance.md) measured 11.094 seconds
for cold scoring and 0.126 seconds with all required null samples already in
memory. Dataset loading took another 7.375 seconds. These are local computation
measurements, not end-to-end UI latency guarantees.

An additional fixed-budget experiment used query
`COc1ccc2[nH]cc(c2c1)C(CC)CN` with the same fingerprint, threshold, seed,
reference pool, and inclusive score aggregation:

| Target | Ligands | E, 999 trials | E, 99,999 trials | Sampling time, 999 / 99,999 |
| --- | ---: | ---: | ---: | ---: |
| CHEMBL5043, Endoplasmic reticulum aminopeptidase 2 | 72 | 47.027 | 43.235 | 0.019 / 1.754 s |
| CHEMBL1946, Melatonin receptor type 1B | 687 | 92.537 | 109.694 | 0.096 / 9.431 s |

Sampling timings share one packed index and exclude loading (7.343 seconds),
index construction (2.412 seconds), queueing, and API/job-database operations.
The budgets share deterministic random streams and are not independent
replicates. The first target had the best screening E-value for this query;
neither example approaches E = 0.05 after refinement. More trials do not imply
smaller E-values or improve the fingerprint's biological discrimination.

The [machine-readable study](precision-study-results.json) also records
9,999-trial measurements, exceedance counts, and pointwise 95% intervals for the
limiting Monte Carlo tail probability. It is not a biological validation study.
Reproduce it without a server or changes to application data:

```bash
pixi run python scripts/study_precision.py \
  --dataset-id e898df53d9bb5240f087 \
  --target-id CHEMBL5043 --target-id CHEMBL1946 \
  --trials 999 9999 99999 \
  --output var/precision-study.json
```

The dataset ID above identifies this installation's immutable ChEMBL build.
Select an available ID and matching target IDs when using another dataset.
The script creates and removes a temporary cache and creates no jobs or keys.

## Proposed default workflow

| Phase | User experience | Statistical meaning |
| --- | --- | --- |
| Structure input | Immediate 2D review | Input depiction only |
| Evidence | Molecular matches appear as soon as raw scoring finishes | Similarity evidence; statistical ranking pending |
| Screening | Available coarse P/E estimates and method labels | Provisional 999-trial Monte Carlo estimates, or the applicable existing analytic method |
| Refinement | Results remain usable while a separate task runs | Predetermined 199,999-trial Monte Carlo computation for positive-score comparisons |
| Completed analysis | Retain/export an immutable final artifact | Actual budget, uncertainty, method, seed, and analysis version are recorded |

For the first release of this workflow, use a fixed **199,999-trial** refinement
profile and retain an optional **999,999-trial** research profile. These are
engineering starting points for a pilot, not validated universal sample-size
requirements. Zero raw scores retain their exact P = 1 and need no refinement.
Existing eligible analytic evaluations retain their method label and model
limitations; they are not presented as empirical Monte Carlo estimates.

Prioritize visible candidates, but eventually refine every positive-score
comparison before declaring the whole search complete. A top-20 display limit
does not establish that unrefined targets are unimportant. A cancelled or
budget-limited refinement retains a partial/unresolved status and its completed
screening artifact. Additional precision can be requested as a new analysis.

Reserve foreground capacity for new searches; background precision tasks must
not monopolize the current two-worker pool. Preload the selected reference data
and precompute commonly used single-query size pairs outside interactive
requests. Use the existing process pool and persistent queue as a starting
point, adding explicit priorities/worker roles rather than one long job that
blocks other users. A warm first useful display within two seconds is a
proposed acceptance target, not an observed guarantee. Cold-load readiness and
queue waiting must be reported separately.

Preserve row selection and scroll position as statistics arrive. Mark updated
rankings and let users apply them without unexpectedly moving the selected
target. Label examples include `Screening`, `Refining`, `Complete`,
`Resolution-limited`, and `Uncertain`. Completion means the declared numerical
procedure finished, not that a target has been experimentally confirmed.

## Precision and reporting rules

For fixed-budget Monte Carlo, retain `p_mc = (k + 1) / (B + 1)` and
`E_mc = N * p_mc`, where `k` counts null scores at least as large as the observed
score, `B` is the actual trial count, and `N` is the original search target count.

| Trials B | Minimum reported P | Minimum reported E, N = 1,517 | Zero-exceedance pointwise 95% one-sided E-tail upper bound |
| ---: | ---: | ---: | ---: |
| 999 | 0.001 | 1.517 | 4.5423 |
| 99,999 | 0.00001 | 0.01517 | 0.04545 |
| 199,999 | 0.000005 | 0.007585 | 0.02272 |
| 999,999 | 0.000001 | 0.001517 | 0.004545 |

The upper-bound column is `N * (1 - 0.05**(1/B))` for zero exceedances. It bounds
the limiting null tail expectation at one predetermined sample size for one
comparison. It is not a simultaneous confidence statement over all targets,
an anytime-valid bound, or a binding probability. A fine numerical floor does
not establish small relative error near a decision threshold.

Show the point estimate, actual trial/exceedance counts, resolution, and a
clearly identified Monte Carlo uncertainty interval. The interval describes
simulation uncertainty under the chosen null; it excludes uncertainty about
whether that null represents real target nonbinding. Interval overlap with an
E threshold should be visible as uncertainty about a stable classification.
An interval can remain ambiguous at the largest budget.

Never replace zero exceedances by P = 0 or by the unsupported statement
`P < 1/(B+1)`. Never reduce N from 1,517 to the number of returned or refined
targets. Keep per-query multiplicity explicit; an analysis claiming significance
across a batch also needs an explicitly defined batch-wide testing family.
E-values remain expected counts, not false-discovery rates or binding
probabilities.

Ordinary pointwise intervals cannot justify repeatedly looking and stopping
when a preferred value appears. An ordinary `(k+1)/(B+1)` value at such a chosen
stopping time must not silently become a final fixed-budget P-value. A later
adaptive implementation should use a published anytime-valid Monte Carlo
method or confidence-sequence/resampling-risk method, with its sampling
assumptions checked and simulation-error guarantees allocated across the
declared testing family. The sequential-testing literature's betting
**e-values** are a different quantity from SEA **E-values**.

## Provenance, cache, and compatibility

Keep reference data, fingerprints, threshold fitting, and statistical precision
profiles separate. The existing `sea calibrate` path calls `save_dataset`, which
re-standardizes/fingerprints the pool and refits calibration. Increasing a
sampling budget alone should not require that expensive reference rebuild.

Use a versioned statistics-profile ID containing the base dataset/model IDs,
null definition, implementation version, seed/stream policy, and planned
budget. Record the effective profile, actual budget, exceedance counts,
interval definition, completion state, and source screening-job ID in API,
Python, JSON, and CSV results. Preserve existing completed results as immutable
artifacts; a refined result is a linked analysis, not an unexplained rewrite.

The current cache stores one fixed-length sorted array per size pair and
checks its shape against the dataset's trial count. New precision profiles
need distinct paths and metadata. True resumable refinement also needs saved
RNG state/sample count, an atomic append/checkpoint protocol, and coordinated
work claims. Sorting alone loses the order needed to reconstruct a stream
checkpoint. Separate cached sampling blocks can support different fixed final
budgets without regenerating all earlier samples. Do not simply change the
trial count associated with an existing cache file.

Legacy API/Python calls should retain their documented behavior while a new
precision profile is introduced explicitly. The new workflow needs an explicit
way to wait for screening versus final results. Agents and exports must be able
to require a completed precision artifact and detect unresolved results.

## Validation before a new default

1. Reproduce fixed-budget scalar results, RNG streams, threshold boundaries,
   cache reloads, cancellation/resumption, and all reported interval formulas.
2. Check null P-value calibration and Monte Carlo variability with independent
   simulations. Validate any sequential procedure over complete stopping paths
   and the full testing family, rather than checking each iteration separately.
3. Evaluate biological retrieval with global scaffold holdouts, removing held-out
   molecules from all target sets and fitting the background only on training
   data. Stratify by query fingerprint density, target ligand count, and family.
   Use measured inactive data for false-positive/precision claims; unmeasured
   ChEMBL pairs must not be treated as inactive. The existing holdout command
   reports known-target retrieval and does not establish those error rates.
4. Verify that the precision profile improves stability near the intended E
   thresholds across multiple seeds. A satisfactory global Gumbel KS statistic
   alone does not validate extrapolation to an extreme single-query tail.
5. Measure warm/cold/uncached p50 and p95 end-to-end latency under concurrent
   users, including worker starvation, time to first evidence, refinement
   duration, mobile display, and interruption/recovery. Keep final provenance
   and unresolved states correct when the precision task outlasts a UI session.

Statistical calibration and biological validation should determine which
interpretive thresholds the product presents. A smaller displayed number is
not a success criterion by itself. Query-conditioned backgrounds, different
fingerprints, tail fitting, importance sampling, or FDR procedures can be
separate research experiments; none should silently replace the current null
or be selected merely because it produces smaller P/E values.

## Primary methodological sources

- Phipson and Smyth (2010), [Permutation P-values Should Never Be Zero](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf): finite randomized-test P-values and multiple-testing consequences. Formal validity remains conditional on the null/sampling assumptions.
- Stoepker and Castro (2024), [Inference with Sequential Monte-Carlo Computation of p-values: Fast and Valid Approaches](https://arxiv.org/abs/2409.18908): distinguishes valid decisions, P-value estimation, sequential stopping, and reporting.
- Fischer and Ramdas (2025), [Sequential Monte Carlo testing by betting](https://doi.org/10.1093/jrsssb/qkaf014): anytime-valid Monte Carlo tests under exchangeability; a possible later alternative, not implemented here.
