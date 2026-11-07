# Validation record

Environment: Linux x86_64, Pixi 0.75.0, Python 3.12.16, and RDKit 2026.09.2. Dependencies are pinned in `pixi.lock`.

Recorded result: `pixi run test` completed with **10 passed in 2.67s**. `pixi install --locked --offline` also succeeded.

## Completed checks

- SEA calculations: structure standardization, salt removal, deduplication, similarity-threshold boundaries, and set-comparison symmetry.
- Equations from the original Supplementary Methods: Z-scores, very small p-values, and log-value precision.
- Monte Carlo evaluation: determinism with a fixed seed, p-value resolution, and set-size limits.
- SQLite fixture matching the ChEMBL schema: quality filters, median potency, target inclusion, provenance, release validation, and deterministic dataset IDs.
- API: authentication, job ownership isolation, input errors, CSV export, cookie sessions, Origin checks, and request-body limits.
- Jobs: exclusive claiming, cancellation of queued and running jobs, and recovery of interrupted jobs after a supervisor restart.
- Headless Chrome: sign-in, separate target searches for two molecules, and display of six result rows using the actual API and a two-process worker. Screenshots at desktop width (1280 px) and mobile width (390 px) were visually inspected.

To repeat the browser smoke check, install Node and Chrome as development tools and run:

```bash
pixi run python scripts/smoke.py
```

Example measurements for the synthetic demo (30 molecules, three sets, two separate molecule queries, and 999 Monte Carlo trials): 0.064 seconds of job computation and 88,016 KiB worker peak RSS. The smoke check, including demo construction, server/worker/Chrome startup, the search, and screenshot capture, took approximately 4.6 seconds. These measurements describe the small demo, not full-ChEMBL performance.

## Checks not yet performed

The full ChEMBL download and reference-database build have not been run. Biological prediction accuracy on real data has not been evaluated with the scaffold holdout. The build, calibration, and evaluation CLI commands are implemented; run `sea build` and `sea evaluate` on real data before making scientific performance claims. The synthetic demo's target names describe chemical groups, not actual proteins.

In the development sandbox, TestClient's internal event-loop communication stalled. API and browser checks were therefore run in an execution environment that allowed local communication.


## ChEMBL query performance follow-up

A full ChEMBL 37 reference dataset was subsequently built with 1,517 targets and 441,860 unique molecules. Investigation of delayed predictions found that Monte Carlo sampling was also being performed for zero-score targets. Since null raw scores are always nonnegative, these comparisons have the exact p-value 1 and now return `method=exact_zero` without sampling. Workers also reuse deterministic null distributions across jobs for the same dataset and size pair.

The follow-up regression suite passed all 12 tests in 2.87 seconds, including the exact zero-score result and reuse of the same Monte Carlo samples across jobs. On one observed single-molecule ChEMBL 37 query, scoring and significance evaluation across 1,517 targets took 49.64 seconds, excluding dataset loading: 1,469 comparisons used the exact zero-score result and 48 used Monte Carlo sampling. Latency depends on the query, ligand-set sizes, and worker cache state; this is not a general performance guarantee or a biological accuracy benchmark.

The web interface now displays waiting/running status and percentage progress, and disables submission buttons while a job is active to prevent accidental duplicate jobs.


## Target-family and species annotations

The full ChEMBL 37 reference dataset was annotated from its original SQLite source after verifying the source release and SHA256. All 1,517 targets received an organism and classification record. Checked examples include CHEMBL233 (Mu-type opioid receptor, GPCR), CHEMBL203 (Epidermal growth factor receptor, Kinase), and CHEMBL240 (KCNH2, Ion channel). All targets in this human-filtered reference database have organism Homo sapiens and taxonomy ID 9606.

The regression suite passed 15 tests in 3.18 seconds. Tests cover ancestry-based family labels, multiple classification paths, missing metadata, species from different organisms, immutable dataset preservation, source checksum validation, annotation-cache refresh, historical result enrichment, new job metadata, CSV columns, and target search. The existing seven ChEMBL jobs (140 returned target rows) were validated against the new response schema with no missing species or family fields; their stored scores and result records were not rewritten. Browser smoke checks also verified target ID, GPCR/Kinase/Ion channel labels, and species/taxonomy columns using a temporary rendering fixture that was not saved to any dataset.
