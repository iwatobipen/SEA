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
