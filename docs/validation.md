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

At the initial validation, the full ChEMBL download and reference-database build had not been run; the later ChEMBL 37 build is recorded below. Biological prediction accuracy on real data has not been evaluated with the scaffold holdout. The build, calibration, and evaluation CLI commands are implemented; run `sea build` and `sea evaluate` on real data before making scientific performance claims. The synthetic demo's target names describe chemical groups, not actual proteins.

In the development sandbox, TestClient's internal event-loop communication stalled. API and browser checks were therefore run in an execution environment that allowed local communication.


## ChEMBL query performance follow-up

A full ChEMBL 37 reference dataset was subsequently built with 1,517 targets and 441,860 unique molecules. Investigation of delayed predictions found that Monte Carlo sampling was also being performed for zero-score targets. Since null raw scores are always nonnegative, these comparisons have the exact p-value 1 and now return `method=exact_zero` without sampling. Workers also reuse deterministic null distributions across jobs for the same dataset and size pair.

The follow-up regression suite passed all 12 tests in 2.87 seconds, including the exact zero-score result and reuse of the same Monte Carlo samples across jobs. On one observed single-molecule ChEMBL 37 query, scoring and significance evaluation across 1,517 targets took 49.64 seconds, excluding dataset loading: 1,469 comparisons used the exact zero-score result and 48 used Monte Carlo sampling. Latency depends on the query, ligand-set sizes, and worker cache state; this is not a general performance guarantee or a biological accuracy benchmark.

The web interface now displays waiting/running status and percentage progress, and disables submission buttons while a job is active to prevent accidental duplicate jobs.


## Target-family and species annotations

The full ChEMBL 37 reference dataset was annotated from its original SQLite source after verifying the source release and SHA256. All 1,517 targets received an organism and classification record. Checked examples include CHEMBL233 (Mu-type opioid receptor, GPCR), CHEMBL203 (Epidermal growth factor receptor, Kinase), and CHEMBL240 (KCNH2, Ion channel). All targets in this human-filtered reference database have organism Homo sapiens and taxonomy ID 9606.

The regression suite passed 15 tests in 3.18 seconds. Tests cover ancestry-based family labels, multiple classification paths, missing metadata, species from different organisms, immutable dataset preservation, source checksum validation, annotation-cache refresh, historical result enrichment, new job metadata, CSV columns, and target search. The existing seven ChEMBL jobs (140 returned target rows) were validated against the new response schema with no missing species or family fields; their stored scores and result records were not rewritten. Browser smoke checks also verified target ID, GPCR/Kinase/Ion channel labels, and species/taxonomy columns using a temporary rendering fixture that was not saved to any dataset.

## Drawing and structure-file input

The regression suite passed 27 tests in 3.18 seconds after adding the structure-import API and locally bundled Ketcher 3.18.0 editor. New checks cover multi-record SDF, V2000/V3000 MOL, isomeric SMILES, salt fragments and charges, molecule names, repeated records, quoted CSV with embedded newlines, empty input, invalid record numbers, rejection of partial imports, record limits, rejection of query structures and mislabeled multi-record MOL files, authentication, browser-session request headers, and a prediction submitted from converted SDF structures.

The browser smoke check imported SDF and MOL into target prediction and SMI/CSV into the two comparison inputs, verified that an invalid SDF left the previous input intact, rejected an empty drawing, drew a carbon-carbon bond through Ketcher's actual toolbar and browser pointer events, and appended it exactly once despite a double click. Editor structures were also added to both comparison inputs. Two imported molecules were submitted through the real API and worker, producing six result rows with the existing target metadata columns. The check took 8.07 seconds and recorded no external HTTP requests. Desktop and mobile screenshots include both the editor dialog and completed results; the mobile dialog was checked for horizontal overflow, with scrolling contained inside the editor viewport.

A Python wheel build was checked to contain the structure-import module and all 21 bundled Ketcher assets and license notices. The wheel was approximately 9.18 MB; no frontend build or runtime Node installation is needed.

## Background service scripts

The regression suite passed 31 tests in 5.54 seconds after adding the start/stop controller. Service checks use real API and worker processes with an isolated data directory and port. They cover readiness, repeated startup without duplicate processes, repeated shutdown, termination of worker children, an occupied port, stale PID metadata pointing to an unrelated live process, and cleanup after API startup failure. Legacy integer PID files are also checked before upgrading their metadata.

Both executable shell wrappers were exercised through the locked Pixi environment from outside the repository, using a data-directory path containing spaces and a custom port with one worker. Start, repeated start, stop, and repeated stop succeeded. Shell syntax checks and the Pixi `start`/`stop` task help commands also passed.

## Original-paper comparison

The [comparison report](paper-comparison.md) and [machine-readable results](paper-comparison-results.json) record an executed fingerprint diagnosis and statistical replay of all 60,516 official supplementary rows. It distinguishes independent aggregation/tail checks, replay from supplied raw scores, and end-to-end biological validation. The latter remains unproven. The current ChEMBL dataset and service model were not modified.

The full regression suite passed 34 tests in 5.26 seconds. The added offline checks cover the reported path/Morgan score difference, independent bit-set aggregation, high-precision tail evaluation, and replay ranking with self-matches excluded. The documented optional Pixi environment was installed from the lock file; both cached-input and fresh official-download comparison runs completed and produced identical JSON results.

## English result explorer and UI refinement

The full regression suite passed 34 tests in 5.56 seconds. Application interface copy, new validation errors, worker errors, and new demo labels are English. Existing built-in demo target IDs receive English display labels without rewriting old datasets or analysis results. The calculation algorithms and ChEMBL reference data are unchanged.

The redesigned white-and-navy interface groups results by query and displays ranked target names, IDs, family badges, species, E-values, and highest saved contributing-pair Tanimoto values. The selected-target panel contains statistics, taxonomy, method/resolution information, molecular pairs, expandable SMILES, and copy controls. The browser check covers query switching, family and text filters, retained SEA ranks when display order changes, keyboard activation, saved-pair expansion, extreme-value display from log values, exact-zero results without fabricated similarity, empty bundles and filters, and authenticated JSON/CSV exports. Three actual demo jobs exercise separate-molecule prediction, two-set comparison, ensemble prediction, and history reopening. Rendering fixtures are temporary and are never stored as scientific results.

The updated smoke runner isolates its dataset, API key, jobs, and Chrome profile in a temporary workspace. SDF/MOL/SMI/CSV conversion and actual Ketcher toolbar/pointer drawing still pass, including rejection of invalid records and empty drawings and protection against double submission. No external web resources or browser runtime exceptions were observed. The run took 13.07 seconds, including a read-only visual review of a previously completed ChEMBL 37 result containing 20 returned targets. This UI review does not measure prediction accuracy or new ChEMBL query latency.

Desktop (1280 px) and mobile (390 px) screenshots were inspected, including long ChEMBL target names, GPCR/Kinase family badges, species, molecular structures, and the editor dialog. Desktop ChEMBL results fit without horizontal table scrolling; mobile results use stacked rows and a separate selected-target panel, with no document overflow. Screenshots are saved under `var/smoke-*.png` and are not committed. Supply `SEA_SMOKE_PREVIEW=/absolute/path/to/result.json` to the smoke runner for an optional visual review of an existing result bundle; it is rendered only in the isolated browser.

## Python client and usage guide

The full regression suite passed 41 tests in 5.98 seconds after adding `SEAClient`. Seven client checks cover single/batch/ensemble request grouping, generator inputs, authentication headers and proxy path prefixes, connection cleanup, polling timeouts and resumption without resubmission or cancellation, failed/cancelled jobs, structured HTTP errors, rejection of invalid local arguments before submission, and integration with the actual API and scoring engine. The integration check also covers UTF-8 BOM CSV input, target search, comparison, history, cancellation, and CSV export.

An additional check started an isolated API and one worker with a temporary synthetic dataset and API key. Separate-molecule prediction returned four rows, ensemble prediction returned two rows with a two-ligand query, set comparison and CSV export succeeded, and a two-record SDF was converted. Both `examples/python_client.py` and the existing raw-HTTP `examples/agent.py` completed successfully. The documented local two-set calculation also ran without an HTTP request. These checks took 4.77 seconds before service shutdown; all temporary services were then stopped. Existing application services, ChEMBL datasets, and jobs were not changed.

The [Python API guide](python-api.md) documents configuration, dataset selection, prediction and comparison, file input, job lifecycle, errors, result metadata, exports, and local calculation. Importing `SEAClient` was separately checked not to import RDKit. The client delegates normalization and scoring to the existing server; no calculation algorithm or statistical model changed.

## Worker batching, persistent background cache, and input previews

The full regression suite passed 59 tests in 6.05 seconds. New numerical checks compare batched fingerprint intersections with RDKit's scalar reference, including empty/full fingerprints, multiple query molecules, exact threshold equality and the next representable float, deterministic sampling for five size pairs, swapped set sizes, and exact cache reloads. Additional checks cover invalid/truncated cache regeneration, unwritable-cache fallback, cancellation before sampling and between intersection chunks, and authentication and salt/charge preservation in the input-depiction API. Existing prediction, batch, ensemble, comparison, cancellation, ownership, and service lifecycle checks still pass.

The executed [worker benchmark](performance.md) used ChEMBL 37 dataset `e898df53d9bb5240f087` (1,517 targets, 441,860 molecules). For `COc1ccc2[nH]cc(c2c1)C(CC)CN`, scalar cold scoring took 41.928 seconds and optimized cold scoring took 11.094 seconds, including index construction: a 3.78-fold improvement. Memory- and disk-cache paths took 0.126 and 0.132 seconds, respectively; these demonstrate reuse, rather than a speedup over the previous warm worker. All 1,517 target score/statistic/evidence records and all 46,953 sorted Monte Carlo scores were identical. Dataset loading took 7.375 seconds, separately from these scoring measurements. The script does not create jobs, API keys, or persistent application cache files. [Machine-readable results](worker-performance-results.json) record the configuration and timings.

The final browser smoke check passed in 16.02 seconds with isolated services and no external resources or JavaScript runtime exceptions. It verified typed-SMILES previews before any job was submitted, previews after loading an example and adding a drawing, full-width single-molecule cards, pagination in molecule order, invalid-SMILES feedback without input loss, clearing during pending rendering, file-import previews for all three input areas, and desktop/mobile overflow. Existing Ketcher drawing, structure-file conversion, result exploration, exports, comparison, ensemble, and history checks also pass. Desktop and mobile input screenshots (`var/smoke-input-desktop.png` and `var/smoke-input-mobile.png`) were visually inspected. Input depictions retain salts and charges; search normalization still runs separately.

The existing application API and workers were restarted after confirming there were no queued or running jobs. All 28 completed jobs remained in place. Live health, the preview API schema, updated UI assets, and depiction authentication were verified after startup. ChEMBL source data, fingerprint settings, calibration, trial counts, and completed results were not rewritten.

## Statistical precision and UX design study

The [design study](statistics-and-ux.md) separates proposed future behavior from executed measurements. `scripts/study_precision.py` evaluated two targets at fixed budgets of 999, 9,999, and 99,999 trials using the existing ChEMBL 37 reference and example query. For CHEMBL5043, sampling increased from 0.019 to 1.754 seconds and E changed from 47.027 to 43.235. For CHEMBL1946, sampling increased from 0.096 to 9.431 seconds and E changed from 92.537 to 109.694. The shared packed index took 2.412 seconds to construct after 7.343 seconds of data loading. The measurements show the cost/precision tradeoff, not improved biological prediction accuracy or independent replicates across trial budgets.

[Recorded results](precision-study-results.json) include exceedance counts and fixed-budget pointwise 95% Clopper–Pearson intervals for the limiting null tail probability. The interval arithmetic and zero-exceedance precision table were checked independently. These intervals are not sequential or simultaneous guarantees. The study used temporary caches and did not create jobs, issue keys, recalibrate data, change application defaults, or restart services. No application implementation changed in this documentation study, so the existing 59-test regression result remains the latest application validation.

## Implemented fixed-budget precision and interactive results

On November 8, 2026, the complete regression suite passed **84 tests in 6.95 seconds**. New checks cover deterministic prefix extension and shorter-budget reads against the unchanged scalar reference, interrupted-block RNG recovery, atomic cache reloads, corruption/checksum rejection, coordinated concurrent sampling, and screening reads while a refinement writer holds the stream lock. Interval tests independently invert binomial tails, including zero/all exceedances. Production profile budgets and versioned identities are asserted separately from the smaller lifecycle-test budgets.

Service tests verify evidence with pending P/E values, immutable 999-trial screening artifacts, automatic linked refinement, full comparison coverage despite top-K output limits, the original per-query target count in E, authenticated previews, idempotent refinement requests, reserved foreground claims, cancellation/resumption, legacy-database migration, and bounded reference preloading. Client checks cover evidence/screening/final waits, a single polling deadline across linked jobs, terminal errors with retained previews, invalid arguments before submission, and compatibility with legacy calls. OpenAPI describes the phase/profile/statistical fields; JSON and CSV retain the actual method, budget, counts, uncertainty and lineage.

The isolated real-process browser check completed Balanced **199,999-trial** searches, comparison, and ensemble workflows. It verified pause/resume, usable evidence during refinement, stable query/selection/order before updates, explicit application of final statistics, final export readiness, history reopening, English copy, desktop/mobile layout, structure-file inputs and local Ketcher drawing. The same runner completed a real **999,999-trial Research** refinement through the Python client and checked source immutability, idempotent requests and the documented `sea warm` command. These checks use production budgets and a synthetic 30-molecule reference; they do not measure biological prediction accuracy.

An isolated ChEMBL 37 run used the unchanged 1,517-target, 441,860-molecule reference and the previously benchmarked query. API preload readiness took 8.395 seconds. Cold time from submission to evidence was **0.733 seconds**, and the 999-trial screening artifact was ready after **12.409 seconds**. While background refinement was actively sampling beyond its initial prefix, a repeated foreground search returned evidence after **0.739 seconds** and screening after **0.804 seconds**. Screening artifacts remained unchanged. Refinement was deliberately paused with a largest saved prefix of 56,423 trials; the full cold ChEMBL refinement was **not** timed to completion. These are single local observations, not p50/p95 latency guarantees or evidence of biological calibration. The [machine-readable report](phased-performance-results.json) records the limits and settings.

Reproduce the latency check with a local immutable reference:

```bash
pixi run python scripts/benchmark_phases.py \
  --dataset-id e898df53d9bb5240f087 \
  --output var/phased-performance.json
```

The script uses temporary keys, jobs, services, and cache files. It hardlinks the existing reference files and leaves the lab database and reference unchanged. It needs an available local port (default 19031). The browser/SDK production-budget check is `pixi run python scripts/smoke.py` and additionally requires Chrome and Node.

Increasing the fixed budget addresses Monte Carlo resolution only. Independent null calibration across seeds, tail validation, biological scaffold-holdout evaluation, measured inactive data, and concurrent multi-user percentile benchmarks remain separate research work; no new sensitivity, false-positive, FDR, or binding-probability claim is made.

After validation, the lab API and two-worker supervisor were restarted while no jobs were queued/running. Hashes of all **30 existing completed result artifacts**, API-key count, and reference-file identities/sizes/timestamps were preserved. Live health, the typed precision schema, updated UI and preview authentication passed. Existing jobs retain the legacy profile; no lab jobs or keys were created for the measurements.


## ChEMBL 37 parameter selection and independent validation

Recorded on November 8, 2026. **104 regression tests passed** in the locked Pixi environment. New tests exercise explicit path/Morgan fingerprints, empirical Z including zero variance, immutable parameter variants, exact-cohort recipe reconstruction and rejection, approved-reference preloading, and both representations through legacy/Balanced prediction and comparison workers. Study tests check global scaffold exclusion, exact raw-score aggregation, uniform-set prefix sampling, known-active retrieval metrics, paired adoption criteria, and the beta-binomial null-audit calculation. The browser smoke check also passed the existing Balanced/Research workflow after dataset-defined fingerprint support.

The [parameter study](chembl37-parameter-validation.md) evaluated 30 fixed candidates on 512 selection scaffolds, froze one, and tested on 512 separate scaffolds with global removal from every training target/background. Morgan radius 2 / threshold 0.45 passed predeclared MRR and Recall@10 gates at both 32,767 and 199,999 null trials. The independent audit found anti-conservative ensemble Gumbel tails in the original model, so the approved variant uses empirical P/E for all positive scores. No excess rejection was detected in 24 declared empirical-null checks. This does not prove biological inactive calibration, universal tail validity, or family-wide improvement. Maximum Tanimoto outperformed the selected SEA model for single-query retrieval; small CYP450 and Transporter strata decreased.

The source dataset remains unchanged. New dataset `9d52be7ff211a45d2709` and model `59cb366620706b4e3870` retain the same source structures, targets, filters and activity provenance. The published recipe checks every regenerated fingerprint, the evidence digest and exact new content identity. Full measurements, protocols, null audits, family strata and limitations are linked from the study report. An isolated two-worker benchmark returned cold evidence/screening in 0.707/5.989 s; full cold refinement was paused and has no measured completion latency. The 20 most frequent actual target-size nulls were prepared at 199,999 trials outside interactive requests.

The idle lab services were restarted with the approved reference. All **30 existing completed result hashes**, all **6 API-key records**, and identities/sizes/timestamps of all **5 original dataset/manifest/annotation files** were preserved. Live health, preview/dataset authentication, OpenAPI and UI version 7 were checked. Application validation and study scripts created no lab jobs or keys.

The final isolated Chrome smoke run passed in 29.81 s, including approved-reference preference, retaining an explicit legacy choice, the long parameter label at mobile width, 199,999-trial Balanced completion, pause/resume, immutable screening artifacts, file import/drawing and 999,999-trial Research completion. Manifest fixtures are used only for display selection and never submitted as scientific jobs.
