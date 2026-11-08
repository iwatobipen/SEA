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
