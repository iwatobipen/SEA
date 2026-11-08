# SEA Open

An open-source implementation of the Similarity Ensemble Approach (SEA), which relates protein targets through the chemical similarity of their ligand sets. Built with Pixi, RDKit, and FastAPI, it provides an English web interface, API key authentication, persistent jobs, and separate worker processes.

Original paper: Keiser et al., *Relating protein pharmacology by ligand chemistry*, Nature Biotechnology 25, 197–206 (2007), https://doi.org/10.1038/nbt1284

## Getting started

Requirements: Linux x86_64 and Pixi.

```bash
pixi install --locked
pixi run demo
pixi run sea key alice
```

Save the API key printed by the last command. Only its SHA256 hash is stored in the database. Start the API server and workers in the background:

```bash
./scripts/start.sh
```

Open http://127.0.0.1:8000 and sign in with your API key. Interactive API documentation is available at `/docs`, and the machine-readable specification is at `/openapi.json`. The demo contains synthetic alcohol, carboxylic acid, and aromatic compound sets; it cannot predict biological targets.

Reference datasets and the service database are stored in `SEA_HOME` (default: `var/`). Run `pixi run sea datasets` to list dataset IDs. If your sandbox does not allow writes to `/tmp`, create a temporary directory inside the workspace and use `TMPDIR="$PWD/.cache/tmp" pixi …`.

## Starting and stopping services

Run these commands from the repository directory after `pixi install --locked` and the initial dataset/API-key setup above:

```bash
./scripts/start.sh
./scripts/stop.sh
```

`start.sh` uses the locked Pixi environment, starts the API and a supervisor with two worker processes, and waits for API health and worker readiness before returning. The default URL is http://127.0.0.1:8000. Services continue running after the command or terminal exits. Repeating the start command keeps the existing processes; repeating the stop command succeeds when they are already stopped. The scripts can also be invoked by absolute path from another directory.

`stop.sh` sends SIGTERM to the recorded API and worker supervisor and waits up to 30 seconds for shutdown. It verifies the executable, command, checkout, data directory, and process start time before signalling a PID. Stale PID files are discarded without signalling unrelated processes. Data, datasets, API keys, and completed results remain on disk. An interrupted running job is requeued when the worker supervisor next starts.

Logs are appended to `var/logs/api.log` and `var/logs/worker.log`; PID metadata is stored alongside them. Inspect the logs with:

```bash
tail -f var/logs/api.log var/logs/worker.log
```

Configure the port and worker count on startup, or allow more time for shutdown:

```bash
./scripts/start.sh --port 8080 --workers 4
./scripts/stop.sh --timeout 60
```

Stop services before starting them with different settings. Use `--host` to change the API bind address. For a separate data directory, use the same `SEA_HOME` for both commands; relative paths are resolved against the repository directory:

```bash
SEA_HOME="$PWD/var-lab" ./scripts/start.sh
SEA_HOME="$PWD/var-lab" ./scripts/stop.sh
```

The custom directory must have its own reference dataset and API keys. Logs and PID metadata are stored under its `logs/` directory. Environment settings such as `SEA_SECURE_COOKIE` and query limits are inherited by newly started services. The scripts do not load `.env` files automatically.

The equivalent Pixi tasks are `pixi run start` and `pixi run stop`; they also accept `--home /absolute/data/path`. For foreground development, run `pixi run serve` and `pixi run worker` in separate terminals and stop each with Ctrl+C. Stop foreground instances before switching to the background scripts. An occupied API port causes startup to fail; processes started during a failed startup are stopped, while previously running services are retained.

## Structure input

Target prediction and both sides of set comparison support three input methods:

- Paste SMILES, one molecule per line.
- Click **Draw structure** to open the bundled [Ketcher molecular editor](https://github.com/epam/ketcher). Draw a molecule and click **Add to input** to append it to the selected input. Add separate compounds one at a time; disconnected fragments in one drawing represent one input structure. Closing the editor without adding leaves the input unchanged.
- Select a structure file to replace the selected input after validation. Supported extensions are `.sdf`/`.sd` (multiple records), `.mol` (one molecule, V2000 or V3000), `.smi`/`.smiles`/`.txt` (SMILES with an optional whitespace-separated name), and `.csv` (a `smiles` column, with an optional `name` column). CSV supports quoted fields and embedded newlines. Files must contain UTF-8 text.

Imported structures appear as isomeric SMILES in the input area for review before submitting a search or comparison. Record order and repeated molecules are retained during import. Stereochemistry, charges, and disconnected fragments are preserved in these SMILES; the existing SEA standardization and ensemble deduplication run when the job is submitted. SDF properties other than the molecule name are not used in the search.

A 2D preview appears below each input after you pause typing SMILES. It also updates after loading an example, importing a file, or adding a drawing. Target search and both comparison sets show four molecules per page, with Previous/Next controls for larger inputs. These previews show the structures before search standardization, including salts and charges. Invalid SMILES display a rendering message without changing the input or creating a job. Search validation still checks every submitted molecule.

Use concrete small-molecule structures. Query atoms/bonds, R groups, and reactions cannot be used as SEA ligands. On narrow screens, scroll horizontally inside the Ketcher editor to reach the full toolbar and drawing area.

An invalid structure rejects the entire import and displays its one-based record number. Previously entered structures remain intact. CSV record numbers count data rows after the header; SMILES text uses line numbers; SDF uses molecule record numbers. Imports allow up to 1,000 records by default (`SEA_MAX_SET_MOLECULES`). Prediction and comparison limits still apply when submitting the resulting molecules. The browser allows files up to 10 MB, and the API request limit includes JSON encoding overhead, so files close to this size may need to be split.

The editor and all its assets are hosted locally, with no CDN or external structure-conversion service. The Ketcher version, archive checksum, and Apache license are recorded in [the bundled editor notice](src/sea/web/vendor/ketcher/README.md). No additional Pixi dependencies or frontend build step are required.

## Exploring search results

The English interface uses a white background with dark-blue controls. **Target search** supports separate searches for each input molecule or a single ensemble search. Choose 20, 50, or 100 returned targets per query. **Load example** fills the input with an example appropriate to the selected library. ChEMBL libraries appear before synthetic demos, with newer release numbers first.

Completed results show the reference library, fingerprint, query structure, and computation time. Use **Viewing query** to move between queries. The ranked list puts the target name, ID, protein-family badges, and species together, alongside the E-value and highest Tanimoto among saved contributing pairs. Select a target by clicking its name or focusing it and pressing Enter to inspect its statistics and molecular matches. Taxonomy IDs and full classification paths are available in the selected result and classification tooltips.

Filter the returned results by target name, ID, family, or species. Display sorting can use SEA rank, highest saved Tanimoto, or target name; the original SEA rank remains visible. These controls operate on the current query's returned results. CSV and JSON downloads contain the complete stored analysis across all queries, independent of display filters.

The detail panel shows query/reference structures side by side, four-decimal Tanimoto values, expandable SMILES with copy buttons, and up to three saved pairs initially. **Show all saved pairs** reveals the remaining saved evidence. Saved pairs may be a subset of all contributing pairs; both counts are displayed. When no pair reaches the calibrated threshold, the UI leaves the best-Tanimoto field empty instead of inferring a similarity from missing evidence.

**Calculation details & provenance** contains the dataset manifest, preprocessing records, threshold, execution information, and job ID. Extremely small p-values and E-values are displayed using their stored natural logarithms when numeric values underflow to zero. Monte Carlo resolution and raw-score tie breaking are explained alongside the selected result. E-values and Tanimoto similarities are not binding probabilities.

At narrow screen widths, the list and selected detail stack vertically, and target rows become compact cards. **Back to results** returns from a selected target to the list. **Job history** reopens previous analyses using the same result explorer. Interface copy, validation errors, and new synthetic-demo names are English; original datasets and stored analyses are preserved. Older built-in demo names are translated for display using their stable demo IDs.

## Statistical precision without blocking exploration

The web UI defaults to **Balanced**. It publishes molecular evidence first (P/E are pending), saves a **999-trial screening** artifact, and automatically queues a separate **199,999-trial refinement**. **Research** uses 999,999 final trials. Eligible calibrated Gumbel evaluations keep their analytic method; zero raw scores keep exact P = 1 and require no sampling. This changes numerical precision without rebuilding the reference or changing its fingerprints, threshold, or statistical null.

You can start another search while refinement runs. The current query, rows, and selection remain in place until you choose **Apply latest statistics** or **Apply final results**. Use **Pause refinement** and **Resume refinement** to interrupt and continue the declared budget. Closing the browser leaves server work running; reopen the refinement from **Job history**. With two or more workers, one worker is reserved for new searches/screening, and all other workers prefer queued foreground jobs before refinements. One-worker installations process jobs serially.

The result banner distinguishes **Evidence**, **Screening**, **Refinement**, and **Final**. Each row records its own phase, actual trial count, exceedance count, and pointwise 95% Clopper–Pearson interval for the limiting null-tail probability and its E-value scaling. **Resolution-limited** means zero exceedances; the estimate is `(k + 1) / (B + 1)`, never zero or an upper bound. **Uncertain** means the interval crosses the E = 0.05 reference threshold. This threshold is an interpretation aid, not a biologically validated cutoff. Final means the fixed computation finished; uncertainty and model limitations may remain.

Refinement covers every positive-score Monte Carlo comparison, including targets outside the displayed top K. E always uses the original target count per query; a batch does not silently become a batch-wide test. For 1,517 targets, Balanced has a minimum reported E of 0.007585, versus 1.517 at 999 trials. More trials need not lower E or improve biological discrimination. Cold refinement across many target sizes can take minutes or longer; evidence remains usable during that work. See [statistical precision and UX](docs/statistics-and-ux.md) for assumptions and limits.

Completed screening and final JSON/CSV artifacts are separate and immutable. Exports identify the profile/version, planned and actual budgets, counts, uncertainty, source/refinement job IDs, and completion state. Downloads are available for stored completed artifacts; previews are available through the preview API. Existing API/Python callers retain the **legacy** default; opt in with `statistics_profile="balanced"` or `"research"` and explicitly choose which phase to wait for.

The API and worker supervisor preload one reference before readiness by default (newest ChEMBL library, otherwise newest available). Set `SEA_PRELOAD_DATASETS=DATASET_ID` or a comma-separated list of at most two IDs to choose references; `SEA_PRELOAD_DATASETS=none` disables preloading. Startup includes this load rather than charging it to the first selected-library request. Other references still load on demand.

Precompute frequent target-size distributions outside interactive requests:

```bash
pixi run sea warm DATASET_ID --profile balanced --query-size 1 --limit 20
```

This prepares the 20 most frequent actual target sizes, retains atomic sampling checkpoints, and does not create a dataset or a job. It can be interrupted and rerun. Different profiles reuse the same deterministic unsorted prefix and saved RNG state. The `cache/null-v2/` cache is separate from the legacy cache, coordinates writers, and allows screening to read completed prefixes while refinement extends them. Set-size pairs, dataset, null definition, implementation version, threshold and seed identify streams. Worker sample-array memory is capped at 64 MiB per reference cache. Run warming before interactive use when practical; cold work still consumes CPU even with a reserved foreground worker.

## Building a ChEMBL reference database

Choose a release from the [official distribution](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/). The `--sha256` argument must be the SHA256 checksum of the **compressed SQLite archive `chembl_<release>_sqlite.tar.gz`**, exactly as downloaded. It is not the checksum of the extracted `.db` file, `checksums.txt`, or another ChEMBL download such as the SDF, MySQL, or PostgreSQL archive.

In the release directory `chembl_<release>/`, open `checksums.txt` and find the row whose filename is `chembl_<release>_sqlite.tar.gz`. Copy the 64-character hexadecimal value from that row's `SHA256` column.

| Placeholder | Meaning |
| --- | --- |
| `RELEASE_NUMBER` | Numeric ChEMBL release, such as `36`. |
| `SQLITE_ARCHIVE_SHA256` | SHA256 for `chembl_RELEASE_NUMBER_sqlite.tar.gz` from that release's `checksums.txt`. |

```bash
pixi run sea download --release RELEASE_NUMBER --sha256 SQLITE_ARCHIVE_SHA256
pixi run sea build --sqlite var/downloads/chembl_RELEASE_NUMBER.db --release ChEMBL_RELEASE_NUMBER
```

For example, the [ChEMBL 36 checksum file](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/checksums.txt) lists the following entry for the [SQLite archive](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/chembl_36_sqlite.tar.gz):

```text
b25820eef0f0481ad7712bdf4bac3b45f354e3cbacb76be1fdbf4205d6b48fb9  chembl_36_sqlite.tar.gz
```

The corresponding commands are:

```bash
pixi run sea download --release 36 --sha256 b25820eef0f0481ad7712bdf4bac3b45f354e3cbacb76be1fdbf4205d6b48fb9
pixi run sea build --sqlite var/downloads/chembl_36.db --release ChEMBL_36
```

The download command fetches `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/chembl_36_sqlite.tar.gz`, verifies the archive's SHA256, and then extracts the SQLite database into `var/downloads/`. It prints the extracted database path; use that path for `--sqlite`. With a custom `SEA_HOME`, downloads are stored under `$SEA_HOME/downloads/` instead.

This example pins release 36; it does not select the latest release. For another release, change both the release number and its matching SQLite archive checksum. Allow several gigabytes for the download, plus extraction space.

You can also pass an existing official SQLite database directly to `sea build --sqlite`. The builder checks the requested release against the database's version table. Its provenance manifest records the SHA256 of the extracted SQLite input file, which differs from the archive checksum used by `sea download`.

Default filters are defined in `config/chembl.json`: human single-protein targets, assay confidence 9, binding assays, Ki/Kd/IC50 measurements with an exact equality relation, and valid pChEMBL values. Measurements with quality flags or marked as potential duplicates are excluded. Parent structures are standardized, and median potency is calculated separately for each target, standardized structure, and measurement type. A ligand is considered active if at least one measurement type has a median pChEMBL ≥ 6. Targets must have at least 10 active ligands. Measurement types are not averaged together.

Each immutable JSON dataset records targets, ligands, measurement provenance, the source SHA256 checksum, filters, exclusion counts, the RDKit version, and calibration details. Changes to data, settings, or calibration produce a new dataset ID. Rebuilding with identical inputs and calibration produces the same ID.

```bash
pixi run sea calibrate DATASET_ID --sizes 1,10,30,100,300,1000 --repeats 128
pixi run sea evaluate DATASET_ID --output var/evaluation.json
```

Calibration for larger sets takes longer. The initial build calibrates set sizes up to 100 molecules and uses Monte Carlo evaluation for larger sets. The evaluation command splits compounds by Murcko scaffold across all targets and recalibrates the background model using only the training partition. It compares Recall@1/5/10 and MRR for known-target retrieval against a maximum-Tanimoto baseline. Unmeasured pairs are not treated as negatives. Evaluation uses at most 100 queries by default.

## Target families and species

ChEMBL search results display the preferred target name, ChEMBL target ID, protein family, organism name, and NCBI taxonomy ID in the web interface. The same information is available in JSON results, CSV exports, and `/api/v1/targets`. Family labels come from ChEMBL's protein classification hierarchy, linked through target components. GPCR, Kinase, Ion channel, Nuclear receptor, and Transporter labels summarize their corresponding ancestors; other targets retain a ChEMBL category such as Protease or Hydrolase. All assigned classification paths are preserved, including multiple classifications. Families are not inferred from target names.

New `sea build` runs include this metadata automatically. To add annotations to an existing dataset without rebuilding fingerprints or recalibrating scores:

```bash
pixi run sea datasets
pixi run sea annotate DATASET_ID --sqlite var/downloads/chembl_37/chembl_37_sqlite/chembl_37.db
```

Replace `DATASET_ID` with the 20-character ID of the matching ChEMBL dataset printed by the first command. Use the same source SQLite file used to build that dataset. The annotation command verifies both the release and the source SHA256, then writes an atomic annotation file under `$SEA_HOME/annotations/`. The original dataset ID, fingerprints, calibration, stored scores, and completed jobs remain intact. Both existing and future job results are enriched when retrieved. Classification-file updates are picked up without restarting the service; refresh the browser and reopen an existing result to see them.

The JSON `target` object includes `id`, `name`, `family`, `families`, `protein_classes` (an array of classification paths), `organism`, and `tax_id`. Results enriched by a separate annotation file include its provenance in `target_annotations`. CSV exports add `target_family`, `target_families`, `target_protein_classes`, `target_organism`, and `target_tax_id`. Target search also matches family names and species. Missing classifications are shown as Unclassified; unavailable species information is left empty. Synthetic demo sets do not receive biological family or species assignments.

The default reference build filters to human proteins, so its annotated species is normally Homo sapiens (9606). Annotations still use each target's actual ChEMBL organism and taxonomy ID rather than assuming the configured species.

## Python API

The [Python API guide](docs/python-api.md) covers target searches, batches and ensembles, set comparisons, SDF/MOL/SMILES/CSV input, job resumption, exports, errors, target metadata, and local calculations without the server. The client uses the same server and workers as the web application and is available in the existing Pixi environment:

```python
import os
from sea import SEAClient

with SEAClient() as client:  # Reads SEA_API_KEY and optional SEA_URL.
    result = client.predict(
        "COc1ccc2[nH]cc(c2c1)C(CC)CN",
        dataset_id=os.environ["SEA_DATASET_ID"],
        top_k=10,
    )

for row in result["results"]:
    print(row["target"]["name"], row["target"]["family"], row["target"]["organism"], row["e_value"])
```

Use an issued API key and the exact dataset ID printed by `pixi run sea datasets`; a release name such as `ChEMBL_37` is not a dataset ID. Run a saved script with `pixi run python your_script.py`. After setting `SEA_API_KEY` and `SEA_DATASET_ID`, the complete runnable example is `pixi run python examples/python_client.py`. `predict()` waits for completion; use `submit_prediction()` and `wait()` separately to retain a job ID for resumption or CSV export.

## Agent and raw HTTP access

Set `SEA_API_KEY` to your API key and `SEA_DATASET_ID` to a dataset ID.

```bash
curl -H "X-API-Key: $SEA_API_KEY" http://127.0.0.1:8000/api/v1/datasets
```

[examples/agent.py](examples/agent.py) submits a request, polls its status, and retrieves the result:

```bash
pixi run python examples/agent.py
```

For `POST /api/v1/predictions`, `queries` is an array of ligand sets. `[["CCO"],["CCCO"]]` searches each molecule separately; `[["CCO","CCCO"]]` searches one ensemble. `POST /api/v1/comparisons` accepts `left` and `right` ligand sets. Both endpoints require `dataset_id` and return HTTP 202 with a `job_id`.

Agents can convert structure files with authenticated `POST /api/v1/structures/import`. Send JSON with `format` (`sdf`, `mol`, `smiles`, or `csv`) and `content` (the file text). The response contains `format`, `count`, and `molecules`, each with `record`, `name`, and `smiles`. Submit these SMILES through the existing prediction or comparison endpoints:

```python
import os
from pathlib import Path
import httpx

headers = {"X-API-Key": os.environ["SEA_API_KEY"]}
base = "http://127.0.0.1:8000/api/v1"
converted = httpx.post(
    f"{base}/structures/import", headers=headers,
    json={"format": "sdf", "content": Path("molecules.sdf").read_text(encoding="utf-8")},
)
converted.raise_for_status()
smiles = [m["smiles"] for m in converted.json()["molecules"]]
job = httpx.post(
    f"{base}/predictions", headers=headers,
    json={"dataset_id": os.environ["SEA_DATASET_ID"], "queries": [[s] for s in smiles]},
)
job.raise_for_status()
print(job.json()["job_id"])
```

Use `"queries": [smiles]` for a single ensemble. Invalid imports return HTTP 422 with structured error details, including `records` when individual structures are invalid; they never return a partial molecule list or create a job.

- `GET /api/v1/jobs/{id}`: job status and progress.
- `GET /api/v1/jobs/{id}/result`: JSON results; add `?format=csv` for CSV.
- `GET /api/v1/jobs/{id}/preview`: latest evidence/partial statistics; HTTP 409 until a snapshot exists.
- `POST /api/v1/jobs/{id}/refine`: request/retrieve a linked refinement with `{"statistics_profile":"balanced"}` or `"research"`; one job per source/profile.
- `POST /api/v1/jobs/{id}/resume`: requeue a cancelled or failed refinement, retaining checkpoints.
- `POST /api/v1/jobs/{id}/cancel`: cancel a search or pause a refinement.
- `GET /api/v1/jobs`: your most recent 100 jobs.
- `GET /api/v1/targets?dataset_id=...&q=...`: search reference targets.
- `GET /api/v1/structure?smiles=...`: authenticated SVG depiction after standardization; add `standardize=false` to depict the entered structure, including its salts and charges.

Default limits are 100 molecules per prediction request, 1,000 molecules per set comparison side, and a 10 MB request body. Configure these with `SEA_MAX_QUERY_MOLECULES`, `SEA_MAX_SET_MOLECULES`, and `SEA_MAX_UPLOAD_BYTES`. Invalid inputs return HTTP 422, authentication failures return 401, another user's jobs return 404, and requests for unfinished results return 409.

## Sharing within a lab

Run the API and `pixi run worker` on the same host. The worker supervisor starts two processes by default; use `pixi run sea worker --count 4` to change this. Restarting the supervisor requeues interrupted jobs. Computation is deterministic, and retrying a job has no external write effects. Back up the service database and datasets. Jobs are not deleted automatically.

Issue a separate API key for each user or agent. Browser sign-in uses an eight-hour HttpOnly, SameSite=Strict session cookie; API keys are not saved in browser storage. For shared deployments over TLS, set `SEA_SECURE_COOKIE=1` and configure a TLS reverse proxy using [docs/Caddyfile](docs/Caddyfile) as a starting point. Key issuance and database construction are administrator CLI operations. Bind the API to 127.0.0.1 and expose only the proxy to the lab network.

## Validation and limitations

```bash
pixi run test
```

See [docs/method.md](docs/method.md) for calculation details and [docs/validation.md](docs/validation.md) for the recorded checks. RDKit and ChEMBL replace the original Daylight fingerprints and MDDR dataset, so numerical agreement with the original paper is not claimed. E-values are not binding probabilities. Legacy calls and preliminary screening use 999 Monte Carlo trials (minimum P = 0.001). The web UI automatically follows screening with fixed-budget refinement; this increases numerical resolution without establishing biological accuracy.

An executed [comparison with the original paper](docs/paper-comparison.md) diagnoses fingerprint-dependent Tanimoto differences and replays all 60,516 published comparisons from their raw scores. It verifies statistical ordering, with explicit numerical discrepancies and limits; it does not establish equivalent biological prediction accuracy. Reproduce it in the optional Pixi environment:

```bash
pixi install --locked -e paper-validation
pixi run --locked -e paper-validation python scripts/compare_paper.py --download --output var/paper-comparison.json
```

Zero-score comparisons return the exact p-value 1 without Monte Carlo sampling (`method=exact_zero`). Workers use batched fingerprint intersections for Monte Carlo evaluation and retain deterministic background samples in memory and under `$SEA_HOME/cache/null-v1/`. Files are isolated by the immutable dataset ID and set sizes, shared by workers in the same deployment, and reusable after a restart. The web interface displays queue state and progress and prevents repeated submissions while a job is active.

See [worker performance](docs/performance.md) for the measured speedup, cache behavior, memory use, and a reproducible benchmark. After updating a running installation, restart its API and workers to load the new code:

```bash
./scripts/stop.sh
./scripts/start.sh
```

The [statistical precision and UX design study](docs/statistics-and-ux.md) examines fixed-budget refinement, simulation uncertainty, cached backgrounds, and progressive results. It includes an executed precision/cost experiment and proposes a future workflow; it does not change the current search defaults.

Change the Monte Carlo trial count with `sea build --trials 99999` or `sea calibrate DATASET_ID --trials 99999`. More trials improve p-value resolution and increase computation time. The result's `execution` field records job duration and worker peak RSS in KiB.

Source code is licensed under MIT. The synthetic demo data is CC0. ChEMBL-derived data uses CC BY-SA 3.0, separately from the software license. Include the manifest's attribution and release information when redistributing a built dataset.
