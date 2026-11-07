# SEA Open

An open-source implementation of the Similarity Ensemble Approach (SEA), which relates protein targets through the chemical similarity of their ligand sets. Built with Pixi, RDKit, and FastAPI, it provides a Japanese web interface, API key authentication, persistent jobs, and separate worker processes.

Original paper: Keiser et al., *Relating protein pharmacology by ligand chemistry*, Nature Biotechnology 25, 197–206 (2007), https://doi.org/10.1038/nbt1284

## Getting started

Requirements: Linux x86_64 and Pixi.

```bash
pixi install --locked
pixi run demo
pixi run sea key alice
```

Save the API key printed by the last command. Only its SHA256 hash is stored in the database. Start the worker and API server in separate terminals:

```bash
pixi run worker
pixi run serve
```

Open http://127.0.0.1:8000 and sign in with your API key. Interactive API documentation is available at `/docs`, and the machine-readable specification is at `/openapi.json`. The demo contains synthetic alcohol, carboxylic acid, and aromatic compound sets; it cannot predict biological targets.

Reference datasets and the service database are stored in `SEA_HOME` (default: `var/`). Run `pixi run sea datasets` to list dataset IDs. If your sandbox does not allow writes to `/tmp`, create a temporary directory inside the workspace and use `TMPDIR="$PWD/.cache/tmp" pixi …`.

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

## Agent access

Set `SEA_API_KEY` to your API key and `SEA_DATASET_ID` to a dataset ID.

```bash
curl -H "X-API-Key: $SEA_API_KEY" http://127.0.0.1:8000/api/v1/datasets
```

[examples/agent.py](examples/agent.py) submits a request, polls its status, and retrieves the result:

```bash
pixi run python examples/agent.py
```

For `POST /api/v1/predictions`, `queries` is an array of ligand sets. `[["CCO"],["CCCO"]]` searches each molecule separately; `[["CCO","CCCO"]]` searches one ensemble. `POST /api/v1/comparisons` accepts `left` and `right` ligand sets. Both endpoints require `dataset_id` and return HTTP 202 with a `job_id`.

- `GET /api/v1/jobs/{id}`: job status and progress.
- `GET /api/v1/jobs/{id}/result`: JSON results; add `?format=csv` for CSV.
- `POST /api/v1/jobs/{id}/cancel`: request cancellation.
- `GET /api/v1/jobs`: your most recent 100 jobs.
- `GET /api/v1/targets?dataset_id=...&q=...`: search reference targets.

Default limits are 100 molecules per prediction request, 1,000 molecules per set comparison side, and a 10 MB request body. Configure these with `SEA_MAX_QUERY_MOLECULES`, `SEA_MAX_SET_MOLECULES`, and `SEA_MAX_UPLOAD_BYTES`. Invalid inputs return HTTP 422, authentication failures return 401, another user's jobs return 404, and requests for unfinished results return 409.

## Sharing within a lab

Run the API and `pixi run worker` on the same host. The worker supervisor starts two processes by default; use `pixi run sea worker --count 4` to change this. Restarting the supervisor requeues interrupted jobs. Computation is deterministic, and retrying a job has no external write effects. Back up the service database and datasets. Jobs are not deleted automatically.

Issue a separate API key for each user or agent. Browser sign-in uses an eight-hour HttpOnly, SameSite=Strict session cookie; API keys are not saved in browser storage. For shared deployments over TLS, set `SEA_SECURE_COOKIE=1` and configure a TLS reverse proxy using [docs/Caddyfile](docs/Caddyfile) as a starting point. Key issuance and database construction are administrator CLI operations. Bind the API to 127.0.0.1 and expose only the proxy to the lab network.

## Validation and limitations

```bash
pixi run test
```

See [docs/method.md](docs/method.md) for calculation details and [docs/validation.md](docs/validation.md) for the recorded checks. RDKit and ChEMBL replace the original Daylight fingerprints and MDDR dataset, so numerical agreement with the original paper is not claimed. E-values are not binding probabilities. With the default 999 Monte Carlo trials, the smallest possible p-value is 0.001; this resolution may be insufficient when searching many targets.

Zero-score comparisons return the exact p-value 1 without Monte Carlo sampling (`method=exact_zero`). Workers reuse deterministic background samples across jobs for the same dataset and set sizes. The web interface displays queue state and progress and prevents repeated submissions while a job is active.

Change the Monte Carlo trial count with `sea build --trials 99999` or `sea calibrate DATASET_ID --trials 99999`. More trials improve p-value resolution and increase computation time. The result's `execution` field records job duration and worker peak RSS in KiB.

Source code is licensed under MIT. The synthetic demo data is CC0. ChEMBL-derived data uses CC BY-SA 3.0, separately from the software license. Include the manifest's attribution and release information when redistributing a built dataset.
