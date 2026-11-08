# Using SEA Open from Python

`SEAClient` provides a synchronous Python interface to the same HTTP API used by
the web application. It submits jobs, waits for completion, and returns ordinary
Python dictionaries. Structure normalization, fingerprints, target annotations,
calibration, and scoring are performed by the configured server and its workers.
The client does not load a reference database or RDKit locally.

## Setup

From the repository checkout, use the existing locked Pixi environment:

```bash
pixi install --locked
./scripts/start.sh
pixi run sea datasets
```

Use your existing API key. If you need a separate key for a script, an
administrator can issue one with `pixi run sea key python-client`. Keys belong
to a user/agent identity; that identity determines which jobs it can access.

Set the following variables in the terminal that will run Python. Replace the
two placeholder values with your issued API key and an ID printed by
`pixi run sea datasets`:

```bash
export SEA_API_KEY="YOUR_ISSUED_API_KEY"
export SEA_DATASET_ID="YOUR_20_CHARACTER_DATASET_ID"
export SEA_URL="http://127.0.0.1:8000"
```

`SEA_DATASET_ID` is a content-addressed identifier such as the `id` field in a
dataset manifest. It is **not** a release name such as `ChEMBL_37`. Choose a
ChEMBL reference dataset for biological target hypotheses; the synthetic demo
contains chemical groups rather than biological targets. Retain the exact
dataset ID with an analysis so its data and calibration can be identified later.

Run a saved Python script with `pixi run python your_script.py`, or launch
`pixi run python` for an interactive session. No additional client dependencies
or package installation are needed inside this checkout's Pixi environment.
For an existing Python 3.12+ environment outside Pixi, `python -m pip install -e .`
from the checkout installs the project and its declared dependencies. A remote
server's workers must still be running; importing the client does not start them.

## Search for targets in a few lines

```python
import os
from sea import SEAClient

with SEAClient() as client:
    result = client.predict(
        "COc1ccc2[nH]cc(c2c1)C(CC)CN",
        dataset_id=os.environ["SEA_DATASET_ID"],
        top_k=10,
    )

for row in result["results"]:
    target = row["target"]
    print(
        target["id"], target["name"], target["family"], target["organism"],
        f"E={row['e_value']:.4g}", f"raw={row['raw_score']:.3f}",
    )
```

The context manager closes the HTTP connection pool. `SEAClient()` reads
`SEA_API_KEY` and `SEA_URL`; the URL defaults to `http://127.0.0.1:8000`. You can
also pass them explicitly:

```python
client = SEAClient(
    api_key=os.environ["SEA_API_KEY"],
    base_url="http://127.0.0.1:8000",
    timeout=60,
)
try:
    for dataset in client.list_datasets():
        print(dataset["id"], dataset["source"], dataset["release"], dataset["target_count"])
finally:
    client.close()
```

`base_url` is the server root, **without `/api/v1`**. A path prefix is supported,
for example `https://example.org/lab/sea/` when a proxy serves SEA there.
`timeout` configures HTTPX's network-operation timeout; it is independent of
the job-completion polling timeout.

A runnable example is provided in [examples/python_client.py](../examples/python_client.py):

```bash
pixi run python examples/python_client.py
```

It prints the complete result as JSON. Set `SEA_QUERY` to use another SMILES.
The existing [examples/agent.py](../examples/agent.py) demonstrates the raw HTTP
workflow using only `httpx`, without the Python client wrapper.

## Separate molecules and ensembles

By default, a flat list runs one search per molecule. Set `ensemble=True` to
treat the same list as one query ligand set:

```python
import os
from sea import SEAClient

smiles = ["COc1ccc2[nH]cc(c2c1)C(CC)CN", "COc1ccc2[nH]cc(CCN)c2c1"]
dataset_id = os.environ["SEA_DATASET_ID"]

with SEAClient() as client:
    separate = client.predict(smiles, dataset_id=dataset_id, top_k=10)
    ensemble = client.predict(smiles, dataset_id=dataset_id, ensemble=True, top_k=10)
```

In `separate`, `query_index` is 0 or 1. In `ensemble`, all results have
`query_index=0`; server-side deduplication determines the query's unique ligand
count. The client also accepts iterables of SMILES, including generators.

For several ensembles in one job, use the explicit submission interface:

```python
with SEAClient() as client:
    job_id = client.submit_prediction(
        [["CCO", "CCCO"], ["CCN", "CCCN"]],
        dataset_id=dataset_id,
        top_k=20,
    )
    result = client.wait(job_id)
```

Each inner list is a separate query set. `submit_prediction()` rejects a flat
list to prevent accidental changes in query grouping; use `predict()` for that
convenience form. The default server limit is 100 input molecules in total per
prediction request. `top_k` applies to each query and may be 1–1,000. Large
batches are not split automatically.

## Compare two ligand sets

```python
with SEAClient() as client:
    comparison = client.compare(
        ["CCO", "CCCO"], ["CCCCO", "CCO"],
        dataset_id=dataset_id,
    )

score = comparison["results"][0]
print(score["raw_score"], score["p_value"], score["e_value"])
```

Each side can be a single SMILES or an iterable. Set comparisons have
`query_index=None`, no target annotation, and one statistical comparison.
The default server limit is 1,000 molecules on each side. The specified reference
dataset supplies the fingerprint definition and calibrated background model.

## Import SDF, MOL, SMILES or CSV

File conversion uses the same authenticated endpoint as the web application:

```python
with SEAClient() as client:
    imported = client.import_file("molecules.sdf")
    smiles = [m["smiles"] for m in imported["molecules"]]
    result = client.predict(smiles, dataset_id=dataset_id, top_k=10)
```

Supported extensions are `.sdf`, `.sd`, `.mol`, `.smi`, `.smiles`, `.txt`, and
`.csv`, with case-insensitive detection. Files must be UTF-8 text. CSV requires
a `smiles` column and may contain a `name` column. To override an extension, use
`client.import_file("compounds.data", format="sdf")`.

For content already in memory, use
`client.import_structures(content, format="sdf")`; the format must be `sdf`,
`mol`, `smiles`, or `csv`. The result contains `count` and `molecules`, with
one-based record numbers, names, and isomeric SMILES. Conversion retains repeated
molecules and validates the whole file before returning. A malformed record
raises `httpx.HTTPStatusError` with structured record errors; no partial import
or prediction job is created. Server submission performs standardization and
deduplication afterward. Default upload and record limits are the same as the UI.

## Submit now, inspect or resume later

Use the submission methods when you need to retain a job ID, export CSV, or
continue an analysis from another Python process using the same API key:

```python
from pathlib import Path
import json

with SEAClient() as client:
    job_id = client.submit_prediction([["CCO"]], dataset_id=dataset_id, top_k=10)
    print("Retain this job ID:", job_id)
    print(client.get_job(job_id))

    result = client.wait(job_id, timeout=3600, poll_interval=1)
    Path("result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    Path("result.csv").write_text(client.get_csv(job_id), encoding="utf-8")
```

`submit_comparison(left, right, dataset_id=...)` provides the same job-ID workflow
for set comparisons. `get_result(job_id)` fetches a completed result without
polling; it raises HTTP 409 if the job is unfinished. `list_jobs()` returns the
current API-key owner's most recent 100 jobs. `cancel(job_id)` requests
cancellation; a running job may take time to reach a cancellation check.

`predict()` and `compare()` default to `wait_timeout=3600` and `poll_interval=1`.
`wait()` calls its equivalent parameter `timeout`. These are polling budgets;
an in-flight HTTP operation uses the separate constructor timeout and may finish
after the polling period. The helpers are synchronous and block the calling
thread while waiting. They do not automatically retry job-creation requests.

## Handle job errors and timeouts

```python
import httpx
from sea import SEAJobError, SEAJobTimeout

with SEAClient() as client:
    try:
        job_id = client.submit_prediction([["CCO"]], dataset_id=dataset_id)
        result = client.wait(job_id, timeout=60)
    except SEAJobTimeout as exc:
        print("Still available on the server:", exc.job_id)
        # Later: client.wait(exc.job_id, timeout=3600)
    except SEAJobError as exc:
        print(exc.job_id, exc.status, exc.error)
    except httpx.HTTPStatusError as exc:
        print(exc.response.status_code, exc.response.json())
    except httpx.RequestError as exc:
        print("The HTTP request could not complete:", exc)
```

A polling timeout **does not cancel the server job**. Resume with its ID instead
of resubmitting the same analysis. Blocking helpers also expose the ID through
`SEAJobTimeout` and `SEAJobError`. HTTP errors retain the server's status and
JSON body, including input record details. Authentication failures are HTTP 401;
invalid inputs are normally HTTP 422. For durable workflows, retain the ID before
waiting so it is also available if a network request fails.

## Results and target metadata

The returned bundle uses the same schema as the web/API result:

| Field | Meaning |
| --- | --- |
| `dataset` | Reference manifest, ID, model ID, fingerprint, RDKit version, and attribution |
| `threshold` | Calibrated Tanimoto cutoff used in raw-score aggregation |
| `normalization` | Input preprocessing and deduplication records |
| `results` | Ranked per-query scores and evidence |
| `execution` | Server computation duration and worker peak memory |

Each prediction row includes `query_index`, `target`, `raw_score`,
`contributing_pairs`, `left_size`, `right_size`, `similar_ligands`, `z_score`,
`p_value`, `e_value`, `log_p_value`, `log_e_value`, `method`, `comparisons`,
and `p_resolution`. `target` contains `id`, `name`, `family`, `families`,
`protein_classes`, `organism`, and `tax_id`. Each saved ligand pair contains
query/reference SMILES and Tanimoto. Saved evidence can contain fewer pairs than
`contributing_pairs`; no saved pair means no pair reached the threshold, rather
than an independently measured maximum Tanimoto of zero.

Use `log_p_value` and `log_e_value` for extremely small values that underflow to
zero. Logs are natural logarithms. E-values are expected random hits across all
searched targets, not binding probabilities; they can exceed 1. Monte Carlo
resolution can create tied E-values. See [method.md](method.md) and the
[original-paper comparison](paper-comparison.md) for interpretation.

To search target annotations without running a prediction:

```python
with SEAClient() as client:
    matches = client.list_targets(dataset_id, query="GPCR", offset=0, limit=50)
```

The search matches target name, ID, family and species. Results include the
reference ligand count. Pagination uses `offset` and `limit` (1–1,000).

If pandas is installed, a prediction bundle can be flattened for notebook work:

```python
import pandas as pd

table = pd.json_normalize(result["results"])
columns = ["query_index", "target.name", "target.family", "target.organism", "e_value", "raw_score"]
print(table.reindex(columns=columns))
```

## Local calculation without the HTTP service

The underlying Python calculation functions are also available in the Pixi
environment. For example, a local two-set calculation can reuse a built dataset:

```python
import os
from sea.chemistry import normalize, raw_score
from sea.data import load_dataset
from sea.statistics import significance

data, molecules = load_dataset(os.environ["SEA_DATASET_ID"])
left, left_records = normalize(["CCO", "CCCO"])
right, right_records = normalize(["CCCCO", "CCO"])
raw, pair_count, pairs = raw_score(left, right, data["model"]["threshold"])
statistics = significance(
    raw, len(left), len(right), list(molecules.values()), data["model"], comparisons=1,
)
print(raw, pair_count, statistics["e_value"])
```

This loads the reference dataset into the local process and uses local CPU/memory.
Set `SEA_HOME` before starting Python if the data are in another directory.
The default is the checkout's `var/` under Pixi. No API key, server, worker,
or queued job is involved. The locked RDKit version must match the dataset.
Use `SEAClient` when you want server-side target searches, persistent jobs,
ownership, shared workers, and the same API result bundle as the web application.
