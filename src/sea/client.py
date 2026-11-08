"""Synchronous Python client for the SEA Open HTTP API.

Only httpx and the Python standard library are imported here. Calculations run
on the configured SEA server and its workers, using the requested dataset.
"""
from collections.abc import Iterable
import math
import os
from pathlib import Path
import time
from typing import Any, Literal
from urllib.parse import quote

import httpx


StructureFormat = Literal["sdf", "mol", "smiles", "csv"]
_FORMATS = {".sdf": "sdf", ".sd": "sdf", ".mol": "mol", ".smi": "smiles",
            ".smiles": "smiles", ".txt": "smiles", ".csv": "csv"}


class SEAJobError(RuntimeError):
    """A server job failed or was cancelled; inspect job_id, status and error."""

    def __init__(self, job: dict[str, Any]):
        self.job_id = job["id"]
        self.status = job["status"]
        self.error = job.get("error")
        super().__init__(f"SEA job {self.job_id} {self.status}" + (f": {self.error}" if self.error else "."))


class SEAJobTimeout(TimeoutError):
    """Polling expired. The server job continues and can be resumed by job_id."""

    def __init__(self, job_id: str):
        self.job_id = job_id
        super().__init__(f"Timed out waiting for SEA job {job_id}. The job was not cancelled; call wait() again to resume.")


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number.")
    return value


def _smiles(values: str | Iterable[str]) -> list[str]:
    result = [values] if isinstance(values, str) else list(values)
    if not result or any(not isinstance(s, str) or not s.strip() for s in result):
        raise ValueError("Provide at least one nonempty SMILES string.")
    return result


def _job_path(job_id: str) -> str:
    if not isinstance(job_id, str) or not job_id:
        raise ValueError("Provide a job ID.")
    return "jobs/" + quote(job_id, safe="")


class SEAClient:
    """Use the existing API with an API key and explicit reference dataset IDs.

    api_key defaults to SEA_API_KEY. base_url defaults to SEA_URL, then
    http://127.0.0.1:8000. Supply the server root, without /api/v1.
    timeout is the httpx timeout in seconds, separate from job wait timeouts.
    Use a context manager, or call close() when finished.
    """

    def __init__(self, api_key: str | None = None, *, base_url: str | None = None,
                 timeout: float = 30.0, transport: httpx.BaseTransport | None = None):
        key = os.getenv("SEA_API_KEY") if api_key is None else api_key
        if not key or not key.strip():
            raise ValueError("Provide api_key or set the SEA_API_KEY environment variable.")
        root = base_url if base_url is not None else os.getenv("SEA_URL", "http://127.0.0.1:8000")
        self._http = httpx.Client(base_url=root.rstrip("/") + "/", headers={"X-API-Key": key},
                                  timeout=_positive(timeout, "timeout"), transport=transport)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def close(self) -> None:
        self._http.close()

    def _request(self, method: str, path: str, **kwargs) -> httpx.Response:
        response = self._http.request(method, "api/v1/" + path, **kwargs)
        response.raise_for_status()
        return response

    def list_datasets(self) -> list[dict[str, Any]]:
        """List available reference manifests; select and retain an explicit ID."""
        return self._request("GET", "datasets").json()

    def list_targets(self, dataset_id: str, *, query: str = "", offset: int = 0,
                     limit: int = 50) -> list[dict[str, Any]]:
        return self._request("GET", "targets", params={
            "dataset_id": dataset_id, "q": query, "offset": offset, "limit": limit}).json()

    def list_jobs(self) -> list[dict[str, Any]]:
        """Return the authenticated owner's most recent 100 jobs."""
        return self._request("GET", "jobs").json()

    def get_job(self, job_id: str) -> dict[str, Any]:
        return self._request("GET", _job_path(job_id)).json()

    def get_result(self, job_id: str) -> dict[str, Any]:
        """Fetch a completed job. The server returns HTTP 409 if unfinished."""
        return self._request("GET", _job_path(job_id) + "/result").json()

    def get_csv(self, job_id: str) -> str:
        return self._request("GET", _job_path(job_id) + "/result", params={"format": "csv"}).text

    def cancel(self, job_id: str) -> dict[str, Any]:
        """Request cancellation; a running job may need time to stop."""
        return self._request("POST", _job_path(job_id) + "/cancel").json()

    def submit_prediction(self, queries: Iterable[Iterable[str]], *, dataset_id: str,
                          top_k: int = 20) -> str:
        """Submit explicit ensembles: [['CCO'], ['CCN']] runs two searches."""
        if isinstance(queries, str):
            raise ValueError("queries must contain ligand sets, for example [['CCO']].")
        sets = []
        for values in queries:
            if isinstance(values, str):
                raise ValueError("Each query must be a ligand set. Use predict() for a flat SMILES list.")
            sets.append(_smiles(values))
        if not sets:
            raise ValueError("Provide at least one query ligand set.")
        response = self._request("POST", "predictions", json={
            "dataset_id": dataset_id, "queries": sets, "top_k": top_k})
        return response.json()["job_id"]

    def submit_comparison(self, left: str | Iterable[str], right: str | Iterable[str], *,
                          dataset_id: str) -> str:
        response = self._request("POST", "comparisons", json={
            "dataset_id": dataset_id, "left": _smiles(left), "right": _smiles(right)})
        return response.json()["job_id"]

    def wait(self, job_id: str, *, timeout: float = 3600.0,
             poll_interval: float = 1.0) -> dict[str, Any]:
        """Poll until complete and return the result bundle; never cancel on timeout.

        timeout bounds the polling period. In-flight HTTP operations use the
        constructor's separate httpx timeout and may finish after that period.
        """
        deadline = time.monotonic() + _positive(timeout, "timeout")
        interval = _positive(poll_interval, "poll_interval")
        while time.monotonic() < deadline:
            job = self.get_job(job_id)
            if job["status"] == "completed":
                return self.get_result(job_id)
            if job["status"] in ("failed", "cancelled"):
                raise SEAJobError(job)
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            time.sleep(min(interval, remaining))
        raise SEAJobTimeout(job_id)

    def predict(self, smiles: str | Iterable[str], *, dataset_id: str, ensemble: bool = False,
                top_k: int = 20, wait_timeout: float = 3600.0,
                poll_interval: float = 1.0) -> dict[str, Any]:
        """Search one SMILES or a batch, wait for completion, and return results.

        By default, each molecule is a separate query. Set ensemble=True to
        combine all molecules into one query set. Use submit_prediction() to
        obtain the job ID before waiting, for durable workflows or CSV export.
        """
        _positive(wait_timeout, "wait_timeout")
        _positive(poll_interval, "poll_interval")
        structures = _smiles(smiles)
        job_id = self.submit_prediction([structures] if ensemble else [[s] for s in structures],
                                        dataset_id=dataset_id, top_k=top_k)
        return self.wait(job_id, timeout=wait_timeout, poll_interval=poll_interval)

    def compare(self, left: str | Iterable[str], right: str | Iterable[str], *, dataset_id: str,
                wait_timeout: float = 3600.0, poll_interval: float = 1.0) -> dict[str, Any]:
        """Compare two ligand sets and wait for the result bundle."""
        _positive(wait_timeout, "wait_timeout")
        _positive(poll_interval, "poll_interval")
        job_id = self.submit_comparison(left, right, dataset_id=dataset_id)
        return self.wait(job_id, timeout=wait_timeout, poll_interval=poll_interval)

    def import_structures(self, content: str, *, format: StructureFormat) -> dict[str, Any]:
        """Convert and validate all records on the server before returning any."""
        return self._request("POST", "structures/import", json={"format": format, "content": content}).json()

    def import_file(self, path: str | Path, *, format: StructureFormat | None = None) -> dict[str, Any]:
        """Read a UTF-8 structure file; infer supported formats from its extension."""
        path = Path(path)
        selected = format if format is not None else _FORMATS.get(path.suffix.lower())
        if selected not in ("sdf", "mol", "smiles", "csv"):
            raise ValueError("Use SDF, MOL, SMI, SMILES, TXT or CSV, or specify format explicitly.")
        return self.import_structures(path.read_text(encoding="utf-8-sig"), format=selected)
