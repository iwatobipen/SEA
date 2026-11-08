"""Python client request semantics, lifecycle, errors, and real API integration."""
import json

import httpx
import pytest
from fastapi.testclient import TestClient

from sea import SEAClient, SEAJobError, SEAJobTimeout
import sea.client as client_module
from sea.api import app
from sea.cli import demo
from sea.jobs import claim, create_key, process


@pytest.mark.parametrize('wait_for', ['evidence', 'screening', 'final'])
def test_precision_wait_follows_declared_phase(wait_for):
    seen = []
    screen = {'analysis': {'phase': 'screening', 'complete': False, 'refinement_job_id': 'child'}, 'results': []}
    final = {'analysis': {'phase': 'final', 'complete': True}, 'results': []}
    def handler(request):
        seen.append(request.url.path)
        if request.url.path.endswith('/preview'):
            return httpx.Response(200, json=screen)
        if request.url.path.endswith('/source/result'):
            return httpx.Response(200, json=screen)
        if request.url.path.endswith('/child/result'):
            return httpx.Response(200, json=final)
        return httpx.Response(200, json={'id': request.url.path.split('/')[-1], 'status': 'completed', 'refinement_job_id': 'child'})
    with SEAClient('key', transport=httpx.MockTransport(handler)) as client:
        result = client.wait('source', wait_for=wait_for)
    assert result == (final if wait_for == 'final' else screen)
    assert any('/child/' in path for path in seen) == (wait_for == 'final')


def test_precision_timeout_keeps_one_deadline_and_reports_child(monkeypatch):
    clock = [0.]
    monkeypatch.setattr(client_module.time, 'monotonic', lambda: clock[0])
    monkeypatch.setattr(client_module.time, 'sleep', lambda value: clock.__setitem__(0, clock[0] + value))
    def handler(request):
        clock[0] += .02
        if request.url.path.endswith('/result'):
            return httpx.Response(200, json={'analysis': {'complete': False, 'refinement_job_id': 'child'}})
        return httpx.Response(200, json={'id': 'child', 'status': 'completed' if request.url.path.endswith('/source') else 'running', 'refinement_job_id': 'child'})
    with SEAClient('key', transport=httpx.MockTransport(handler)) as client:
        with pytest.raises(SEAJobTimeout) as caught:
            client.wait('source', timeout=.1, poll_interval=.01)
    assert caught.value.job_id == 'child' and clock[0] <= .12


def test_prediction_shapes_and_proxy_base_path(monkeypatch):
    monkeypatch.setenv("SEA_API_KEY", "test-key")
    monkeypatch.setenv("SEA_URL", "https://example.test/lab/sea/")
    submitted = []

    def handler(request):
        assert request.headers["X-API-Key"] == "test-key"
        if request.method == "POST":
            assert request.url.path == "/lab/sea/api/v1/predictions"
            submitted.append(json.loads(request.content))
            return httpx.Response(202, json={"job_id": "job-one"})
        if request.url.path.endswith("/result"):
            return httpx.Response(200, json={"results": [{"raw_score": 1}]})
        return httpx.Response(200, json={"id": "job-one", "status": "completed"})

    with SEAClient(transport=httpx.MockTransport(handler)) as client:
        assert client.predict("CCO", dataset_id="pinned-dataset")["results"]
        client.predict(["CCO", "CCN"], dataset_id="pinned-dataset", top_k=50)
        client.predict((s for s in ["CCO", "CCN"]), dataset_id="pinned-dataset", ensemble=True)
        client.submit_prediction([["CCO", "CCCO"], ["CCN"]], dataset_id="pinned-dataset")
    assert [body["queries"] for body in submitted] == [
        [["CCO"]], [["CCO"], ["CCN"]], [["CCO", "CCN"]], [["CCO", "CCCO"], ["CCN"]]]
    assert submitted[1]["top_k"] == 50
    assert all(body["dataset_id"] == "pinned-dataset" for body in submitted)
    assert client._http.is_closed


def test_polling_timeout_can_resume_without_resubmitting_or_cancelling(monkeypatch):
    clock = [0.]
    monkeypatch.setattr(client_module.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(client_module.time, "sleep", lambda seconds: clock.__setitem__(0, clock[0] + seconds))
    state = {"completed": False, "requests": []}

    def handler(request):
        state["requests"].append((request.method, request.url.path))
        if request.url.path.endswith("/result"):
            return httpx.Response(200, json={"results": []})
        return httpx.Response(200, json={"id": "existing-job", "status": "completed" if state["completed"] else "running"})

    with SEAClient("key", transport=httpx.MockTransport(handler)) as client:
        with pytest.raises(SEAJobTimeout) as caught:
            client.wait("existing-job", timeout=.1, poll_interval=.06)
        assert caught.value.job_id == "existing-job"
        assert all(method == "GET" for method, _ in state["requests"])
        state["completed"] = True
        assert client.wait(caught.value.job_id) == {"results": []}
        assert state["requests"][-1][1].endswith("/existing-job/result")


@pytest.mark.parametrize("status", ["failed", "cancelled"])
@pytest.mark.parametrize("wait_for", ["evidence", "final"])
def test_terminal_jobs_raise_with_the_job_id(status, wait_for):
    def handler(request):
        return httpx.Response(200, json={"id": "terminal-job", "status": status, "error": "Example error"})

    with SEAClient("key", transport=httpx.MockTransport(handler)) as client:
        with pytest.raises(SEAJobError) as caught:
            client.wait("terminal-job", wait_for=wait_for)
        assert caught.value.job_id == "terminal-job"
        assert caught.value.status == status
        assert caught.value.error == "Example error"


def test_request_errors_preserve_http_status_and_structured_detail():
    def handler(request):
        return httpx.Response(422, json={"detail": {"code": "invalid_structures", "records": [{"record": 2}]}})

    with SEAClient("key", transport=httpx.MockTransport(handler)) as client:
        with pytest.raises(httpx.HTTPStatusError) as caught:
            client.import_structures("invalid", format="sdf")
        assert caught.value.response.status_code == 422
        assert caught.value.response.json()["detail"]["records"][0]["record"] == 2


def test_invalid_local_arguments_do_not_submit_jobs(monkeypatch):
    monkeypatch.delenv("SEA_API_KEY", raising=False)
    with pytest.raises(ValueError, match="SEA_API_KEY"):
        SEAClient()
    requests = []
    transport = httpx.MockTransport(lambda request: requests.append(request))
    with SEAClient("key", transport=transport) as client:
        with pytest.raises(ValueError):
            client.predict('CCO', dataset_id='dataset', statistics_profile='adaptive')
        with pytest.raises(ValueError):
            client.predict('CCO', dataset_id='dataset', wait_for='unknown')
        for values in ([], "", ["CCO", ""], ["CCO", 123]):
            with pytest.raises(ValueError):
                client.predict(values, dataset_id="dataset")
        with pytest.raises(ValueError, match="ligand set"):
            client.submit_prediction(["CCO"], dataset_id="dataset")
        for invalid in (0, -1, float("inf"), float("nan")):
            with pytest.raises(ValueError):
                client.predict("CCO", dataset_id="dataset", wait_timeout=invalid)
            with pytest.raises(ValueError):
                client.compare(["CCO"], ["CCN"], dataset_id="dataset", poll_interval=invalid)
    assert not requests


def test_client_with_real_api_scoring_import_history_and_exports(tmp_path, monkeypatch):
    monkeypatch.setenv("SEA_HOME", str(tmp_path / "service"))
    dataset = demo()
    key = create_key("python-client-test")
    csv_file = tmp_path / "structures.CSV"
    csv_file.write_text('\ufeffsmiles,name\nCCO,ethanol\nCCCO,propanol\n', encoding="utf-8")

    with TestClient(app) as api_client:
        def handler(request):
            # Drive the actual job engine when the SDK polls, without a network
            # listener or a long-lived worker for this integration test.
            if request.method == "GET" and request.url.path.startswith("/api/v1/jobs/") and not request.url.path.endswith("/result"):
                job = claim("python-client-test")
                if job:
                    process(job)
            response = api_client.request(request.method, str(request.url), headers=request.headers, content=request.content)
            return httpx.Response(response.status_code, headers=response.headers, content=response.content)

        with SEAClient(key, transport=httpx.MockTransport(handler)) as client:
            assert client.list_datasets()[0]["id"] == dataset
            assert len(client.list_targets(dataset, query="Alcohol")) == 1
            converted = client.import_file(csv_file)
            smiles = [m["smiles"] for m in converted["molecules"]]
            assert smiles == ["CCO", "CCCO"]
            result = client.predict(smiles, dataset_id=dataset, top_k=1)
            assert len(result["results"]) == 2
            assert [r["query_index"] for r in result["results"]] == [0, 1]
            comparison = client.compare(smiles, ["CCO"], dataset_id=dataset)
            assert comparison["results"][0]["comparisons"] == 1
            assert comparison["results"][0]["query_index"] is None
            job_id = client.submit_prediction([smiles], dataset_id=dataset)
            result = client.wait(job_id)
            assert result["results"][0]["left_size"] == 2
            assert client.get_result(job_id) == result
            assert "target_family" in client.get_csv(job_id).splitlines()[0]
            assert {j["id"] for j in client.list_jobs()} >= {job_id}
            cancelled = client.submit_comparison("CCO", "CCN", dataset_id=dataset)
            assert client.cancel(cancelled)["status"] == "cancelled"
            with pytest.raises(SEAJobError) as caught:
                client.wait(cancelled)
            assert caught.value.status == "cancelled"
            with pytest.raises(httpx.HTTPStatusError) as caught:
                client.import_structures("CCO valid\ninvalid bad\n", format="smiles")
            assert caught.value.response.json()["detail"]["records"][0]["record"] == 2
