"""Exercise actual detached processes in an isolated SEA_HOME."""
import importlib.util
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
from types import SimpleNamespace

import httpx
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/services.py"
spec = importlib.util.spec_from_file_location("sea_services", SCRIPT)
services_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(services_module)


def command(home, action, *arguments):
    return subprocess.run(
        [sys.executable, str(SCRIPT), action, "--home", str(home), *arguments],
        cwd=home.parent, capture_output=True, text=True, timeout=40,
    )


def free_port():
    with socket.create_server(("127.0.0.1", 0)) as server:
        return server.getsockname()[1]


def test_start_duplicate_start_stop_and_duplicate_stop(tmp_path):
    home = tmp_path / "lab data"
    port = free_port()
    controller = services_module.Services(home)
    records, workers = {}, []
    try:
        result = command(home, "start", "--port", str(port))
        assert result.returncode == 0, result.stderr + result.stdout
        records = {name: json.loads(controller.pidfile(name).read_text()) for name in ("api", "worker")}
        pid = records["worker"]["pid"]
        workers = [{"pid": int(p), "start_time": services_module.identity(int(p))}
                   for p in Path(f"/proc/{pid}/task/{pid}/children").read_text().split()]
        assert len(workers) == 2
        with httpx.Client(trust_env=False) as client:
            assert client.get(f"http://127.0.0.1:{port}/health").json() == {"status": "ok"}
        # Earlier local launches recorded just the PID, without start-time metadata.
        for name, info in records.items():
            controller.pidfile(name).write_text(str(info["pid"]) + "\n")
        result = command(home, "start", "--port", str(port))
        assert result.returncode == 0 and "already running" in result.stdout
        assert records == {name: json.loads(controller.pidfile(name).read_text()) for name in records}
        result = command(home, "stop")
        assert result.returncode == 0, result.stderr
        assert all(not services_module.alive(info) for info in [*records.values(), *workers])
        assert not list(controller.logs.glob("*.pid"))
        result = command(home, "stop")
        assert result.returncode == 0 and "already stopped" in result.stdout
    finally:
        command(home, "stop")


def test_occupied_port_does_not_start_workers(tmp_path):
    home = tmp_path / "service"
    with socket.create_server(("127.0.0.1", 0)) as server:
        result = command(home, "start", "--port", str(server.getsockname()[1]))
    assert result.returncode == 1
    assert not list((home / "logs").glob("*.pid"))


def test_stale_pid_does_not_stop_unrelated_process(tmp_path):
    home = tmp_path / "service"
    controller = services_module.Services(home)
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"], cwd=ROOT,
                               env={**os.environ, "SEA_HOME": str(home)})
    try:
        controller.remember("api", {"pid": sleeper.pid, "start_time": services_module.identity(sleeper.pid)})
        result = command(home, "stop")
        assert result.returncode == 0
        assert "Ignoring stale api" in result.stdout
        assert sleeper.poll() is None
    finally:
        sleeper.terminate()
        sleeper.wait(timeout=5)


def test_failed_api_start_rolls_back_new_worker(tmp_path, monkeypatch):
    controller = services_module.Services(tmp_path / "service")
    launch = controller.launch
    launched = []
    def fail_api(name, arguments):
        child, info = launch(name, ["-c", "raise SystemExit(1)"] if name == "api" else arguments)
        launched.append(info)
        return child, info
    monkeypatch.setattr(controller, "launch", fail_api)
    with pytest.raises(RuntimeError, match="api exited"):
        controller.start(SimpleNamespace(host="127.0.0.1", port=free_port(), workers=2, timeout=10))
    assert all(not services_module.alive(info) for info in launched)
    assert not list(controller.logs.glob("*.pid"))
