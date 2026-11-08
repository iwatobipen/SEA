"""Start and stop this checkout's API and worker on Linux."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time

import httpx

ROOT = Path(__file__).resolve().parents[1]


def identity(pid):
    """A PID plus Linux start time distinguishes reused process IDs."""
    try:
        fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
        return fields[19] if fields[0] not in ("Z", "X") else None
    except (FileNotFoundError, ProcessLookupError):
        return None


def alive(info):
    return info["start_time"] is not None and identity(info["pid"]) == info["start_time"]


class Services:
    def __init__(self, home):
        self.home = Path(home).expanduser().resolve()
        self.logs = self.home / "logs"
        self.logs.mkdir(parents=True, exist_ok=True)

    def pidfile(self, name):
        return self.logs / f"{name}.pid"

    def remember(self, name, info):
        target = self.pidfile(name)
        temporary = target.with_suffix(f".pid.{os.getpid()}.tmp")
        temporary.write_text(json.dumps(info) + "\n")
        temporary.replace(target)

    def find(self, name):
        path = self.pidfile(name)
        if not path.exists():
            return None
        try:
            recorded = json.loads(path.read_text())
            # Accept the integer PID files used by earlier local launches.
            pid = recorded if type(recorded) is int else recorded["pid"]
            if type(pid) is not int or pid <= 1:
                raise ValueError("Invalid PID")
            start_time = identity(pid)
            if start_time is not None:
                proc = Path(f"/proc/{pid}")
                argv = [part.decode() for part in (proc / "cmdline").read_bytes().split(b"\0") if part]
                prefix = ["-m", "uvicorn", "sea.api:app"] if name == "api" else ["-m", "sea", "worker"]
                env = dict(part.decode().split("=", 1) for part in (proc / "environ").read_bytes().split(b"\0") if b"=" in part)
                process_home = Path(env.get("SEA_HOME", "var"))
                if not process_home.is_absolute():
                    process_home = ROOT / process_home
                matches = (
                    argv[1:4] == prefix
                    and Path(argv[0]).resolve() == Path(sys.executable).resolve()
                    and (proc / "cwd").resolve() == ROOT
                    and process_home.resolve() == self.home
                    and (type(recorded) is int or recorded.get("start_time") == start_time)
                )
                if matches:
                    return {"pid": pid, "start_time": start_time, "argv": argv}
        except (FileNotFoundError, ProcessLookupError, ValueError, TypeError, KeyError, IndexError):
            pass
        print(f"Ignoring stale {name} PID file; no process was signalled.")
        path.unlink(missing_ok=True)
        return None

    def launch(self, name, arguments):
        env = {**os.environ, "SEA_HOME": str(self.home)}
        env["PYTHONPATH"] = str(ROOT / "src") + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        with (self.logs / f"{name}.log").open("ab") as log:
            process = subprocess.Popen(
                [sys.executable, *arguments], cwd=ROOT, env=env,
                stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                start_new_session=True,
            )
        info = {"pid": process.pid, "start_time": identity(process.pid)}
        self.remember(name, info)
        return process, info

    def stop(self, processes, timeout):
        for name, info in processes:
            if alive(info):
                try:
                    os.kill(info["pid"], signal.SIGTERM)
                except ProcessLookupError:
                    pass
        deadline = time.monotonic() + timeout
        remaining = list(processes)
        while remaining:
            for name, info in remaining[:]:
                if not alive(info):
                    self.pidfile(name).unlink(missing_ok=True)
                    remaining.remove((name, info))
                    print(f"Stopped {name}.")
            if not remaining:
                return
            if time.monotonic() >= deadline:
                names = ", ".join(name for name, _ in remaining)
                raise RuntimeError(f"Still stopping after {timeout:g}s: {names}. PID files retained; inspect {self.logs} and retry with a longer --timeout.")
            time.sleep(.1)

    def start(self, args):
        existing = {name: self.find(name) for name in ("api", "worker")}
        expected = {"api": {"--host": args.host, "--port": str(args.port)}, "worker": {"--count": str(args.workers)}}
        defaults = {"--host": "127.0.0.1", "--port": "8000", "--count": "2"}
        for name, info in existing.items():
            if not info:
                continue
            for option, value in expected[name].items():
                argv = info["argv"]
                actual = argv[argv.index(option) + 1] if option in argv else defaults[option]
                if actual != value:
                    raise RuntimeError(f"{name} already runs with {option} {actual}. Stop it before changing this setting.")
        if not existing["api"]:
            # Refuse an occupied port before starting any new worker.
            with socket.create_server((args.host, args.port), family=socket.AF_INET6 if ":" in args.host else socket.AF_INET):
                pass
        launched = []
        children = []
        try:
            for name, command in (
                ("worker", ["-m", "sea", "worker", "--count", str(args.workers)]),
                ("api", ["-m", "uvicorn", "sea.api:app", "--host", args.host, "--port", str(args.port)]),
            ):
                if not existing[name]:
                    child, info = self.launch(name, command)
                    children.append((name, child))
                    launched.append((name, info))
                    existing[name] = info
            host = {"0.0.0.0": "127.0.0.1", "::": "::1"}.get(args.host, args.host)
            host = f"[{host}]" if ":" in host else host
            url = f"http://{host}:{args.port}"
            deadline = time.monotonic() + args.timeout
            with httpx.Client(timeout=1, trust_env=False) as client:
                while time.monotonic() < deadline:
                    for name, child in children:
                        if child.poll() is not None:
                            raise RuntimeError(f"{name} exited; inspect {self.logs / (name + '.log')}.")
                    try:
                        response = client.get(url + "/health")
                        healthy = response.status_code == 200 and response.json() == {"status": "ok"}
                    except (httpx.HTTPError, ValueError):
                        healthy = False
                    # Python 3.12's Linux worker supervisor forks these children
                    # only after acquiring its lock and initializing the queue.
                    pid = existing["worker"]["pid"]
                    try:
                        worker_pids = Path(f"/proc/{pid}/task/{pid}/children").read_text().split()
                        worker_ready = len(worker_pids) == args.workers and all(identity(int(child)) for child in worker_pids)
                    except FileNotFoundError:
                        worker_ready = False
                    if healthy and worker_ready and all(alive(info) for info in existing.values()):
                        break
                    time.sleep(.1)
                else:
                    raise RuntimeError(f"Services were not ready within {args.timeout:g}s; inspect {self.logs}.")
            for name, info in existing.items():
                self.remember(name, {key: info[key] for key in ("pid", "start_time")})
            print("SEA services are ready." if launched else "SEA services are already running.")
            print(f"Web app: {url}\nLogs: {self.logs}")
        except BaseException:
            self.stop(launched, 15)
            raise
        finally:
            for _, child in children:
                if child.poll() is not None:
                    child.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("start", "stop"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--home", default=os.getenv("SEA_HOME", str(ROOT / "var")), help="Data directory; use the same value for start and stop.")
        sub.add_argument("--timeout", type=float, default=30, help="Seconds to wait for readiness or shutdown (default: 30).")
        if command == "start":
            sub.add_argument("--host", default="127.0.0.1")
            sub.add_argument("--port", type=int, default=8000)
            sub.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.timeout <= 0 or (args.command == "start" and (not 1 <= args.port <= 65535 or not 1 <= args.workers <= 32)):
        parser.error("Use a positive timeout, a port from 1 to 65535, and 1 to 32 workers.")
    services = Services(args.home)
    try:
        with (services.logs / "services.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if args.command == "start":
                services.start(args)
            else:
                processes = [(name, info) for name in ("api", "worker") if (info := services.find(name))]
                if processes:
                    services.stop(processes, args.timeout)
                else:
                    print("SEA services are already stopped.")
    except (OSError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
