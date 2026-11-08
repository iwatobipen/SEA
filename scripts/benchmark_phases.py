"""Measure isolated ChEMBL evidence/screening latency during precision work.

The full cold refinement is deliberately paused, not claimed complete. Lab
datasets are hardlinked read-only by convention; keys/jobs/cache are temporary.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import httpx

from sea import SEAClient
from sea.data import home
from sea.jobs import create_key


def until(get, accept, timeout=120):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = get()
        if accept(value):
            return value
        time.sleep(.05)
    raise TimeoutError('Timed out measuring a phase transition.')


def measure(dataset_id, smiles, port):
    original_home = home().resolve()
    Path('.cache').mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='sea-phases-', dir='.cache') as directory:
        root = Path(directory).resolve()
        for folder in ('datasets', 'manifests', 'annotations'):
            source = original_home / folder / f'{dataset_id}.json'
            if source.exists():
                (root / folder).mkdir(exist_ok=True)
                os.link(source, root / folder / source.name)
        old_home = os.environ.get('SEA_HOME')
        os.environ['SEA_HOME'] = str(root)
        key = create_key('isolated-phase-benchmark')
        env = {**os.environ, 'SEA_PRELOAD_DATASETS': dataset_id}
        processes = []
        started = time.perf_counter()
        try:
            with (root / 'services.log').open('w') as log:
                processes.append(subprocess.Popen([sys.executable, '-m', 'uvicorn', 'sea.api:app', '--host', '127.0.0.1', '--port', str(port)], env=env, stdout=log, stderr=log))
                processes.append(subprocess.Popen([sys.executable, '-m', 'sea', 'worker', '--count', '2'], env=env, stdout=log, stderr=log))
                base = f'http://127.0.0.1:{port}'
                with httpx.Client(base_url=base, trust_env=False, timeout=5) as http:
                    def health():
                        if any(p.poll() is not None for p in processes):
                            raise RuntimeError((root / 'services.log').read_text()[-2000:])
                        try:
                            return http.get('/health').status_code == 200
                        except httpx.HTTPError:
                            return False
                    until(health, bool)
                readiness = time.perf_counter() - started
                with SEAClient(key, base_url=base) as client:
                    request_started = time.perf_counter()
                    source = client.submit_prediction([[smiles]], dataset_id=dataset_id, statistics_profile='balanced')
                    evidence = client.wait(source, wait_for='evidence', timeout=120, poll_interval=.05)
                    first_evidence = time.perf_counter() - request_started
                    screen = client.wait(source, wait_for='screening', timeout=120, poll_interval=.05)
                    screen_seconds = time.perf_counter() - request_started
                    child = screen['analysis']['refinement_job_id']
                    until(lambda: client.get_job(child), lambda job: job['status'] == 'running' and job['progress'] > 0)
                    print('Cold screening complete; measuring foreground search during refinement.', flush=True)
                    second_started = time.perf_counter()
                    second = client.submit_prediction([[smiles]], dataset_id=dataset_id, statistics_profile='balanced')
                    second_evidence = client.wait(second, wait_for='evidence', timeout=120, poll_interval=.05)
                    second_evidence_seconds = time.perf_counter() - second_started
                    second_screen = client.wait(second, wait_for='screening', timeout=120, poll_interval=.05)
                    second_screen_seconds = time.perf_counter() - second_started
                    for cid in (child, second_screen['analysis']['refinement_job_id']):
                        client.cancel(cid)
                        until(lambda: client.get_job(cid), lambda job: job['status'] in ('cancelled', 'completed'))
                    immutable = client.get_result(source) == screen
                    checkpoint_counts = []
                    import numpy as np
                    for path in (root / 'cache' / 'null-v2' / dataset_id).glob('*.npz'):
                        with np.load(path) as archive:
                            checkpoint_counts.append(len(archive['values']))
                    assert immutable
                    assert evidence['analysis']['phase'] == 'evidence'
                    assert second_evidence['analysis']['phase'] in ('evidence', 'screening')
                    return {
                        'dataset_id': dataset_id, 'model_id': screen['dataset']['model_id'],
                        'query': smiles, 'targets': screen['dataset']['target_count'],
                        'reference_molecules': screen['dataset']['molecule_count'],
                        'workers': 2, 'foreground_reserved': 1,
                        'api_preload_readiness_seconds': readiness,
                        'cold_evidence_seconds': first_evidence,
                        'cold_screening_seconds': screen_seconds,
                        'evidence_during_refinement_seconds': second_evidence_seconds,
                        'screening_during_refinement_seconds': second_screen_seconds,
                        'source_immutable': immutable,
                        'refinement_status': client.get_job(child)['status'],
                        'refinement_complete': client.get_job(child)['status'] == 'completed',
                        'planned_trials': screen['analysis']['profile']['planned_trials'],
                        'cached_size_pairs': len(checkpoint_counts),
                        'largest_saved_prefix': max(checkpoint_counts, default=0),
                        'scope': 'Single local run using isolated API, two workers, database and empty cache. Preload readiness is separate. Full cold ChEMBL refinement deliberately paused; no final-latency or biological-accuracy claim. Repeated query reuses completed sampling prefixes during background extension.',
                    }
        finally:
            for process in processes:
                process.terminate()
            for process in processes:
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            if old_home is None:
                os.environ.pop('SEA_HOME', None)
            else:
                os.environ['SEA_HOME'] = old_home


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True)
    parser.add_argument('--smiles', default='COc1ccc2[nH]cc(c2c1)C(CC)CN')
    parser.add_argument('--port', type=int, default=19031)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = measure(args.dataset_id, args.smiles, args.port)
    text = json.dumps(report, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.write_text(text, encoding='utf-8')
    print(text)
