"""Optional real-process + browser smoke check; requires Chrome and Node."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import httpx
from rdkit import Chem
from sea.cli import demo
from sea.jobs import create_key
from sea import SEAClient

started = time.perf_counter()
Path('.cache').mkdir(exist_ok=True)
Path('var').mkdir(exist_ok=True)
# Keep smoke keys, jobs, datasets and browser sessions out of the lab workspace.
with tempfile.TemporaryDirectory(prefix='sea-browser-', dir='.cache') as temporary:
    root = Path(temporary).resolve()
    os.environ['SEA_HOME'] = str(root / 'home')
    dataset = demo()
    key = create_key('browser-smoke')
    processes = []
    try:
        processes.append(subprocess.Popen([sys.executable, '-m', 'uvicorn', 'sea.api:app', '--port', '19021'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))
        processes.append(subprocess.Popen([sys.executable, '-m', 'sea', 'worker'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))
        processes.append(subprocess.Popen(['google-chrome', '--headless', '--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage', '--no-first-run', f'--user-data-dir={root / "chrome"}', '--remote-debugging-port=19022', 'about:blank'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))
        for _ in range(60):
            try:
                if httpx.get('http://127.0.0.1:19021/health').status_code == 200 and httpx.get('http://127.0.0.1:19022/json').status_code == 200:
                    break
            except httpx.HTTPError:
                pass
            time.sleep(.5)
        else:
            raise RuntimeError('Services failed to start')
        blocks = [Chem.MolToMolBlock(Chem.MolFromSmiles(smiles)) for smiles in ('CCO', 'CCCO')]
        subprocess.run(['node', 'scripts/browser-smoke.mjs'], check=True, env={
            **os.environ, 'SEA_SMOKE_KEY': key, 'SEA_SMOKE_MOL': blocks[0],
            'SEA_SMOKE_SDF': ''.join(block + '$$$$\n' for block in blocks),
        }, timeout=180)
        with SEAClient(key, base_url='http://127.0.0.1:19021') as client:
            source = min((job for job in client.list_jobs() if job['kind'] == 'prediction'), key=lambda job: job['created'])['id']
            original = client.get_result(source)
            research_started = time.perf_counter()
            research_id = client.refine(source, statistics_profile='research')
            assert client.refine(source, statistics_profile='research') == research_id
            final = client.wait(research_id, timeout=120, poll_interval=.1)
            assert final['analysis']['complete'] and final['analysis']['profile']['planned_trials'] == 999999
            assert all(row['statistics']['phase'] == 'final' and
                       (row['method'] != 'monte_carlo' or row['statistics']['trials'] == 999999)
                       for row in final['results'])
            assert client.get_result(source) == original
            assert client.get_job(source)['refinement_job_id'] == original['analysis']['refinement_job_id']
            assert 'stream_policy' in client.get_csv(research_id).splitlines()[0]
            warmed = subprocess.run([sys.executable, '-m', 'sea', 'warm', dataset, '--profile', 'balanced', '--limit', '1'],
                                    check=True, capture_output=True, text=True, timeout=30)
            assert json.loads(warmed.stdout)['status'] == 'ready'
            print(json.dumps({'sdk_research_final_budget': 999999, 'research_seconds': time.perf_counter() - research_started,
                              'source_immutable': True, 'idempotent_refinement': True, 'warm_cli': 'passed'}))
        print(json.dumps({'elapsed_seconds': time.perf_counter() - started, 'dataset_id': dataset}))
    finally:
        for process in processes:
            process.terminate()
        for process in processes:
            try:
                process.wait(timeout=8)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
