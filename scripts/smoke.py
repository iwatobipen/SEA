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
        }, timeout=90)
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
