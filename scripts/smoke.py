"""Optional real-process + browser smoke check; requires Chrome and Node."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import httpx
from sea.cli import demo
from sea.jobs import create_key

started=time.perf_counter()
dataset=demo()
key=create_key('browser-smoke')
processes=[]
try:
    processes.append(subprocess.Popen([sys.executable,'-m','uvicorn','sea.api:app','--port','19021'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL))
    processes.append(subprocess.Popen([sys.executable,'-m','sea','worker'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL))
    profile=Path('var/chrome-smoke').resolve()
    processes.append(subprocess.Popen(['google-chrome','--headless','--no-sandbox','--disable-gpu','--disable-dev-shm-usage','--no-first-run',f'--user-data-dir={profile}','--remote-debugging-port=19022','about:blank'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL))
    for _ in range(60):
        try:
            if httpx.get('http://127.0.0.1:19021/health').status_code==200 and httpx.get('http://127.0.0.1:19022/json').status_code==200:
                break
        except httpx.HTTPError:
            pass
        time.sleep(.5)
    else:
        raise RuntimeError('Services failed to start')
    subprocess.run(['node','scripts/browser-smoke.mjs'],check=True,env={**os.environ,'SEA_SMOKE_KEY':key},timeout=60)
    print(json.dumps({'elapsed_seconds':time.perf_counter()-started,'dataset_id':dataset}))
finally:
    for process in processes:
        process.terminate()
    for process in processes:
        try:process.wait(timeout=8)
        except subprocess.TimeoutExpired:
            process.kill();process.wait()
