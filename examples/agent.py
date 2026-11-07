"""SEA_API_KEY and SEA_DATASET_ID must be set. No external model dependency."""
import json
import os
import time
import httpx

with httpx.Client(base_url=os.getenv('SEA_URL','http://127.0.0.1:8000'),headers={'X-API-Key':os.environ['SEA_API_KEY']},timeout=30) as client:
    response=client.post('/api/v1/predictions',json={
        'dataset_id':os.environ['SEA_DATASET_ID'], 'queries':[['CCO']], 'top_k':10})
    response.raise_for_status()
    job_id=response.json()['job_id']
    deadline=time.monotonic()+3600
    while time.monotonic()<deadline:
        response=client.get(f'/api/v1/jobs/{job_id}')
        response.raise_for_status()
        status=response.json()
        if status['status']=='completed':
            result=client.get(f'/api/v1/jobs/{job_id}/result')
            result.raise_for_status()
            print(json.dumps(result.json(),ensure_ascii=False,indent=2))
            break
        if status['status'] in ('failed','cancelled'):
            raise RuntimeError(status)
        time.sleep(1)
    else:
        raise TimeoutError(job_id)
