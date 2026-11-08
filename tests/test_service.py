import json
import pytest
from fastapi.testclient import TestClient
from sea.api import app
from sea.cli import demo
from sea.jobs import create_key, claim, process, connection


@pytest.fixture()
def service(tmp_path,monkeypatch):
    monkeypatch.setenv('SEA_HOME',str(tmp_path))
    dataset=demo()
    a,b=create_key('alice'),create_key('bob')
    client=TestClient(app)
    return client,dataset,{'X-API-Key':a},{'X-API-Key':b}


def test_auth_jobs_and_results(service):
    client,dataset,a,b=service
    assert client.get('/api/v1/datasets').status_code==401
    assert client.get('/api/v1/datasets',headers=a).json()[0]['id']==dataset
    payload={'dataset_id':dataset,'left':['CCO','CCCO'],'right':['CCCCO','CCO']}
    response=client.post('/api/v1/comparisons',headers=a,json=payload)
    assert response.status_code==202
    jid=response.json()['job_id']
    assert client.get(f'/api/v1/jobs/{jid}',headers=b).status_code==404
    assert client.get(f'/api/v1/jobs/{jid}/result',headers=a).status_code==409
    job=claim('test')
    assert claim('test-2') is None
    process(job)
    status=client.get(f'/api/v1/jobs/{jid}',headers=a).json()
    assert status['status']=='completed',status
    result=client.get(f'/api/v1/jobs/{jid}/result',headers=a).json()
    assert result['results'][0]['comparisons']==1
    assert result['results'][0]['method']=='monte_carlo'
    assert len(result['normalization']['left'])==2
    csv=client.get(f'/api/v1/jobs/{jid}/result?format=csv',headers=a)
    assert dataset in csv.text
    assert client.get('/api/v1/structure',params={'smiles':'CCO'},headers=a).status_code==200
    assert client.get('/').status_code==200
    assert client.get('/openapi.json').status_code==200


def test_prediction_validation_and_cancel(service):
    client,dataset,a,b=service
    bad=client.post('/api/v1/predictions',headers=a,json={'dataset_id':dataset,'queries':[['CCO','bad']]})
    assert bad.status_code==422 and bad.json()['detail']['rows'][0]['row']==2
    assert client.post('/api/v1/predictions',headers=a,json={'dataset_id':'../bad','queries':[['CCO']]}).status_code==422
    payload={'dataset_id':dataset,'queries':[['CCO']],'top_k':1}
    jid=client.post('/api/v1/predictions',headers=a,json=payload).json()['job_id']
    process(claim('test'))
    result=client.get(f'/api/v1/jobs/{jid}/result',headers=a).json()
    assert len(result['results'])==1 and result['results'][0]['comparisons']==3
    jid=client.post('/api/v1/predictions',headers=a,json=payload).json()['job_id']
    assert client.post(f'/api/v1/jobs/{jid}/cancel',headers=b).status_code==404
    assert client.post(f'/api/v1/jobs/{jid}/cancel',headers=a).json()['status']=='cancelled'
    assert claim('test') is None


def test_session_and_limit(service,monkeypatch):
    client,dataset,a,b=service
    assert client.post('/api/v1/session',json={'api_key':a['X-API-Key']}).status_code==200
    assert client.get('/api/v1/datasets').status_code==200
    payload={'dataset_id':dataset,'queries':[['CCO']]}
    assert client.post('/api/v1/predictions',json=payload).status_code==403
    assert client.post('/api/v1/predictions',headers={'X-SEA-Request':'1','Origin':'https://attacker.example'},json=payload).status_code==403
    assert client.delete('/api/v1/session',headers={'X-SEA-Request':'1'}).status_code==200
    assert client.get('/api/v1/datasets').status_code==401
    monkeypatch.setenv('SEA_MAX_UPLOAD_BYTES','20')
    assert client.post('/api/v1/session',content='a'*21).status_code==413


def test_claim_is_atomic_and_running_cancel(service):
    from concurrent.futures import ThreadPoolExecutor
    from sea.jobs import submit
    client,dataset,a,b=service
    jid=client.post('/api/v1/comparisons',headers=a,json={'dataset_id':dataset,'left':['CCO'],'right':['CCCO']}).json()['job_id']
    with ThreadPoolExecutor(max_workers=2) as executor:
        results=list(executor.map(claim,['first','second']))
    assert sum(r is not None for r in results)==1
    job=next(r for r in results if r)
    client.post(f'/api/v1/jobs/{jid}/cancel',headers=a)
    process(job)
    assert client.get(f'/api/v1/jobs/{jid}',headers=a).json()['status']=='cancelled'


def test_worker_restart_recovers_claimed_job(service):
    import os
    import subprocess
    import sys
    import time
    client,dataset,a,b=service
    jid=client.post('/api/v1/comparisons',headers=a,json={'dataset_id':dataset,'left':['CCO'],'right':['CCCO']}).json()['job_id']
    assert claim('crashed-worker')['id']==jid
    worker=subprocess.Popen([sys.executable,'-m','sea','worker','--count','2'],env=os.environ.copy())
    try:
        for _ in range(100):
            row=client.get(f'/api/v1/jobs/{jid}',headers=a).json()
            if row['status']=='completed':break
            if row['status']=='failed':pytest.fail(row['error'])
            time.sleep(.1)
        assert row['status']=='completed'
    finally:
        worker.terminate()
        worker.wait(timeout=10)


def test_background_samples_are_reused_across_jobs(service,monkeypatch):
    from sea.background import NullCache
    from sea.jobs import background_cache
    client,dataset,a,b=service
    original=NullCache.sample
    calls=[]
    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)
    monkeypatch.setattr(NullCache,'sample',counted)
    payload={'dataset_id':dataset,'left':['CCO'],'right':['CCO']}
    outputs=[]
    for _ in range(2):
        jid=client.post('/api/v1/comparisons',headers=a,json=payload).json()['job_id']
        process(claim('test'))
        outputs.append(client.get(f'/api/v1/jobs/{jid}/result',headers=a).json()['results'])
    assert outputs[0]==outputs[1]
    assert len(calls)==1  # One null distribution, reused by the second job.
    from sea.data import home
    cache=background_cache(str(home().resolve()),dataset)
    assert (1,1) in cache
    assert background_cache(str(home().resolve()),dataset) is cache
    assert background_cache(str(home().resolve()),'different-dataset') is not cache


def test_input_preview_preserves_fragments_and_requires_authentication(service):
    client,dataset,a,b=service
    params={'smiles':'CC(=O)[O-].[Na+]','standardize':'false'}
    assert client.get('/api/v1/structure',params=params).status_code==401
    response=client.get('/api/v1/structure',params=params,headers=a)
    assert response.status_code==200 and response.headers['content-type'].startswith('image/svg+xml')
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    mol=Chem.MolFromSmiles(params['smiles'])
    drawer=rdMolDraw2D.MolDraw2DSVG(300,180)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    assert response.text==drawer.GetDrawingText()
    assert client.get('/api/v1/structure',params={'smiles':params['smiles']},headers=a).text!=response.text
    for invalid in ('', 'bad', 'C'*10001):
        assert client.get('/api/v1/structure',params={'smiles':invalid,'standardize':'false'},headers=a).status_code==422


def test_structure_import_and_prediction(service, monkeypatch):
    from rdkit import Chem
    client,dataset,a,b=service
    mol = Chem.MolFromSmiles('CCO')
    mol.SetProp('_Name','Ethanol')
    sdf = Chem.MolToMolBlock(mol)+'$$$$\n'
    payload={'format':'sdf','content':sdf}
    assert client.post('/api/v1/structures/import',json=payload).status_code==401
    response=client.post('/api/v1/structures/import',headers=a,json=payload)
    assert response.status_code==200,response.text
    converted=response.json()
    assert converted['molecules']==[{'record':1,'name':'Ethanol','smiles':'CCO'}]
    jid=client.post('/api/v1/predictions',headers=a,json={
        'dataset_id':dataset,'queries':[[m['smiles']] for m in converted['molecules']]
    }).json()['job_id']
    process(claim('test'))
    assert client.get(f'/api/v1/jobs/{jid}/result',headers=a).json()['results'][0]['raw_score']>0
    invalid=client.post('/api/v1/structures/import',headers=a,json={'format':'sdf','content':sdf+'invalid\n$$$$\n'})
    assert invalid.status_code==422
    assert invalid.json()['detail']['records'][0]['record']==2
    assert client.post('/api/v1/structures/import',headers=a,json={'format':'pdb','content':'CCO'}).status_code==422
    monkeypatch.setenv('SEA_MAX_SET_MOLECULES','1')
    assert client.post('/api/v1/structures/import',headers=a,json={'format':'smiles','content':'CCO\nCCCO'}).json()['detail']['code']=='too_many_molecules'
    assert client.post('/api/v1/session',json={'api_key':a['X-API-Key']}).status_code==200
    assert client.post('/api/v1/structures/import',json=payload).status_code==403
    assert client.post('/api/v1/structures/import',headers={'X-SEA-Request':'1'},json=payload).status_code==200
