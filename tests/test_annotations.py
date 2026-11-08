import csv
import io
import json
import sqlite3
import pytest
from fastapi.testclient import TestClient
from sea.annotations import (read_target_annotations, annotate_dataset, load_annotations,
                             annotation_path, target_info, enrich_results)
from sea.data import save_dataset, digest_file, load_dataset
from sea.api import app
from sea.jobs import create_key, submit, claim, process, connection


@pytest.fixture
def annotated_dataset(tmp_path,monkeypatch):
    monkeypatch.setenv('SEA_HOME',str(tmp_path/'service'))
    source=tmp_path/'chembl.db'
    con=sqlite3.connect(source)
    con.executescript('''
    CREATE TABLE version(name TEXT); INSERT INTO version VALUES ('ChEMBL_TEST');
    CREATE TABLE target_dictionary(tid INTEGER,chembl_id TEXT,pref_name TEXT,organism TEXT,tax_id INTEGER);
    INSERT INTO target_dictionary VALUES (1,'C1','Opioid receptor','Homo sapiens',9606),
      (2,'C2','Protein kinase','Mus musculus',10090),(3,'C3','Potassium channel','Rattus norvegicus',10116),
      (4,'C4','Unclassified target',NULL,NULL);
    CREATE TABLE protein_classification(protein_class_id INTEGER,parent_id INTEGER,pref_name TEXT,class_level INTEGER);
    INSERT INTO protein_classification VALUES (0,NULL,'Protein class',0),(1,0,'Enzyme',1),
      (6,1,'Kinase',2),(11,0,'Membrane receptor',1),
      (1020,11,'Family A G protein-coupled receptor',2),(2000,1020,'Opioid receptor',3),
      (1005,0,'Ion channel',1),(1019,1005,'Voltage-gated ion channel',2);
    CREATE TABLE target_components(tid INTEGER,component_id INTEGER);
    INSERT INTO target_components VALUES (1,1),(2,2),(3,3);
    CREATE TABLE component_class(component_id INTEGER,protein_class_id INTEGER);
    INSERT INTO component_class VALUES (1,2000),(1,1020),(2,6),(3,1019);
    ''')
    con.commit();con.close()
    targets=[{'id':f'C{i}','name':'Legacy name','smiles':['C'*n+'O' for n in range(1,13)]} for i in range(1,5)]
    dataset=save_dataset(targets,{'source':'ChEMBL','release':'ChEMBL_TEST','source_sha256':digest_file(source)},sizes=(1,10,12),repeats=16)
    return dataset,source


def test_hierarchy_and_missing_metadata(annotated_dataset):
    _,source=annotated_dataset
    with sqlite3.connect(source) as con:
        metadata=read_target_annotations(con,['C1','C2','C3','C4'])
    assert metadata['C1']['family']=='GPCR'
    assert metadata['C1']['families']==['GPCR']
    assert metadata['C1']['protein_classes'][-1][-1]=='Opioid receptor'
    assert metadata['C2']['family']=='Kinase'
    assert metadata['C2']['organism']=='Mus musculus' and metadata['C2']['tax_id']==10090
    assert metadata['C3']['family']=='Ion channel'
    assert metadata['C4']['family']=='Unclassified'
    assert metadata['C4']['organism'] is None
    # A suggestive name alone must never imply a classification.
    assert target_info({'id':'unknown','name':'Kinase-like test'})['family'] is None


def test_annotation_preserves_dataset_and_validates_source(annotated_dataset):
    dataset,source=annotated_dataset
    from sea.data import home
    before=digest_file(home()/'datasets'/f'{dataset}.json')
    report=annotate_dataset(dataset,source)
    assert report['target_count']==4
    assert before==digest_file(home()/'datasets'/f'{dataset}.json')
    old={'dataset':{'id':dataset},'results':[{'target':{'id':'C1','name':'Old name'},'raw_score':1}]}
    updated=enrich_results(old)
    assert updated['results'][0]['target']['family']=='GPCR'
    assert updated['results'][0]['target']['name']=='Opioid receptor'
    assert updated['results'][0]['raw_score']==1
    first_id=report['annotation_id']
    assert annotate_dataset(dataset,source)['annotation_id']==first_id
    # Sidecar replacement invalidates the read cache without restarting the API.
    path=annotation_path(dataset)
    document=json.loads(path.read_text())
    document['targets']['C1']['name']='Updated target name'
    path.write_text(json.dumps(document))
    assert load_annotations(dataset)['targets']['C1']['name']=='Updated target name'
    with sqlite3.connect(source) as con:
        con.execute("UPDATE target_dictionary SET pref_name='Changed source' WHERE tid=1")
    with pytest.raises(ValueError,match='SHA256'):
        annotate_dataset(dataset,source)
    assert load_annotations(dataset)['targets']['C1']['name']=='Updated target name'


def test_result_json_csv_and_target_search_enrich_old_jobs(annotated_dataset):
    dataset,source=annotated_dataset
    client=TestClient(app)
    headers={'X-API-Key':create_key('test')}
    page=client.get('/')
    assert page.headers['cache-control']=='no-cache'
    assert '<html lang="en">' in page.text
    assert 'Target / family / species' in page.text
    assert client.get('/static/app.js?v=2').headers['cache-control']=='no-cache'
    response=client.post('/api/v1/predictions',headers=headers,json={'dataset_id':dataset,'queries':[['CO']],'top_k':4})
    jid=response.json()['job_id']
    process(claim('test'))
    assert client.get(f'/api/v1/jobs/{jid}',headers=headers).json()['status']=='completed'
    # This job completed before annotations existed.
    with connection() as con:
        stored_before=con.execute('SELECT result FROM jobs WHERE id=?',(jid,)).fetchone()[0]
    annotate_dataset(dataset,source)
    bundle=client.get(f'/api/v1/jobs/{jid}/result',headers=headers).json()
    target=next(r['target'] for r in bundle['results'] if r['target']['id']=='C2')
    assert target['name']=='Protein kinase' and target['family']=='Kinase'
    assert target['organism']=='Mus musculus' and target['tax_id']==10090
    assert bundle['target_annotations']['release']=='ChEMBL_TEST'
    exported=client.get(f'/api/v1/jobs/{jid}/result?format=csv',headers=headers)
    row=next(r for r in csv.DictReader(io.StringIO(exported.text)) if r['target_id']=='C3')
    assert row['target_family']=='Ion channel' and row['target_tax_id']=='10116'
    assert json.loads(row['target_protein_classes'])[0][0]=='Ion channel'
    targets=client.get('/api/v1/targets',headers=headers,params={'dataset_id':dataset,'q':'GPCR'}).json()
    assert len(targets)==1 and targets[0]['id']=='C1'
    with connection() as con:
        assert con.execute('SELECT result FROM jobs WHERE id=?',(jid,)).fetchone()[0]==stored_before
    # Future jobs carry metadata in their stored result as well.
    jid=client.post('/api/v1/predictions',headers=headers,json={'dataset_id':dataset,'queries':[['CO']],'top_k':4}).json()['job_id']
    process(claim('test'))
    with connection() as con:
        stored=json.loads(con.execute('SELECT result FROM jobs WHERE id=?',(jid,)).fetchone()[0])
    assert all(r['target']['organism'] or r['target']['id']=='C4' for r in stored['results'])
    assert stored['target_annotations']['release']=='ChEMBL_TEST'
