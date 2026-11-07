import json
import sqlite3
import pytest
from sea.data import build_chembl, load_dataset


def test_chembl_filters_and_deterministic_build(tmp_path,monkeypatch):
    monkeypatch.setenv('SEA_HOME',str(tmp_path/'out'))
    source=tmp_path/'chembl.db'
    con=sqlite3.connect(source)
    con.executescript('''
    CREATE TABLE version(name TEXT); INSERT INTO version VALUES ('ChEMBL_TEST');
    CREATE TABLE activities(molregno INTEGER, assay_id INTEGER, standard_type TEXT,pchembl_value REAL,activity_id INTEGER,doc_id INTEGER,standard_relation TEXT,data_validity_comment TEXT,potential_duplicate INTEGER);
    CREATE TABLE assays(assay_id INTEGER,tid INTEGER,confidence_score INTEGER,assay_type TEXT);
    INSERT INTO assays VALUES(1,1,9,'B'); INSERT INTO assays VALUES(2,1,8,'B');
    CREATE TABLE target_dictionary(tid INTEGER,chembl_id TEXT,pref_name TEXT,tax_id INTEGER,target_type TEXT,organism TEXT);
    INSERT INTO target_dictionary VALUES(1,'CHEMBL_TARGET','Test protein',9606,'SINGLE PROTEIN','Homo sapiens');
    CREATE TABLE protein_classification(protein_class_id INTEGER,parent_id INTEGER,pref_name TEXT,class_level INTEGER);
    INSERT INTO protein_classification VALUES(0,NULL,'Protein class',0),(1,0,'Enzyme',1),(6,1,'Kinase',2);
    CREATE TABLE target_components(tid INTEGER,component_id INTEGER);
    INSERT INTO target_components VALUES(1,1);
    CREATE TABLE component_class(component_id INTEGER,protein_class_id INTEGER);
    INSERT INTO component_class VALUES(1,6);
    CREATE TABLE molecule_dictionary(molregno INTEGER,chembl_id TEXT);
    CREATE TABLE molecule_hierarchy(molregno INTEGER,parent_molregno INTEGER);
    CREATE TABLE compound_structures(molregno INTEGER,canonical_smiles TEXT);
    ''')
    for i in range(1,13):
        con.execute('INSERT INTO molecule_dictionary VALUES (?,?)',(i,f'CHEMBL_MOL{i}'))
        con.execute('INSERT INTO compound_structures VALUES (?,?)',(i,'C'*i+'O'))
        con.execute('INSERT INTO activities VALUES (?,?,?,?,?,?,?,?,?)',(i,1,'Ki',7,i,1,'=',None,0))
    # High potency must not override a median below the threshold.
    con.execute("INSERT INTO activities VALUES (12,1,'Ki',3,30,1,'=',NULL,0)")
    # Invalid-quality records do not change the result.
    con.execute("INSERT INTO activities VALUES (1,1,'Ki',0,31,1,'=',NULL,1)")
    con.execute("INSERT INTO activities VALUES (1,2,'Ki',0,32,1,'=',NULL,0)")
    con.commit();con.close()
    settings={'min_pchembl':6,'min_ligands':10}
    first=build_chembl(source,'ChEMBL_TEST',settings,repeats=16)
    second=build_chembl(source,'ChEMBL_TEST',settings,repeats=16)
    assert first==second
    data,pool=load_dataset(first)
    assert len(data['targets'][0]['smiles'])==11
    assert data['targets'][0]['family']=='Kinase'
    assert data['targets'][0]['organism']=='Homo sapiens'
    assert data['targets'][0]['tax_id']==9606
    assert data['manifest']['counts']['quality_excluded']==1
    assert data['manifest']['data_license']=='CC-BY-SA-3.0'
    assert all(data['targets'][0]['evidence'].values())
    with pytest.raises(ValueError): build_chembl(source,'ChEMBL_WRONG',settings,repeats=16)
