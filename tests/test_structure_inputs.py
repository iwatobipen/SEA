import pytest
from rdkit import Chem

from sea.structure_inputs import import_structures, StructureImportError


def molblock(smiles, name="", v3000=False):
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp("_Name", name)
    return Chem.MolToMolBlock(mol, forceV3000=v3000)


def test_sdf_preserves_records_stereo_fragments_and_names():
    content = molblock("C[C@H](O)C(=O)O", "Chiral molecule") + "$$$$\n"
    content += molblock("CC(=O)[O-].[Na+]", "Salt", v3000=True) + "$$$$\n"
    result = import_structures(content, "sdf")
    assert result["count"] == 2
    first, second = result["molecules"]
    assert first == {"record": 1, "name": "Chiral molecule", "smiles": "C[C@H](O)C(=O)O"}
    assert second == {"record": 2, "name": "Salt", "smiles": "CC(=O)[O-].[Na+]"}
    assert import_structures(molblock("F/C=C/F", v3000=True), "mol")["molecules"][0]["smiles"] == "F/C=C/F"


def test_invalid_record_rejects_whole_import_with_position():
    sdf = molblock("CCO", "Good") + "$$$$\ninvalid molecule\n$$$$\n"
    with pytest.raises(StructureImportError) as caught:
        import_structures(sdf, "sdf")
    detail = caught.value.args[0]
    assert detail["code"] == "invalid_structures"
    assert [row["record"] for row in detail["records"]] == [2]
    with pytest.raises(StructureImportError) as caught:
        import_structures("CCO first\nwrong second\n[Na+] inorganic\n", "smiles")
    assert [row["record"] for row in caught.value.args[0]["records"]] == [2, 3]


def test_csv_and_smiles_text():
    csv = '\ufeffSMILES,name\r\nCCO,"name, with comma"\r\nCCCO,"two\nlines"\r\n'
    result = import_structures(csv, "csv")
    assert result["molecules"] == [
        {"record": 1, "name": "name, with comma", "smiles": "CCO"},
        {"record": 2, "name": "two\nlines", "smiles": "CCCO"},
    ]
    result = import_structures("CCO ethanol\n\nOCC\tduplicate retained\n", "smiles")
    assert result["count"] == 2
    assert [m["smiles"] for m in result["molecules"]] == ["CCO", "CCO"]
    assert result["molecules"][1]["record"] == 3
    with pytest.raises(StructureImportError) as caught:
        import_structures("name,smiles\nmissing,\n", "csv")
    assert caught.value.args[0]["records"][0]["record"] == 1


@pytest.mark.parametrize("content,format,code", [
    ("", "sdf", "empty_or_binary_input"),
    ("a\x00b", "mol", "empty_or_binary_input"),
    ("smiles,name\n", "csv", "empty_input"),
    ("compound\nCCO", "csv", "missing_smiles_column"),
    ('smiles,name\nCCO,"unterminated', "csv", "invalid_csv"),
    ("not a molfile", "mol", "invalid_structures"),
])
def test_empty_or_malformed_input(content, format, code):
    with pytest.raises(StructureImportError) as caught:
        import_structures(content, format)
    assert caught.value.args[0]["code"] == code


def test_record_limit_does_not_silently_truncate():
    with pytest.raises(StructureImportError) as caught:
        import_structures("CCO\nOCC", "smiles", maximum=1)
    assert caught.value.args[0]["code"] == "too_many_molecules"


def test_query_structures_and_mislabeled_sdf_are_rejected():
    query = Chem.MolToMolBlock(Chem.MolFromSmarts('[C,N]CC'))
    with pytest.raises(StructureImportError) as caught:
        import_structures(query, "mol")
    assert caught.value.args[0]["records"][0]["record"] == 1
    with pytest.raises(StructureImportError):
        import_structures("CC*", "smiles")
    with pytest.raises(StructureImportError) as caught:
        import_structures(molblock("CCO")+'$$$$\n'+molblock("CCN")+'$$$$\n', "mol")
    assert caught.value.args[0]["code"] == "multiple_mol_records"
