"""Convert structure files to validated isomeric SMILES without dropping records."""
import csv
import io
from typing import Literal

from rdkit import Chem

from .chemistry import molecule

StructureFormat = Literal["sdf", "mol", "smiles", "csv"]


class StructureImportError(ValueError):
    pass


def _records(content: str, format: StructureFormat):
    if format == "sdf":
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO(content.encode("utf-8")), sanitize=True, removeHs=True,
            strictParsing=True,
        )
        for number, mol in enumerate(supplier, 1):
            yield number, mol.GetProp("_Name") if mol is not None and mol.HasProp("_Name") else "", mol
    elif format == "mol":
        if any(line.strip() == "$$$$" for line in content.splitlines()):
            raise StructureImportError({"code": "multiple_mol_records", "message": "SDFレコードを含むファイルはSDF形式で読み込んでください"})
        mol = Chem.MolFromMolBlock(content, sanitize=True, removeHs=True, strictParsing=True)
        yield 1, mol.GetProp("_Name") if mol is not None and mol.HasProp("_Name") else "", mol
    elif format == "csv":
        reader = csv.reader(io.StringIO(content, newline=""), strict=True)
        header = next(reader, [])
        columns = [value.strip().casefold() for value in header]
        if "smiles" not in columns:
            raise StructureImportError({"code": "missing_smiles_column", "message": "smiles列が必要です"})
        smiles_index = columns.index("smiles")
        name_index = columns.index("name") if "name" in columns else None
        for number, row in enumerate(reader, 1):
            if not any(value.strip() for value in row):
                continue
            smiles = row[smiles_index].strip() if len(row) > smiles_index else ""
            name = row[name_index].strip() if name_index is not None and len(row) > name_index else ""
            yield number, name, Chem.MolFromSmiles(smiles)
    elif format == "smiles":
        for number, line in enumerate(content.splitlines(), 1):
            if not line.strip():
                continue
            parts = line.strip().split(maxsplit=1)
            yield number, parts[1] if len(parts) > 1 else "", Chem.MolFromSmiles(parts[0])
    else:
        raise StructureImportError({"code": "unsupported_format", "message": "未対応の構造形式です"})


def import_structures(content: str, format: StructureFormat, maximum: int = 1000):
    """Validate the entire input before returning any molecules.

    Fragments, charges and stereochemistry remain in the returned SMILES.
    SEA standardization and deduplication happen when a job is submitted.
    """
    content = content.removeprefix("\ufeff")
    if not content.strip() or "\x00" in content:
        raise StructureImportError({"code": "empty_or_binary_input", "message": "空でないテキスト形式の構造ファイルを指定してください"})
    imported, errors, count = [], [], 0
    try:
        for number, name, mol in _records(content, format):
            count += 1
            if count > maximum:
                raise StructureImportError({"code": "too_many_molecules", "limit": maximum, "message": f"読み込み上限は{maximum}分子です"})
            try:
                if mol is None or not mol.GetNumAtoms():
                    raise ValueError("構造を読み取れません")
                if any(atom.HasQuery() or atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()) or any(bond.HasQuery() for bond in mol.GetBonds()):
                    raise ValueError("クエリ原子・結合やR基を含まない、具体的な分子構造を入力してください")
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                molecule(smiles)  # Use the same eligibility rules as job submission.
                imported.append({"record": number, "name": name, "smiles": smiles})
            except Exception as exc:
                errors.append({"record": number, "name": name, "message": str(exc)})
    except csv.Error as exc:
        raise StructureImportError({"code": "invalid_csv", "message": f"CSVを読み取れません: {exc}"}) from exc
    except StructureImportError:
        raise
    except (ValueError, RuntimeError) as exc:
        raise StructureImportError({"code": "invalid_file", "message": f"構造ファイルを読み取れません: {exc}"}) from exc
    if errors:
        raise StructureImportError({"code": "invalid_structures", "message": "不正な構造が含まれています", "records": errors})
    if not imported:
        raise StructureImportError({"code": "empty_input", "message": "構造が含まれていません"})
    return {"format": format, "count": len(imported), "molecules": imported}
