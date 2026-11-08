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
            raise StructureImportError({"code": "multiple_mol_records", "message": "Import files containing SDF records as SDF, not MOL."})
        mol = Chem.MolFromMolBlock(content, sanitize=True, removeHs=True, strictParsing=True)
        yield 1, mol.GetProp("_Name") if mol is not None and mol.HasProp("_Name") else "", mol
    elif format == "csv":
        reader = csv.reader(io.StringIO(content, newline=""), strict=True)
        header = next(reader, [])
        columns = [value.strip().casefold() for value in header]
        if "smiles" not in columns:
            raise StructureImportError({"code": "missing_smiles_column", "message": "A smiles column is required."})
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
        raise StructureImportError({"code": "unsupported_format", "message": "Unsupported structure format."})


def import_structures(content: str, format: StructureFormat, maximum: int = 1000):
    """Validate the entire input before returning any molecules.

    Fragments, charges and stereochemistry remain in the returned SMILES.
    SEA standardization and deduplication happen when a job is submitted.
    """
    content = content.removeprefix("\ufeff")
    if not content.strip() or "\x00" in content:
        raise StructureImportError({"code": "empty_or_binary_input", "message": "Provide a nonempty structure file in text format."})
    imported, errors, count = [], [], 0
    try:
        for number, name, mol in _records(content, format):
            count += 1
            if count > maximum:
                raise StructureImportError({"code": "too_many_molecules", "limit": maximum, "message": f"The import limit is {maximum} molecules."})
            try:
                if mol is None or not mol.GetNumAtoms():
                    raise ValueError("Unable to read the molecular structure.")
                if any(atom.HasQuery() or atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()) or any(bond.HasQuery() for bond in mol.GetBonds()):
                    raise ValueError("Enter a concrete molecule without query atoms, query bonds or R groups.")
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                molecule(smiles)  # Use the same eligibility rules as job submission.
                imported.append({"record": number, "name": name, "smiles": smiles})
            except Exception as exc:
                errors.append({"record": number, "name": name, "message": str(exc)})
    except csv.Error as exc:
        raise StructureImportError({"code": "invalid_csv", "message": f"Unable to read the CSV file: {exc}"}) from exc
    except StructureImportError:
        raise
    except (ValueError, RuntimeError) as exc:
        raise StructureImportError({"code": "invalid_file", "message": f"Unable to read the structure file: {exc}"}) from exc
    if errors:
        raise StructureImportError({"code": "invalid_structures", "message": "The file contains invalid molecular structures.", "records": errors})
    if not imported:
        raise StructureImportError({"code": "empty_input", "message": "No molecular structures were found."})
    return {"format": format, "count": len(imported), "molecules": imported}
