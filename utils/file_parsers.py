"""
File parsing utilities for chemistry file formats
Supports: .smi, .smiles, .sdf, .mol, .csv, .tsv, .txt, .rxn
"""

from rdkit import Chem
import logging
import csv
import io
from typing import List, Dict, Tuple, Optional

logger = logging.getLogger(__name__)


def parse_smiles_file(file_content: str) -> List[Dict[str, str]]:
    """
    Parse SMILES file format (.smi, .smiles, .txt)
    Format: SMILES [ID/Name] (tab or space separated)

    Returns list of dicts with 'smiles' and optional 'name'
    """
    molecules = []
    lines = file_content.strip().split('\n')

    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.startswith('#'):  # Skip empty lines and comments
            continue

        # Try tab-separated first, then space-separated
        parts = line.split('\t') if '\t' in line else line.split(None, 1)

        smiles = parts[0].strip()
        name = parts[1].strip() if len(parts) > 1 else f"Molecule_{line_num}"

        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES on line {line_num}: {smiles}")
            continue

        molecules.append({
            'smiles': smiles,
            'name': name,
            'line_number': line_num
        })

    return molecules


def parse_sdf_file(file_content: str) -> List[Dict[str, str]]:
    """
    Parse SDF file format (.sdf)

    Returns list of dicts with 'smiles', 'name', and optional properties
    """
    molecules = []

    # Use RDKit's SDMolSupplier with a string stream
    supplier = Chem.SDMolSupplier()
    supplier.SetData(file_content)

    for idx, mol in enumerate(supplier):
        if mol is None:
            logger.warning(f"Could not parse molecule {idx + 1} in SDF file")
            continue

        smiles = Chem.MolToSmiles(mol)

        # Try to get molecule name from properties
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Molecule_{idx + 1}"

        # Extract all properties
        properties = {}
        for prop_name in mol.GetPropNames():
            if not prop_name.startswith('_'):  # Skip internal properties
                properties[prop_name] = mol.GetProp(prop_name)

        molecules.append({
            'smiles': smiles,
            'name': name,
            'mol_index': idx + 1,
            'properties': properties
        })

    return molecules


def parse_csv_file(file_content: str, smiles_column: Optional[str] = None,
                   name_column: Optional[str] = None, delimiter: str = ',') -> List[Dict[str, str]]:
    """
    Parse CSV/TSV file format (.csv, .tsv)

    Args:
        file_content: CSV file content as string
        smiles_column: Name of column containing SMILES (auto-detect if None)
        name_column: Name of column containing names (auto-detect if None)
        delimiter: Delimiter character (',' for CSV, '\t' for TSV)

    Returns list of dicts with 'smiles', 'name', and other columns
    """
    molecules = []

    # Parse CSV
    csv_file = io.StringIO(file_content)
    reader = csv.DictReader(csv_file, delimiter=delimiter)

    if not reader.fieldnames:
        raise ValueError("CSV file has no header row")

    # Auto-detect SMILES column
    if smiles_column is None:
        possible_names = ['smiles', 'SMILES', 'Smiles', 'smile', 'SMILE', 'canonical_smiles', 'structure']
        for name in possible_names:
            if name in reader.fieldnames:
                smiles_column = name
                break

        if smiles_column is None:
            # Use first column as default
            smiles_column = reader.fieldnames[0]
            logger.warning(f"Could not auto-detect SMILES column, using first column: {smiles_column}")

    # Auto-detect name column
    if name_column is None:
        possible_names = ['name', 'Name', 'NAME', 'id', 'ID', 'Id', 'compound_id', 'molecule_name']
        for name in possible_names:
            if name in reader.fieldnames:
                name_column = name
                break

    # Parse rows
    for row_num, row in enumerate(reader, 1):
        if smiles_column not in row:
            logger.warning(f"Row {row_num} missing SMILES column '{smiles_column}'")
            continue

        smiles = row[smiles_column].strip()
        if not smiles:
            continue

        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES on row {row_num}: {smiles}")
            continue

        name = row.get(name_column, f"Molecule_{row_num}") if name_column else f"Molecule_{row_num}"

        # Include all other columns as properties
        properties = {k: v for k, v in row.items() if k not in [smiles_column, name_column]}

        molecules.append({
            'smiles': smiles,
            'name': name,
            'row_number': row_num,
            'properties': properties
        })

    return molecules


def parse_mol_file(file_content: str) -> Dict[str, str]:
    """
    Parse single MOL file format (.mol)

    Returns dict with 'smiles' and 'name'
    """
    mol = Chem.MolFromMolBlock(file_content)

    if mol is None:
        raise ValueError("Invalid MOL file format")

    smiles = Chem.MolToSmiles(mol)
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Molecule_1"

    # Extract properties
    properties = {}
    for prop_name in mol.GetPropNames():
        if not prop_name.startswith('_'):
            properties[prop_name] = mol.GetProp(prop_name)

    return {
        'smiles': smiles,
        'name': name,
        'properties': properties
    }


def parse_rxn_file(file_content: str) -> Dict[str, any]:
    """
    Parse RXN file format (.rxn) - MDL Reaction format

    Returns dict with 'reactants', 'products', and reaction metadata
    """
    try:
        rxn = Chem.ReactionFromRxnBlock(file_content)

        if rxn is None:
            raise ValueError("Invalid RXN file format")

        reactants = []
        for i in range(rxn.GetNumReactantTemplates()):
            reactant = rxn.GetReactantTemplate(i)
            if reactant:
                smiles = Chem.MolToSmiles(reactant)
                reactants.append(smiles)

        products = []
        for i in range(rxn.GetNumProductTemplates()):
            product = rxn.GetProductTemplate(i)
            if product:
                smiles = Chem.MolToSmiles(product)
                products.append(smiles)

        # Try to construct reaction SMARTS
        reaction_smarts = '.'.join(reactants) + '>>' + '.'.join(products)

        return {
            'reactants': reactants,
            'products': products,
            'reaction_smarts': reaction_smarts,
            'num_reactants': len(reactants),
            'num_products': len(products)
        }

    except Exception as e:
        raise ValueError(f"Error parsing RXN file: {str(e)}")


def detect_file_format(filename: str, file_content: str) -> str:
    """
    Detect file format from filename and content

    Returns format string: 'smiles', 'sdf', 'csv', 'tsv', 'mol', 'rxn', 'txt'
    """
    filename_lower = filename.lower()

    # Check extension
    if filename_lower.endswith('.smi') or filename_lower.endswith('.smiles'):
        return 'smiles'
    elif filename_lower.endswith('.sdf'):
        return 'sdf'
    elif filename_lower.endswith('.mol'):
        return 'mol'
    elif filename_lower.endswith('.csv'):
        return 'csv'
    elif filename_lower.endswith('.tsv') or filename_lower.endswith('.tab'):
        return 'tsv'
    elif filename_lower.endswith('.rxn'):
        return 'rxn'
    elif filename_lower.endswith('.txt'):
        # Try to detect from content
        if file_content.startswith('$$$$') or 'M  END' in file_content:
            return 'sdf'
        elif '\t' in file_content:
            return 'tsv'
        else:
            return 'smiles'

    # Fallback: try to detect from content
    if file_content.startswith('$$$$') or 'M  END' in file_content:
        return 'sdf'
    elif '$RXN' in file_content[:100]:
        return 'rxn'
    elif ',' in file_content and '\n' in file_content:
        return 'csv'
    elif '\t' in file_content:
        return 'tsv'
    else:
        return 'smiles'


def parse_molecule_file(filename: str, file_content: str,
                       file_format: Optional[str] = None) -> List[Dict[str, str]]:
    """
    Universal molecule file parser

    Args:
        filename: Original filename
        file_content: File content as string
        file_format: Explicit format (auto-detect if None)

    Returns list of molecule dicts with 'smiles', 'name', and optional properties
    """
    if file_format is None:
        file_format = detect_file_format(filename, file_content)

    logger.info(f"Parsing file '{filename}' as format: {file_format}")

    try:
        if file_format == 'smiles' or file_format == 'txt':
            molecules = parse_smiles_file(file_content)
        elif file_format == 'sdf':
            molecules = parse_sdf_file(file_content)
        elif file_format == 'csv':
            molecules = parse_csv_file(file_content, delimiter=',')
        elif file_format == 'tsv':
            molecules = parse_csv_file(file_content, delimiter='\t')
        elif file_format == 'mol':
            mol_data = parse_mol_file(file_content)
            molecules = [mol_data]  # Return as single-item list
        else:
            raise ValueError(f"Unsupported file format: {file_format}")

        logger.info(f"Successfully parsed {len(molecules)} molecules from '{filename}'")
        return molecules

    except Exception as e:
        logger.error(f"Error parsing file '{filename}': {str(e)}")
        raise ValueError(f"Failed to parse file: {str(e)}")


def validate_file_size(file_size_bytes: int, max_size_mb: int = 50) -> bool:
    """
    Validate file size

    Args:
        file_size_bytes: File size in bytes
        max_size_mb: Maximum allowed size in megabytes

    Returns True if valid, raises ValueError if too large
    """
    max_size_bytes = max_size_mb * 1024 * 1024

    if file_size_bytes > max_size_bytes:
        raise ValueError(f"File size ({file_size_bytes / (1024*1024):.1f} MB) exceeds maximum allowed size ({max_size_mb} MB)")

    return True


def validate_molecule_count(count: int, max_molecules: int = 10000) -> bool:
    """
    Validate number of molecules

    Args:
        count: Number of molecules
        max_molecules: Maximum allowed molecules

    Returns True if valid, raises ValueError if too many
    """
    if count > max_molecules:
        raise ValueError(f"Number of molecules ({count}) exceeds maximum allowed ({max_molecules})")

    return True
