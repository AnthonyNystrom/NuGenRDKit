from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdchem, AllChem, Descriptors, rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
import logging

molecular_structure_bp = Blueprint('molecular_structure', __name__)
logger = logging.getLogger(__name__)

@molecular_structure_bp.route('/smiles_to_mol', methods=['POST'])
@molecular_structure_bp.route('/mol', methods=['POST'])
def smiles_to_mol():
    """Convert SMILES string to MOL format"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        mol_block = Chem.MolToMolBlock(mol)
        
        return jsonify({
            'mol_block': mol_block,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error converting SMILES to MOL: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/mol_to_smiles', methods=['POST'])
def mol_to_smiles():
    """Convert MOL format to SMILES string"""
    try:
        data = request.get_json()
        mol_block = data.get('mol_block')
        
        if not mol_block:
            return jsonify({'error': 'MOL block required', 'success': False}), 400
            
        mol = Chem.MolFromMolBlock(mol_block)
        if mol is None:
            return jsonify({'error': 'Invalid MOL block', 'success': False}), 400
            
        smiles = Chem.MolToSmiles(mol)
        canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        
        return jsonify({
            'smiles': smiles,
            'canonical_smiles': canonical_smiles,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error converting MOL to SMILES: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/inchi_to_smiles', methods=['POST'])
def inchi_to_smiles():
    """Convert InChI to SMILES"""
    try:
        data = request.get_json()
        inchi = data.get('inchi')
        
        if not inchi:
            return jsonify({'error': 'InChI string required', 'success': False}), 400
            
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            return jsonify({'error': 'Invalid InChI string', 'success': False}), 400
            
        smiles = Chem.MolToSmiles(mol)
        
        return jsonify({
            'smiles': smiles,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error converting InChI to SMILES: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/smiles_to_inchi', methods=['POST'])
@molecular_structure_bp.route('/inchi', methods=['POST'])
def smiles_to_inchi():
    """Convert SMILES to InChI"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        inchi = Chem.MolToInchi(mol)
        inchi_key = Chem.MolToInchiKey(mol)
        
        return jsonify({
            'inchi': inchi,
            'inchi_key': inchi_key,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error converting SMILES to InChI: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/canonicalize', methods=['POST'])
def canonicalize_smiles():
    """Canonicalize SMILES string"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        canonical_smiles = Chem.MolToSmiles(mol)
        
        return jsonify({
            'original_smiles': smiles,
            'canonical_smiles': canonical_smiles,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error canonicalizing SMILES: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/validate', methods=['POST'])
def validate_structure():
    """Validate molecular structure"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        is_valid = mol is not None
        
        result = {
            'smiles': smiles,
            'is_valid': is_valid,
            'success': True
        }
        
        if is_valid:
            result.update({
                'canonical_smiles': Chem.MolToSmiles(mol),
                'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
                'num_atoms': mol.GetNumAtoms(),
                'num_bonds': mol.GetNumBonds()
            })
        
        return jsonify(result)
        
    except Exception as e:
        logger.error(f"Error validating structure: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/standardize', methods=['POST'])
def standardize_molecule():
    """Standardize molecular structure"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Standardization steps
        Chem.SanitizeMol(mol)
        standardized_smiles = Chem.MolToSmiles(mol)
        
        return jsonify({
            'original_smiles': smiles,
            'standardized_smiles': standardized_smiles,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error standardizing molecule: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/batch_convert', methods=['POST'])
def batch_convert_structures():
    """Convert structures from an uploaded file"""
    try:
        from utils.file_parsers import parse_molecule_file

        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'Empty filename', 'success': False}), 400

        output_format = request.form.get('output_format', 'canonicalize')

        # Parse molecules from file
        file_content = file.read().decode('utf-8')
        molecules = parse_molecule_file(file.filename, file_content)

        if not molecules:
            return jsonify({'error': 'No valid molecules found in file', 'success': False}), 400

        results = []
        num_successful = 0

        for mol_data in molecules:
            try:
                smiles = mol_data['smiles']
                mol = Chem.MolFromSmiles(smiles)

                if mol is None:
                    results.append({
                        'name': mol_data['name'],
                        'input_smiles': smiles,
                        'error': 'Invalid SMILES',
                        'output': None
                    })
                    continue

                # Convert based on format
                output_data = {'name': mol_data['name'], 'input_smiles': smiles}

                if output_format == 'canonicalize':
                    output_data['canonical_smiles'] = Chem.MolToSmiles(mol)
                    output_data['output'] = output_data['canonical_smiles']

                elif output_format == 'inchi':
                    inchi = Chem.MolToInchi(mol)
                    inchi_key = Chem.MolToInchiKey(mol)
                    output_data['inchi'] = inchi
                    output_data['inchi_key'] = inchi_key
                    output_data['output'] = inchi

                elif output_format == 'mol':
                    mol_block = Chem.MolToMolBlock(mol)
                    output_data['mol_block'] = mol_block
                    output_data['output'] = mol_block

                elif output_format == 'validate':
                    output_data['is_valid'] = True
                    output_data['canonical_smiles'] = Chem.MolToSmiles(mol)
                    output_data['output'] = 'Valid'

                output_data['error'] = None
                results.append(output_data)
                num_successful += 1

            except Exception as e:
                logger.warning(f"Error converting molecule {mol_data.get('name', 'unknown')}: {str(e)}")
                results.append({
                    'name': mol_data.get('name', 'unknown'),
                    'input_smiles': mol_data.get('smiles', ''),
                    'error': str(e),
                    'output': None
                })

        return jsonify({
            'success': True,
            'file_name': file.filename,
            'output_format': output_format,
            'num_molecules': len(molecules),
            'num_successful': num_successful,
            'results': results
        })

    except Exception as e:
        logger.error(f"Error in batch convert: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

# ===========================
# HYDROGEN OPERATIONS
# ===========================

@molecular_structure_bp.route('/add_hydrogens', methods=['POST'])
def add_hydrogens():
    """Add explicit hydrogens to molecule"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        add_coords = data.get('add_coords', False)

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Add hydrogens
        mol_with_h = Chem.AddHs(mol, addCoords=add_coords)
        smiles_with_h = Chem.MolToSmiles(mol_with_h)

        return jsonify({
            'original_smiles': smiles,
            'smiles_with_hydrogens': smiles_with_h,
            'num_atoms_original': mol.GetNumAtoms(),
            'num_atoms_with_h': mol_with_h.GetNumAtoms(),
            'num_hydrogens_added': mol_with_h.GetNumAtoms() - mol.GetNumAtoms(),
            'success': True
        })

    except Exception as e:
        logger.error(f"Error adding hydrogens: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/remove_hydrogens', methods=['POST'])
def remove_hydrogens():
    """Remove explicit hydrogens from molecule"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Remove hydrogens
        mol_no_h = Chem.RemoveHs(mol)
        smiles_no_h = Chem.MolToSmiles(mol_no_h)

        return jsonify({
            'original_smiles': smiles,
            'smiles_without_hydrogens': smiles_no_h,
            'num_atoms_original': mol.GetNumAtoms(),
            'num_atoms_without_h': mol_no_h.GetNumAtoms(),
            'num_hydrogens_removed': mol.GetNumAtoms() - mol_no_h.GetNumAtoms(),
            'success': True
        })

    except Exception as e:
        logger.error(f"Error removing hydrogens: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

# ===========================
# STEREOCHEMISTRY OPERATIONS
# ===========================

@molecular_structure_bp.route('/stereochemistry', methods=['POST'])
def get_stereochemistry():
    """Get stereochemistry information"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Assign stereochemistry
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Find chiral centers
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        # Get stereo info
        stereo_info = {
            'chiral_centers': [
                {
                    'atom_idx': idx,
                    'chirality': chirality
                } for idx, chirality in chiral_centers
            ],
            'num_chiral_centers': len(chiral_centers),
            'has_stereochemistry': len(chiral_centers) > 0,
            'smiles_with_stereo': Chem.MolToSmiles(mol, isomericSmiles=True),
            'smiles_without_stereo': Chem.MolToSmiles(mol, isomericSmiles=False)
        }

        return jsonify({
            'smiles': smiles,
            'stereochemistry': stereo_info,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error analyzing stereochemistry: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/enumerate_stereoisomers', methods=['POST'])
def enumerate_stereoisomers():
    """Enumerate all possible stereoisomers"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        max_isomers = data.get('max_isomers', 32)

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Enumerate stereoisomers
        from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

        opts = StereoEnumerationOptions(maxIsomers=max_isomers)
        isomers = list(EnumerateStereoisomers(mol, options=opts))

        stereoisomers = [
            {
                'smiles': Chem.MolToSmiles(isomer, isomericSmiles=True),
                'index': idx
            } for idx, isomer in enumerate(isomers)
        ]

        return jsonify({
            'original_smiles': smiles,
            'num_stereoisomers': len(stereoisomers),
            'stereoisomers': stereoisomers,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error enumerating stereoisomers: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

# ===========================
# KEKULIZATION & AROMATICITY
# ===========================

@molecular_structure_bp.route('/kekulize', methods=['POST'])
def kekulize_molecule():
    """Convert aromatic representation to Kekule form"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Kekulize
        Chem.Kekulize(mol, clearAromaticFlags=True)
        kekule_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)

        return jsonify({
            'original_smiles': smiles,
            'kekule_smiles': kekule_smiles,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error kekulizing molecule: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/aromaticity', methods=['POST'])
def set_aromaticity():
    """Set aromaticity for molecule"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Set aromaticity
        Chem.SetAromaticity(mol)
        aromatic_smiles = Chem.MolToSmiles(mol)

        # Count aromatic atoms and bonds
        num_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        num_aromatic_bonds = sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic())

        return jsonify({
            'original_smiles': smiles,
            'aromatic_smiles': aromatic_smiles,
            'num_aromatic_atoms': num_aromatic_atoms,
            'num_aromatic_bonds': num_aromatic_bonds,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error setting aromaticity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

# ===========================
# FRAGMENT OPERATIONS
# ===========================

@molecular_structure_bp.route('/fragment', methods=['POST'])
def fragment_molecule():
    """Fragment molecule into disconnected pieces"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Get fragments
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)

        fragments = [
            {
                'smiles': Chem.MolToSmiles(frag),
                'num_atoms': frag.GetNumAtoms(),
                'molecular_weight': Descriptors.MolWt(frag),
                'index': idx
            } for idx, frag in enumerate(frags)
        ]

        return jsonify({
            'original_smiles': smiles,
            'num_fragments': len(fragments),
            'fragments': fragments,
            'is_single_fragment': len(fragments) == 1,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error fragmenting molecule: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/murcko_scaffold', methods=['POST'])
def get_murcko_scaffold():
    """Get Murcko scaffold (core structure)"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        generic = data.get('generic', False)

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Get scaffold
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold)

        result = {
            'original_smiles': smiles,
            'scaffold_smiles': scaffold_smiles,
            'scaffold_num_atoms': scaffold.GetNumAtoms(),
            'success': True
        }

        # Get generic scaffold if requested
        if generic:
            generic_scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
            result['generic_scaffold_smiles'] = Chem.MolToSmiles(generic_scaffold)

        return jsonify(result)

    except Exception as e:
        logger.error(f"Error getting Murcko scaffold: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

# ===========================
# RING ANALYSIS
# ===========================

@molecular_structure_bp.route('/ring_info', methods=['POST'])
def get_ring_info():
    """Analyze ring systems in molecule"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Get ring info
        ring_info = mol.GetRingInfo()

        # Atom rings
        atom_rings = ring_info.AtomRings()

        # Bond rings
        bond_rings = ring_info.BondRings()

        # Analyze ring sizes
        ring_sizes = [len(ring) for ring in atom_rings]

        result = {
            'smiles': smiles,
            'num_rings': ring_info.NumRings(),
            'ring_sizes': ring_sizes,
            'atom_rings': [list(ring) for ring in atom_rings],
            'num_aromatic_rings': Descriptors.NumAromaticRings(mol),
            'num_aliphatic_rings': Descriptors.NumAliphaticRings(mol),
            'num_saturated_rings': Descriptors.NumSaturatedRings(mol),
            'success': True
        }

        return jsonify(result)

    except Exception as e:
        logger.error(f"Error analyzing rings: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

# ===========================
# ADDITIONAL FORMAT SUPPORT
# ===========================

@molecular_structure_bp.route('/sdf_to_smiles', methods=['POST'])
def sdf_to_smiles():
    """Convert SDF to SMILES"""
    try:
        data = request.get_json()
        sdf_block = data.get('sdf_block')

        if not sdf_block:
            return jsonify({'error': 'SDF block required', 'success': False}), 400

        mol = Chem.MolFromMolBlock(sdf_block)
        if mol is None:
            return jsonify({'error': 'Invalid SDF block', 'success': False}), 400

        smiles = Chem.MolToSmiles(mol)

        return jsonify({
            'smiles': smiles,
            'canonical_smiles': Chem.MolToSmiles(Chem.MolFromSmiles(smiles)),
            'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
            'success': True
        })

    except Exception as e:
        logger.error(f"Error converting SDF: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/smiles_to_sdf', methods=['POST'])
def smiles_to_sdf():
    """Convert SMILES to SDF"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Add 2D coordinates for SDF
        AllChem.Compute2DCoords(mol)

        sdf_block = Chem.MolToMolBlock(mol)

        return jsonify({
            'sdf_block': sdf_block,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error converting to SDF: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/pdb_to_smiles', methods=['POST'])
def pdb_to_smiles():
    """Convert PDB to SMILES"""
    try:
        data = request.get_json()
        pdb_block = data.get('pdb_block')

        if not pdb_block:
            return jsonify({'error': 'PDB block required', 'success': False}), 400

        mol = Chem.MolFromPDBBlock(pdb_block)
        if mol is None:
            return jsonify({'error': 'Invalid PDB block', 'success': False}), 400

        smiles = Chem.MolToSmiles(mol)

        return jsonify({
            'smiles': smiles,
            'canonical_smiles': Chem.MolToSmiles(Chem.MolFromSmiles(smiles)),
            'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
            'success': True
        })

    except Exception as e:
        logger.error(f"Error converting PDB: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/smiles_to_pdb', methods=['POST'])
def smiles_to_pdb():
    """Convert SMILES to PDB"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Add 3D coordinates for PDB
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        pdb_block = Chem.MolToPDBBlock(mol)

        return jsonify({
            'pdb_block': pdb_block,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error converting to PDB: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

# ===========================
# ADVANCED CLEANUP OPERATIONS
# ===========================

@molecular_structure_bp.route('/neutralize', methods=['POST'])
def neutralize_charges():
    """Neutralize charges in molecule"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Neutralize
        uncharger = rdMolStandardize.Uncharger()
        neutral_mol = uncharger.uncharge(mol)
        neutral_smiles = Chem.MolToSmiles(neutral_mol)

        return jsonify({
            'original_smiles': smiles,
            'neutralized_smiles': neutral_smiles,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error neutralizing molecule: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/remove_fragments', methods=['POST'])
def remove_fragments():
    """Remove salt/solvent fragments and keep largest fragment"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Remove fragments (keep largest)
        chooser = rdMolStandardize.LargestFragmentChooser()
        cleaned_mol = chooser.choose(mol)
        cleaned_smiles = Chem.MolToSmiles(cleaned_mol)

        return jsonify({
            'original_smiles': smiles,
            'cleaned_smiles': cleaned_smiles,
            'original_num_atoms': mol.GetNumAtoms(),
            'cleaned_num_atoms': cleaned_mol.GetNumAtoms(),
            'success': True
        })

    except Exception as e:
        logger.error(f"Error removing fragments: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@molecular_structure_bp.route('/cleanup', methods=['POST'])
def cleanup_molecule():
    """Full molecular cleanup (standardize, neutralize, remove fragments)"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Full cleanup pipeline
        # 1. Remove fragments
        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)

        # 2. Neutralize charges
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)

        # 3. Standardize
        Chem.SanitizeMol(mol)

        cleaned_smiles = Chem.MolToSmiles(mol)

        return jsonify({
            'original_smiles': smiles,
            'cleaned_smiles': cleaned_smiles,
            'canonical_smiles': Chem.MolToSmiles(Chem.MolFromSmiles(cleaned_smiles)),
            'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
            'success': True
        })

    except Exception as e:
        logger.error(f"Error cleaning up molecule: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500