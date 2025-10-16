from flask import Blueprint, request, jsonify, send_file
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, GraphDescriptors, rdMolDescriptors
from rdkit.Chem.EState import EState_VSA, EState
import logging
import sys
import os
import csv
import io

# Add utils to path for file parsing
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.file_parsers import parse_molecule_file, validate_file_size, validate_molecule_count

descriptors_bp = Blueprint('descriptors', __name__)
logger = logging.getLogger(__name__)

@descriptors_bp.route('/basic', methods=['POST'])
def basic_descriptors():
    """Calculate basic molecular descriptors"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        descriptors = {
            'molecular_weight': Descriptors.MolWt(mol),
            'exact_molecular_weight': Descriptors.ExactMolWt(mol),
            'heavy_atom_count': Descriptors.HeavyAtomCount(mol),
            'num_heteroatoms': Descriptors.NumHeteroatoms(mol),
            'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'num_h_acceptors': Descriptors.NumHAcceptors(mol),
            'num_h_donors': Descriptors.NumHDonors(mol),
            'num_rings': Descriptors.RingCount(mol),
            'num_aromatic_rings': Descriptors.NumAromaticRings(mol),
            'num_saturated_rings': Descriptors.NumSaturatedRings(mol),
            'num_aliphatic_rings': Descriptors.NumAliphaticRings(mol),
            'balaban_j': Descriptors.BalabanJ(mol) if hasattr(Descriptors, 'BalabanJ') else None,
            'bertz_ct': Descriptors.BertzCT(mol) if hasattr(Descriptors, 'BertzCT') else None
        }
        
        return jsonify({
            'smiles': smiles,
            'descriptors': descriptors,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating basic descriptors: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@descriptors_bp.route('/lipinski', methods=['POST'])
def lipinski_descriptors():
    """Calculate Lipinski Rule of 5 descriptors"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        descriptors = {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'num_h_donors': Descriptors.NumHDonors(mol),
            'num_h_acceptors': Descriptors.NumHAcceptors(mol),
            'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol)
        }
        
        # Check Lipinski violations
        violations = 0
        violations += 1 if descriptors['molecular_weight'] > 500 else 0
        violations += 1 if descriptors['logp'] > 5 else 0
        violations += 1 if descriptors['num_h_donors'] > 5 else 0
        violations += 1 if descriptors['num_h_acceptors'] > 10 else 0
        
        descriptors['lipinski_violations'] = violations
        descriptors['passes_lipinski'] = violations <= 1
        
        return jsonify({
            'smiles': smiles,
            'descriptors': descriptors,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating Lipinski descriptors: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@descriptors_bp.route('/logp', methods=['POST'])
def logp_descriptors():
    """Calculate various LogP descriptors"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        descriptors = {
            'wildman_crippen_logp': Descriptors.MolLogP(mol),
            'molar_refractivity': Descriptors.MolMR(mol),
            'crippen_logp': Crippen.MolLogP(mol),
            'crippen_mr': Crippen.MolMR(mol)
        }
        
        return jsonify({
            'smiles': smiles,
            'descriptors': descriptors,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating LogP descriptors: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@descriptors_bp.route('/topological', methods=['POST'])
def topological_descriptors():
    """Calculate topological descriptors"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        descriptors = {
            'balaban_j': GraphDescriptors.BalabanJ(mol),
            'bertz_ct': GraphDescriptors.BertzCT(mol),
            'chi0': GraphDescriptors.Chi0(mol),
            'chi1': GraphDescriptors.Chi1(mol),
            'chi0n': GraphDescriptors.Chi0n(mol),
            'chi1n': GraphDescriptors.Chi1n(mol),
            'chi2n': GraphDescriptors.Chi2n(mol) if hasattr(GraphDescriptors, 'Chi2n') else None,
            'chi3n': GraphDescriptors.Chi3n(mol) if hasattr(GraphDescriptors, 'Chi3n') else None,
            'chi4n': GraphDescriptors.Chi4n(mol) if hasattr(GraphDescriptors, 'Chi4n') else None,
            'chi0v': GraphDescriptors.Chi0v(mol),
            'chi1v': GraphDescriptors.Chi1v(mol),
            'chi2v': GraphDescriptors.Chi2v(mol) if hasattr(GraphDescriptors, 'Chi2v') else None,
            'chi3v': GraphDescriptors.Chi3v(mol) if hasattr(GraphDescriptors, 'Chi3v') else None,
            'chi4v': GraphDescriptors.Chi4v(mol) if hasattr(GraphDescriptors, 'Chi4v') else None,
            'hallkieralpha': GraphDescriptors.HallKierAlpha(mol) if hasattr(GraphDescriptors, 'HallKierAlpha') else None,
            'kappa1': GraphDescriptors.Kappa1(mol) if hasattr(GraphDescriptors, 'Kappa1') else None,
            'kappa2': GraphDescriptors.Kappa2(mol) if hasattr(GraphDescriptors, 'Kappa2') else None,
            'kappa3': GraphDescriptors.Kappa3(mol) if hasattr(GraphDescriptors, 'Kappa3') else None
        }
        
        return jsonify({
            'smiles': smiles,
            'descriptors': descriptors,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating topological descriptors: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@descriptors_bp.route('/all', methods=['POST'])
def all_descriptors():
    """Calculate all available descriptors"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        descriptors = {}
        
        # Get all descriptor functions
        descriptor_functions = [
            (name, func) for name, func in Descriptors.descList
        ]
        
        for name, func in descriptor_functions:
            try:
                descriptors[name] = func(mol)
            except Exception as e:
                logger.warning(f"Failed to calculate descriptor {name}: {str(e)}")
                descriptors[name] = None
        
        return jsonify({
            'smiles': smiles,
            'descriptors': descriptors,
            'num_descriptors': len(descriptors),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating all descriptors: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@descriptors_bp.route('/list', methods=['GET'])
def list_descriptors():
    """List all available descriptors"""
    try:
        descriptor_names = [name for name, func in Descriptors.descList]
        
        return jsonify({
            'descriptors': descriptor_names,
            'count': len(descriptor_names),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error listing descriptors: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@descriptors_bp.route('/vsa', methods=['POST'])
def vsa_descriptors():
    """Calculate VSA (Van der Waals Surface Area) descriptors"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        descriptors = {
            'tpsa': Descriptors.TPSA(mol),
            'labute_asa': Descriptors.LabuteASA(mol) if hasattr(Descriptors, 'LabuteASA') else None
        }
        
        # Add available VSA descriptors dynamically
        vsa_descriptor_names = [name for name in dir(Descriptors) 
                               if any(vsa in name for vsa in ['VSA', 'SlogP', 'SMR', 'PEOE', 'EState'])]
        
        for desc_name in vsa_descriptor_names:
            try:
                desc_func = getattr(Descriptors, desc_name)
                if callable(desc_func):
                    descriptors[desc_name] = desc_func(mol)
            except Exception as e:
                logger.warning(f"Failed to calculate VSA descriptor {desc_name}: {str(e)}")
                descriptors[desc_name] = None
        
        return jsonify({
            'smiles': smiles,
            'descriptors': descriptors,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating VSA descriptors: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@descriptors_bp.route('/batch_file', methods=['POST'])
def batch_descriptors_file():
    """Calculate descriptors for multiple molecules from uploaded file"""
    try:
        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected', 'success': False}), 400

        # Get descriptor type
        descriptor_type = request.form.get('descriptor_type', 'basic')  # basic, lipinski, all, etc.
        output_format = request.form.get('output_format', 'json')  # json or csv

        # Read file content
        file_content = file.read().decode('utf-8')
        validate_file_size(len(file_content.encode('utf-8')), max_size_mb=50)

        # Parse molecule file
        molecules = parse_molecule_file(file.filename, file_content)

        if not molecules:
            return jsonify({'error': 'No valid molecules found in file', 'success': False}), 400

        validate_molecule_count(len(molecules), max_molecules=10000)

        # Calculate descriptors for each molecule
        results = []

        for i, mol_data in enumerate(molecules):
            try:
                smiles = mol_data['smiles']
                mol = Chem.MolFromSmiles(smiles)

                if mol is None:
                    results.append({
                        'index': i,
                        'name': mol_data.get('name', f'Molecule_{i+1}'),
                        'smiles': smiles,
                        'error': 'Invalid SMILES',
                        'descriptors': None
                    })
                    continue

                # Calculate descriptors based on type
                descriptors = {}

                if descriptor_type == 'basic' or descriptor_type == 'all':
                    descriptors.update({
                        'molecular_weight': Descriptors.MolWt(mol),
                        'exact_molecular_weight': Descriptors.ExactMolWt(mol),
                        'heavy_atom_count': Descriptors.HeavyAtomCount(mol),
                        'num_heteroatoms': Descriptors.NumHeteroatoms(mol),
                        'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                        'num_h_acceptors': Descriptors.NumHAcceptors(mol),
                        'num_h_donors': Descriptors.NumHDonors(mol),
                        'num_rings': Descriptors.RingCount(mol),
                        'num_aromatic_rings': Descriptors.NumAromaticRings(mol),
                        'num_saturated_rings': Descriptors.NumSaturatedRings(mol),
                        'num_aliphatic_rings': Descriptors.NumAliphaticRings(mol)
                    })

                if descriptor_type == 'lipinski' or descriptor_type == 'all':
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    hbd = Descriptors.NumHDonors(mol)
                    hba = Descriptors.NumHAcceptors(mol)

                    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

                    descriptors.update({
                        'logp': logp,
                        'lipinski_violations': violations,
                        'passes_lipinski': violations <= 1
                    })

                if descriptor_type == 'topological' or descriptor_type == 'all':
                    descriptors.update({
                        'balaban_j': GraphDescriptors.BalabanJ(mol),
                        'bertz_ct': GraphDescriptors.BertzCT(mol),
                        'chi0': GraphDescriptors.Chi0(mol),
                        'chi1': GraphDescriptors.Chi1(mol)
                    })

                if descriptor_type == 'all':
                    # Add all available descriptors
                    for name, func in Descriptors.descList:
                        try:
                            if name not in descriptors:
                                descriptors[name] = func(mol)
                        except Exception as e:
                            logger.warning(f"Failed to calculate descriptor {name}: {str(e)}")

                results.append({
                    'index': i,
                    'name': mol_data.get('name', f'Molecule_{i+1}'),
                    'smiles': smiles,
                    'error': None,
                    'descriptors': descriptors
                })

            except Exception as e:
                results.append({
                    'index': i,
                    'name': mol_data.get('name', f'Molecule_{i+1}'),
                    'smiles': mol_data['smiles'],
                    'error': str(e),
                    'descriptors': None
                })

        # Return results in requested format
        if output_format == 'csv':
            # Generate CSV
            output = io.StringIO()

            # Get all descriptor names
            all_descriptor_names = set()
            for result in results:
                if result['descriptors']:
                    all_descriptor_names.update(result['descriptors'].keys())

            descriptor_names = sorted(list(all_descriptor_names))

            # Write CSV
            fieldnames = ['index', 'name', 'smiles'] + descriptor_names + ['error']
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()

            for result in results:
                row = {
                    'index': result['index'],
                    'name': result['name'],
                    'smiles': result['smiles'],
                    'error': result['error'] or ''
                }

                if result['descriptors']:
                    for desc_name in descriptor_names:
                        row[desc_name] = result['descriptors'].get(desc_name, '')

                writer.writerow(row)

            # Create response
            output.seek(0)
            return send_file(
                io.BytesIO(output.getvalue().encode('utf-8')),
                mimetype='text/csv',
                as_attachment=True,
                download_name=f'descriptors_{file.filename.rsplit(".", 1)[0]}.csv'
            )

        else:  # JSON format
            return jsonify({
                'file_name': file.filename,
                'descriptor_type': descriptor_type,
                'num_molecules': len(molecules),
                'num_successful': len([r for r in results if r['error'] is None]),
                'results': results,
                'success': True
            })

    except Exception as e:
        logger.error(f"Error in batch descriptor calculation: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500