from flask import Blueprint, request, jsonify, send_file
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Avalon import pyAvalonTools
import numpy as np
import logging
import sys
import os
import csv
import io

# Add utils to path for file parsing
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.file_parsers import parse_molecule_file, validate_file_size, validate_molecule_count

fingerprints_bp = Blueprint('fingerprints', __name__)
logger = logging.getLogger(__name__)

def bitvect_to_list(bitvect):
    """Convert RDKit bit vector to list"""
    return [int(x) for x in bitvect.ToBitString()]

def bitvect_to_array(bitvect):
    """Convert RDKit bit vector to numpy array"""
    return np.array(bitvect_to_list(bitvect))

@fingerprints_bp.route('/morgan', methods=['POST'])
def morgan_fingerprint():
    """Generate Morgan (ECFP) fingerprints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        radius = data.get('radius', 2)
        n_bits = data.get('n_bits', 2048)
        use_features = data.get('use_features', False)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate Morgan fingerprint
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, radius, nBits=n_bits, useFeatures=use_features
        )
        
        fp_array = bitvect_to_list(fp)
        
        return jsonify({
            'smiles': smiles,
            'fingerprint': fp_array,
            'radius': radius,
            'n_bits': n_bits,
            'use_features': use_features,
            'on_bits': fp.GetNumOnBits(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating Morgan fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/rdkit', methods=['POST'])
def rdkit_fingerprint():
    """Generate RDKit topological fingerprints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        min_path = data.get('min_path', 1)
        max_path = data.get('max_path', 7)
        fp_size = data.get('fp_size', 2048)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate RDKit fingerprint
        fp = Chem.RDKFingerprint(mol, minPath=min_path, maxPath=max_path, fpSize=fp_size)
        fp_array = bitvect_to_list(fp)
        
        return jsonify({
            'smiles': smiles,
            'fingerprint': fp_array,
            'min_path': min_path,
            'max_path': max_path,
            'fp_size': fp_size,
            'on_bits': fp.GetNumOnBits(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating RDKit fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/maccs', methods=['POST'])
def maccs_fingerprint():
    """Generate MACCS keys fingerprints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate MACCS keys
        fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
        fp_array = bitvect_to_list(fp)
        
        return jsonify({
            'smiles': smiles,
            'fingerprint': fp_array,
            'fp_size': 167,  # MACCS keys are always 167 bits
            'on_bits': fp.GetNumOnBits(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating MACCS fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/avalon', methods=['POST'])
def avalon_fingerprint():
    """Generate Avalon fingerprints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        n_bits = data.get('n_bits', 512)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate Avalon fingerprint
        fp = pyAvalonTools.GetAvalonFP(mol, nBits=n_bits)
        fp_array = bitvect_to_list(fp)
        
        return jsonify({
            'smiles': smiles,
            'fingerprint': fp_array,
            'n_bits': n_bits,
            'on_bits': fp.GetNumOnBits(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating Avalon fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/atom_pairs', methods=['POST'])
def atom_pairs_fingerprint():
    """Generate Atom Pairs fingerprints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        n_bits = data.get('n_bits', 2048)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate Atom Pairs fingerprint
        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits)
        fp_array = bitvect_to_list(fp)
        
        return jsonify({
            'smiles': smiles,
            'fingerprint': fp_array,
            'n_bits': n_bits,
            'on_bits': fp.GetNumOnBits(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating Atom Pairs fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/topological_torsions', methods=['POST'])
def topological_torsions_fingerprint():
    """Generate Topological Torsions fingerprints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        n_bits = data.get('n_bits', 2048)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate Topological Torsions fingerprint
        fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=n_bits)
        fp_array = bitvect_to_list(fp)
        
        return jsonify({
            'smiles': smiles,
            'fingerprint': fp_array,
            'n_bits': n_bits,
            'on_bits': fp.GetNumOnBits(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating Topological Torsions fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/pattern', methods=['POST'])
def pattern_fingerprint():
    """Generate pattern-based fingerprints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        fp_size = data.get('fp_size', 2048)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate pattern fingerprint
        fp = Chem.PatternFingerprint(mol, fpSize=fp_size)
        fp_array = bitvect_to_list(fp)
        
        return jsonify({
            'smiles': smiles,
            'fingerprint': fp_array,
            'fp_size': fp_size,
            'on_bits': fp.GetNumOnBits(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating pattern fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/layered', methods=['POST'])
def layered_fingerprint():
    """Generate layered fingerprints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        fp_size = data.get('fp_size', 2048)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate layered fingerprint
        fp = Chem.LayeredFingerprint(mol, fpSize=fp_size)
        fp_array = bitvect_to_list(fp)
        
        return jsonify({
            'smiles': smiles,
            'fingerprint': fp_array,
            'fp_size': fp_size,
            'on_bits': fp.GetNumOnBits(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating layered fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/compare', methods=['POST'])
def compare_fingerprints():
    """Compare two fingerprints using various similarity metrics"""
    try:
        data = request.get_json()
        smiles1 = data.get('smiles1')
        smiles2 = data.get('smiles2')
        fp_type = data.get('fp_type', 'morgan')
        radius = data.get('radius', 2)
        n_bits = data.get('n_bits', 2048)
        
        if not smiles1 or not smiles2:
            return jsonify({'error': 'Both SMILES strings required', 'success': False}), 400
            
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 is None or mol2 is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400
            
        # Generate fingerprints based on type
        if fp_type == 'morgan':
            fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1, radius, nBits=n_bits)
            fp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2, radius, nBits=n_bits)
        elif fp_type == 'rdkit':
            fp1 = Chem.RDKFingerprint(mol1, fpSize=n_bits)
            fp2 = Chem.RDKFingerprint(mol2, fpSize=n_bits)
        elif fp_type == 'maccs':
            fp1 = rdMolDescriptors.GetMACCSKeysFingerprint(mol1)
            fp2 = rdMolDescriptors.GetMACCSKeysFingerprint(mol2)
        elif fp_type == 'avalon':
            fp1 = pyAvalonTools.GetAvalonFP(mol1, nBits=n_bits)
            fp2 = pyAvalonTools.GetAvalonFP(mol2, nBits=n_bits)
        elif fp_type == 'atom_pairs':
            fp1 = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol1, nBits=n_bits)
            fp2 = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol2, nBits=n_bits)
        elif fp_type == 'topological_torsions':
            fp1 = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol1, nBits=n_bits)
            fp2 = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol2, nBits=n_bits)
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400
            
        # Calculate similarity metrics
        tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)
        dice = DataStructs.DiceSimilarity(fp1, fp2)
        cosine = DataStructs.CosineSimilarity(fp1, fp2)
        sokal = DataStructs.SokalSimilarity(fp1, fp2)
        
        return jsonify({
            'smiles1': smiles1,
            'smiles2': smiles2,
            'fingerprint_type': fp_type,
            'similarities': {
                'tanimoto': tanimoto,
                'dice': dice,
                'cosine': cosine,
                'sokal': sokal
            },
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error comparing fingerprints: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@fingerprints_bp.route('/batch_file', methods=['POST'])
def batch_fingerprints_file():
    """Generate fingerprints for multiple molecules from uploaded file"""
    try:
        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected', 'success': False}), 400

        # Get fingerprint parameters
        fp_type = request.form.get('fp_type', 'morgan')
        radius = int(request.form.get('radius', 2))
        n_bits = int(request.form.get('n_bits', 2048))
        output_format = request.form.get('output_format', 'json')  # json, csv, or numpy

        # Read file content
        file_content = file.read().decode('utf-8')
        validate_file_size(len(file_content.encode('utf-8')), max_size_mb=50)

        # Parse molecule file
        molecules = parse_molecule_file(file.filename, file_content)

        if not molecules:
            return jsonify({'error': 'No valid molecules found in file', 'success': False}), 400

        validate_molecule_count(len(molecules), max_molecules=10000)

        # Generate fingerprints for each molecule
        results = []
        fingerprint_matrix = []

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
                        'fingerprint': None
                    })
                    if output_format == 'numpy':
                        fingerprint_matrix.append([0] * n_bits)
                    continue

                # Generate fingerprint based on type
                if fp_type == 'morgan':
                    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
                elif fp_type == 'rdkit':
                    fp = Chem.RDKFingerprint(mol, fpSize=n_bits)
                elif fp_type == 'maccs':
                    fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
                    n_bits = 167  # MACCS is always 167 bits
                elif fp_type == 'avalon':
                    fp = pyAvalonTools.GetAvalonFP(mol, nBits=n_bits)
                elif fp_type == 'atom_pairs':
                    fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits)
                elif fp_type == 'topological_torsions':
                    fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=n_bits)
                else:
                    return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400

                fp_array = bitvect_to_list(fp)
                fingerprint_matrix.append(fp_array)

                results.append({
                    'index': i,
                    'name': mol_data.get('name', f'Molecule_{i+1}'),
                    'smiles': smiles,
                    'error': None,
                    'fingerprint': fp_array if output_format == 'json' else None,
                    'on_bits': fp.GetNumOnBits()
                })

            except Exception as e:
                results.append({
                    'index': i,
                    'name': mol_data.get('name', f'Molecule_{i+1}'),
                    'smiles': mol_data['smiles'],
                    'error': str(e),
                    'fingerprint': None
                })
                if output_format == 'numpy':
                    fingerprint_matrix.append([0] * n_bits)

        # Return results in requested format
        if output_format == 'csv':
            # Generate CSV with fingerprint bits as columns
            output = io.StringIO()

            fieldnames = ['index', 'name', 'smiles'] + [f'bit_{i}' for i in range(n_bits)] + ['error']
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()

            for i, result in enumerate(results):
                row = {
                    'index': result['index'],
                    'name': result['name'],
                    'smiles': result['smiles'],
                    'error': result['error'] or ''
                }

                if result['fingerprint']:
                    for bit_idx, bit_val in enumerate(result['fingerprint']):
                        row[f'bit_{bit_idx}'] = bit_val
                else:
                    for bit_idx in range(n_bits):
                        row[f'bit_{bit_idx}'] = 0

                writer.writerow(row)

            output.seek(0)
            return send_file(
                io.BytesIO(output.getvalue().encode('utf-8')),
                mimetype='text/csv',
                as_attachment=True,
                download_name=f'fingerprints_{fp_type}_{file.filename.rsplit(".", 1)[0]}.csv'
            )

        elif output_format == 'numpy':
            # Return numpy array as binary
            fp_array = np.array(fingerprint_matrix, dtype=np.int8)

            output = io.BytesIO()
            np.save(output, fp_array)
            output.seek(0)

            return send_file(
                output,
                mimetype='application/octet-stream',
                as_attachment=True,
                download_name=f'fingerprints_{fp_type}_{file.filename.rsplit(".", 1)[0]}.npy'
            )

        else:  # JSON format
            return jsonify({
                'file_name': file.filename,
                'fingerprint_type': fp_type,
                'radius': radius if fp_type == 'morgan' else None,
                'n_bits': n_bits,
                'num_molecules': len(molecules),
                'num_successful': len([r for r in results if r['error'] is None]),
                'results': results,
                'success': True
            })

    except Exception as e:
        logger.error(f"Error in batch fingerprint generation: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500