from flask import Blueprint, request, jsonify
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors, Descriptors
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Avalon import pyAvalonTools
import numpy as np
import logging
import sys
import os

# Add utils to path for file parsing
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.file_parsers import parse_molecule_file, validate_file_size, validate_molecule_count

similarity_bp = Blueprint('similarity', __name__)
logger = logging.getLogger(__name__)

# Add route aliases for other metrics that frontend might call
@similarity_bp.route('/dice', methods=['POST'])
def dice_similarity():
    """Calculate Dice similarity - route to tanimoto with dice metric"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')
        
        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400
            
        # Use the tanimoto function but override the metric
        request_data = {
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fp_type': fp_type
        }
        
        # Temporarily override request data
        original_data = request.get_json
        request.get_json = lambda: request_data
        
        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)
        
        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400
            
        # Generate fingerprints using same logic as tanimoto
        if fp_type == 'morgan':
            query_fp = GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
            target_fp = GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)
        elif fp_type == 'rdkit':
            query_fp = Chem.RDKFingerprint(query_mol, fpSize=2048)
            target_fp = Chem.RDKFingerprint(target_mol, fpSize=2048)
        elif fp_type == 'maccs':
            query_fp = rdMolDescriptors.GetMACCSKeysFingerprint(query_mol)
            target_fp = rdMolDescriptors.GetMACCSKeysFingerprint(target_mol)
        elif fp_type == 'avalon':
            query_fp = pyAvalonTools.GetAvalonFP(query_mol, nBits=2048)
            target_fp = pyAvalonTools.GetAvalonFP(target_mol, nBits=2048)
        elif fp_type == 'atom_pairs':
            query_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(query_mol, nBits=2048)
            target_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(target_mol, nBits=2048)
        elif fp_type == 'topological_torsions':
            query_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(query_mol, nBits=2048)
            target_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(target_mol, nBits=2048)
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400
            
        dice_sim = DataStructs.DiceSimilarity(query_fp, target_fp)
        
        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'dice_similarity': dice_sim,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating Dice similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/cosine', methods=['POST'])  
def cosine_similarity():
    """Calculate Cosine similarity - route to tanimoto with cosine metric"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')
        
        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400
            
        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)
        
        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400
            
        # Generate fingerprints
        if fp_type == 'morgan':
            query_fp = GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
            target_fp = GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)
        elif fp_type == 'rdkit':
            query_fp = Chem.RDKFingerprint(query_mol, fpSize=2048)
            target_fp = Chem.RDKFingerprint(target_mol, fpSize=2048)
        elif fp_type == 'maccs':
            query_fp = rdMolDescriptors.GetMACCSKeysFingerprint(query_mol)
            target_fp = rdMolDescriptors.GetMACCSKeysFingerprint(target_mol)
        elif fp_type == 'avalon':
            query_fp = pyAvalonTools.GetAvalonFP(query_mol, nBits=2048)
            target_fp = pyAvalonTools.GetAvalonFP(target_mol, nBits=2048)
        elif fp_type == 'atom_pairs':
            query_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(query_mol, nBits=2048)
            target_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(target_mol, nBits=2048)
        elif fp_type == 'topological_torsions':
            query_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(query_mol, nBits=2048)
            target_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(target_mol, nBits=2048)
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400
            
        cosine_sim = DataStructs.CosineSimilarity(query_fp, target_fp)
        
        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'cosine_similarity': cosine_sim,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating Cosine similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/sokal', methods=['POST'])
def sokal_similarity():
    """Calculate Sokal similarity - route to tanimoto with sokal metric"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')
        
        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400
            
        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)
        
        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400
            
        # Generate fingerprints
        if fp_type == 'morgan':
            query_fp = GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
            target_fp = GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)
        elif fp_type == 'rdkit':
            query_fp = Chem.RDKFingerprint(query_mol, fpSize=2048)
            target_fp = Chem.RDKFingerprint(target_mol, fpSize=2048)
        elif fp_type == 'maccs':
            query_fp = rdMolDescriptors.GetMACCSKeysFingerprint(query_mol)
            target_fp = rdMolDescriptors.GetMACCSKeysFingerprint(target_mol)
        elif fp_type == 'avalon':
            query_fp = pyAvalonTools.GetAvalonFP(query_mol, nBits=2048)
            target_fp = pyAvalonTools.GetAvalonFP(target_mol, nBits=2048)
        elif fp_type == 'atom_pairs':
            query_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(query_mol, nBits=2048)
            target_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(target_mol, nBits=2048)
        elif fp_type == 'topological_torsions':
            query_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(query_mol, nBits=2048)
            target_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(target_mol, nBits=2048)
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400
            
        sokal_sim = DataStructs.SokalSimilarity(query_fp, target_fp)
        
        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'sokal_similarity': sokal_sim,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating Sokal similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/tanimoto', methods=['POST'])
def tanimoto_similarity():
    """Calculate Tanimoto similarity between molecules"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')
        radius = data.get('radius', 2)
        n_bits = data.get('n_bits', 2048)
        
        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400
            
        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)
        
        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400
            
        # Generate fingerprints
        if fp_type == 'morgan':
            query_fp = GetMorganFingerprintAsBitVect(query_mol, radius, nBits=n_bits)
            target_fp = GetMorganFingerprintAsBitVect(target_mol, radius, nBits=n_bits)
        elif fp_type == 'rdkit':
            query_fp = Chem.RDKFingerprint(query_mol, fpSize=n_bits)
            target_fp = Chem.RDKFingerprint(target_mol, fpSize=n_bits)
        elif fp_type == 'maccs':
            query_fp = rdMolDescriptors.GetMACCSKeysFingerprint(query_mol)
            target_fp = rdMolDescriptors.GetMACCSKeysFingerprint(target_mol)
        elif fp_type == 'avalon':
            query_fp = pyAvalonTools.GetAvalonFP(query_mol, nBits=n_bits)
            target_fp = pyAvalonTools.GetAvalonFP(target_mol, nBits=n_bits)
        elif fp_type == 'atom_pairs':
            query_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(query_mol, nBits=n_bits)
            target_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(target_mol, nBits=n_bits)
        elif fp_type == 'topological_torsions':
            query_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(query_mol, nBits=n_bits)
            target_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(target_mol, nBits=n_bits)
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400
            
        tanimoto = DataStructs.TanimotoSimilarity(query_fp, target_fp)
        
        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'tanimoto_similarity': tanimoto,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating Tanimoto similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/bulk_similarity', methods=['POST'])
def bulk_similarity():
    """Calculate similarity between a query molecule and multiple targets"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles_list = data.get('target_smiles_list', [])
        fp_type = data.get('fp_type', 'morgan')
        radius = data.get('radius', 2)
        n_bits = data.get('n_bits', 2048)
        similarity_metric = data.get('similarity_metric', 'tanimoto')
        threshold = data.get('threshold', 0.0)
        
        if not query_smiles:
            return jsonify({'error': 'Query SMILES required', 'success': False}), 400
            
        if not target_smiles_list:
            return jsonify({'error': 'Target SMILES list required', 'success': False}), 400
            
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return jsonify({'error': 'Invalid query SMILES', 'success': False}), 400
            
        # Generate query fingerprint
        if fp_type == 'morgan':
            query_fp = GetMorganFingerprintAsBitVect(query_mol, radius, nBits=n_bits)
        elif fp_type == 'rdkit':
            query_fp = Chem.RDKFingerprint(query_mol, fpSize=n_bits)
        elif fp_type == 'maccs':
            query_fp = rdMolDescriptors.GetMACCSKeysFingerprint(query_mol)
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400
            
        results = []
        
        for i, target_smiles in enumerate(target_smiles_list):
            try:
                target_mol = Chem.MolFromSmiles(target_smiles)
                if target_mol is None:
                    results.append({
                        'index': i,
                        'target_smiles': target_smiles,
                        'similarity': None,
                        'error': 'Invalid SMILES'
                    })
                    continue
                    
                # Generate target fingerprint
                if fp_type == 'morgan':
                    target_fp = GetMorganFingerprintAsBitVect(target_mol, radius, nBits=n_bits)
                elif fp_type == 'rdkit':
                    target_fp = Chem.RDKFingerprint(target_mol, fpSize=n_bits)
                elif fp_type == 'maccs':
                    target_fp = rdMolDescriptors.GetMACCSKeysFingerprint(target_mol)
                    
                # Calculate similarity
                if similarity_metric == 'tanimoto':
                    similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
                elif similarity_metric == 'dice':
                    similarity = DataStructs.DiceSimilarity(query_fp, target_fp)
                elif similarity_metric == 'cosine':
                    similarity = DataStructs.CosineSimilarity(query_fp, target_fp)
                else:
                    similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
                    
                if similarity >= threshold:
                    results.append({
                        'index': i,
                        'target_smiles': target_smiles,
                        'similarity': similarity,
                        'error': None
                    })
                    
            except Exception as e:
                results.append({
                    'index': i,
                    'target_smiles': target_smiles,
                    'similarity': None,
                    'error': str(e)
                })
        
        # Sort by similarity descending
        results.sort(key=lambda x: x['similarity'] if x['similarity'] is not None else -1, reverse=True)
        
        return jsonify({
            'query_smiles': query_smiles,
            'fingerprint_type': fp_type,
            'similarity_metric': similarity_metric,
            'threshold': threshold,
            'num_targets': len(target_smiles_list),
            'num_hits': len([r for r in results if r['similarity'] is not None]),
            'results': results,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error in bulk similarity search: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/substructure_search', methods=['POST'])
@similarity_bp.route('/substructure', methods=['POST'])
def substructure_search():
    """Perform substructure search"""
    try:
        data = request.get_json()
        pattern_smiles = data.get('pattern_smiles') or data.get('pattern')
        target_smiles_list = data.get('target_smiles_list') or data.get('molecules', [])
        use_chirality = data.get('use_chirality', False)
        
        if not pattern_smiles:
            return jsonify({'error': 'Pattern SMILES required', 'success': False}), 400
            
        if not target_smiles_list:
            return jsonify({'error': 'Target SMILES list required', 'success': False}), 400
            
        pattern_mol = Chem.MolFromSmiles(pattern_smiles)
        if pattern_mol is None:
            return jsonify({'error': 'Invalid pattern SMILES', 'success': False}), 400
            
        matches = []
        
        for i, target_smiles in enumerate(target_smiles_list):
            try:
                target_mol = Chem.MolFromSmiles(target_smiles)
                if target_mol is None:
                    continue
                    
                # Perform substructure match
                has_match = target_mol.HasSubstructMatch(pattern_mol, useChirality=use_chirality)
                
                if has_match:
                    # Get match details
                    match_atoms = target_mol.GetSubstructMatch(pattern_mol)
                    matches.append({
                        'index': i,
                        'target_smiles': target_smiles,
                        'match_atoms': list(match_atoms) if match_atoms else []
                    })
                    
            except Exception as e:
                logger.warning(f"Error matching target {i}: {str(e)}")
                continue
        
        return jsonify({
            'pattern_smiles': pattern_smiles,
            'use_chirality': use_chirality,
            'num_targets': len(target_smiles_list),
            'num_matches': len(matches),
            'matches': matches,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error in substructure search: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/maximum_common_substructure', methods=['POST'])
def maximum_common_substructure():
    """Find maximum common substructure between molecules"""
    try:
        data = request.get_json()
        smiles1 = data.get('smiles1')
        smiles2 = data.get('smiles2')
        
        if not smiles1 or not smiles2:
            return jsonify({'error': 'Both SMILES strings required', 'success': False}), 400
            
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 is None or mol2 is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400
            
        # Find MCS
        from rdkit.Chem import rdFMCS
        
        mcs = rdFMCS.FindMCS([mol1, mol2])
        
        result = {
            'smiles1': smiles1,
            'smiles2': smiles2,
            'mcs_smarts': mcs.smartsString,
            'num_atoms': mcs.numAtoms,
            'num_bonds': mcs.numBonds,
            'canceled': mcs.canceled,
            'success': True
        }
        
        # Try to get MCS as molecule
        try:
            if mcs.smartsString:
                mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                if mcs_mol:
                    result['mcs_smiles'] = Chem.MolToSmiles(mcs_mol)
        except:
            pass
            
        return jsonify(result)
        
    except Exception as e:
        logger.error(f"Error finding MCS: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/diverse_subset', methods=['POST'])
def diverse_subset():
    """Select diverse subset from molecule list"""
    try:
        data = request.get_json()
        smiles_list = data.get('smiles_list', [])
        num_compounds = data.get('num_compounds', 10)
        fp_type = data.get('fp_type', 'morgan')
        radius = data.get('radius', 2)
        n_bits = data.get('n_bits', 2048)
        
        if not smiles_list:
            return jsonify({'error': 'SMILES list required', 'success': False}), 400
            
        if num_compounds >= len(smiles_list):
            return jsonify({'error': 'Number of compounds must be less than input size', 'success': False}), 400
            
        # Generate molecules and fingerprints
        mols = []
        fps = []
        valid_indices = []
        
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mols.append(mol)
                valid_indices.append(i)
                
                if fp_type == 'morgan':
                    fp = GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
                elif fp_type == 'rdkit':
                    fp = Chem.RDKFingerprint(mol, fpSize=n_bits)
                elif fp_type == 'maccs':
                    fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
                else:
                    return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400
                    
                fps.append(fp)
        
        if len(fps) < num_compounds:
            return jsonify({'error': 'Not enough valid molecules for selection', 'success': False}), 400
        
        # Perform diversity selection using MaxMin picker
        from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
        
        picker = MaxMinPicker()
        selected_indices = picker.LazyPick(fps, len(fps), num_compounds)
        
        selected_compounds = []
        for idx in selected_indices:
            original_idx = valid_indices[idx]
            selected_compounds.append({
                'original_index': original_idx,
                'smiles': smiles_list[original_idx]
            })
        
        return jsonify({
            'input_size': len(smiles_list),
            'valid_molecules': len(mols),
            'requested_compounds': num_compounds,
            'fingerprint_type': fp_type,
            'selected_compounds': selected_compounds,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error in diverse subset selection: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/bulk_similarity_file', methods=['POST'])
def bulk_similarity_file():
    """Calculate similarity between a query molecule and molecules from uploaded file"""
    try:
        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected', 'success': False}), 400

        # Get query SMILES from form data
        query_smiles = request.form.get('query_smiles')
        if not query_smiles:
            return jsonify({'error': 'Query SMILES required', 'success': False}), 400

        # Validate query SMILES
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return jsonify({'error': 'Invalid query SMILES', 'success': False}), 400

        # Get optional parameters
        fp_type = request.form.get('fp_type', 'morgan')
        radius = int(request.form.get('radius', 2))
        n_bits = int(request.form.get('n_bits', 2048))
        similarity_metric = request.form.get('similarity_metric', 'tanimoto')
        threshold = float(request.form.get('threshold', 0.0))

        # Read file content
        file_content = file.read().decode('utf-8')

        # Validate file size (50 MB max)
        validate_file_size(len(file_content.encode('utf-8')), max_size_mb=50)

        # Parse molecule file
        molecules = parse_molecule_file(file.filename, file_content)

        if not molecules:
            return jsonify({'error': 'No valid molecules found in file', 'success': False}), 400

        # Validate molecule count
        validate_molecule_count(len(molecules), max_molecules=10000)

        # Extract SMILES list
        target_smiles_list = [mol['smiles'] for mol in molecules]

        # Generate query fingerprint
        if fp_type == 'morgan':
            query_fp = GetMorganFingerprintAsBitVect(query_mol, radius, nBits=n_bits)
        elif fp_type == 'rdkit':
            query_fp = Chem.RDKFingerprint(query_mol, fpSize=n_bits)
        elif fp_type == 'maccs':
            query_fp = rdMolDescriptors.GetMACCSKeysFingerprint(query_mol)
        elif fp_type == 'avalon':
            query_fp = pyAvalonTools.GetAvalonFP(query_mol, nBits=n_bits)
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400

        results = []

        for i, mol_data in enumerate(molecules):
            try:
                target_smiles = mol_data['smiles']
                target_mol = Chem.MolFromSmiles(target_smiles)

                if target_mol is None:
                    results.append({
                        'index': i,
                        'name': mol_data.get('name', f'Molecule_{i+1}'),
                        'target_smiles': target_smiles,
                        'similarity': None,
                        'error': 'Invalid SMILES'
                    })
                    continue

                # Generate target fingerprint
                if fp_type == 'morgan':
                    target_fp = GetMorganFingerprintAsBitVect(target_mol, radius, nBits=n_bits)
                elif fp_type == 'rdkit':
                    target_fp = Chem.RDKFingerprint(target_mol, fpSize=n_bits)
                elif fp_type == 'maccs':
                    target_fp = rdMolDescriptors.GetMACCSKeysFingerprint(target_mol)
                elif fp_type == 'avalon':
                    target_fp = pyAvalonTools.GetAvalonFP(target_mol, nBits=n_bits)

                # Calculate similarity
                if similarity_metric == 'tanimoto':
                    similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
                elif similarity_metric == 'dice':
                    similarity = DataStructs.DiceSimilarity(query_fp, target_fp)
                elif similarity_metric == 'cosine':
                    similarity = DataStructs.CosineSimilarity(query_fp, target_fp)
                else:
                    similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)

                if similarity >= threshold:
                    results.append({
                        'index': i,
                        'name': mol_data.get('name', f'Molecule_{i+1}'),
                        'target_smiles': target_smiles,
                        'similarity': similarity,
                        'error': None
                    })

            except Exception as e:
                results.append({
                    'index': i,
                    'name': mol_data.get('name', f'Molecule_{i+1}'),
                    'target_smiles': mol_data['smiles'],
                    'similarity': None,
                    'error': str(e)
                })

        # Sort by similarity descending
        results.sort(key=lambda x: x['similarity'] if x['similarity'] is not None else -1, reverse=True)

        return jsonify({
            'query_smiles': query_smiles,
            'fingerprint_type': fp_type,
            'similarity_metric': similarity_metric,
            'threshold': threshold,
            'num_targets': len(molecules),
            'num_hits': len([r for r in results if r['similarity'] is not None]),
            'file_name': file.filename,
            'results': results,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error in bulk similarity file upload: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/substructure_search_file', methods=['POST'])
def substructure_search_file():
    """Perform substructure search on molecules from uploaded file"""
    try:
        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected', 'success': False}), 400

        # Get pattern SMILES from form data
        pattern_smiles = request.form.get('pattern_smiles') or request.form.get('pattern')
        if not pattern_smiles:
            return jsonify({'error': 'Pattern SMILES required', 'success': False}), 400

        # Validate pattern SMILES
        pattern_mol = Chem.MolFromSmiles(pattern_smiles)
        if pattern_mol is None:
            return jsonify({'error': 'Invalid pattern SMILES', 'success': False}), 400

        # Get optional parameters
        use_chirality = request.form.get('use_chirality', 'false').lower() == 'true'

        # Read file content
        file_content = file.read().decode('utf-8')

        # Validate file size
        validate_file_size(len(file_content.encode('utf-8')), max_size_mb=50)

        # Parse molecule file
        molecules = parse_molecule_file(file.filename, file_content)

        if not molecules:
            return jsonify({'error': 'No valid molecules found in file', 'success': False}), 400

        # Validate molecule count
        validate_molecule_count(len(molecules), max_molecules=10000)

        matches = []

        for i, mol_data in enumerate(molecules):
            try:
                target_smiles = mol_data['smiles']
                target_mol = Chem.MolFromSmiles(target_smiles)

                if target_mol is None:
                    continue

                # Perform substructure match
                has_match = target_mol.HasSubstructMatch(pattern_mol, useChirality=use_chirality)

                if has_match:
                    # Get match details
                    match_atoms = target_mol.GetSubstructMatch(pattern_mol)
                    matches.append({
                        'index': i,
                        'name': mol_data.get('name', f'Molecule_{i+1}'),
                        'target_smiles': target_smiles,
                        'match_atoms': list(match_atoms) if match_atoms else []
                    })

            except Exception as e:
                logger.warning(f"Error matching target {i}: {str(e)}")
                continue

        return jsonify({
            'pattern_smiles': pattern_smiles,
            'use_chirality': use_chirality,
            'num_targets': len(molecules),
            'num_matches': len(matches),
            'file_name': file.filename,
            'matches': matches,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error in substructure search file upload: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/diverse_subset_file', methods=['POST'])
def diverse_subset_file():
    """Select diverse subset from molecules in uploaded file"""
    try:
        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected', 'success': False}), 400

        # Get parameters from form data
        num_compounds = int(request.form.get('num_compounds', 10))
        fp_type = request.form.get('fp_type', 'morgan')
        radius = int(request.form.get('radius', 2))
        n_bits = int(request.form.get('n_bits', 2048))

        # Read file content
        file_content = file.read().decode('utf-8')

        # Validate file size
        validate_file_size(len(file_content.encode('utf-8')), max_size_mb=50)

        # Parse molecule file
        molecules = parse_molecule_file(file.filename, file_content)

        if not molecules:
            return jsonify({'error': 'No valid molecules found in file', 'success': False}), 400

        # Validate molecule count
        validate_molecule_count(len(molecules), max_molecules=10000)

        if num_compounds >= len(molecules):
            return jsonify({'error': 'Number of compounds must be less than input size', 'success': False}), 400

        # Generate molecules and fingerprints
        mols = []
        fps = []
        valid_indices = []

        for i, mol_data in enumerate(molecules):
            smiles = mol_data['smiles']
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mols.append(mol)
                valid_indices.append(i)

                if fp_type == 'morgan':
                    fp = GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
                elif fp_type == 'rdkit':
                    fp = Chem.RDKFingerprint(mol, fpSize=n_bits)
                elif fp_type == 'maccs':
                    fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
                else:
                    return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400

                fps.append(fp)

        if len(fps) < num_compounds:
            return jsonify({'error': 'Not enough valid molecules for selection', 'success': False}), 400

        # Perform diversity selection using MaxMin picker
        from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

        picker = MaxMinPicker()
        # Convert RDKit LazyGenerator to list for iteration
        selected_indices = []
        lazy_pick = picker.LazyPick(fps, len(fps), num_compounds)
        for idx in range(num_compounds):
            selected_indices.append(lazy_pick[idx])

        selected_compounds = []
        for idx in selected_indices:
            original_idx = valid_indices[idx]
            mol_data = molecules[original_idx]
            selected_compounds.append({
                'original_index': original_idx,
                'name': mol_data.get('name', f'Molecule_{original_idx+1}'),
                'smiles': mol_data['smiles']
            })

        return jsonify({
            'input_size': len(molecules),
            'valid_molecules': len(mols),
            'requested_compounds': num_compounds,
            'fingerprint_type': fp_type,
            'file_name': file.filename,
            'selected_compounds': selected_compounds,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error in diverse subset file upload: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

# Additional similarity metrics
def get_fingerprint_for_similarity(mol, fp_type='morgan', radius=2, n_bits=2048):
    """Helper function to generate fingerprints for similarity calculations"""
    if fp_type == 'morgan':
        return GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    elif fp_type == 'rdkit':
        return Chem.RDKFingerprint(mol, fpSize=n_bits)
    elif fp_type == 'maccs':
        from rdkit.Chem import MACCSkeys
        return MACCSkeys.GenMACCSKeys(mol)
    elif fp_type == 'avalon':
        return pyAvalonTools.GetAvalonFP(mol, nBits=n_bits)
    elif fp_type == 'atom_pairs':
        return rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits)
    elif fp_type == 'topological_torsions':
        return rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=n_bits)
    else:
        return GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

@similarity_bp.route('/braunblanquet', methods=['POST'])
def braunblanquet_similarity():
    """Calculate BraunBlanquet similarity coefficient"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.BraunBlanquetSimilarity(query_fp, target_fp)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'BraunBlanquet',
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating BraunBlanquet similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/kulczynski', methods=['POST'])
def kulczynski_similarity():
    """Calculate Kulczynski similarity coefficient"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.KulczynskiSimilarity(query_fp, target_fp)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'Kulczynski',
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating Kulczynski similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/mcconnaughey', methods=['POST'])
def mcconnaughey_similarity():
    """Calculate McConnaughey similarity coefficient"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.McConnaugheySimilarity(query_fp, target_fp)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'McConnaughey',
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating McConnaughey similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/rogotgoldberg', methods=['POST'])
def rogotgoldberg_similarity():
    """Calculate RogotGoldberg similarity coefficient"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.RogotGoldbergSimilarity(query_fp, target_fp)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'RogotGoldberg',
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating RogotGoldberg similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/russel', methods=['POST'])
def russel_similarity():
    """Calculate Russel similarity coefficient"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.RusselSimilarity(query_fp, target_fp)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'Russel',
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating Russel similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/tversky', methods=['POST'])
def tversky_similarity():
    """Calculate Tversky similarity coefficient with alpha and beta parameters"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')
        alpha = data.get('alpha', 0.5)
        beta = data.get('beta', 0.5)

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.TverskySimilarity(query_fp, target_fp, alpha, beta)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'Tversky',
            'alpha': alpha,
            'beta': beta,
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating Tversky similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/asymmetric', methods=['POST'])
def asymmetric_similarity():
    """Calculate Asymmetric similarity coefficient"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.AsymmetricSimilarity(query_fp, target_fp)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'Asymmetric',
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating Asymmetric similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/allbit', methods=['POST'])
def allbit_similarity():
    """Calculate AllBit similarity coefficient"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.AllBitSimilarity(query_fp, target_fp)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'AllBit',
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating AllBit similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@similarity_bp.route('/onbit', methods=['POST'])
def onbit_similarity():
    """Calculate OnBit similarity coefficient"""
    try:
        data = request.get_json()
        query_smiles = data.get('query_smiles') or data.get('smiles1')
        target_smiles = data.get('target_smiles') or data.get('smiles2')
        fp_type = data.get('fp_type') or data.get('fingerprint_type', 'morgan')

        if not query_smiles or not target_smiles:
            return jsonify({'error': 'Both query and target SMILES required', 'success': False}), 400

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        query_fp = get_fingerprint_for_similarity(query_mol, fp_type)
        target_fp = get_fingerprint_for_similarity(target_mol, fp_type)

        similarity = DataStructs.OnBitSimilarity(query_fp, target_fp)

        return jsonify({
            'query_smiles': query_smiles,
            'target_smiles': target_smiles,
            'fingerprint_type': fp_type,
            'similarity': similarity,
            'metric': 'OnBit',
            'success': True
        })
    except Exception as e:
        logger.error(f"Error calculating OnBit similarity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500