from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdMolAlign
import logging
import json
import sys
import os

# Add utils to path for file parsing
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.file_parsers import parse_molecule_file, parse_sdf_file, validate_file_size

coordinates_bp = Blueprint('coordinates', __name__)
logger = logging.getLogger(__name__)

def mol_to_coordinates(mol, conf_id=0):
    """Extract coordinates from molecule conformer"""
    if mol.GetNumConformers() == 0:
        return None
        
    conf = mol.GetConformer(conf_id)
    coordinates = []
    
    for atom_idx in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        coordinates.append({
            'atom_idx': atom_idx,
            'symbol': atom.GetSymbol(),
            'x': pos.x,
            'y': pos.y,
            'z': pos.z
        })
        
    return coordinates

@coordinates_bp.route('/generate_2d', methods=['POST'])
@coordinates_bp.route('/2d', methods=['POST'])
def generate_2d_coordinates():
    """Generate 2D coordinates for a molecule"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate 2D coordinates
        rdDepictor.Compute2DCoords(mol)
        
        coordinates = mol_to_coordinates(mol)
        
        return jsonify({
            'smiles': smiles,
            'coordinates': coordinates,
            'num_atoms': mol.GetNumAtoms(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating 2D coordinates: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/generate_3d', methods=['POST'])
@coordinates_bp.route('/3d', methods=['POST'])
def generate_3d_coordinates():
    """Generate 3D coordinates for a molecule"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        num_conformers = data.get('num_conformers', 1)
        optimize = data.get('optimize', True)
        random_seed = data.get('random_seed', -1)
        max_iterations = data.get('max_iterations', 200)
        force_field = data.get('force_field', 'UFF')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D conformers
        conf_ids = AllChem.EmbedMultipleConfs(
            mol, 
            numConfs=num_conformers,
            randomSeed=random_seed,
            maxAttempts=50
        )
        
        if not conf_ids:
            return jsonify({'error': 'Could not generate 3D coordinates', 'success': False}), 400
            
        energies = []
        conformers = []
        
        for conf_id in conf_ids:
            if optimize:
                # Optimize geometry based on force field
                if force_field.upper() == 'UFF':
                    result = AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=max_iterations)
                    if result == 0:  # Optimization successful
                        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                        energy = ff.CalcEnergy() if ff else None
                    else:
                        energy = None
                elif force_field.upper() in ['MMFF', 'MMFF94', 'MMFF94S']:
                    mmffVariant = 'MMFF94s' if force_field.upper() == 'MMFF94S' else 'MMFF94'
                    result = AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=max_iterations, mmffVariant=mmffVariant)
                    if result == 0:  # Optimization successful
                        # For MMFF, we need to create the force field with properties
                        mmffProps = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=mmffVariant)
                        if mmffProps is not None:
                            ff = AllChem.MMFFGetMoleculeForceField(mol, mmffProps, confId=conf_id)
                            energy = ff.CalcEnergy() if ff else None
                        else:
                            energy = None
                    else:
                        energy = None
                else:
                    result = AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=max_iterations)  # Fallback to UFF
                    if result == 0:  # Optimization successful
                        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                        energy = ff.CalcEnergy() if ff else None
                    else:
                        energy = None
                energies.append(energy)
            else:
                energies.append(None)
                
            coordinates = mol_to_coordinates(mol, conf_id)
            conformers.append({
                'conformer_id': conf_id,
                'energy': energies[-1],
                'coordinates': coordinates
            })
        
        # Sort by energy if optimized and add ranking/relative energy info
        if optimize and any(e is not None for e in energies):
            # Sort conformers by energy
            valid_conformers = [c for c in conformers if c['energy'] is not None]
            valid_conformers.sort(key=lambda x: x['energy'])
            
            # Add relative energy and rank information
            min_energy = valid_conformers[0]['energy'] if valid_conformers else 0
            for i, conformer in enumerate(valid_conformers):
                # Convert energy difference to kcal/mol (if energy is in hartree, multiply by 627.509)
                relative_energy = (conformer['energy'] - min_energy) * 627.509 if conformer['energy'] is not None else 0
                conformer['relative_energy_kcal_mol'] = relative_energy
                conformer['rank'] = i + 1
                conformer['energy_hartree'] = conformer['energy']  # Keep original energy
            
            conformers = valid_conformers
        else:
            # For non-optimized conformers, set default values
            for i, conformer in enumerate(conformers):
                conformer['relative_energy_kcal_mol'] = 0.0
                conformer['rank'] = i + 1
                conformer['energy_hartree'] = conformer.get('energy', None)
        
        return jsonify({
            'smiles': smiles,
            'num_conformers_requested': num_conformers,
            'num_conformers_generated': len(conformers),
            'optimized': optimize,
            'force_field': force_field if optimize else None,
            'conformers': conformers,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating 3D coordinates: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/optimize_geometry', methods=['POST'])
@coordinates_bp.route('/optimize', methods=['POST'])
def optimize_geometry():
    """Optimize molecular geometry using force field"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        force_field = data.get('force_field', 'UFF')
        max_iterations = data.get('max_iterations', 200)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Optimize based on force field
        if force_field.upper() == 'UFF':
            result = AllChem.UFFOptimizeMolecule(mol, maxIters=max_iterations)
            if result == 0:  # Optimization successful
                ff = AllChem.UFFGetMoleculeForceField(mol)
                energy = ff.CalcEnergy() if ff else None
            else:
                energy = None
        elif force_field.upper() in ['MMFF', 'MMFF94', 'MMFF94S']:
            # RDKit's MMFFOptimizeMolecule supports both MMFF94 and MMFF94s variants
            # The 's' variant uses a static parameterization
            mmffVariant = 'MMFF94s' if force_field.upper() == 'MMFF94S' else 'MMFF94'
            result = AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iterations, mmffVariant=mmffVariant)
            if result == 0:  # Optimization successful
                # For MMFF, we need to create the force field with properties
                mmffProps = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=mmffVariant)
                if mmffProps is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, mmffProps)
                    energy = ff.CalcEnergy() if ff else None
                else:
                    energy = None
            else:
                energy = None
        else:
            return jsonify({'error': f'Unsupported force field: {force_field}', 'success': False}), 400
            
        coordinates = mol_to_coordinates(mol)
        
        return jsonify({
            'smiles': smiles,
            'force_field': force_field,
            'optimization_energy': energy,
            'max_iterations': max_iterations,
            'coordinates': coordinates,
            'num_atoms': mol.GetNumAtoms(),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error optimizing geometry: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/align_molecules', methods=['POST'])
def align_molecules():
    """Align two molecules and calculate RMSD"""
    try:
        data = request.get_json()
        reference_smiles = data.get('reference_smiles')
        probe_smiles = data.get('probe_smiles')
        use_substructure = data.get('use_substructure', False)
        
        if not reference_smiles or not probe_smiles:
            return jsonify({'error': 'Both reference and probe SMILES required', 'success': False}), 400
            
        ref_mol = Chem.MolFromSmiles(reference_smiles)
        probe_mol = Chem.MolFromSmiles(probe_smiles)
        
        if ref_mol is None or probe_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400
            
        # Add hydrogens and generate 3D coordinates
        ref_mol = Chem.AddHs(ref_mol)
        probe_mol = Chem.AddHs(probe_mol)
        
        AllChem.EmbedMolecule(ref_mol, randomSeed=42)
        AllChem.EmbedMolecule(probe_mol, randomSeed=42)
        
        AllChem.UFFOptimizeMolecule(ref_mol)
        AllChem.UFFOptimizeMolecule(probe_mol)
        
        # Perform alignment
        if use_substructure:
            # Find common substructure and align
            match = probe_mol.GetSubstructMatch(ref_mol)
            if not match:
                # Try reverse match
                match = ref_mol.GetSubstructMatch(probe_mol)
                if not match:
                    return jsonify({'error': 'No common substructure found for alignment', 'success': False}), 400
                else:
                    # Swap molecules for alignment
                    ref_mol, probe_mol = probe_mol, ref_mol
                    reference_smiles, probe_smiles = probe_smiles, reference_smiles
            
            rmsd = rdMolAlign.AlignMol(probe_mol, ref_mol, atomMap=[(i, match[i]) for i in range(len(match))])
        else:
            # Align by optimal assignment
            rmsd = rdMolAlign.AlignMol(probe_mol, ref_mol)
        
        # Get aligned coordinates
        ref_coords = mol_to_coordinates(ref_mol)
        probe_coords = mol_to_coordinates(probe_mol)
        
        return jsonify({
            'reference_smiles': reference_smiles,
            'probe_smiles': probe_smiles,
            'rmsd': rmsd,
            'use_substructure': use_substructure,
            'reference_coordinates': ref_coords,
            'probe_coordinates': probe_coords,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error aligning molecules: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/conformer_search', methods=['POST'])
@coordinates_bp.route('/multiple', methods=['POST'])
def conformer_search():
    """Perform systematic conformer search"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        num_conformers = data.get('num_conformers', 50)
        energy_window = data.get('energy_window', 10.0)  # kcal/mol
        rmsd_threshold = data.get('rmsd_threshold', 1.0)  # Increased default to get more diverse conformers
        optimize = data.get('optimize', True)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate multiple conformers with better diversity
        conf_ids = AllChem.EmbedMultipleConfs(
            mol, 
            numConfs=num_conformers,
            maxAttempts=300,
            randomSeed=data.get('random_seed', 42),
            pruneRmsThresh=-1,  # Don't prune based on RMSD during generation
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
            ETversion=2,  # Use version 2 of the experimental torsion preferences
            enforceChirality=True
        )
        
        if not conf_ids:
            return jsonify({'error': 'Could not generate conformers', 'success': False}), 400
            
        conformers = []
        energies = []
        
        # Calculate energies
        for conf_id in conf_ids:
            if optimize:
                result = AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                if result != 0:  # Optimization failed
                    continue
                # Get final energy
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                energy = ff.CalcEnergy() if ff else None
            else:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                energy = ff.CalcEnergy() if ff else None
                
            if energy is not None:
                energies.append((conf_id, energy))
        
        # Sort by energy
        energies.sort(key=lambda x: x[1])
        
        if not energies:
            return jsonify({'error': 'No valid conformers generated', 'success': False}), 400
            
        # Filter by energy window (convert hartree to kcal/mol: 1 hartree = 627.509 kcal/mol)
        min_energy = energies[0][1]
        filtered_conformers = [(cid, e) for cid, e in energies if (e - min_energy) * 627.509 <= energy_window]
        
        # Apply RMSD pruning to keep diverse conformers (only if requested)
        if len(filtered_conformers) > 1 and rmsd_threshold > 0 and rmsd_threshold < 2.0:
            diverse_conformers = []
            for conf_id, energy in filtered_conformers:
                is_diverse = True
                for existing_id, _ in diverse_conformers:
                    try:
                        rmsd = AllChem.GetConformerRMS(mol, conf_id, existing_id, prealigned=True)
                        if rmsd < rmsd_threshold:
                            is_diverse = False
                            break
                    except:
                        # If RMSD calculation fails, consider it diverse
                        pass
                if is_diverse:
                    diverse_conformers.append((conf_id, energy))
            filtered_conformers = diverse_conformers
        
        # Build results
        for i, (conf_id, energy) in enumerate(filtered_conformers):
            coordinates = mol_to_coordinates(mol, conf_id)
            relative_energy = (energy - min_energy) * 627.509  # Convert to kcal/mol
            
            conformers.append({
                'rank': i + 1,
                'conformer_id': conf_id,
                'energy_hartree': energy,
                'relative_energy_kcal_mol': relative_energy,
                'coordinates': coordinates
            })
        
        return jsonify({
            'smiles': smiles,
            'num_conformers_requested': num_conformers,
            'num_conformers_generated': len(conf_ids),
            'num_conformers_optimized': len(energies),
            'num_conformers_in_window': len(conformers),
            'energy_window_kcal_mol': energy_window,
            'optimized': optimize,
            'conformers': conformers,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error in conformer search: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/export/<format>', methods=['POST'])
def export_coordinates(format):
    """Export molecular coordinates in various formats"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        coordinate_data = data.get('coordinate_data')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        if not coordinate_data:
            return jsonify({'error': 'Coordinate data required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens for 3D formats
        if format.lower() in ['sdf', 'mol', 'xyz']:
            mol = Chem.AddHs(mol)
            
        # Handle different coordinate data structures
        coordinates = None
        
        if 'conformers' in coordinate_data and coordinate_data['conformers']:
            # Multiple conformer data - use the first (lowest energy) conformer
            conformer_data = coordinate_data['conformers'][0]
            coordinates = conformer_data.get('coordinates', [])
        elif 'coordinates' in coordinate_data:
            # Single conformer data
            coordinates = coordinate_data['coordinates']
        elif 'mol_block' in coordinate_data:
            # Already have mol block - can return directly for SDF/MOL
            if format.lower() in ['sdf', 'mol']:
                return jsonify({
                    'smiles': smiles,
                    'format': format,
                    'file_content': coordinate_data['mol_block'],
                    'success': True
                })
        
        if not coordinates:
            return jsonify({'error': 'No coordinates found in data', 'success': False}), 400
            
        # Create conformer with coordinates
        conf = Chem.Conformer(mol.GetNumAtoms())
        for coord in coordinates:
            atom_idx = coord.get('atom_idx', coord.get('index', 0))
            x = coord.get('x', 0.0)
            y = coord.get('y', 0.0) 
            z = coord.get('z', 0.0)
            
            # Ensure atom index is within bounds
            if atom_idx < mol.GetNumAtoms():
                conf.SetAtomPosition(atom_idx, (float(x), float(y), float(z)))
            
        mol.AddConformer(conf)
        
        # Export in requested format
        if format.lower() == 'sdf':
            file_content = Chem.MolToMolBlock(mol)
            if not file_content:
                return jsonify({'error': 'Failed to generate SDF format', 'success': False}), 500
                
        elif format.lower() == 'mol':
            file_content = Chem.MolToMolBlock(mol)
            if not file_content:
                return jsonify({'error': 'Failed to generate MOL format', 'success': False}), 500
                
        elif format.lower() == 'xyz':
            # Generate XYZ format manually
            lines = [str(len(coordinates))]
            lines.append(f"Generated from SMILES: {smiles}")
            
            for coord in coordinates:
                atom_idx = coord.get('atom_idx', coord.get('index', 0))
                if atom_idx < mol.GetNumAtoms():
                    atom = mol.GetAtomWithIdx(atom_idx)
                    symbol = atom.GetSymbol()
                else:
                    # Fallback to symbol from coordinate data
                    symbol = coord.get('symbol', 'C')
                    
                x = coord.get('x', 0.0)
                y = coord.get('y', 0.0)
                z = coord.get('z', 0.0)
                lines.append(f"{symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
                    
            file_content = '\n'.join(lines)
            
        else:
            return jsonify({'error': f'Unsupported format: {format}', 'success': False}), 400
            
        return jsonify({
            'smiles': smiles,
            'format': format,
            'file_content': file_content,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error exporting coordinates as {format}: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/import_file', methods=['POST'])
def import_coordinates_file():
    """Import molecular coordinates from uploaded file (SDF, MOL, etc.)"""
    try:
        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected', 'success': False}), 400

        # Read file content
        file_content = file.read().decode('utf-8')
        validate_file_size(len(file_content.encode('utf-8')), max_size_mb=50)

        # Detect format from filename
        filename_lower = file.filename.lower()

        if filename_lower.endswith('.sdf'):
            # Parse SDF file
            molecules = parse_sdf_file(file_content)

            results = []
            for mol_data in molecules:
                try:
                    smiles = mol_data['smiles']
                    mol = Chem.MolFromSmiles(smiles)

                    if mol is None:
                        continue

                    # Try to get 3D coordinates if available from original SDF
                    # Note: RDKit's parse_sdf_file converts to SMILES, losing coordinates
                    # For full 3D support, we need to re-read the SDF
                    results.append({
                        'name': mol_data['name'],
                        'smiles': smiles,
                        'properties': mol_data.get('properties', {}),
                        'has_3d_coords': False  # SDF through SMILES loses coords
                    })

                except Exception as e:
                    logger.warning(f"Error processing molecule: {str(e)}")
                    continue

            return jsonify({
                'file_name': file.filename,
                'format': 'sdf',
                'num_molecules': len(results),
                'molecules': results,
                'success': True,
                'note': 'For 3D coordinate preservation, use direct SDF processing'
            })

        elif filename_lower.endswith('.mol'):
            # Parse single MOL file
            mol = Chem.MolFromMolBlock(file_content)

            if mol is None:
                return jsonify({'error': 'Invalid MOL file format', 'success': False}), 400

            # Extract coordinates if available
            coordinates = None
            if mol.GetNumConformers() > 0:
                coordinates = mol_to_coordinates(mol)

            smiles = Chem.MolToSmiles(mol)

            return jsonify({
                'file_name': file.filename,
                'format': 'mol',
                'smiles': smiles,
                'num_atoms': mol.GetNumAtoms(),
                'num_bonds': mol.GetNumBonds(),
                'has_3d_coords': coordinates is not None,
                'coordinates': coordinates,
                'success': True
            })

        else:
            # Try parsing as SMILES file
            molecules = parse_molecule_file(file.filename, file_content)

            if not molecules:
                return jsonify({'error': 'Could not parse file format', 'success': False}), 400

            results = []
            for mol_data in molecules:
                results.append({
                    'name': mol_data['name'],
                    'smiles': mol_data['smiles'],
                    'has_3d_coords': False
                })

            return jsonify({
                'file_name': file.filename,
                'format': 'smiles',
                'num_molecules': len(results),
                'molecules': results,
                'success': True,
                'note': 'SMILES format does not include 3D coordinates'
            })

    except Exception as e:
        logger.error(f"Error importing coordinates from file: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/o3a_align', methods=['POST'])
def o3a_alignment():
    """Align molecules using Open3DAlign (shape-based alignment)"""
    try:
        data = request.get_json()
        # Support both parameter naming conventions
        reference_smiles = data.get('reference_smiles') or data.get('ref_smiles')
        probe_smiles = data.get('probe_smiles')
        use_crippen = data.get('use_crippen', False)

        if not reference_smiles or not probe_smiles:
            return jsonify({'error': 'Both reference and probe SMILES required', 'success': False}), 400

        ref_mol = Chem.MolFromSmiles(reference_smiles)
        probe_mol = Chem.MolFromSmiles(probe_smiles)

        if ref_mol is None or probe_mol is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400

        # Add hydrogens and generate 3D coordinates
        ref_mol = Chem.AddHs(ref_mol)
        probe_mol = Chem.AddHs(probe_mol)

        AllChem.EmbedMolecule(ref_mol, randomSeed=42)
        AllChem.EmbedMolecule(probe_mol, randomSeed=42)

        AllChem.UFFOptimizeMolecule(ref_mol)
        AllChem.UFFOptimizeMolecule(probe_mol)

        # Perform O3A alignment
        if use_crippen:
            o3a = rdMolAlign.GetCrippenO3A(probe_mol, ref_mol)
        else:
            o3a = rdMolAlign.GetO3A(probe_mol, ref_mol)

        o3a.Align()
        score = o3a.Score()

        # Get aligned coordinates
        ref_coords = mol_to_coordinates(ref_mol)
        probe_coords = mol_to_coordinates(probe_mol)

        return jsonify({
            'reference_smiles': reference_smiles,
            'probe_smiles': probe_smiles,
            'alignment_method': 'Crippen O3A' if use_crippen else 'O3A',
            'score': score,
            'reference_coordinates': ref_coords,
            'probe_coordinates': probe_coords,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error performing O3A alignment: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/constrained_embed', methods=['POST'])
def constrained_embedding():
    """Generate 3D coordinates with distance constraints"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        constraints = data.get('constraints', [])  # List of {atom1, atom2, min_dist, max_dist}

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Use constrained embedding if constraints provided
        if constraints:
            from rdkit.Chem import rdDistGeom

            # Create bounds matrix
            bounds_matrix = rdDistGeom.GetMoleculeBoundsMatrix(mol)

            # Apply custom constraints
            for constraint in constraints:
                atom1 = constraint.get('atom1', 0)
                atom2 = constraint.get('atom2', 1)
                min_dist = constraint.get('min_dist', 1.0)
                max_dist = constraint.get('max_dist', 2.0)

                if atom1 < mol.GetNumAtoms() and atom2 < mol.GetNumAtoms():
                    bounds_matrix[atom1, atom2] = min_dist
                    bounds_matrix[atom2, atom1] = max_dist

            # Embed with constraints
            AllChem.EmbedMolecule(mol, useRandomCoords=False, boundsMatrix=bounds_matrix)
        else:
            # Standard embedding
            AllChem.EmbedMolecule(mol, randomSeed=42)

        # Optimize
        AllChem.UFFOptimizeMolecule(mol)

        coordinates = mol_to_coordinates(mol)

        return jsonify({
            'smiles': smiles,
            'num_constraints': len(constraints),
            'coordinates': coordinates,
            'num_atoms': mol.GetNumAtoms(),
            'success': True
        })

    except Exception as e:
        logger.error(f"Error with constrained embedding: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@coordinates_bp.route('/geometry_transforms', methods=['POST'])
def geometry_transforms():
    """Perform geometry transformations (bond lengths, angles, dihedrals)"""
    try:
        from rdkit.Chem import rdMolTransforms

        data = request.get_json()
        smiles = data.get('smiles')
        transform_type = data.get('transform_type')  # 'bond', 'angle', 'dihedral', 'centroid'

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)

        result = {'smiles': smiles}

        if transform_type == 'bond' and mol.GetNumBonds() > 0:
            # Get bond length for first bond
            atom1_idx = data.get('atom1', 0)
            atom2_idx = data.get('atom2', 1)

            if atom1_idx < mol.GetNumAtoms() and atom2_idx < mol.GetNumAtoms():
                bond_length = rdMolTransforms.GetBondLength(mol.GetConformer(), atom1_idx, atom2_idx)
                result['bond_length'] = bond_length

                # Set new bond length if provided
                new_length = data.get('new_length')
                if new_length:
                    rdMolTransforms.SetBondLength(mol.GetConformer(), atom1_idx, atom2_idx, new_length)
                    result['new_bond_length'] = new_length

        elif transform_type == 'angle' and mol.GetNumAtoms() >= 3:
            # Get angle for three atoms
            atom1_idx = data.get('atom1', 0)
            atom2_idx = data.get('atom2', 1)
            atom3_idx = data.get('atom3', 2)

            if all(idx < mol.GetNumAtoms() for idx in [atom1_idx, atom2_idx, atom3_idx]):
                angle_deg = rdMolTransforms.GetAngleDeg(mol.GetConformer(), atom1_idx, atom2_idx, atom3_idx)
                result['angle_degrees'] = angle_deg

                # Set new angle if provided
                new_angle = data.get('new_angle')
                if new_angle:
                    rdMolTransforms.SetAngleDeg(mol.GetConformer(), atom1_idx, atom2_idx, atom3_idx, new_angle)
                    result['new_angle_degrees'] = new_angle

        elif transform_type == 'dihedral' and mol.GetNumAtoms() >= 4:
            # Get dihedral angle for four atoms
            atom1_idx = data.get('atom1', 0)
            atom2_idx = data.get('atom2', 1)
            atom3_idx = data.get('atom3', 2)
            atom4_idx = data.get('atom4', 3)

            if all(idx < mol.GetNumAtoms() for idx in [atom1_idx, atom2_idx, atom3_idx, atom4_idx]):
                dihedral_deg = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), atom1_idx, atom2_idx, atom3_idx, atom4_idx)
                result['dihedral_degrees'] = dihedral_deg

                # Set new dihedral if provided
                new_dihedral = data.get('new_dihedral')
                if new_dihedral:
                    rdMolTransforms.SetDihedralDeg(mol.GetConformer(), atom1_idx, atom2_idx, atom3_idx, atom4_idx, new_dihedral)
                    result['new_dihedral_degrees'] = new_dihedral

        elif transform_type == 'centroid':
            # Compute centroid
            centroid = rdMolTransforms.ComputeCentroid(mol.GetConformer())
            result['centroid'] = {'x': centroid.x, 'y': centroid.y, 'z': centroid.z}

        elif transform_type == 'principal_axes':
            # Compute principal axes and moments
            axes, moments = rdMolTransforms.ComputePrincipalAxesAndMoments(mol.GetConformer())
            result['principal_moments'] = [moments[0], moments[1], moments[2]]

        result['coordinates'] = mol_to_coordinates(mol)
        result['success'] = True

        return jsonify(result)

    except Exception as e:
        logger.error(f"Error with geometry transforms: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500