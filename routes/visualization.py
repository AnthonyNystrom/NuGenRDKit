from flask import Blueprint, request, jsonify, send_file
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor, AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.Chem import Scaffolds, DataStructs
import io
import base64
import tempfile
import os
import logging
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

visualization_bp = Blueprint('visualization', __name__)
logger = logging.getLogger(__name__)

def mol_to_svg(mol, size=(300, 300), highlight_atoms=None, highlight_bonds=None, show_atom_labels=False, show_atom_indices=False):
    """Convert molecule to SVG string"""
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    
    # Set drawing options
    drawer.drawOptions().addAtomIndices = show_atom_indices
    drawer.drawOptions().addStereoAnnotation = False
    
    # For atom labels, we need to modify how atoms are displayed
    if show_atom_labels:
        # Force all atoms to show their symbols (not just heteroatoms)
        for atom in mol.GetAtoms():
            atom.SetProp('atomLabel', atom.GetSymbol())
    
    if highlight_atoms or highlight_bonds:
        drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms or [], highlightBonds=highlight_bonds or [])
    else:
        drawer.DrawMolecule(mol)
    
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def mol_to_png_base64(mol, size=(300, 300), highlight_atoms=None, highlight_bonds=None, show_atom_labels=False, show_atom_indices=False):
    """Convert molecule to base64 encoded PNG"""
    if show_atom_labels or show_atom_indices:
        # Use RDKit drawer for more control when labels/indices are needed
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        drawer.drawOptions().addAtomIndices = show_atom_indices
        drawer.drawOptions().addStereoAnnotation = False
        
        # For atom labels, force all atoms to show their symbols
        if show_atom_labels:
            for atom in mol.GetAtoms():
                atom.SetProp('atomLabel', atom.GetSymbol())
        
        if highlight_atoms or highlight_bonds:
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms or [], highlightBonds=highlight_bonds or [])
        else:
            drawer.DrawMolecule(mol)
        
        drawer.FinishDrawing()
        img_data = drawer.GetDrawingText()
        img_base64 = base64.b64encode(img_data).decode('utf-8')
        return img_base64
    else:
        # Use standard RDKit image generation
        img = Draw.MolToImage(mol, size=size, highlightAtoms=highlight_atoms, highlightBonds=highlight_bonds)
        
        buffer = io.BytesIO()
        img.save(buffer, format='PNG')
        buffer.seek(0)
        
        img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
        return img_base64

@visualization_bp.route('/draw_svg', methods=['POST'])
def draw_molecule_svg():
    """Draw molecule as SVG"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        width = data.get('width', 300)
        height = data.get('height', 300)
        highlight_atoms = data.get('highlight_atoms', [])
        highlight_bonds = data.get('highlight_bonds', [])
        show_atom_labels = data.get('show_atom_labels', False)
        show_atom_indices = data.get('show_atom_indices', False)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate 2D coordinates if needed
        if mol.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(mol)
            
        svg_text = mol_to_svg(mol, (width, height), highlight_atoms, highlight_bonds, show_atom_labels, show_atom_indices)
        
        return jsonify({
            'smiles': smiles,
            'svg': svg_text,
            'width': width,
            'height': height,
            'show_atom_labels': show_atom_labels,
            'show_atom_indices': show_atom_indices,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error drawing SVG: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/draw_png', methods=['POST'])
def draw_molecule_png():
    """Draw molecule as PNG (base64 encoded)"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        width = data.get('width', 300)
        height = data.get('height', 300)
        highlight_atoms = data.get('highlight_atoms', [])
        highlight_bonds = data.get('highlight_bonds', [])
        show_atom_labels = data.get('show_atom_labels', False)
        show_atom_indices = data.get('show_atom_indices', False)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Generate 2D coordinates if needed
        if mol.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(mol)
            
        png_base64 = mol_to_png_base64(mol, (width, height), highlight_atoms, highlight_bonds, show_atom_labels, show_atom_indices)
        
        return jsonify({
            'smiles': smiles,
            'png_base64': png_base64,
            'width': width,
            'height': height,
            'show_atom_labels': show_atom_labels,
            'show_atom_indices': show_atom_indices,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error drawing PNG: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/draw_grid', methods=['POST'])
def draw_molecule_grid():
    """Draw multiple molecules in a grid"""
    try:
        data = request.get_json()
        smiles_list = data.get('smiles_list', [])
        legends = data.get('legends', [])
        mols_per_row = data.get('mols_per_row', 4)
        mol_size = data.get('mol_size', (200, 200))
        format_type = data.get('format', 'png')  # 'png' or 'svg'
        
        if not smiles_list:
            return jsonify({'error': 'SMILES list required', 'success': False}), 400
            
        # Convert SMILES to molecules
        mols = []
        valid_legends = []
        
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                # Generate 2D coordinates
                rdDepictor.Compute2DCoords(mol)
                mols.append(mol)
                
                # Add legend
                if i < len(legends):
                    valid_legends.append(legends[i])
                else:
                    valid_legends.append(smiles[:20] + '...' if len(smiles) > 20 else smiles)
            else:
                logger.warning(f"Invalid SMILES: {smiles}")
                
        if not mols:
            return jsonify({'error': 'No valid molecules found', 'success': False}), 400
            
        if format_type == 'svg':
            # Create SVG grid
            img = Draw.MolsToGridImage(
                mols, 
                molsPerRow=mols_per_row, 
                subImgSize=mol_size,
                legends=valid_legends,
                useSVG=True
            )
            
            return jsonify({
                'num_molecules': len(mols),
                'mols_per_row': mols_per_row,
                'format': 'svg',
                'svg': img,
                'success': True
            })
        else:
            # Create PNG grid
            img = Draw.MolsToGridImage(
                mols, 
                molsPerRow=mols_per_row, 
                subImgSize=mol_size,
                legends=valid_legends
            )
            
            buffer = io.BytesIO()
            img.save(buffer, format='PNG')
            buffer.seek(0)
            
            img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
            
            return jsonify({
                'num_molecules': len(mols),
                'mols_per_row': mols_per_row,
                'format': 'png',
                'png_base64': img_base64,
                'success': True
            })
            
    except Exception as e:
        logger.error(f"Error drawing grid: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/draw_reaction', methods=['POST'])
def draw_reaction():
    """Draw chemical reaction"""
    try:
        data = request.get_json()
        reaction_smarts = data.get('reaction_smarts')
        width = data.get('width', 800)
        height = data.get('height', 300)
        format_type = data.get('format', 'svg')
        
        if not reaction_smarts:
            return jsonify({'error': 'Reaction SMARTS required', 'success': False}), 400
            
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return jsonify({'error': 'Invalid reaction SMARTS', 'success': False}), 400
            
        if format_type == 'svg':
            drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            
            return jsonify({
                'reaction_smarts': reaction_smarts,
                'format': 'svg',
                'svg': drawer.GetDrawingText(),
                'width': width,
                'height': height,
                'success': True
            })
        else:
            # PNG format
            img = Draw.ReactionToImage(rxn, subImgSize=(width//4, height))
            
            buffer = io.BytesIO()
            img.save(buffer, format='PNG')
            buffer.seek(0)
            
            img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
            
            return jsonify({
                'reaction_smarts': reaction_smarts,
                'format': 'png',
                'png_base64': img_base64,
                'width': width,
                'height': height,
                'success': True
            })
            
    except Exception as e:
        logger.error(f"Error drawing reaction: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/highlight_substructure', methods=['POST'])
def highlight_substructure():
    """Draw molecule with highlighted substructure"""
    try:
        data = request.get_json()
        molecule_smiles = data.get('molecule_smiles')
        pattern_smarts = data.get('pattern_smarts')
        width = data.get('width', 300)
        height = data.get('height', 300)
        format_type = data.get('format', 'svg')
        
        if not molecule_smiles or not pattern_smarts:
            return jsonify({'error': 'Both molecule SMILES and pattern SMARTS required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(molecule_smiles)
        if mol is None:
            return jsonify({'error': 'Invalid molecule SMILES', 'success': False}), 400
            
        pattern = Chem.MolFromSmarts(pattern_smarts)
        if pattern is None:
            return jsonify({'error': 'Invalid pattern SMARTS', 'success': False}), 400
            
        # Find matches
        matches = mol.GetSubstructMatches(pattern)
        
        if not matches:
            return jsonify({'error': 'Pattern not found in molecule', 'success': False}), 400
            
        # Use first match for highlighting
        highlight_atoms = list(matches[0])
        
        # Find bonds within the match
        highlight_bonds = []
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in highlight_atoms and bond.GetEndAtomIdx() in highlight_atoms:
                highlight_bonds.append(bond.GetIdx())
                
        # Generate 2D coordinates if needed
        if mol.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(mol)
            
        if format_type == 'svg':
            svg_text = mol_to_svg(mol, (width, height), highlight_atoms, highlight_bonds, 
                                show_atom_labels=data.get('show_atom_labels', False),
                                show_atom_indices=data.get('show_atom_indices', False))
            
            return jsonify({
                'molecule_smiles': molecule_smiles,
                'pattern_smarts': pattern_smarts,
                'num_matches': len(matches),
                'highlighted_atoms': highlight_atoms,
                'highlighted_bonds': highlight_bonds,
                'format': 'svg',
                'svg': svg_text,
                'success': True
            })
        else:
            png_base64 = mol_to_png_base64(mol, (width, height), highlight_atoms, highlight_bonds,
                                         show_atom_labels=data.get('show_atom_labels', False),
                                         show_atom_indices=data.get('show_atom_indices', False))
            
            return jsonify({
                'molecule_smiles': molecule_smiles,
                'pattern_smarts': pattern_smarts,
                'num_matches': len(matches),
                'highlighted_atoms': highlight_atoms,
                'highlighted_bonds': highlight_bonds,
                'format': 'png',
                'png_base64': png_base64,
                'success': True
            })
            
    except Exception as e:
        logger.error(f"Error highlighting substructure: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/draw_3d', methods=['POST'])
@visualization_bp.route('/3d', methods=['POST'])
def draw_3d_molecule():
    """Generate 3D molecular structure (returns coordinates for 3D visualization)"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        optimize = data.get('optimize', True)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        if optimize:
            AllChem.UFFOptimizeMolecule(mol)
            
        # Extract 3D coordinates
        conf = mol.GetConformer()
        atoms = []
        bonds = []
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atoms.append({
                'index': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'x': pos.x,
                'y': pos.y,
                'z': pos.z,
                'atomic_number': atom.GetAtomicNum()
            })
            
        for bond in mol.GetBonds():
            bonds.append({
                'begin_atom': bond.GetBeginAtomIdx(),
                'end_atom': bond.GetEndAtomIdx(),
                'bond_type': str(bond.GetBondType()),
                'bond_order': bond.GetBondTypeAsDouble()
            })
            
        return jsonify({
            'smiles': smiles,
            'optimized': optimize,
            'num_atoms': len(atoms),
            'num_bonds': len(bonds),
            'atoms': atoms,
            'bonds': bonds,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating 3D structure: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/surface', methods=['POST'])
def molecular_surface():
    """Generate molecular surface data (simplified approach)"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        
        # Calculate molecular surface properties
        # Note: This is a simplified calculation - real surface calculations require more complex algorithms
        tpsa = Descriptors.TPSA(mol)
        logp = Descriptors.MolLogP(mol)
        vol = AllChem.ComputeMolVolume(mol)
        
        # Estimate surface area (simplified approximation)
        surface_area = tpsa * 1.5  # Very rough approximation
        
        # Get atom data for surface representation
        conf = mol.GetConformer()
        atoms = []
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atoms.append({
                'index': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'x': pos.x,
                'y': pos.y,
                'z': pos.z,
                'vdw_radius': Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum())
            })
            
        return jsonify({
            'smiles': smiles,
            'surface_area': surface_area,
            'volume': vol,
            'tpsa': tpsa,
            'logp': logp,
            'atoms': atoms,
            'num_atoms': len(atoms),
            'surface_type': 'van_der_waals',
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error generating molecular surface: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/pharmacophore', methods=['POST'])
def pharmacophore_analysis():
    """Generate pharmacophore features"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        
        # Identify pharmacophore features (simplified approach)
        features = []
        
        # H-bond donors (OH, NH)
        donor_patterns = [
            Chem.MolFromSmarts('[OH]'),
            Chem.MolFromSmarts('[NH]'),
            Chem.MolFromSmarts('[NH2]')
        ]
        
        for pattern in donor_patterns:
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    features.append({
                        'type': 'hydrogen_donor',
                        'atoms': list(match),
                        'description': 'H-bond donor'
                    })
        
        # H-bond acceptors (N, O with lone pairs)
        acceptor_patterns = [
            Chem.MolFromSmarts('[O]'),
            Chem.MolFromSmarts('[N]')
        ]
        
        for pattern in acceptor_patterns:
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    # Skip if already identified as donor
                    if not any(atom in sum([f['atoms'] for f in features if f['type'] == 'hydrogen_donor'], []) for atom in match):
                        features.append({
                            'type': 'hydrogen_acceptor',
                            'atoms': list(match),
                            'description': 'H-bond acceptor'
                        })
        
        # Aromatic rings
        aromatic_pattern = Chem.MolFromSmarts('c1ccccc1')
        if aromatic_pattern:
            matches = mol.GetSubstructMatches(aromatic_pattern)
            for match in matches:
                features.append({
                    'type': 'aromatic_ring',
                    'atoms': list(match),
                    'description': 'Aromatic ring'
                })
        
        # Hydrophobic regions (aliphatic carbons)
        hydrophobic_pattern = Chem.MolFromSmarts('[CH3,CH2,CH]')
        if hydrophobic_pattern:
            matches = mol.GetSubstructMatches(hydrophobic_pattern)
            for match in matches:
                features.append({
                    'type': 'hydrophobic',
                    'atoms': list(match),
                    'description': 'Hydrophobic region'
                })
        
        # Get 3D coordinates for features
        conf = mol.GetConformer()
        for feature in features:
            coords = []
            for atom_idx in feature['atoms']:
                pos = conf.GetAtomPosition(atom_idx)
                coords.append({'x': pos.x, 'y': pos.y, 'z': pos.z})
            feature['coordinates'] = coords
            
            # Calculate centroid
            if coords:
                feature['centroid'] = {
                    'x': sum(c['x'] for c in coords) / len(coords),
                    'y': sum(c['y'] for c in coords) / len(coords),
                    'z': sum(c['z'] for c in coords) / len(coords)
                }
        
        # Get atom coordinates for 3D visualization
        conf = mol.GetConformer()
        atoms = []
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atoms.append({
                'index': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'x': pos.x,
                'y': pos.y,
                'z': pos.z
            })
        
        return jsonify({
            'smiles': smiles,
            'features': features,
            'feature_count': len(features),
            'feature_types': list(set(f['type'] for f in features)),
            'atoms': atoms,
            'num_atoms': len(atoms),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error analyzing pharmacophore: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/scaffold', methods=['POST'])
def scaffold_analysis():
    """Analyze molecular scaffold"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens and generate 3D coordinates for visualization
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        
        # Generate different scaffold representations
        try:
            # Murcko scaffold (most common)
            murcko_scaffold = Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
            murcko_smiles = Chem.MolToSmiles(murcko_scaffold) if murcko_scaffold else None
            
            # Generic scaffold (remove side chains)
            generic_scaffold = Scaffolds.MurckoScaffold.MakeScaffoldGeneric(murcko_scaffold) if murcko_scaffold else None
            generic_smiles = Chem.MolToSmiles(generic_scaffold) if generic_scaffold else None
            
        except Exception as e:
            logger.warning(f"Error generating Murcko scaffold: {str(e)}")
            murcko_smiles = None
            generic_smiles = None
        
        # Simple ring system analysis
        ring_info = mol.GetRingInfo()
        ring_systems = []
        
        for ring in ring_info.AtomRings():
            ring_mol = Chem.PathToSubmol(mol, ring)
            if ring_mol:
                ring_smiles = Chem.MolToSmiles(ring_mol)
                ring_systems.append({
                    'smiles': ring_smiles,
                    'size': len(ring),
                    'atoms': list(ring),
                    'is_aromatic': all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
                })
        
        # Side chain analysis (atoms not in rings)
        ring_atoms = set()
        for ring in ring_info.AtomRings():
            ring_atoms.update(ring)
        
        side_chain_atoms = []
        for atom in mol.GetAtoms():
            if atom.GetIdx() not in ring_atoms:
                side_chain_atoms.append({
                    'index': atom.GetIdx(),
                    'symbol': atom.GetSymbol(),
                    'neighbors': [n.GetIdx() for n in atom.GetNeighbors()]
                })
        
        # Get atom coordinates for 3D visualization
        conf = mol.GetConformer()
        atoms = []
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atoms.append({
                'index': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'x': pos.x,
                'y': pos.y,
                'z': pos.z
            })
        
        return jsonify({
            'smiles': smiles,
            'scaffold_smiles': murcko_smiles,
            'generic_scaffold_smiles': generic_smiles,
            'ring_systems': ring_systems,
            'num_rings': len(ring_systems),
            'side_chain_atoms': side_chain_atoms,
            'num_side_chain_atoms': len(side_chain_atoms),
            'scaffold_complexity': len(ring_systems) + len(side_chain_atoms) / 10,
            'atoms': atoms,
            'num_atoms': len(atoms),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error analyzing scaffold: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/batch_file', methods=['POST'])
def batch_file_visualization():
    """Generate visualizations for molecules from an uploaded file"""
    try:
        from utils.file_parsers import parse_molecule_file
        from flask import send_file
        import io
        import zipfile
        import tempfile
        import os

        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'Empty filename', 'success': False}), 400

        viz_type = request.form.get('viz_type', '2d')
        output_format = request.form.get('output_format', 'grid')
        width = int(request.form.get('width', 300))
        height = int(request.form.get('height', 300))
        drawing_style = request.form.get('drawing_style', 'default')

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
                        'smiles': smiles,
                        'error': 'Invalid SMILES',
                        'svg': None
                    })
                    continue

                # Generate visualization based on type
                if viz_type == '2d':
                    from rdkit.Chem import Draw
                    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()
                    svg = drawer.GetDrawingText()

                    results.append({
                        'name': mol_data['name'],
                        'smiles': smiles,
                        'svg': svg,
                        'error': None
                    })
                    num_successful += 1

                elif viz_type == '3d':
                    # For 3D, we'll generate a 2D depiction as fallback
                    from rdkit.Chem import Draw
                    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()
                    svg = drawer.GetDrawingText()

                    results.append({
                        'name': mol_data['name'],
                        'smiles': smiles,
                        'svg': svg,
                        'error': None
                    })
                    num_successful += 1

                else:
                    # Default to 2D
                    from rdkit.Chem import Draw
                    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()
                    svg = drawer.GetDrawingText()

                    results.append({
                        'name': mol_data['name'],
                        'smiles': smiles,
                        'svg': svg,
                        'error': None
                    })
                    num_successful += 1

            except Exception as e:
                logger.warning(f"Error visualizing molecule {mol_data.get('name', 'unknown')}: {str(e)}")
                results.append({
                    'name': mol_data.get('name', 'unknown'),
                    'smiles': mol_data.get('smiles', ''),
                    'error': str(e),
                    'svg': None
                })

        # Return results based on format
        if output_format == 'zip':
            # Create ZIP file with all images
            memory_file = io.BytesIO()
            with zipfile.ZipFile(memory_file, 'w', zipfile.ZIP_DEFLATED) as zf:
                for idx, result in enumerate(results):
                    if result.get('svg') and not result.get('error'):
                        filename = f"{result['name']}.svg".replace(' ', '_')
                        zf.writestr(filename, result['svg'])

            memory_file.seek(0)
            return send_file(
                memory_file,
                mimetype='application/zip',
                as_attachment=True,
                download_name=f'visualizations_{file.filename}.zip'
            )
        else:
            # Return JSON with SVGs for grid display
            return jsonify({
                'success': True,
                'file_name': file.filename,
                'viz_type': viz_type,
                'num_molecules': len(molecules),
                'num_successful': num_successful,
                'results': results
            })

    except Exception as e:
        logger.error(f"Error in batch file visualization: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/similarity_map', methods=['POST'])
def similarity_map():
    """Generate similarity map showing atomic contributions to a property or fingerprint similarity"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        ref_smiles = data.get('ref_smiles')  # Reference molecule for similarity
        fp_type = data.get('fp_type', 'morgan')
        map_type = data.get('map_type', 'similarity')  # 'similarity' or 'weights'
        weights = data.get('weights')  # Custom atomic weights

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Generate 2D coordinates if needed
        if mol.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(mol)

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(6, 6))

        if map_type == 'similarity' and ref_smiles:
            # Similarity map between two molecules
            ref_mol = Chem.MolFromSmiles(ref_smiles)
            if ref_mol is None:
                return jsonify({'error': 'Invalid reference SMILES', 'success': False}), 400

            if fp_type == 'morgan':
                weights = SimilarityMaps.GetAtomicWeightsForFingerprint(ref_mol, mol, SimilarityMaps.GetMorganFingerprint)
                SimilarityMaps.GetSimilarityMapFromWeights(mol, weights, draw2d=rdMolDraw2D.MolDraw2DCairo(400, 400), contourLines=10, alpha=0.5)
            elif fp_type == 'ap':
                weights = SimilarityMaps.GetAtomicWeightsForFingerprint(ref_mol, mol, SimilarityMaps.GetAPFingerprint)
                SimilarityMaps.GetSimilarityMapFromWeights(mol, weights, draw2d=rdMolDraw2D.MolDraw2DCairo(400, 400), contourLines=10, alpha=0.5)
            elif fp_type == 'tt':
                weights = SimilarityMaps.GetAtomicWeightsForFingerprint(ref_mol, mol, SimilarityMaps.GetTTFingerprint)
                SimilarityMaps.GetSimilarityMapFromWeights(mol, weights, draw2d=rdMolDraw2D.MolDraw2DCairo(400, 400), contourLines=10, alpha=0.5)
        elif weights:
            # Custom weights provided
            SimilarityMaps.GetSimilarityMapFromWeights(mol, weights, draw2d=rdMolDraw2D.MolDraw2DCairo(400, 400), contourLines=10, alpha=0.5)
        else:
            # Generate default Morgan fingerprint similarity map
            SimilarityMaps.GetSimilarityMapFromWeights(
                mol,
                [1.0] * mol.GetNumAtoms(),  # Equal weights
                draw2d=rdMolDraw2D.MolDraw2DCairo(400, 400),
                contourLines=10,
                alpha=0.5
            )

        # Save to buffer
        buffer = io.BytesIO()
        plt.savefig(buffer, format='png', bbox_inches='tight', dpi=150)
        buffer.seek(0)
        plt.close(fig)

        # Convert to base64
        img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')

        return jsonify({
            'success': True,
            'smiles': smiles,
            'ref_smiles': ref_smiles,
            'fp_type': fp_type,
            'map_type': map_type,
            'image': f'data:image/png;base64,{img_base64}'
        })

    except Exception as e:
        logger.error(f"Error generating similarity map: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/fingerprint_bit', methods=['POST'])
def fingerprint_bit():
    """Visualize specific fingerprint bit and the atoms/bonds that activate it"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        fp_type = data.get('fp_type', 'morgan')
        bit_id = data.get('bit_id', 0)
        radius = data.get('radius', 2)

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Generate the fingerprint bit image
        if fp_type == 'morgan':
            # Get fingerprint with bit info
            bi = {}
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, 2048, bitInfo=bi)
            if bit_id in bi:
                # DrawMorganBit with useSVG=False returns PIL Image
                img = Draw.DrawMorganBit(mol, bit_id, bi, useSVG=False)
                # If it's bytes, wrap in BytesIO
                if isinstance(img, bytes):
                    img_base64 = base64.b64encode(img).decode('utf-8')
                else:
                    buffer = io.BytesIO()
                    img.save(buffer, format='PNG')
                    buffer.seek(0)
                    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
            else:
                # Provide helpful suggestions
                available_bits = sorted(list(bi.keys()))[:10]
                return jsonify({
                    'error': f'Bit {bit_id} not set in fingerprint',
                    'suggestion': f'Try one of these bits: {", ".join(map(str, available_bits))}',
                    'available_bits': available_bits,
                    'success': False
                }), 400
        elif fp_type == 'rdkit':
            bi = {}
            fp = Chem.RDKFingerprint(mol, fpSize=2048, bitInfo=bi)
            if bit_id in bi:
                # DrawRDKitBit with useSVG=False returns PIL Image
                img = Draw.DrawRDKitBit(mol, bit_id, bi, useSVG=False)
                # If it's bytes, wrap in BytesIO
                if isinstance(img, bytes):
                    img_base64 = base64.b64encode(img).decode('utf-8')
                else:
                    buffer = io.BytesIO()
                    img.save(buffer, format='PNG')
                    buffer.seek(0)
                    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
            else:
                # Provide helpful suggestions
                available_bits = sorted(list(bi.keys()))[:10]
                return jsonify({
                    'error': f'Bit {bit_id} not set in fingerprint',
                    'suggestion': f'Try one of these bits: {", ".join(map(str, available_bits))}',
                    'available_bits': available_bits,
                    'success': False
                }), 400
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400

        return jsonify({
            'success': True,
            'smiles': smiles,
            'fp_type': fp_type,
            'bit_id': bit_id,
            'image': f'data:image/png;base64,{img_base64}'
        })

    except Exception as e:
        logger.error(f"Error visualizing fingerprint bit: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/fingerprint_env', methods=['POST'])
def fingerprint_env():
    """Visualize atom environment that contributes to fingerprint"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        fp_type = data.get('fp_type', 'morgan')
        radius = data.get('radius', 2)
        atom_id = data.get('atom_id', 0)

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        if atom_id >= mol.GetNumAtoms():
            return jsonify({'error': f'Invalid atom_id: {atom_id}', 'success': False}), 400

        # Generate the environment image
        # DrawMorganEnv draws atom environments for Morgan fingerprints
        # For RDKit fingerprints, we highlight the atom and its neighbors
        if fp_type == 'morgan':
            img = Draw.DrawMorganEnv(mol, atomId=atom_id, radius=radius, molSize=(300, 300), useSVG=False)
        elif fp_type == 'rdkit':
            # For RDKit fingerprints, highlight the atom and its local environment
            # Since DrawRDKitEnv requires bondPath which we don't have without bit info,
            # we'll draw the molecule with the atom highlighted
            from rdkit.Chem import rdDepictor
            rdDepictor.Compute2DCoords(mol)

            # Highlight the atom and its neighbors
            highlight_atoms = [atom_id]
            neighbor_atoms = [x.GetIdx() for x in mol.GetAtomWithIdx(atom_id).GetNeighbors()]
            highlight_atoms.extend(neighbor_atoms)

            drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
            drawer.FinishDrawing()
            img = drawer.GetDrawingText()
        else:
            return jsonify({'error': f'Unsupported fingerprint type: {fp_type}', 'success': False}), 400

        # Convert to base64
        if isinstance(img, bytes):
            img_base64 = base64.b64encode(img).decode('utf-8')
        else:
            buffer = io.BytesIO()
            img.save(buffer, format='PNG')
            buffer.seek(0)
            img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')

        return jsonify({
            'success': True,
            'smiles': smiles,
            'fp_type': fp_type,
            'atom_id': atom_id,
            'radius': radius,
            'image': f'data:image/png;base64,{img_base64}'
        })

    except Exception as e:
        logger.error(f"Error visualizing fingerprint environment: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@visualization_bp.route('/matrix_grid', methods=['POST'])
def matrix_grid():
    """Generate matrix grid layout for molecule comparison"""
    try:
        data = request.get_json()
        smiles_matrix = data.get('smiles_matrix', [])  # 2D array of SMILES
        labels = data.get('labels', [])
        size = data.get('size', 200)
        highlight_substructure = data.get('highlight_substructure')

        if not smiles_matrix or not isinstance(smiles_matrix, list):
            return jsonify({'error': 'SMILES matrix required (2D array)', 'success': False}), 400

        # Convert SMILES to mol objects
        mol_matrix = []
        for row in smiles_matrix:
            mol_row = []
            for smiles in row:
                if smiles:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        mol_row.append(mol)
                    else:
                        mol_row.append(None)
                else:
                    mol_row.append(None)
            mol_matrix.append(mol_row)

        # Get substructure if provided
        highlight_atoms_list = []
        if highlight_substructure:
            subst_mol = Chem.MolFromSmarts(highlight_substructure)
            if subst_mol:
                for row in mol_matrix:
                    row_highlights = []
                    for mol in row:
                        if mol and mol.HasSubstructMatch(subst_mol):
                            matches = mol.GetSubstructMatches(subst_mol)
                            row_highlights.append(matches[0] if matches else [])
                        else:
                            row_highlights.append([])
                    highlight_atoms_list.append(row_highlights)

        # Generate matrix image
        # MolsMatrixToGridImage takes a 2D list directly
        if highlight_atoms_list:
            img = Draw.MolsMatrixToGridImage(
                mol_matrix,
                subImgSize=(size, size),
                highlightAtomListsMatrix=highlight_atoms_list
            )
        else:
            img = Draw.MolsMatrixToGridImage(
                mol_matrix,
                subImgSize=(size, size)
            )

        # Convert to base64
        buffer = io.BytesIO()
        img.save(buffer, format='PNG')
        buffer.seek(0)
        img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')

        return jsonify({
            'success': True,
            'num_rows': len(mol_matrix),
            'num_cols': len(mol_matrix[0]) if mol_matrix else 0,
            'image': f'data:image/png;base64,{img_base64}'
        })

    except Exception as e:
        logger.error(f"Error generating matrix grid: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500