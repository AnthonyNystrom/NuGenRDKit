from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions, Descriptors
import logging
import time
import itertools
import sys
import os

# Add utils to path for file parsing
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.file_parsers import parse_molecule_file, parse_rxn_file, validate_file_size, validate_molecule_count

reactions_bp = Blueprint('reactions', __name__)
logger = logging.getLogger(__name__)

@reactions_bp.route('/process', methods=['POST'])
def process_reaction():
    """Unified endpoint to process different types of chemical reactions"""
    try:
        start_time = time.time()
        data = request.get_json()
        reaction_type = data.get('reaction_type', 'smarts')
        max_products = data.get('max_products', 50)
        mode = data.get('mode', 'all')
        
        if reaction_type == 'enumeration':
            return process_enumeration(data, max_products, mode, start_time)
        else:
            return process_smarts_reaction(data, max_products, mode, start_time)
            
    except Exception as e:
        logger.error(f"Error processing reaction: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

def process_enumeration(data, max_products, mode, start_time):
    """Process product enumeration with building blocks and core structure"""
    building_blocks = data.get('building_blocks', [])
    core_structure = data.get('core_structure', '')
    
    if not building_blocks or not core_structure:
        return jsonify({'error': 'Building blocks and core structure required', 'success': False}), 400
    
    products = []
    unique_smiles = set()
    
    # Parse core structure to find attachment points
    core_mol = Chem.MolFromSmiles(core_structure)
    if core_mol is None:
        return jsonify({'error': 'Invalid core structure SMILES', 'success': False}), 400
    
    # Check if core structure has attachment points
    has_attachment_point = '[*' in core_structure or 'R' in core_structure
    
    for i, building_block in enumerate(building_blocks):
        if len(products) >= max_products:
            break
            
        try:
            bb_mol = Chem.MolFromSmiles(building_block)
            if bb_mol is None:
                continue
            
            if has_attachment_point:
                # Use proper attachment point substitution
                try:
                    # Create reaction SMARTS for attachment
                    reaction_smarts = f"{core_structure}.{building_block}>>[*:1]"
                    
                    # Simple substitution for now - more sophisticated methods could be used
                    product_smiles = core_structure
                    
                    # Handle different attachment point formats
                    if '[*:1]' in product_smiles:
                        product_smiles = product_smiles.replace('[*:1]', building_block)
                    elif '[*]' in product_smiles:
                        product_smiles = product_smiles.replace('[*]', building_block)
                    elif 'R1' in product_smiles:
                        product_smiles = product_smiles.replace('R1', building_block)
                        
                    # Clean up the SMILES and validate
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol is not None:
                        # Sanitize and canonicalize
                        Chem.SanitizeMol(product_mol)
                        canonical_smiles = Chem.MolToSmiles(product_mol)
                        
                        if mode == 'unique' and canonical_smiles in unique_smiles:
                            continue
                            
                        unique_smiles.add(canonical_smiles)
                        mw = Descriptors.MolWt(product_mol)
                        
                        products.append({
                            'smiles': canonical_smiles,
                            'name': f'Product {len(products) + 1}',
                            'molecular_weight': mw,
                            'building_block_index': i,
                            'building_block': building_block,
                            'core_structure': core_structure,
                            'yield': None,
                            'reaction_center': f'Attachment point substitution with {building_block}'
                        })
                        
                except Exception as e:
                    logger.warning(f"Attachment substitution failed for {building_block}: {str(e)}")
                    continue
            else:
                # No explicit attachment points - try simple concatenation
                try:
                    # Try different connection approaches
                    connection_attempts = [
                        f"{core_structure}.{building_block}",  # Separate molecules
                        f"{core_structure}{building_block}",   # Direct concatenation
                        f"{core_structure}C{building_block}",  # Carbon bridge
                        f"{core_structure}O{building_block}",  # Oxygen bridge
                        f"{core_structure}N{building_block}"   # Nitrogen bridge
                    ]
                    
                    for attempt in connection_attempts:
                        try:
                            product_mol = Chem.MolFromSmiles(attempt)
                            if product_mol is not None:
                                Chem.SanitizeMol(product_mol)
                                canonical_smiles = Chem.MolToSmiles(product_mol)
                                
                                if mode == 'unique' and canonical_smiles in unique_smiles:
                                    continue
                                    
                                unique_smiles.add(canonical_smiles)
                                mw = Descriptors.MolWt(product_mol)
                                
                                products.append({
                                    'smiles': canonical_smiles,
                                    'name': f'Product {len(products) + 1}',
                                    'molecular_weight': mw,
                                    'building_block_index': i,
                                    'building_block': building_block,
                                    'core_structure': core_structure,
                                    'connection_type': attempt.replace(core_structure, '').replace(building_block, '') or 'direct',
                                    'yield': None,
                                    'reaction_center': f'Direct connection via {attempt.replace(core_structure, "").replace(building_block, "") or "bond"} with {building_block}'
                                })
                                break  # Use first successful connection
                                
                        except Exception:
                            continue
                            
                except Exception as e:
                    logger.warning(f"Connection attempt failed for {building_block}: {str(e)}")
                    continue
                
        except Exception as e:
            logger.warning(f"Failed to process building block {building_block}: {str(e)}")
            continue
    
    processing_time = round((time.time() - start_time) * 1000, 1)
    
    return jsonify({
        'reaction_name': 'Product Enumeration',
        'reaction_type': 'enumeration',
        'core_structure': core_structure,
        'building_blocks': building_blocks,
        'products': products,
        'unique_products': len(unique_smiles),
        'processing_time': processing_time,
        'mode': mode,
        'success': True
    })

def process_smarts_reaction(data, max_products, mode, start_time):
    """Process SMARTS-based chemical reactions"""
    reaction_smarts = data.get('reaction_smarts', '')
    reactants = data.get('reactants', [])
    
    if not reaction_smarts or not reactants:
        return jsonify({'error': 'Reaction SMARTS and reactants required', 'success': False}), 400
    
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return jsonify({'error': 'Invalid reaction SMARTS', 'success': False}), 400
        
        # Convert reactant SMILES to molecules
        reactant_mols = []
        for smiles in reactants:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return jsonify({'error': f'Invalid reactant SMILES: {smiles}', 'success': False}), 400
            reactant_mols.append(mol)
        
        products = []
        unique_smiles = set()
        
        # Handle different numbers of reactants
        num_reactant_templates = rxn.GetNumReactantTemplates()
        
        if num_reactant_templates == 1:
            # Single reactant - apply reaction to each reactant
            for i, reactant_mol in enumerate(reactant_mols):
                if len(products) >= max_products:
                    break
                    
                try:
                    product_sets = rxn.RunReactants((reactant_mol,))
                    
                    for product_set in product_sets:
                        if len(products) >= max_products:
                            break
                            
                        valid_products = []
                        for product_mol in product_set:
                            if product_mol:
                                try:
                                    Chem.SanitizeMol(product_mol)
                                    smiles = Chem.MolToSmiles(product_mol)
                                    
                                    if mode == 'unique' and smiles in unique_smiles:
                                        continue
                                    
                                    unique_smiles.add(smiles)
                                    mw = Descriptors.MolWt(product_mol)
                                    
                                    products.append({
                                        'smiles': smiles,
                                        'name': f'Product {len(products) + 1}',
                                        'molecular_weight': mw,
                                        'reactant_index': i,
                                        'reactant': reactants[i],
                                        'yield': None,
                                        'reaction_center': 'Single reactant transformation'
                                    })
                                    
                                except Exception as e:
                                    logger.warning(f"Failed to sanitize product: {str(e)}")
                                    continue
                                    
                except Exception as e:
                    logger.warning(f"Failed to process reactant {reactants[i]}: {str(e)}")
                    continue
                    
        else:
            # Multiple reactants - intelligently match reactants to templates
            if len(reactant_mols) >= num_reactant_templates:
                # First, try to match reactants to reaction templates
                reactant_groups = [[] for _ in range(num_reactant_templates)]
                
                # For each reactant, check which template(s) it matches
                for mol_idx, mol in enumerate(reactant_mols):
                    for template_idx in range(num_reactant_templates):
                        # Check if this molecule matches this reactant template
                        test_reactants = [None] * num_reactant_templates
                        test_reactants[template_idx] = mol
                        
                        # Try a partial match to see if this molecule fits this position
                        matches_template = False
                        try:
                            # Create a test tuple with only this molecule in the right position
                            for combo in itertools.combinations(reactant_mols, num_reactant_templates):
                                if mol in combo:
                                    try:
                                        # Test if reaction works with this molecule in this position
                                        test_products = rxn.RunReactants(combo)
                                        if test_products:
                                            matches_template = True
                                            break
                                    except:
                                        pass
                        except:
                            pass
                        
                        if matches_template or True:  # For now, add all molecules to all groups
                            reactant_groups[template_idx].append((mol_idx, mol))
                
                # Generate all valid combinations
                valid_combinations = []
                
                # If we have molecules for each template position, create combinations
                if all(reactant_groups):
                    # For 2-reactant reactions, try all pairs
                    if num_reactant_templates == 2:
                        for mol1_idx, mol1 in reactant_groups[0]:
                            for mol2_idx, mol2 in reactant_groups[1]:
                                if mol1_idx != mol2_idx:  # Don't use same molecule twice
                                    valid_combinations.append(((mol1_idx, mol2_idx), (mol1, mol2)))
                    else:
                        # For reactions with more templates, use combinatorial approach
                        for reactant_combo in itertools.combinations(enumerate(reactant_mols), num_reactant_templates):
                            indices, mols = zip(*reactant_combo)
                            valid_combinations.append((indices, mols))
                else:
                    # Fallback to trying all combinations
                    for reactant_combo in itertools.combinations(enumerate(reactant_mols), num_reactant_templates):
                        indices, mols = zip(*reactant_combo)
                        valid_combinations.append((indices, mols))
                
                # Process valid combinations
                for combo_indices, reactant_combo in valid_combinations:
                    if len(products) >= max_products:
                        break
                        
                    try:
                        product_sets = rxn.RunReactants(reactant_combo)
                        
                        for product_set in product_sets:
                            if len(products) >= max_products:
                                break
                                
                            for product_mol in product_set:
                                if product_mol:
                                    try:
                                        Chem.SanitizeMol(product_mol)
                                        smiles = Chem.MolToSmiles(product_mol)
                                        
                                        if mode == 'unique' and smiles in unique_smiles:
                                            continue
                                        
                                        unique_smiles.add(smiles)
                                        mw = Descriptors.MolWt(product_mol)
                                        
                                        products.append({
                                            'smiles': smiles,
                                            'name': f'Product {len(products) + 1}',
                                            'molecular_weight': mw,
                                            'reactants': [Chem.MolToSmiles(mol) for mol in reactant_combo],
                                            'yield': None,
                                            'reaction_center': 'Multi-reactant combination'
                                        })
                                        
                                    except Exception as e:
                                        logger.warning(f"Failed to sanitize product: {str(e)}")
                                        continue
                                        
                    except Exception as e:
                        logger.warning(f"Failed to process reactant combination: {str(e)}")
                        continue
        
        processing_time = round((time.time() - start_time) * 1000, 1)
        
        # Determine reaction name from SMARTS pattern
        reaction_name = get_reaction_name(reaction_smarts)
        
        return jsonify({
            'reaction_name': reaction_name,
            'reaction_type': 'smarts',
            'reaction_smarts': reaction_smarts,
            'reactants': reactants,
            'products': products,
            'unique_products': len(unique_smiles),
            'processing_time': processing_time,
            'mode': mode,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error in SMARTS reaction processing: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

def get_reaction_name(reaction_smarts):
    """Get a descriptive name for common reaction SMARTS patterns"""
    reaction_patterns = {
        '[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])[Cl:3]': 'Acid Chloride Formation',
        '[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4]': 'Amide Formation',
        '[C:1][OH:2].[C:3](=[O:4])[OH:5]>>[C:1][O:2][C:3](=[O:4])': 'Esterification',
        '[c:1][N+:2](=[O:3])[O-:4]>>[c:1][N:2]': 'Nitro Reduction',
        '[C:1]=[C:2]>>[C:1][C:2]': 'Alkene Reduction'
    }
    
    return reaction_patterns.get(reaction_smarts, 'Custom Reaction')

@reactions_bp.route('/parse_smarts', methods=['POST'])
def parse_reaction_smarts():
    """Parse reaction SMARTS string"""
    try:
        data = request.get_json()
        reaction_smarts = data.get('reaction_smarts')
        
        if not reaction_smarts:
            return jsonify({'error': 'Reaction SMARTS required', 'success': False}), 400
            
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return jsonify({'error': 'Invalid reaction SMARTS', 'success': False}), 400
            
        # Get reaction information
        num_reactants = rxn.GetNumReactantTemplates()
        num_products = rxn.GetNumProductTemplates()
        
        reactants = []
        for i in range(num_reactants):
            reactant = rxn.GetReactantTemplate(i)
            reactants.append({
                'index': i,
                'smarts': Chem.MolToSmarts(reactant),
                'num_atoms': reactant.GetNumAtoms()
            })
            
        products = []
        for i in range(num_products):
            product = rxn.GetProductTemplate(i)
            products.append({
                'index': i,
                'smarts': Chem.MolToSmarts(product),
                'num_atoms': product.GetNumAtoms()
            })
        
        return jsonify({
            'reaction_smarts': reaction_smarts,
            'num_reactants': num_reactants,
            'num_products': num_products,
            'reactants': reactants,
            'products': products,
            'is_valid': True,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error parsing reaction SMARTS: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/run_reaction', methods=['POST'])
def run_reaction():
    """Run a chemical reaction on reactant molecules"""
    try:
        data = request.get_json()
        reaction_smarts = data.get('reaction_smarts')
        reactant_smiles = data.get('reactant_smiles', [])
        max_products = data.get('max_products', 100)
        
        if not reaction_smarts:
            return jsonify({'error': 'Reaction SMARTS required', 'success': False}), 400
            
        if not reactant_smiles:
            return jsonify({'error': 'Reactant SMILES required', 'success': False}), 400
            
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return jsonify({'error': 'Invalid reaction SMARTS', 'success': False}), 400
            
        # Convert SMILES to molecules
        reactants = []
        for smiles in reactant_smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return jsonify({'error': f'Invalid reactant SMILES: {smiles}', 'success': False}), 400
            reactants.append(mol)
            
        # Check if number of reactants matches reaction requirements
        if len(reactants) != rxn.GetNumReactantTemplates():
            return jsonify({
                'error': f'Reaction requires {rxn.GetNumReactantTemplates()} reactants, got {len(reactants)}',
                'success': False
            }), 400
            
        # Run the reaction
        products = rxn.RunReactants(reactants)
        
        product_results = []
        for i, product_set in enumerate(products[:max_products]):
            product_smiles = []
            for product in product_set:
                if product:
                    # Clean up the product
                    Chem.SanitizeMol(product)
                    product_smiles.append(Chem.MolToSmiles(product))
                    
            if product_smiles:
                product_results.append({
                    'product_set': i,
                    'products': product_smiles
                })
        
        return jsonify({
            'reaction_smarts': reaction_smarts,
            'reactants': reactant_smiles,
            'num_product_sets': len(product_results),
            'products': product_results,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error running reaction: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/validate_reaction', methods=['POST'])
def validate_reaction():
    """Validate a chemical reaction"""
    try:
        data = request.get_json()
        reactant_smiles = data.get('reactant_smiles', [])
        product_smiles = data.get('product_smiles', [])
        
        if not reactant_smiles or not product_smiles:
            return jsonify({'error': 'Both reactant and product SMILES required', 'success': False}), 400
            
        # Convert to molecules
        reactants = []
        for smiles in reactant_smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return jsonify({'error': f'Invalid reactant SMILES: {smiles}', 'success': False}), 400
            reactants.append(mol)
            
        products = []
        for smiles in product_smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return jsonify({'error': f'Invalid product SMILES: {smiles}', 'success': False}), 400
            products.append(mol)
            
        # Calculate mass balance
        def get_atom_count(mols):
            atom_count = {}
            for mol in mols:
                for atom in mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    atom_count[symbol] = atom_count.get(symbol, 0) + 1
            return atom_count
            
        reactant_atoms = get_atom_count(reactants)
        product_atoms = get_atom_count(products)
        
        # Check mass balance
        balanced = reactant_atoms == product_atoms
        
        # Calculate molecular weights
        reactant_mw = sum(Chem.Descriptors.MolWt(mol) for mol in reactants)
        product_mw = sum(Chem.Descriptors.MolWt(mol) for mol in products)
        mass_difference = abs(reactant_mw - product_mw)
        
        result = {
            'reactants': reactant_smiles,
            'products': product_smiles,
            'mass_balanced': balanced,
            'reactant_atoms': reactant_atoms,
            'product_atoms': product_atoms,
            'reactant_molecular_weight': reactant_mw,
            'product_molecular_weight': product_mw,
            'mass_difference': mass_difference,
            'success': True
        }
        
        if not balanced:
            # Find missing/extra atoms
            all_atoms = set(list(reactant_atoms.keys()) + list(product_atoms.keys()))
            missing_atoms = {}
            extra_atoms = {}
            
            for atom in all_atoms:
                reactant_count = reactant_atoms.get(atom, 0)
                product_count = product_atoms.get(atom, 0)
                difference = product_count - reactant_count
                
                if difference > 0:
                    extra_atoms[atom] = difference
                elif difference < 0:
                    missing_atoms[atom] = -difference
                    
            result['missing_atoms'] = missing_atoms
            result['extra_atoms'] = extra_atoms
        
        return jsonify(result)
        
    except Exception as e:
        logger.error(f"Error validating reaction: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/reaction_center', methods=['POST'])
def find_reaction_center():
    """Find the reaction center (atoms that change during reaction)"""
    try:
        data = request.get_json()
        reactant_smiles = data.get('reactant_smiles', [])
        product_smiles = data.get('product_smiles', [])
        
        if not reactant_smiles or not product_smiles:
            return jsonify({'error': 'Both reactant and product SMILES required', 'success': False}), 400
            
        if len(reactant_smiles) != 1 or len(product_smiles) != 1:
            return jsonify({'error': 'Single reactant and product required for reaction center analysis', 'success': False}), 400
            
        reactant = Chem.MolFromSmiles(reactant_smiles[0])
        product = Chem.MolFromSmiles(product_smiles[0])
        
        if reactant is None or product is None:
            return jsonify({'error': 'Invalid SMILES string(s)', 'success': False}), 400
            
        # Find maximum common substructure to identify unchanged parts
        from rdkit.Chem import rdFMCS
        
        mcs = rdFMCS.FindMCS([reactant, product], completeRingsOnly=True)
        
        if mcs.numAtoms == 0:
            return jsonify({'error': 'No common substructure found', 'success': False}), 400
            
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        
        # Find atoms that are part of the reaction center
        reactant_match = reactant.GetSubstructMatch(mcs_mol)
        product_match = product.GetSubstructMatch(mcs_mol)
        
        # Identify atoms not in the common substructure
        reactant_center = []
        product_center = []
        
        for atom_idx in range(reactant.GetNumAtoms()):
            if atom_idx not in reactant_match:
                atom = reactant.GetAtomWithIdx(atom_idx)
                reactant_center.append({
                    'atom_idx': atom_idx,
                    'symbol': atom.GetSymbol(),
                    'hybridization': str(atom.GetHybridization()),
                    'formal_charge': atom.GetFormalCharge()
                })
                
        for atom_idx in range(product.GetNumAtoms()):
            if atom_idx not in product_match:
                atom = product.GetAtomWithIdx(atom_idx)
                product_center.append({
                    'atom_idx': atom_idx,
                    'symbol': atom.GetSymbol(),
                    'hybridization': str(atom.GetHybridization()),
                    'formal_charge': atom.GetFormalCharge()
                })
        
        return jsonify({
            'reactant': reactant_smiles[0],
            'product': product_smiles[0],
            'mcs_smarts': mcs.smartsString,
            'mcs_atoms': mcs.numAtoms,
            'reactant_center': reactant_center,
            'product_center': product_center,
            'reactant_unchanged_atoms': len(reactant_match),
            'product_unchanged_atoms': len(product_match),
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error finding reaction center: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/enumerate_library', methods=['POST'])
def enumerate_library():
    """Enumerate a combinatorial library using a reaction"""
    try:
        data = request.get_json()
        reaction_smarts = data.get('reaction_smarts')
        reactant_lists = data.get('reactant_lists', [])
        max_products = data.get('max_products', 1000)
        
        if not reaction_smarts:
            return jsonify({'error': 'Reaction SMARTS required', 'success': False}), 400
            
        if not reactant_lists:
            return jsonify({'error': 'Reactant lists required', 'success': False}), 400
            
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return jsonify({'error': 'Invalid reaction SMARTS', 'success': False}), 400
            
        # Validate reactant lists
        if len(reactant_lists) != rxn.GetNumReactantTemplates():
            return jsonify({
                'error': f'Reaction requires {rxn.GetNumReactantTemplates()} reactant lists, got {len(reactant_lists)}',
                'success': False
            }), 400
            
        # Convert SMILES to molecules
        reactant_mols = []
        for i, smiles_list in enumerate(reactant_lists):
            mol_list = []
            for smiles in smiles_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return jsonify({'error': f'Invalid SMILES in reactant list {i}: {smiles}', 'success': False}), 400
                mol_list.append(mol)
            reactant_mols.append(mol_list)
            
        # Generate all combinations
        import itertools
        combinations = list(itertools.product(*reactant_mols))
        
        if len(combinations) > max_products:
            combinations = combinations[:max_products]
            
        # Run reactions
        library = []
        for combo_idx, reactant_combo in enumerate(combinations):
            try:
                products = rxn.RunReactants(reactant_combo)
                
                for product_set in products:
                    product_smiles = []
                    valid_products = True
                    
                    for product in product_set:
                        if product:
                            try:
                                Chem.SanitizeMol(product)
                                product_smiles.append(Chem.MolToSmiles(product))
                            except:
                                valid_products = False
                                break
                        else:
                            valid_products = False
                            break
                            
                    if valid_products and product_smiles:
                        library.append({
                            'combination_index': combo_idx,
                            'reactants': [Chem.MolToSmiles(mol) for mol in reactant_combo],
                            'products': product_smiles
                        })
                        
            except Exception as e:
                logger.warning(f"Failed to process combination {combo_idx}: {str(e)}")
                continue
        
        return jsonify({
            'reaction_smarts': reaction_smarts,
            'num_reactant_lists': len(reactant_lists),
            'reactant_list_sizes': [len(lst) for lst in reactant_lists],
            'total_combinations': len(combinations),
            'successful_reactions': len(library),
            'library': library,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error enumerating library: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/enumerate_library_file', methods=['POST'])
def enumerate_library_file():
    """Enumerate a combinatorial library using a reaction with reactants from uploaded files"""
    try:
        # Get reaction SMARTS from form data
        reaction_smarts = request.form.get('reaction_smarts')
        if not reaction_smarts:
            return jsonify({'error': 'Reaction SMARTS required', 'success': False}), 400

        # Validate reaction SMARTS
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return jsonify({'error': 'Invalid reaction SMARTS', 'success': False}), 400

        num_reactant_templates = rxn.GetNumReactantTemplates()

        # Get files for each reactant position
        reactant_lists = []
        for i in range(num_reactant_templates):
            file_key = f'reactants_{i}'
            if file_key not in request.files:
                return jsonify({'error': f'File for reactant position {i} required (key: {file_key})', 'success': False}), 400

            file = request.files[file_key]
            if file.filename == '':
                return jsonify({'error': f'No file selected for reactant position {i}', 'success': False}), 400

            # Read and parse file
            file_content = file.read().decode('utf-8')
            validate_file_size(len(file_content.encode('utf-8')), max_size_mb=20)

            molecules = parse_molecule_file(file.filename, file_content)
            if not molecules:
                return jsonify({'error': f'No valid molecules in file for reactant position {i}', 'success': False}), 400

            # Extract SMILES
            smiles_list = [mol['smiles'] for mol in molecules]
            reactant_lists.append(smiles_list)

        # Get max products
        max_products = int(request.form.get('max_products', 1000))

        # Convert SMILES to molecules
        reactant_mols = []
        for i, smiles_list in enumerate(reactant_lists):
            mol_list = []
            for smiles in smiles_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return jsonify({'error': f'Invalid SMILES in reactant list {i}: {smiles}', 'success': False}), 400
                mol_list.append(mol)
            reactant_mols.append(mol_list)

        # Generate all combinations
        combinations = list(itertools.product(*reactant_mols))

        if len(combinations) > max_products:
            combinations = combinations[:max_products]

        # Run reactions
        library = []
        for combo_idx, reactant_combo in enumerate(combinations):
            try:
                products = rxn.RunReactants(reactant_combo)

                for product_set in products:
                    product_smiles = []
                    valid_products = True

                    for product in product_set:
                        if product:
                            try:
                                Chem.SanitizeMol(product)
                                product_smiles.append(Chem.MolToSmiles(product))
                            except:
                                valid_products = False
                                break
                        else:
                            valid_products = False
                            break

                    if valid_products and product_smiles:
                        library.append({
                            'combination_index': combo_idx,
                            'reactants': [Chem.MolToSmiles(mol) for mol in reactant_combo],
                            'products': product_smiles
                        })

            except Exception as e:
                logger.warning(f"Failed to process combination {combo_idx}: {str(e)}")
                continue

        return jsonify({
            'reaction_smarts': reaction_smarts,
            'num_reactant_lists': len(reactant_lists),
            'reactant_list_sizes': [len(lst) for lst in reactant_lists],
            'total_combinations': len(combinations),
            'successful_reactions': len(library),
            'library': library,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error enumerating library from files: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/process_file', methods=['POST'])
def process_reaction_file():
    """Process reaction enumeration with building blocks or reactants from file"""
    try:
        # Get reaction type
        reaction_type = request.form.get('reaction_type', 'enumeration')
        max_products = int(request.form.get('max_products', 50))
        mode = request.form.get('mode', 'all')

        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected', 'success': False}), 400

        # Read file content
        file_content = file.read().decode('utf-8')
        validate_file_size(len(file_content.encode('utf-8')), max_size_mb=50)

        if reaction_type == 'enumeration':
            # Get core structure from form data
            core_structure = request.form.get('core_structure')
            if not core_structure:
                return jsonify({'error': 'Core structure required for enumeration', 'success': False}), 400

            # Parse building blocks from file
            molecules = parse_molecule_file(file.filename, file_content)
            if not molecules:
                return jsonify({'error': 'No valid molecules found in file', 'success': False}), 400

            validate_molecule_count(len(molecules), max_molecules=1000)

            building_blocks = [mol['smiles'] for mol in molecules]

            # Call existing enumeration logic
            start_time = time.time()
            result_data = {
                'building_blocks': building_blocks,
                'core_structure': core_structure
            }

            result = process_enumeration(result_data, max_products, mode, start_time)
            return result

        else:  # SMARTS reaction
            reaction_smarts = request.form.get('reaction_smarts')
            if not reaction_smarts:
                return jsonify({'error': 'Reaction SMARTS required', 'success': False}), 400

            # Parse reactants from file
            molecules = parse_molecule_file(file.filename, file_content)
            if not molecules:
                return jsonify({'error': 'No valid molecules found in file', 'success': False}), 400

            validate_molecule_count(len(molecules), max_molecules=1000)

            reactants = [mol['smiles'] for mol in molecules]

            # Call existing SMARTS reaction logic
            start_time = time.time()
            result_data = {
                'reaction_smarts': reaction_smarts,
                'reactants': reactants
            }

            result = process_smarts_reaction(result_data, max_products, mode, start_time)
            return result

    except Exception as e:
        logger.error(f"Error processing reaction from file: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/upload_rxn', methods=['POST'])
def upload_rxn_file():
    """Upload and parse RXN file format"""
    try:
        # Check if file is present
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected', 'success': False}), 400

        if not file.filename.lower().endswith('.rxn'):
            return jsonify({'error': 'File must be .rxn format', 'success': False}), 400

        # Read file content
        file_content = file.read().decode('utf-8')
        validate_file_size(len(file_content.encode('utf-8')), max_size_mb=10)

        # Parse RXN file
        rxn_data = parse_rxn_file(file_content)

        return jsonify({
            'file_name': file.filename,
            'reactants': rxn_data['reactants'],
            'products': rxn_data['products'],
            'reaction_smarts': rxn_data['reaction_smarts'],
            'num_reactants': rxn_data['num_reactants'],
            'num_products': rxn_data['num_products'],
            'success': True
        })

    except Exception as e:
        logger.error(f"Error uploading RXN file: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/reaction_fingerprint', methods=['POST'])
def reaction_fingerprint():
    """Generate reaction fingerprint"""
    try:
        data = request.get_json()
        reaction_smarts = data.get('reaction_smarts')
        fp_type = data.get('fp_type', 'structural')  # 'structural' or 'difference'
        fp_size = data.get('fp_size', 2048)

        if not reaction_smarts:
            return jsonify({'error': 'Reaction SMARTS required', 'success': False}), 400

        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return jsonify({'error': 'Invalid reaction SMARTS', 'success': False}), 400

        # Generate reaction fingerprint
        if fp_type == 'difference':
            fp = AllChem.CreateDifferenceFingerprintForReaction(rxn)
        else:
            fp = AllChem.CreateStructuralFingerprintForReaction(rxn)

        # Convert fingerprint to bit string
        fp_bits = fp.ToBitString()
        on_bits = [i for i, bit in enumerate(fp_bits) if bit == '1']

        return jsonify({
            'reaction_smarts': reaction_smarts,
            'fp_type': fp_type,
            'fp_size': len(fp_bits),
            'fingerprint': {
                'num_bits': len(fp_bits),
                'num_on_bits': len(on_bits)
            },
            'num_on_bits': len(on_bits),
            'on_bits': on_bits[:100],  # Return first 100 on bits
            'success': True
        })

    except Exception as e:
        logger.error(f"Error generating reaction fingerprint: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/brics_react', methods=['POST'])
def brics_reactions():
    """Generate products using BRICS fragmentation and recombination"""
    try:
        from rdkit.Chem import BRICS

        data = request.get_json()
        # Support both 'molecules' (from UI) and 'smiles_list' (legacy)
        smiles_list = data.get('molecules') or data.get('smiles_list', [])
        max_products = data.get('max_products', 50)
        min_fragment_size = data.get('min_fragment_size', 1)

        if not smiles_list or len(smiles_list) < 1:
            return jsonify({'error': 'At least one SMILES required', 'success': False}), 400

        # Fragment all molecules
        all_fragments = set()
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fragments = BRICS.BRICSDecompose(mol, minFragmentSize=min_fragment_size)
                all_fragments.update(fragments)

        # Build new molecules from fragments
        if len(all_fragments) < 2:
            return jsonify({'error': 'Not enough fragments for recombination', 'success': False}), 400

        fragment_list = list(all_fragments)
        products = []
        build_error = None

        # Use BRICS.BRICSBuild to recombine fragments
        # Note: This may fail if fragments are too complex or incompatible
        try:
            # Convert fragment SMILES to mol objects for BRICSBuild
            fragment_mols = []
            for frag_smiles in fragment_list:
                frag_mol = Chem.MolFromSmiles(frag_smiles)
                if frag_mol:
                    fragment_mols.append(frag_mol)

            if len(fragment_mols) >= 2:
                # BRICSBuild requires mol objects
                built_mols = BRICS.BRICSBuild(fragment_mols)
                count = 0
                for built_mol in built_mols:
                    if count >= max_products:
                        break

                    if built_mol:
                        try:
                            Chem.SanitizeMol(built_mol)
                            product_smiles = Chem.MolToSmiles(built_mol)
                            products.append({
                                'smiles': product_smiles,
                                'molecular_weight': Descriptors.MolWt(built_mol),
                                'num_atoms': built_mol.GetNumAtoms()
                            })
                            count += 1
                        except Exception:
                            continue
        except Exception as e:
            build_error = str(e)
            logger.warning(f"BRICS build error: {build_error}")

        return jsonify({
            'input_molecules': len(smiles_list),
            'num_fragments': len(all_fragments),
            'fragments': list(all_fragments),
            'num_products': len(products),
            'products': products,
            'build_attempted': True,
            'build_error': build_error if build_error and len(products) == 0 else None,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error with BRICS reactions: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@reactions_bp.route('/recap_react', methods=['POST'])
def recap_reactions():
    """Analyze molecule using RECAP decomposition"""
    try:
        from rdkit.Chem import Recap

        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Perform RECAP decomposition
        recap_tree = Recap.RecapDecompose(mol)

        # Get all leaf nodes (final fragments)
        fragments = []
        fragment_smiles_set = set()

        def get_all_nodes(node, level=0):
            if hasattr(node, 'mol') and node.mol:
                frag_smiles = Chem.MolToSmiles(node.mol)
                if frag_smiles not in fragment_smiles_set:
                    fragment_smiles_set.add(frag_smiles)
                    fragments.append(frag_smiles)

            # Children is a dictionary, not a list
            if hasattr(node, 'children') and node.children:
                for child_key, child_node in node.children.items():
                    get_all_nodes(child_node, level + 1)

        get_all_nodes(recap_tree)

        return jsonify({
            'smiles': smiles,
            'num_fragments': len(fragments),
            'fragments': fragments,
            'success': True
        })

    except Exception as e:
        logger.error(f"Error with RECAP reactions: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500