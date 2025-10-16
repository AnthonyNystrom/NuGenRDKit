from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Fragments
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
import logging

properties_bp = Blueprint('properties', __name__)
logger = logging.getLogger(__name__)

@properties_bp.route('/physicochemical', methods=['POST'])
def physicochemical_properties():
    """Calculate physicochemical properties"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        properties = {
            'molecular_weight': Descriptors.MolWt(mol),
            'exact_molecular_weight': Descriptors.ExactMolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'logs': 'N/A',  # MolLogS not available in this RDKit version
            'tpsa': Descriptors.TPSA(mol),
            'num_h_donors': Descriptors.NumHDonors(mol),
            'num_h_acceptors': Descriptors.NumHAcceptors(mol),
            'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'num_heavy_atoms': Descriptors.HeavyAtomCount(mol),
            'heavy_atom_count': Descriptors.HeavyAtomCount(mol),  # Alias for frontend
            'fraction_csp3': getattr(Descriptors, 'FractionCsp3', lambda m: 0.0)(mol),  # Fallback if not available
            'molar_refractivity': Descriptors.MolMR(mol),
            'num_aromatic_rings': Descriptors.NumAromaticRings(mol),
            'num_saturated_rings': Descriptors.NumSaturatedRings(mol),
            'num_aliphatic_rings': Descriptors.NumAliphaticRings(mol),
            'ring_count': Descriptors.RingCount(mol),
            'num_rings': Descriptors.RingCount(mol),  # Alias for frontend
            'formal_charge': Chem.rdmolops.GetFormalCharge(mol),
            'molecular_formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
            # Additional properties for comprehensive analysis
            'num_heterocycles': len([ring for ring in mol.GetRingInfo().AtomRings() if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring)]),
            'num_carbocycles': len([ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)]),
            'num_bridgehead_atoms': len([atom for atom in mol.GetAtoms() if atom.IsInRing() and len([neighbor for neighbor in atom.GetNeighbors() if neighbor.IsInRing()]) > 2]),
            'num_aromatic_atoms': len([atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]),
            'num_aromatic_bonds': len([bond for bond in mol.GetBonds() if bond.GetIsAromatic()]),
            'num_stereocenters': len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            'num_unspecified_stereocenters': len([center for center in Chem.FindMolChiralCenters(mol, includeUnassigned=True) if center[1] == '?'])
        }
        
        return jsonify({
            'smiles': smiles,
            'properties': properties,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating physicochemical properties: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/drug_likeness', methods=['POST'])
def drug_likeness():
    """Calculate drug-likeness properties and rules"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Calculate basic properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        
        # Lipinski Rule of 5
        lipinski_violations = 0
        lipinski_violations += 1 if mw > 500 else 0
        lipinski_violations += 1 if logp > 5 else 0
        lipinski_violations += 1 if hbd > 5 else 0
        lipinski_violations += 1 if hba > 10 else 0
        
        # Veber Rules (oral bioavailability)
        veber_compliant = (tpsa <= 140 or hba <= 12) and rotb <= 10
        
        # Egan Rules (BBB permeability)
        egan_compliant = tpsa <= 131.6 and logp <= 5.88
        
        # Muegge Rules
        muegge_violations = 0
        muegge_violations += 1 if not (200 <= mw <= 600) else 0
        muegge_violations += 1 if not (-2 <= logp <= 5) else 0
        muegge_violations += 1 if tpsa > 150 else 0
        muegge_violations += 1 if mol.GetNumAtoms() > 70 else 0
        muegge_violations += 1 if rotb > 15 else 0
        muegge_violations += 1 if Descriptors.RingCount(mol) > 7 else 0
        
        # PAINS filters (simplified)
        pains_alerts = []
        # Add basic PAINS patterns
        pains_smarts = [
            'C1=CC=CC=C1', # Simple benzene check (placeholder)
        ]
        
        for smarts in pains_smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                pains_alerts.append(smarts)
        
        properties = {
            'molecular_weight': mw,
            'logp': logp,
            'num_h_donors': hbd,
            'num_h_acceptors': hba,
            'tpsa': tpsa,
            'num_rotatable_bonds': rotb,
            'lipinski_violations': lipinski_violations,
            'passes_lipinski': lipinski_violations <= 1,
            'veber_compliant': veber_compliant,
            'egan_compliant': egan_compliant,
            'muegge_violations': muegge_violations,
            'passes_muegge': muegge_violations == 0,
            'pains_alerts': len(pains_alerts),
            'drug_like_score': (
                (1 if lipinski_violations <= 1 else 0) +
                (1 if veber_compliant else 0) +
                (1 if egan_compliant else 0) +
                (1 if muegge_violations == 0 else 0)
            ) / 4.0
        }
        
        return jsonify({
            'smiles': smiles,
            'properties': properties,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating drug-likeness: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/fragments', methods=['POST'])
def molecular_fragments():
    """Calculate molecular fragment counts"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Calculate fragment counts
        fragments = {}
        
        # Get all fragment descriptor functions
        if hasattr(Fragments, 'descList'):
            fragment_functions = [
                (name, func) for name, func in Fragments.descList
            ]
        else:
            # Fallback for newer RDKit versions
            fragment_functions = [
                ('fr_Al_COO', getattr(Fragments, 'fr_Al_COO', lambda m: 0)),
                ('fr_Al_OH', getattr(Fragments, 'fr_Al_OH', lambda m: 0)),
                ('fr_ArN', getattr(Fragments, 'fr_ArN', lambda m: 0)),
                ('fr_Ar_COO', getattr(Fragments, 'fr_Ar_COO', lambda m: 0)),
                ('fr_benzene', getattr(Fragments, 'fr_benzene', lambda m: 0)),
                ('fr_NH2', getattr(Fragments, 'fr_NH2', lambda m: 0))
            ]
        
        for name, func in fragment_functions:
            try:
                fragments[name] = func(mol)
            except Exception as e:
                logger.warning(f"Failed to calculate fragment {name}: {str(e)}")
                fragments[name] = None
        
        # Calculate additional ring properties for fragment analysis
        ring_properties = {
            'num_aromatic_rings': Descriptors.NumAromaticRings(mol),
            'num_saturated_rings': Descriptors.NumSaturatedRings(mol),
            'num_aliphatic_rings': Descriptors.NumAliphaticRings(mol),
            'num_heterocycles': len([ring for ring in mol.GetRingInfo().AtomRings() if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring)]),
            'num_carbocycles': len([ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)]),
            'num_bridgehead_atoms': len([atom for atom in mol.GetAtoms() if atom.IsInRing() and len([neighbor for neighbor in atom.GetNeighbors() if neighbor.IsInRing()]) > 2])
        }
        
        return jsonify({
            'smiles': smiles,
            'properties': {
                'fragments': fragments,
                'total_fragments': sum(v for v in fragments.values() if v is not None),
                **ring_properties
            },
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating molecular fragments: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/scaffold', methods=['POST'])
def molecular_scaffold():
    """Calculate molecular scaffold (Murcko scaffold)"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Get Murcko scaffold
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold) if scaffold else None
        
        # Get generic scaffold
        generic_scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold) if scaffold else None
        generic_scaffold_smiles = Chem.MolToSmiles(generic_scaffold) if generic_scaffold else None
        
        result = {
            'smiles': smiles,
            'scaffold_smiles': scaffold_smiles,
            'generic_scaffold_smiles': generic_scaffold_smiles,
            'success': True
        }
        
        if scaffold:
            result['scaffold_molecular_weight'] = Descriptors.MolWt(scaffold)
            result['scaffold_heavy_atoms'] = scaffold.GetNumHeavyAtoms()
            result['scaffold_rings'] = Descriptors.RingCount(scaffold)
            
        return jsonify(result)
        
    except Exception as e:
        logger.error(f"Error calculating molecular scaffold: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/aromaticity', methods=['POST'])
def aromaticity_analysis():
    """Analyze aromaticity properties"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Count aromatic atoms and bonds
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_bonds = sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic())
        
        # Get aromatic ring information
        ring_info = mol.GetRingInfo()
        aromatic_rings = []
        
        for ring in ring_info.AtomRings():
            if len(ring) >= 5:  # Consider rings with 5+ atoms
                ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                if all(atom.GetIsAromatic() for atom in ring_atoms):
                    aromatic_rings.append({
                        'size': len(ring),
                        'atoms': list(ring),
                        'atom_symbols': [atom.GetSymbol() for atom in ring_atoms]
                    })
        
        properties = {
            'num_aromatic_atoms': aromatic_atoms,
            'num_aromatic_bonds': aromatic_bonds,
            'num_aromatic_rings': len(aromatic_rings),
            'aromatic_rings': aromatic_rings,
            'fraction_aromatic_atoms': aromatic_atoms / mol.GetNumAtoms() if mol.GetNumAtoms() > 0 else 0,
            'num_benzene_rings': Descriptors.NumAromaticRings(mol),
            'num_saturated_rings': Descriptors.NumSaturatedRings(mol),
            'num_aliphatic_rings': Descriptors.NumAliphaticRings(mol)
        }
        
        return jsonify({
            'smiles': smiles,
            'properties': properties,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error analyzing aromaticity: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/charges', methods=['POST'])
def charge_analysis():
    """Analyze charge distribution"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        method = data.get('method', 'gasteiger')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Add hydrogens for charge calculation
        mol = Chem.AddHs(mol)
        
        # Calculate partial charges
        if method.lower() == 'gasteiger':
            from rdkit.Chem import rdPartialCharges
            rdPartialCharges.ComputeGasteigerCharges(mol)
            
            charges = []
            total_charge = 0
            
            for atom in mol.GetAtoms():
                charge = float(atom.GetProp('_GasteigerCharge'))
                charges.append({
                    'atom_idx': atom.GetIdx(),
                    'symbol': atom.GetSymbol(),
                    'charge': charge
                })
                total_charge += charge
                
            formal_charge = Chem.rdmolops.GetFormalCharge(mol)
            
            return jsonify({
                'smiles': smiles,
                'method': method,
                'formal_charge': formal_charge,
                'total_partial_charge': total_charge,
                'atom_charges': charges,
                'num_atoms_with_charges': len(charges),
                'success': True
            })
        else:
            return jsonify({'error': f'Unsupported charge method: {method}', 'success': False}), 400
        
    except Exception as e:
        logger.error(f"Error analyzing charges: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/stereochemistry', methods=['POST'])
def stereochemistry_analysis():
    """Analyze stereochemistry"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Assign stereochemistry
        Chem.rdmolops.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        
        # Find stereocenters
        stereocenters = []
        for atom in mol.GetAtoms():
            if atom.HasProp('_CIPCode'):
                stereocenters.append({
                    'atom_idx': atom.GetIdx(),
                    'symbol': atom.GetSymbol(),
                    'cip_code': atom.GetProp('_CIPCode'),
                    'chiral_tag': str(atom.GetChiralTag())
                })
        
        # Find stereo bonds
        stereo_bonds = []
        for bond in mol.GetBonds():
            if bond.GetStereo() != Chem.rdchem.BondStereo.STEREONONE:
                stereo_bonds.append({
                    'bond_idx': bond.GetIdx(),
                    'begin_atom': bond.GetBeginAtomIdx(),
                    'end_atom': bond.GetEndAtomIdx(),
                    'stereo': str(bond.GetStereo())
                })
        
        properties = {
            'num_stereocenters': len(stereocenters),
            'num_stereo_bonds': len(stereo_bonds),
            'stereocenters': stereocenters,
            'stereo_bonds': stereo_bonds,
            'is_chiral': len(stereocenters) > 0,
            'num_undefined_stereocenters': len([sc for sc in stereocenters if sc['cip_code'] == '?'])
        }
        
        return jsonify({
            'smiles': smiles,
            'properties': properties,
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error analyzing stereochemistry: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/admet', methods=['POST'])
def admet_predictions():
    """ADMET predictions using simple rules-based approach"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
            
        # Calculate basic properties for ADMET predictions
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        
        # ADMET predictions based on empirical rules
        predictions = {
            'absorption': {
                'human_intestinal_absorption': 'High' if tpsa <= 140 and mw <= 500 else 'Low',
                'caco2_permeability': 'High' if logp > 0 and tpsa <= 60 else 'Medium' if tpsa <= 90 else 'Low',
                'bioavailability_score': min(1.0, (500 - mw) / 500 * (140 - tpsa) / 140) if mw <= 500 and tpsa <= 140 else 0.0
            },
            'distribution': {
                'bbb_permeation': 'High' if logp <= 3 and tpsa <= 90 and mw <= 450 else 'Low',
                'cns_permeation': 'High' if logp <= 3 and tpsa <= 76 and mw <= 360 else 'Low',
                'plasma_protein_binding': 'High' if logp > 1.5 else 'Low'
            },
            'metabolism': {
                'cyp2d6_substrate': 'Likely' if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()) else 'Unlikely',
                'cyp3a4_substrate': 'Likely' if mw > 300 and logp > 2 else 'Unlikely',
                'cyp_inhibition_risk': 'High' if logp > 3 and mw > 400 else 'Low'
            },
            'excretion': {
                'renal_clearance': 'High' if tpsa > 75 and logp < 2 else 'Low',
                'biliary_excretion': 'High' if mw > 400 else 'Low'
            },
            'toxicity': {
                'hepatotoxicity': 'Risk' if logp > 3 or mw > 500 else 'Low Risk',
                'cardiotoxicity': 'Risk' if rotb > 7 and logp > 4 else 'Low Risk',
                'ames_mutagenicity': 'Non-mutagenic' if not any([
                    mol.HasSubstructMatch(Chem.MolFromSmarts('[N+](=O)[O-]')),  # Nitro groups
                    mol.HasSubstructMatch(Chem.MolFromSmarts('N=N')),  # Azo groups
                ]) else 'Alert'
            }
        }
        
        # Overall ADMET score
        positive_scores = 0
        total_scores = 0
        
        for category in predictions.values():
            for key, value in category.items():
                total_scores += 1
                if value in ['High', 'Good', 'Non-mutagenic', 'Low Risk']:
                    positive_scores += 1
                elif value in ['Medium', 'Likely']:
                    positive_scores += 0.5
        
        admet_score = positive_scores / total_scores if total_scores > 0 else 0.0
        
        return jsonify({
            'smiles': smiles,
            'properties': {
                'molecular_weight': mw,
                'logp': logp,
                'tpsa': tpsa,
                'num_h_donors': hbd,
                'num_h_acceptors': hba,
                'num_rotatable_bonds': rotb,
                'admet_predictions': predictions,
                'admet_score': admet_score
            },
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error calculating ADMET predictions: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/all', methods=['POST'])
def comprehensive_analysis():
    """Comprehensive analysis combining all property types"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400
        
        # Calculate all properties directly
        all_properties = {}
        
        # Physicochemical properties
        try:
            all_properties['molecular_weight'] = Descriptors.MolWt(mol)
            all_properties['exact_molecular_weight'] = Descriptors.ExactMolWt(mol)
            all_properties['logp'] = Descriptors.MolLogP(mol)
            all_properties['tpsa'] = Descriptors.TPSA(mol)
            all_properties['num_h_donors'] = Descriptors.NumHDonors(mol)
            all_properties['num_h_acceptors'] = Descriptors.NumHAcceptors(mol)
            all_properties['num_rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
            all_properties['num_heavy_atoms'] = Descriptors.HeavyAtomCount(mol)
            all_properties['heavy_atom_count'] = Descriptors.HeavyAtomCount(mol)
            all_properties['molar_refractivity'] = Descriptors.MolMR(mol)
            all_properties['num_aromatic_rings'] = Descriptors.NumAromaticRings(mol)
            all_properties['num_saturated_rings'] = Descriptors.NumSaturatedRings(mol)
            all_properties['num_aliphatic_rings'] = Descriptors.NumAliphaticRings(mol)
            all_properties['ring_count'] = Descriptors.RingCount(mol)
            all_properties['num_rings'] = Descriptors.RingCount(mol)
            all_properties['formal_charge'] = Chem.rdmolops.GetFormalCharge(mol)
            all_properties['molecular_formula'] = Chem.rdMolDescriptors.CalcMolFormula(mol)
        except Exception as e:
            logger.warning(f"Failed to calculate physicochemical properties: {str(e)}")
            
        # Drug-likeness properties
        try:
            mw = all_properties.get('molecular_weight', 0)
            logp = all_properties.get('logp', 0)
            hbd = all_properties.get('num_h_donors', 0)
            hba = all_properties.get('num_h_acceptors', 0)
            tpsa = all_properties.get('tpsa', 0)
            rotb = all_properties.get('num_rotatable_bonds', 0)
            
            # Lipinski Rule violations
            lipinski_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
            all_properties['lipinski_violations'] = lipinski_violations
            all_properties['passes_lipinski'] = lipinski_violations <= 1
            all_properties['veber_compliant'] = (tpsa <= 140 or hba <= 12) and rotb <= 10
            all_properties['egan_compliant'] = tpsa <= 131.6 and logp <= 5.88
        except Exception as e:
            logger.warning(f"Failed to calculate drug-likeness properties: {str(e)}")
            
        # Fragment analysis
        try:
            fragments = {}
            if hasattr(Fragments, 'descList'):
                fragment_functions = [(name, func) for name, func in Fragments.descList]
            else:
                fragment_functions = [
                    ('fr_Al_COO', getattr(Fragments, 'fr_Al_COO', lambda m: 0)),
                    ('fr_Al_OH', getattr(Fragments, 'fr_Al_OH', lambda m: 0)),
                    ('fr_ArN', getattr(Fragments, 'fr_ArN', lambda m: 0)),
                    ('fr_benzene', getattr(Fragments, 'fr_benzene', lambda m: 0)),
                    ('fr_NH2', getattr(Fragments, 'fr_NH2', lambda m: 0))
                ]
            
            for name, func in fragment_functions:
                try:
                    fragments[name] = func(mol)
                except Exception as e:
                    logger.warning(f"Failed to calculate fragment {name}: {str(e)}")
                    fragments[name] = None
            
            all_properties['fragments'] = fragments
            all_properties['total_fragments'] = sum(v for v in fragments.values() if v is not None)
        except Exception as e:
            logger.warning(f"Failed to calculate fragment analysis: {str(e)}")
            
        # Aromaticity analysis
        try:
            aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
            aromatic_bonds = sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic())
            all_properties['num_aromatic_atoms'] = aromatic_atoms
            all_properties['num_aromatic_bonds'] = aromatic_bonds
            all_properties['fraction_aromatic_atoms'] = aromatic_atoms / mol.GetNumAtoms() if mol.GetNumAtoms() > 0 else 0
        except Exception as e:
            logger.warning(f"Failed to calculate aromaticity analysis: {str(e)}")
            
        # Stereochemistry analysis
        try:
            Chem.rdmolops.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
            stereocenters = [atom for atom in mol.GetAtoms() if atom.HasProp('_CIPCode')]
            all_properties['num_stereocenters'] = len(stereocenters)
            all_properties['is_chiral'] = len(stereocenters) > 0
        except Exception as e:
            logger.warning(f"Failed to calculate stereochemistry analysis: {str(e)}")
            
        # ADMET predictions
        try:
            mw = all_properties.get('molecular_weight', 0)
            logp = all_properties.get('logp', 0)
            tpsa = all_properties.get('tpsa', 0)
            hbd = all_properties.get('num_h_donors', 0)
            hba = all_properties.get('num_h_acceptors', 0)
            rotb = all_properties.get('num_rotatable_bonds', 0)
            
            admet_predictions = {
                'absorption': {
                    'human_intestinal_absorption': 'High' if tpsa <= 140 and mw <= 500 else 'Low',
                    'bioavailability_score': min(1.0, (500 - mw) / 500 * (140 - tpsa) / 140) if mw <= 500 and tpsa <= 140 else 0.0
                },
                'distribution': {
                    'bbb_permeation': 'High' if logp <= 3 and tpsa <= 90 and mw <= 450 else 'Low',
                    'cns_permeation': 'High' if logp <= 3 and tpsa <= 76 and mw <= 360 else 'Low'
                },
                'metabolism': {
                    'cyp3a4_substrate': 'Likely' if mw > 300 and logp > 2 else 'Unlikely'
                },
                'toxicity': {
                    'hepatotoxicity': 'Risk' if logp > 3 or mw > 500 else 'Low Risk'
                }
            }
            
            all_properties['admet_predictions'] = admet_predictions
            
            # Calculate ADMET score
            positive_scores = sum(1 for category in admet_predictions.values() 
                                for value in category.values() 
                                if value in ['High', 'Good', 'Non-mutagenic', 'Low Risk', 'Unlikely'])
            total_scores = sum(len(category) for category in admet_predictions.values())
            all_properties['admet_score'] = positive_scores / total_scores if total_scores > 0 else 0.0
        except Exception as e:
            logger.warning(f"Failed to calculate ADMET predictions: {str(e)}")
        
        return jsonify({
            'smiles': smiles,
            'properties': all_properties,
            'analysis_type': 'comprehensive',
            'success': True
        })
        
    except Exception as e:
        logger.error(f"Error in comprehensive analysis: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/batch_file', methods=['POST'])
def batch_file_properties():
    """Analyze properties for molecules from an uploaded file"""
    try:
        from utils.file_parsers import parse_molecule_file
        from flask import send_file
        import io
        import csv

        if 'file' not in request.files:
            return jsonify({'error': 'No file provided', 'success': False}), 400

        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'Empty filename', 'success': False}), 400

        analysis_type = request.form.get('analysis_type', 'all')
        output_format = request.form.get('output_format', 'json')

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
                        'properties': None
                    })
                    continue

                # Calculate properties based on analysis type
                properties = {}

                if analysis_type in ['physicochemical', 'all']:
                    properties.update({
                        'molecular_weight': Descriptors.MolWt(mol),
                        'logp': Descriptors.MolLogP(mol),
                        'tpsa': Descriptors.TPSA(mol),
                        'num_h_donors': Descriptors.NumHDonors(mol),
                        'num_h_acceptors': Descriptors.NumHAcceptors(mol),
                        'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                        'heavy_atom_count': Descriptors.HeavyAtomCount(mol),
                        'num_rings': Descriptors.RingCount(mol),
                        'molecular_formula': Chem.rdMolDescriptors.CalcMolFormula(mol)
                    })

                if analysis_type in ['drug_likeness', 'all']:
                    mw = properties.get('molecular_weight', Descriptors.MolWt(mol))
                    logp = properties.get('logp', Descriptors.MolLogP(mol))
                    hbd = properties.get('num_h_donors', Descriptors.NumHDonors(mol))
                    hba = properties.get('num_h_acceptors', Descriptors.NumHAcceptors(mol))

                    lipinski_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
                    properties['lipinski_violations'] = lipinski_violations
                    properties['passes_lipinski'] = lipinski_violations <= 1

                if analysis_type in ['aromaticity', 'all']:
                    properties['num_aromatic_atoms'] = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
                    properties['num_aromatic_rings'] = Descriptors.NumAromaticRings(mol)

                results.append({
                    'name': mol_data['name'],
                    'smiles': smiles,
                    'properties': properties,
                    'error': None
                })
                num_successful += 1

            except Exception as e:
                logger.warning(f"Error processing molecule {mol_data.get('name', 'unknown')}: {str(e)}")
                results.append({
                    'name': mol_data.get('name', 'unknown'),
                    'smiles': mol_data.get('smiles', ''),
                    'error': str(e),
                    'properties': None
                })

        # Return results based on format
        if output_format == 'csv':
            # Create CSV
            output = io.StringIO()
            if results:
                # Get all property keys
                property_keys = set()
                for r in results:
                    if r['properties']:
                        property_keys.update(r['properties'].keys())

                fieldnames = ['name', 'smiles'] + sorted(list(property_keys)) + ['error']
                writer = csv.DictWriter(output, fieldnames=fieldnames)
                writer.writeheader()

                for r in results:
                    row = {'name': r['name'], 'smiles': r['smiles'], 'error': r.get('error', '')}
                    if r['properties']:
                        row.update(r['properties'])
                    writer.writerow(row)

            output.seek(0)
            return send_file(
                io.BytesIO(output.getvalue().encode()),
                mimetype='text/csv',
                as_attachment=True,
                download_name=f'properties_{file.filename}.csv'
            )
        else:
            # Return JSON
            return jsonify({
                'success': True,
                'file_name': file.filename,
                'analysis_type': analysis_type,
                'num_molecules': len(molecules),
                'num_successful': num_successful,
                'results': results
            })

    except Exception as e:
        logger.error(f"Error in batch file properties: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/qed', methods=['POST'])
def qed_drug_likeness():
    """Calculate QED (Quantitative Estimate of Drug-likeness)"""
    try:
        from rdkit.Chem import QED

        data = request.get_json()
        smiles = data.get('smiles')

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Calculate QED score
        qed_score = QED.qed(mol)

        # Get detailed QED properties
        qed_properties = QED.properties(mol)

        return jsonify({
            'smiles': smiles,
            'qed_score': qed_score,
            'properties': {
                'qed': {
                    'qed': qed_score,
                    'MW': qed_properties.MW,
                    'ALOGP': qed_properties.ALOGP,
                    'HBA': qed_properties.HBA,
                    'HBD': qed_properties.HBD,
                    'PSA': qed_properties.PSA,
                    'ROTB': qed_properties.ROTB,
                    'AROM': qed_properties.AROM,
                    'ALERTS': qed_properties.ALERTS
                },
                'qed_score': qed_score,
                'molecular_weight': qed_properties.MW,
                'logp': qed_properties.ALOGP,
                'num_h_acceptors': int(qed_properties.HBA),
                'num_h_donors': int(qed_properties.HBD),
                'tpsa': qed_properties.PSA,
                'num_rotatable_bonds': int(qed_properties.ROTB),
                'num_aromatic_rings': int(qed_properties.AROM)
            },
            'interpretation': 'High' if qed_score >= 0.7 else 'Medium' if qed_score >= 0.4 else 'Low',
            'success': True
        })

    except Exception as e:
        logger.error(f"Error calculating QED: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/brics', methods=['POST'])
def brics_decomposition():
    """Decompose molecule using BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures)"""
    try:
        from rdkit.Chem import BRICS

        data = request.get_json()
        smiles = data.get('smiles')
        min_fragment_size = data.get('min_fragment_size', 1)

        if not smiles:
            return jsonify({'error': 'SMILES string required', 'success': False}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'success': False}), 400

        # Perform BRICS decomposition
        fragments = list(BRICS.BRICSDecompose(mol, minFragmentSize=min_fragment_size))

        # Get BRICS bonds
        bonds = list(BRICS.FindBRICSBonds(mol))

        fragment_info = []
        for frag_smiles in fragments:
            frag_mol = Chem.MolFromSmiles(frag_smiles)
            if frag_mol:
                fragment_info.append({
                    'smiles': frag_smiles,
                    'num_atoms': frag_mol.GetNumAtoms(),
                    'molecular_weight': Descriptors.MolWt(frag_mol)
                })

        return jsonify({
            'smiles': smiles,
            'num_fragments': len(fragments),
            'fragments': fragments,  # Simple list of SMILES strings for UI
            'fragment_details': fragment_info,  # Detailed info with MW, etc
            'num_brics_bonds': len(bonds),
            'success': True
        })

    except Exception as e:
        logger.error(f"Error with BRICS decomposition: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500

@properties_bp.route('/recap', methods=['POST'])
def recap_decomposition():
    """Decompose molecule using RECAP (Retrosynthetic Combinatorial Analysis Procedure)"""
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

        # Get all unique fragments from the RECAP tree
        fragment_smiles_set = set()
        fragments_list = []

        def get_all_nodes(node):
            if hasattr(node, 'mol') and node.mol:
                frag_smiles = Chem.MolToSmiles(node.mol)
                if frag_smiles not in fragment_smiles_set:
                    fragment_smiles_set.add(frag_smiles)
                    fragments_list.append(frag_smiles)

            # Children is a dictionary, not a list
            if hasattr(node, 'children') and node.children:
                for child_key, child_node in node.children.items():
                    get_all_nodes(child_node)

        get_all_nodes(recap_tree)

        fragment_info = []
        for frag_smiles in fragments_list:
            frag_mol = Chem.MolFromSmiles(frag_smiles)
            if frag_mol:
                fragment_info.append({
                    'smiles': frag_smiles,
                    'num_atoms': frag_mol.GetNumAtoms(),
                    'molecular_weight': Descriptors.MolWt(frag_mol)
                    })

        return jsonify({
            'smiles': smiles,
            'num_fragments': len(fragments_list),
            'fragments': fragments_list,  # Simple list of SMILES strings for UI
            'fragment_details': fragment_info,  # Detailed info with MW, etc
            'success': True
        })

    except Exception as e:
        logger.error(f"Error with RECAP decomposition: {str(e)}")
        return jsonify({'error': str(e), 'success': False}), 500