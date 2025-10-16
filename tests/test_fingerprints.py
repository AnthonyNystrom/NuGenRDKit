import pytest
import json
from app import app

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

class TestFingerprints:
    
    def test_morgan_fingerprint(self, client):
        """Test Morgan fingerprint generation"""
        response = client.post('/api/v1/fingerprints/morgan',
                              json={'smiles': 'CCO', 'radius': 2, 'n_bits': 1024})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'fingerprint' in data
        assert len(data['fingerprint']) == 1024
        assert 'on_bits' in data
    
    def test_rdkit_fingerprint(self, client):
        """Test RDKit fingerprint generation"""
        response = client.post('/api/v1/fingerprints/rdkit',
                              json={'smiles': 'CCO', 'fp_size': 2048})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'fingerprint' in data
        assert len(data['fingerprint']) == 2048
        assert 'on_bits' in data
    
    def test_maccs_fingerprint(self, client):
        """Test MACCS keys fingerprint"""
        response = client.post('/api/v1/fingerprints/maccs',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'fingerprint' in data
        assert len(data['fingerprint']) == 167  # MACCS keys are always 167 bits
        assert data['fp_size'] == 167
    
    def test_avalon_fingerprint(self, client):
        """Test Avalon fingerprint generation"""
        response = client.post('/api/v1/fingerprints/avalon',
                              json={'smiles': 'CCO', 'n_bits': 512})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'fingerprint' in data
        assert len(data['fingerprint']) == 512
    
    def test_atom_pairs_fingerprint(self, client):
        """Test Atom Pairs fingerprint"""
        response = client.post('/api/v1/fingerprints/atom_pairs',
                              json={'smiles': 'CCO', 'n_bits': 1024})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'fingerprint' in data
        assert len(data['fingerprint']) == 1024
    
    def test_topological_torsions_fingerprint(self, client):
        """Test Topological Torsions fingerprint"""
        response = client.post('/api/v1/fingerprints/topological_torsions',
                              json={'smiles': 'CCCC', 'n_bits': 2048})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'fingerprint' in data
        assert len(data['fingerprint']) == 2048
    
    def test_pattern_fingerprint(self, client):
        """Test pattern fingerprint"""
        response = client.post('/api/v1/fingerprints/pattern',
                              json={'smiles': 'CCO', 'fp_size': 1024})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'fingerprint' in data
        assert len(data['fingerprint']) == 1024
    
    def test_layered_fingerprint(self, client):
        """Test layered fingerprint"""
        response = client.post('/api/v1/fingerprints/layered',
                              json={'smiles': 'CCO', 'fp_size': 1024})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'fingerprint' in data
        assert len(data['fingerprint']) == 1024
    
    def test_compare_fingerprints(self, client):
        """Test fingerprint comparison"""
        response = client.post('/api/v1/fingerprints/compare',
                              json={
                                  'smiles1': 'CCO',
                                  'smiles2': 'CCC',
                                  'fp_type': 'morgan',
                                  'radius': 2
                              })
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'similarities' in data
        
        similarities = data['similarities']
        assert 'tanimoto' in similarities
        assert 'dice' in similarities
        assert 'cosine' in similarities
        assert 'sokal' in similarities
        
        # All similarities should be between 0 and 1
        for sim_type, sim_value in similarities.items():
            assert 0 <= sim_value <= 1
    
    def test_compare_identical_molecules(self, client):
        """Test fingerprint comparison of identical molecules"""
        response = client.post('/api/v1/fingerprints/compare',
                              json={
                                  'smiles1': 'CCO',
                                  'smiles2': 'CCO',
                                  'fp_type': 'morgan'
                              })
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        
        # Identical molecules should have similarity of 1.0
        assert data['similarities']['tanimoto'] == 1.0
    
    def test_invalid_smiles(self, client):
        """Test error handling for invalid SMILES"""
        response = client.post('/api/v1/fingerprints/morgan',
                              json={'smiles': 'invalid_smiles'})
        assert response.status_code == 400
        data = json.loads(response.data)
        assert data['success'] == False
    
    def test_missing_smiles(self, client):
        """Test error handling for missing SMILES"""
        response = client.post('/api/v1/fingerprints/morgan',
                              json={})
        assert response.status_code == 400
        data = json.loads(response.data)
        assert data['success'] == False
    
    def test_unsupported_fingerprint_type(self, client):
        """Test error handling for unsupported fingerprint type in comparison"""
        response = client.post('/api/v1/fingerprints/compare',
                              json={
                                  'smiles1': 'CCO',
                                  'smiles2': 'CCC',
                                  'fp_type': 'unsupported_type'
                              })
        assert response.status_code == 400
        data = json.loads(response.data)
        assert data['success'] == False