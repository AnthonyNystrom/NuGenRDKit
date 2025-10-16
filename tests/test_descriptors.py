import pytest
import json
from app import app

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

class TestDescriptors:
    
    def test_basic_descriptors(self, client):
        """Test basic descriptor calculation"""
        response = client.post('/api/v1/descriptors/basic',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'descriptors' in data
        
        descriptors = data['descriptors']
        assert 'molecular_weight' in descriptors
        assert 'heavy_atom_count' in descriptors
        assert 'num_rotatable_bonds' in descriptors
        
        # Test specific values for ethanol (CCO)
        assert descriptors['heavy_atom_count'] == 3
        assert abs(descriptors['molecular_weight'] - 46.069) < 0.01
    
    def test_lipinski_descriptors(self, client):
        """Test Lipinski rule descriptors"""
        response = client.post('/api/v1/descriptors/lipinski',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'descriptors' in data
        
        descriptors = data['descriptors']
        assert 'lipinski_violations' in descriptors
        assert 'passes_lipinski' in descriptors
        assert 'num_h_donors' in descriptors
        assert 'num_h_acceptors' in descriptors
        
        # Ethanol should pass Lipinski
        assert descriptors['passes_lipinski'] == True
    
    def test_logp_descriptors(self, client):
        """Test LogP descriptors"""
        response = client.post('/api/v1/descriptors/logp',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'descriptors' in data
        
        descriptors = data['descriptors']
        assert 'wildman_crippen_logp' in descriptors
        assert 'molar_refractivity' in descriptors
    
    def test_topological_descriptors(self, client):
        """Test topological descriptors"""
        response = client.post('/api/v1/descriptors/topological',
                              json={'smiles': 'c1ccccc1'})  # benzene
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'descriptors' in data
        
        descriptors = data['descriptors']
        assert 'balaban_j' in descriptors
        assert 'bertz_ct' in descriptors
        assert 'chi0' in descriptors
    
    def test_all_descriptors(self, client):
        """Test all descriptors calculation"""
        response = client.post('/api/v1/descriptors/all',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'descriptors' in data
        assert 'num_descriptors' in data
        
        # Should have many descriptors (200+)
        assert data['num_descriptors'] > 200
    
    def test_list_descriptors(self, client):
        """Test listing available descriptors"""
        response = client.get('/api/v1/descriptors/list')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'descriptors' in data
        assert 'count' in data
        assert data['count'] > 200
    
    def test_vsa_descriptors(self, client):
        """Test VSA descriptors"""
        response = client.post('/api/v1/descriptors/vsa',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'descriptors' in data
        
        descriptors = data['descriptors']
        assert 'tpsa' in descriptors
        assert 'labute_asa' in descriptors
        
        # Should have VSA descriptors
        vsa_descriptors = [k for k in descriptors.keys() if 'VSA' in k]
        assert len(vsa_descriptors) > 0
    
    def test_invalid_smiles(self, client):
        """Test error handling for invalid SMILES"""
        response = client.post('/api/v1/descriptors/basic',
                              json={'smiles': 'invalid_smiles'})
        assert response.status_code == 400
        data = json.loads(response.data)
        assert data['success'] == False
    
    def test_missing_smiles(self, client):
        """Test error handling for missing SMILES"""
        response = client.post('/api/v1/descriptors/basic',
                              json={})
        assert response.status_code == 400
        data = json.loads(response.data)
        assert data['success'] == False