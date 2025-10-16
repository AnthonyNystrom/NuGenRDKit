import pytest
import json
from app import app

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

class TestMolecularStructure:
    
    def test_smiles_to_mol(self, client):
        """Test SMILES to MOL conversion"""
        response = client.post('/api/v1/structure/smiles_to_mol',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'mol_block' in data
    
    def test_mol_to_smiles(self, client):
        """Test MOL to SMILES conversion"""
        mol_block = """
  Mrv2311 01012023000000000000  

  3  2  0  0  0  0            999 V2000
   -0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
"""
        response = client.post('/api/v1/structure/mol_to_smiles',
                              json={'mol_block': mol_block})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'smiles' in data
        assert 'canonical_smiles' in data
    
    def test_canonicalize_smiles(self, client):
        """Test SMILES canonicalization"""
        response = client.post('/api/v1/structure/canonicalize',
                              json={'smiles': 'OCC'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert data['canonical_smiles'] == 'CCO'
    
    def test_validate_structure_valid(self, client):
        """Test structure validation with valid SMILES"""
        response = client.post('/api/v1/structure/validate',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert data['is_valid'] == True
        assert 'molecular_formula' in data
    
    def test_validate_structure_invalid(self, client):
        """Test structure validation with invalid SMILES"""
        response = client.post('/api/v1/structure/validate',
                              json={'smiles': 'CCC[invalid]'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert data['is_valid'] == False
    
    def test_smiles_to_inchi(self, client):
        """Test SMILES to InChI conversion"""
        response = client.post('/api/v1/structure/smiles_to_inchi',
                              json={'smiles': 'CCO'})
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['success'] == True
        assert 'inchi' in data
        assert 'inchi_key' in data
    
    def test_missing_smiles(self, client):
        """Test error handling for missing SMILES"""
        response = client.post('/api/v1/structure/canonicalize',
                              json={})
        assert response.status_code == 400
        data = json.loads(response.data)
        assert data['success'] == False
        assert 'error' in data
    
    def test_invalid_smiles(self, client):
        """Test error handling for invalid SMILES"""
        response = client.post('/api/v1/structure/canonicalize',
                              json={'smiles': 'invalid_smiles'})
        assert response.status_code == 400
        data = json.loads(response.data)
        assert data['success'] == False