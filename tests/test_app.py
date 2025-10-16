import pytest
import json
from app import app

@pytest.fixture
def client():
    """Create test client"""
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

def test_index(client):
    """Test index endpoint - returns HTML"""
    response = client.get('/')
    assert response.status_code == 200
    assert b'<!DOCTYPE html>' in response.data or b'<html' in response.data

def test_api_info(client):
    """Test API info endpoint - returns JSON"""
    response = client.get('/api')
    assert response.status_code == 200
    data = json.loads(response.data)
    assert 'message' in data
    assert 'endpoints' in data

def test_health(client):
    """Test health endpoint"""
    response = client.get('/health')
    assert response.status_code == 200
    data = json.loads(response.data)
    assert data['status'] == 'healthy'
    assert data['rdkit_working'] == True

def test_invalid_endpoint(client):
    """Test invalid endpoint"""
    response = client.get('/nonexistent')
    assert response.status_code == 404
    data = json.loads(response.data)
    assert data['success'] == False
    assert 'error' in data

def test_method_not_allowed(client):
    """Test method not allowed"""
    response = client.post('/')
    assert response.status_code == 405
    data = json.loads(response.data)
    assert data['success'] == False

def test_invalid_json(client):
    """Test invalid JSON"""
    response = client.post('/api/v1/structure/smiles_to_mol',
                          data='invalid json',
                          content_type='application/json')
    assert response.status_code == 400
    data = json.loads(response.data)
    assert data['success'] == False