"""
Phase-0 error-envelope contract. Sampled endpoints must produce the documented
JSON error shape on invalid input so Phase 1's `json_err(...)` refactor can
guarantee byte-for-byte parity.
"""

from __future__ import annotations

import pytest

# Sample endpoints spanning all 8 blueprints — one per domain is enough to
# catch envelope drift introduced by a refactor.
_SAMPLE_ENDPOINTS = [
    "/api/v1/structure/canonicalize",
    "/api/v1/descriptors/basic",
    "/api/v1/fingerprints/morgan",
    "/api/v1/similarity/tanimoto",
    "/api/v1/coordinates/generate_2d",
    "/api/v1/properties/drug_likeness",
    "/api/v1/visualization/draw_svg",
]


@pytest.mark.parametrize("url", _SAMPLE_ENDPOINTS)
def test_invalid_smiles_returns_error_envelope(client, url):
    response = client.post(url, json={"smiles": "ThisIsNotValidSMILES!!!"})
    assert response.status_code == 400, (
        f"POST {url} with bad SMILES returned {response.status_code}"
    )
    data = response.get_json()
    assert data is not None
    assert data.get("success") is False
    assert "error" in data, f"Missing 'error' key in envelope: {data}"
    assert isinstance(data["error"], str) and data["error"], (
        f"'error' must be a non-empty string: {data}"
    )


@pytest.mark.parametrize("url", _SAMPLE_ENDPOINTS)
def test_missing_payload_returns_error(client, url):
    """
    POST with empty JSON body. The app's before_request validate_json hook
    (app.py:100-109) catches empty-body JSON and returns a 400 envelope; some
    routes accept an empty dict and then fail SMILES parsing downstream — in
    either case the contract is `success=False` + `error`.
    """
    response = client.post(url, json={})
    assert response.status_code == 400
    data = response.get_json()
    assert data is not None
    assert data.get("success") is False
    assert "error" in data


@pytest.mark.parametrize("url", _SAMPLE_ENDPOINTS)
def test_malformed_json_body_returns_error(client, url):
    response = client.post(url, data="{not json}", content_type="application/json")
    assert response.status_code == 400
    data = response.get_json()
    assert data is not None
    assert data.get("success") is False


def test_similarity_missing_query_returns_error(client):
    """Similarity endpoints require a query molecule."""
    response = client.post(
        "/api/v1/similarity/tanimoto",
        json={"target_smiles": "CCO"},
    )
    assert response.status_code == 400
    data = response.get_json()
    assert data.get("success") is False
    assert "error" in data


def test_reactions_run_without_smarts_returns_error(client):
    response = client.post("/api/v1/reactions/run_reaction", json={"reactant_smiles": ["CCO"]})
    assert response.status_code == 400
    data = response.get_json()
    assert data.get("success") is False


def test_health_is_reachable_without_body(client):
    """Regression guard: /health must not require a JSON body."""
    response = client.get("/health")
    assert response.status_code == 200
    data = response.get_json()
    assert data.get("success") is True
