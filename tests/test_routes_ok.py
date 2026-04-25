"""
Phase-0 route smoke: every page URL returns 200 + well-formed HTML, and every
JSON GET route returns JSON. Oracle for Phase 5 page migration — URLs and
response contracts must survive every refactor.
"""

from __future__ import annotations

import pytest

from tests.fixtures.inventory import (
    HTML_DOC_ROUTES,
    JSON_GET_ROUTES,
    OTHER_JSON_GET_ROUTES,
    PAGES,
)


@pytest.mark.parametrize("url", PAGES + HTML_DOC_ROUTES)
def test_page_returns_html(client, url):
    response = client.get(url)
    assert response.status_code == 200, (
        f"GET {url} returned {response.status_code}: {response.data[:200]!r}"
    )
    body = response.data
    assert b"<!DOCTYPE html>" in body or b"<html" in body, f"GET {url} did not return HTML"
    assert b"</html>" in body, f"GET {url} HTML was truncated"


@pytest.mark.parametrize("url", JSON_GET_ROUTES)
def test_json_get_returns_success(client, url):
    response = client.get(url)
    assert response.status_code == 200, (
        f"GET {url} returned {response.status_code}: {response.data[:200]!r}"
    )
    data = response.get_json()
    assert data is not None, f"GET {url} did not return JSON"
    assert data.get("success") is True, f"GET {url} response missing success=True: {data}"


@pytest.mark.parametrize("url", OTHER_JSON_GET_ROUTES)
def test_other_json_get_returns_json(client, url):
    """
    Routes that return JSON in their own format (e.g. OpenAPI). We only
    assert reachability + parseable JSON, not the success envelope.
    """
    response = client.get(url)
    assert response.status_code == 200, (
        f"GET {url} returned {response.status_code}: {response.data[:200]!r}"
    )
    data = response.get_json()
    assert data is not None, f"GET {url} did not return JSON"


def test_unknown_endpoint_returns_404(client):
    response = client.get("/this/path/does/not/exist")
    assert response.status_code == 404
    data = response.get_json()
    assert data is not None, "404 response must be JSON (see app.py:79-82)"
    assert data.get("success") is False


def test_method_not_allowed_returns_405(client):
    """POST on the index page returns the standard 405 envelope."""
    response = client.post("/")
    assert response.status_code == 405
    data = response.get_json()
    assert data is not None
    assert data.get("success") is False
