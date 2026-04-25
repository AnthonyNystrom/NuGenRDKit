"""
Phase-0 API smoke: POST every /api/v1/** endpoint with its inventory payload;
assert success and the documented response keys.

This is the oracle for Phase 1 (backend hardening) and Phase 2 (code-quality
refactor): neither phase may change the response shape of any endpoint
listed in tests/fixtures/inventory.py.
"""

from __future__ import annotations

import io

import pytest

from tests.fixtures.inventory import (
    ALL_ENDPOINTS,
    HTML_DOC_ROUTES,
    JSON_GET_ROUTES,
    OTHER_JSON_GET_ROUTES,
    UPLOAD_ENDPOINTS,
    Endpoint,
    UploadEndpoint,
    fast_endpoints,
    fast_uploads,
    non_alias_endpoints,
)


def _ids_for(endpoints: list[Endpoint]) -> list[str]:
    return [e.url + (" [alias]" if e.alias else "") for e in endpoints]


def _upload_ids_for(uploads: list[UploadEndpoint]) -> list[str]:
    return [u.url for u in uploads]


# --- Fast (default `pytest`) ---------------------------------------------- #


@pytest.mark.parametrize(
    "endpoint",
    fast_endpoints(),
    ids=_ids_for(fast_endpoints()),
)
def test_endpoint_happy_path(client, endpoint):
    response = client.post(endpoint.url, json=endpoint.payload)
    assert response.status_code == 200, (
        f"POST {endpoint.url} returned {response.status_code}: {response.data[:400]!r}"
    )
    data = response.get_json()
    assert data is not None, f"POST {endpoint.url} did not return JSON"
    assert data.get("success") is True, f"POST {endpoint.url} success!=True; response: {data}"
    missing = [k for k in endpoint.expect_keys if k not in data]
    assert not missing, f"POST {endpoint.url} missing required keys {missing}; response: {data}"


@pytest.mark.parametrize(
    "upload",
    fast_uploads(),
    ids=_upload_ids_for(fast_uploads()),
)
def test_upload_happy_path(client, upload):
    data = {upload.file_field: (io.BytesIO(upload.file_body.encode()), upload.filename)}
    data.update(upload.form)
    response = client.post(upload.url, data=data, content_type="multipart/form-data")
    assert response.status_code == 200, (
        f"POST {upload.url} (multipart) returned {response.status_code}: {response.data[:400]!r}"
    )
    body = response.get_json()
    assert body is not None, f"POST {upload.url} did not return JSON"
    assert body.get("success") is True, f"POST {upload.url} success!=True; response: {body}"


# --- Aliases resolve to the same handler --------------------------------- #


@pytest.mark.parametrize(
    "endpoint",
    [e for e in ALL_ENDPOINTS if e.alias],
    ids=[e.url for e in ALL_ENDPOINTS if e.alias],
)
def test_alias_resolves_to_same_handler(client, endpoint):
    """Aliases (e.g. /2d vs /generate_2d) should return 200 with success=True."""
    response = client.post(endpoint.url, json=endpoint.payload)
    assert response.status_code == 200, f"Alias {endpoint.url} returned {response.status_code}"
    assert response.get_json().get("success") is True


# --- Slow paths only under `pytest -m slow` ------------------------------ #

_SLOW = [e for e in non_alias_endpoints() if e.slow and not e.known_broken]
_SLOW_UPLOADS = [u for u in UPLOAD_ENDPOINTS if u.slow and not u.known_broken]


@pytest.mark.slow
@pytest.mark.parametrize("endpoint", _SLOW, ids=_ids_for(_SLOW))
def test_slow_endpoint_happy_path(client, endpoint):
    response = client.post(endpoint.url, json=endpoint.payload)
    assert response.status_code == 200, (
        f"[slow] POST {endpoint.url} -> {response.status_code}: {response.data[:400]!r}"
    )
    assert response.get_json().get("success") is True


@pytest.mark.slow
@pytest.mark.parametrize("upload", _SLOW_UPLOADS, ids=_upload_ids_for(_SLOW_UPLOADS))
def test_slow_upload_happy_path(client, upload):
    data = {upload.file_field: (io.BytesIO(upload.file_body.encode()), upload.filename)}
    data.update(upload.form)
    response = client.post(upload.url, data=data, content_type="multipart/form-data")
    assert response.status_code == 200, (
        f"[slow] POST {upload.url} -> {response.status_code}: {response.data[:400]!r}"
    )
    assert response.get_json().get("success") is True


# --- Known-broken endpoints ---------------------------------------------- #

_KNOWN_BROKEN = [e for e in ALL_ENDPOINTS if e.known_broken]
_KNOWN_BROKEN_UPLOADS = [u for u in UPLOAD_ENDPOINTS if u.known_broken]


# Only register the parametrized test when there is at least one known-broken
# entry; otherwise pytest emits a confusing "got empty parameter set" skip.
if _KNOWN_BROKEN:

    @pytest.mark.parametrize(
        "endpoint",
        _KNOWN_BROKEN,
        ids=[e.url for e in _KNOWN_BROKEN],
    )
    def test_known_broken_endpoint_still_broken(client, endpoint):
        """
        Documented bugs not yet fixed. If these start PASSING, someone fixed
        them without updating the inventory — remove the `known_broken`
        marker and move the entry into the happy-path suite.
        """
        response = client.post(endpoint.url, json=endpoint.payload)
        assert response.status_code >= 400, (
            f"{endpoint.url} no longer returns an error — remove the "
            f"known_broken marker in tests/fixtures/inventory.py. "
            f"(Marker was: {endpoint.known_broken})"
        )


if _KNOWN_BROKEN_UPLOADS:

    @pytest.mark.parametrize(
        "upload",
        _KNOWN_BROKEN_UPLOADS,
        ids=[u.url for u in _KNOWN_BROKEN_UPLOADS],
    )
    def test_known_broken_upload_still_broken(client, upload):
        data = {upload.file_field: (io.BytesIO(upload.file_body.encode()), upload.filename)}
        data.update(upload.form)
        response = client.post(upload.url, data=data, content_type="multipart/form-data")
        assert response.status_code >= 400, (
            f"{upload.url} no longer returns an error — remove the "
            f"known_broken marker. (Marker was: {upload.known_broken})"
        )


# --- Coverage check: every registered /api/v1/** route is in the inventory- #


def test_every_api_route_is_inventoried(client):
    """
    Guard against drift: if someone adds a new /api/v1/** route without
    registering it in tests/fixtures/inventory.py, this test fails and forces
    them to add smoke coverage before merging.
    """
    from app import app as flask_app

    known_urls = {e.url for e in ALL_ENDPOINTS}
    known_urls.update(u.url for u in UPLOAD_ENDPOINTS)
    # GET-only routes are covered by test_routes_ok.py.
    known_urls.update(JSON_GET_ROUTES)
    known_urls.update(OTHER_JSON_GET_ROUTES)
    known_urls.update(HTML_DOC_ROUTES)

    registered = set()
    for rule in flask_app.url_map.iter_rules():
        if not rule.rule.startswith("/api/v1/"):
            continue
        if "<" in rule.rule:
            # Parameterized routes (/api/v1/coordinates/export/<format>)
            # are tested explicitly; skip the auto-discovery.
            continue
        registered.add(rule.rule)

    missing = sorted(registered - known_urls)
    assert not missing, (
        "These routes are registered in the app but missing from "
        "tests/fixtures/inventory.py:\n  " + "\n  ".join(missing)
    )
