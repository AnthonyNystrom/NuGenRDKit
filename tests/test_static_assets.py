"""
Phase-0 static-asset sanity: every file referenced from templates/JS exists
on disk, and the Flask static handler serves it with 200.
"""

from __future__ import annotations

import os

import pytest

from tests.fixtures.inventory import STATIC_ASSETS

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.mark.parametrize("relative_path", STATIC_ASSETS)
def test_static_file_exists_on_disk(relative_path):
    path = os.path.join(_PROJECT_ROOT, relative_path)
    assert os.path.isfile(path), f"Missing static file: {relative_path}"
    assert os.path.getsize(path) > 0, f"Empty static file: {relative_path}"


@pytest.mark.parametrize("relative_path", STATIC_ASSETS)
def test_static_file_served_by_flask(client, relative_path):
    """
    /static/<file> must serve with 200 for every referenced asset. Phase 1
    will add cache-busting `?v=<mtime>`; that query-string is ignored here
    so the test keeps passing.
    """
    assert relative_path.startswith("static/"), relative_path
    url = "/" + relative_path
    response = client.get(url)
    assert response.status_code == 200, f"GET {url} returned {response.status_code}"
    assert response.data, f"GET {url} returned empty body"
