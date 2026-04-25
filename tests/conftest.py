"""
Shared pytest fixtures.

The existing per-file `client` fixtures in test_*.py continue to work; they
simply override this module-level version when present. New tests should
depend on the fixtures defined here.
"""

from __future__ import annotations

import os
import sys

import pytest

# Add project root to sys.path so `from app import app` works regardless of
# where pytest is invoked. Phase 2 will make this unnecessary by fixing the
# package layout.
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from app import app as flask_app  # noqa: E402


def pytest_configure(config):
    """Register custom markers so `pytest --strict-markers` passes."""
    config.addinivalue_line(
        "markers",
        "slow: marks tests that exercise long-running RDKit operations "
        "(3D embeds, conformer search, /all descriptors). "
        "Run with `pytest -m slow`.",
    )


@pytest.fixture
def app():
    """Flask app configured for testing."""
    flask_app.config.update({"TESTING": True})
    yield flask_app


@pytest.fixture
def client(app):
    """Flask test client."""
    return app.test_client()


@pytest.fixture
def aspirin():
    return "CC(=O)OC1=CC=CC=C1C(=O)O"


@pytest.fixture
def caffeine():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def ethanol():
    return "CCO"
