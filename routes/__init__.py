"""
Blueprint registration factory.

`register_blueprints(app)` mounts every /api/v1/** handler. Keeping the list
here (instead of in app.py) means Phase 2 can add new helpers without
touching the top-level app file.
"""

from __future__ import annotations

from flask import Flask

from routes.coordinates import coordinates_bp
from routes.descriptors import descriptors_bp
from routes.fingerprints import fingerprints_bp
from routes.molecular_structure import molecular_structure_bp
from routes.openapi import openapi_bp
from routes.properties import properties_bp
from routes.reactions import reactions_bp
from routes.similarity import similarity_bp
from routes.status import status_bp
from routes.visualization import visualization_bp

_BLUEPRINTS = [
    (molecular_structure_bp, "/api/v1/structure"),
    (descriptors_bp, "/api/v1/descriptors"),
    (fingerprints_bp, "/api/v1/fingerprints"),
    (similarity_bp, "/api/v1/similarity"),
    (coordinates_bp, "/api/v1/coordinates"),
    (properties_bp, "/api/v1/properties"),
    (reactions_bp, "/api/v1/reactions"),
    (visualization_bp, "/api/v1/visualization"),
    # Phase 7 polish: human-readable diagnostic probes. Mounted under
    # /api/status (NOT under /api/v1) so the inventory drift detector
    # for /api/v1/** ignores them.
    (status_bp, "/api/status"),
    # OpenAPI 3.0 spec + Swagger UI (Phase 7). Lives under /api/v1/ so
    # the schema discovery is one nested level alongside the endpoints
    # it documents. The drift detector exempts these in test_api_smoke.
    (openapi_bp, "/api/v1"),
]


def register_blueprints(app: Flask) -> None:
    for bp, prefix in _BLUEPRINTS:
        app.register_blueprint(bp, url_prefix=prefix)
