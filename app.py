"""
NuGenRDKit — Flask application entry point.

Cross-cutting concerns (logging, CORS, rate limits, security headers, CSP,
cache-busting, request IDs) live in `utils.config.configure_app`. This file
only wires the app together: page routes, blueprint registration, error
handlers, and a meaningful health check.
"""

from __future__ import annotations

import logging
import os
import platform
import subprocess
import warnings

# Silence RDKit's import-time deprecation chatter before importing rdkit.
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=".*MorganGenerator.*")
warnings.filterwarnings("ignore", message=".*DEPRECATION.*")

import rdkit  # noqa: E402
from flask import Flask, jsonify, render_template, request  # noqa: E402
from rdkit import RDLogger  # noqa: E402

from routes import register_blueprints  # noqa: E402
from utils.config import configure_app  # noqa: E402
from utils.errors import json_err  # noqa: E402

# Quiet RDKit's native logger unless the caller explicitly asked for verbose.
RDLogger.logger().setLevel(RDLogger.ERROR)


app = Flask(__name__)
configure_app(app)
register_blueprints(app)

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# JSON body validation (POSTs to /api/v1/** with application/json)
# --------------------------------------------------------------------------- #


@app.before_request
def _validate_json_body():
    if request.method != "POST":
        return None
    if (request.content_type or "").startswith("application/json"):
        try:
            if not request.get_json(silent=True):
                return json_err("Invalid or empty JSON", status=400)
        except Exception as exc:
            logger.warning("JSON parsing error: %s", exc)
            return json_err("Invalid JSON format", status=400)
    return None


# --------------------------------------------------------------------------- #
# Error handlers — unified envelope
# --------------------------------------------------------------------------- #


@app.errorhandler(400)
def _bad_request(error):
    logger.warning("Bad request: %s", error)
    return json_err("Bad request", status=400)


@app.errorhandler(404)
def _not_found(error):
    logger.warning("Not found: %s", error)
    return json_err("Endpoint not found", status=404)


@app.errorhandler(405)
def _method_not_allowed(error):
    logger.warning("Method not allowed: %s", error)
    return json_err("Method not allowed", status=405)


@app.errorhandler(413)
def _payload_too_large(error):
    logger.warning("Payload too large: %s", error)
    max_mb = app.config.get("MAX_CONTENT_LENGTH", 0) // (1024 * 1024)
    return json_err(
        f"Payload exceeds maximum size ({max_mb} MB)",
        status=413,
        code="payload_too_large",
    )


@app.errorhandler(429)
def _ratelimit(e):
    logger.warning("Rate limit exceeded: %s", e)
    retry_after = getattr(e, "retry_after", None)
    extra = {"retry_after": str(retry_after)} if retry_after else {}
    return json_err("Rate limit exceeded", status=429, code="rate_limited", **extra)


@app.errorhandler(500)
def _server_error(error):
    logger.exception("Internal server error: %s", error)
    return json_err("Internal server error", status=500, code="internal_error")


# --------------------------------------------------------------------------- #
# HTML page routes
# --------------------------------------------------------------------------- #

# Endpoint names must stay stable — templates call url_for('index'),
# url_for('structure'), etc. Breaking these breaks every page.
_PAGES = (
    ("index", "/", "index.html"),
    ("structure", "/structure", "structure.html"),
    ("descriptors", "/descriptors", "descriptors.html"),
    ("fingerprints", "/fingerprints", "fingerprints.html"),
    ("similarity", "/similarity", "similarity.html"),
    ("coordinates", "/coordinates", "coordinates.html"),
    ("properties", "/properties", "properties.html"),
    ("reactions", "/reactions", "reactions.html"),
    ("visualization", "/visualization", "visualization.html"),
)


def _make_page_view(endpoint: str, template_name: str):
    def _view():
        try:
            return render_template(template_name)
        except Exception as exc:
            logger.exception("Error rendering template %s: %s", template_name, exc)
            return json_err(f"Template error: {exc}", status=500)

    _view.__name__ = endpoint
    return _view


for _endpoint, _path, _template in _PAGES:
    app.add_url_rule(_path, endpoint=_endpoint, view_func=_make_page_view(_endpoint, _template))


# --------------------------------------------------------------------------- #
# API index + health
# --------------------------------------------------------------------------- #


@app.route("/api")
def api_info():
    return jsonify(
        {
            "success": True,
            "message": "RDKit Flask API",
            "version": "1.0.0",
            "endpoints": {
                "structure": "/api/v1/structure",
                "descriptors": "/api/v1/descriptors",
                "fingerprints": "/api/v1/fingerprints",
                "similarity": "/api/v1/similarity",
                "coordinates": "/api/v1/coordinates",
                "properties": "/api/v1/properties",
                "reactions": "/api/v1/reactions",
                "visualization": "/api/v1/visualization",
            },
        }
    )


def _git_sha() -> str:
    """Best-effort short SHA; empty string when unavailable.

    Resolution order:
      1. `GIT_SHA` env var — explicit override; used by Docker images
         where the build pipeline bakes the deploying commit at build
         time (.dockerignore strips .git from the runtime context).
      2. `git rev-parse --short HEAD` — local checkout fallback.
      3. Empty string when neither is available.
    """
    env = os.environ.get("GIT_SHA", "").strip()
    if env:
        # Match the local-checkout format: 7-char short sha when a full
        # 40-char hash was passed in (e.g. ${{ github.sha }} from CI).
        return env[:7] if len(env) >= 40 else env
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
            cwd=os.path.dirname(os.path.abspath(__file__)),
            timeout=1,
        )
        return out.decode().strip()
    except Exception:
        return ""


@app.route("/health")
def health():
    """Liveness + dependency versions. Phase 6 deploy docs reference this."""
    from rdkit import Chem

    try:
        mol = Chem.MolFromSmiles("C")
        rdkit_ok = mol is not None
    except Exception:
        rdkit_ok = False

    try:
        import flask as _flask

        flask_version = _flask.__version__
    except Exception:
        flask_version = "unknown"

    return jsonify(
        {
            "success": True,
            "status": "healthy" if rdkit_ok else "degraded",
            "rdkit_working": rdkit_ok,
            "versions": {
                "rdkit": rdkit.__version__,
                "flask": flask_version,
                "python": platform.python_version(),
            },
            "git_sha": _git_sha(),
        }
    ), (200 if rdkit_ok else 503)


if __name__ == "__main__":
    host = os.environ.get("HOST", "0.0.0.0")
    port = int(os.environ.get("PORT", "8000"))
    debug = os.environ.get("FLASK_ENV", "").lower() == "development"
    app.run(host=host, port=port, debug=debug)
