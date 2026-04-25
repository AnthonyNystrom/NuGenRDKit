"""
OpenAPI spec + Swagger UI for /api/v1/**.

Phase 7 polish.

Why hand-rolled instead of flask-smorest?
-----------------------------------------
flask-smorest gives you marshmallow-validated request/response schemas
for free, but only after you wrap each route with `@blp.arguments(…)`
+ `@blp.response(…)` decorators. With ~100 endpoints across 8
blueprints, that's 100+ files of work and a real risk of breaking
the JSON shapes the Phase-0 oracle pins.

This module instead reflects on `app.url_map` at request time and
emits a spec with:
  • path + HTTP method (auto-discovered)
  • summary + description (from the view function's docstring)
  • a generic `application/json` body schema
  • response 200/400/500 with the standard envelope

It's less type-safe than schema-validated specs, but it ships in
one file, can't break existing routes, and gives clients a
browsable contract today. A future PR can incrementally add
schemas per route.
"""

from __future__ import annotations

import re
from typing import Any

from flask import Blueprint, current_app, jsonify, render_template_string

openapi_bp = Blueprint("openapi", __name__)

_INFO = {
    "title": "NuGenRDKit API",
    "description": (
        "Comprehensive cheminformatics REST API powered by RDKit. "
        "Every endpoint accepts JSON via POST and returns the standard "
        "envelope `{success, ...}` on success or "
        "`{success: false, error, request_id, code}` on failure."
    ),
    "version": "1.0.0",
    "license": {"name": "MIT", "url": "https://opensource.org/licenses/MIT"},
    "contact": {"name": "NuGenRDKit", "url": "https://github.com/AnthonyNystrom/NuGenRDKit"},
}

# Tag descriptions show up in Swagger UI's collapsible sections.
_TAGS = [
    {
        "name": "structure",
        "description": "SMILES ↔ MOL ↔ InChI ↔ SDF ↔ PDB conversions, validation, standardization.",
    },
    {"name": "descriptors", "description": "Lipinski, LogP, topological, VSA, 217+ descriptors."},
    {
        "name": "fingerprints",
        "description": "Morgan, RDKit, MACCS, Avalon, atom pairs, torsions, pattern, layered.",
    },
    {
        "name": "similarity",
        "description": "13 metrics + bulk search + substructure + diversity selection.",
    },
    {
        "name": "coordinates",
        "description": "2D/3D generation, force-field optimization, conformer search, alignment.",
    },
    {"name": "properties", "description": "Drug-likeness, ADMET, QED, BRICS/RECAP, scaffolds."},
    {
        "name": "reactions",
        "description": "SMARTS reactions, library enumeration, reaction fingerprints.",
    },
    {
        "name": "visualization",
        "description": "SVG/PNG drawing, similarity maps, fingerprint bit visualization.",
    },
    {"name": "diagnostics", "description": "Health + status probes (not part of /api/v1)."},
]


# --------------------------------------------------------------------------- #
# Spec generation
# --------------------------------------------------------------------------- #


def _docstring_summary(view) -> tuple[str, str]:
    """Split a docstring into a one-line summary + a multi-line description."""
    if not view or not view.__doc__:
        return ("", "")
    doc = view.__doc__.strip()
    parts = doc.split("\n", 1)
    summary = parts[0].strip()
    description = parts[1].strip() if len(parts) > 1 else ""
    return summary, description


def _tag_for_path(path: str) -> str:
    """Map a URL like /api/v1/structure/canonicalize → 'structure'."""
    m = re.match(r"^/api/v1/([^/]+)/?", path)
    if m:
        return m.group(1)
    if path.startswith("/api/status"):
        return "diagnostics"
    return "misc"


_GENERIC_BODY = {
    "type": "object",
    "description": (
        'Endpoint-specific JSON body. Most endpoints take {"smiles": '
        '"<SMILES string>"} plus optional parameters. See the '
        "view-function docstring for details."
    ),
    "additionalProperties": True,
    "example": {"smiles": "CCO"},
}

_OK_ENVELOPE = {
    "description": "Operation succeeded.",
    "content": {
        "application/json": {
            "schema": {
                "type": "object",
                "required": ["success"],
                "properties": {
                    "success": {"type": "boolean", "example": True},
                },
                "additionalProperties": True,
            }
        }
    },
}

_ERR_ENVELOPE = {
    "description": "Validation or processing error.",
    "content": {
        "application/json": {
            "schema": {
                "type": "object",
                "required": ["success", "error"],
                "properties": {
                    "success": {"type": "boolean", "example": False},
                    "error": {"type": "string"},
                    "request_id": {"type": "string"},
                    "code": {"type": "string"},
                },
            }
        }
    },
}

_RATE_LIMIT_ENVELOPE = {
    "description": "Rate limit exceeded. Retry after the indicated interval.",
    "headers": {
        "Retry-After": {"schema": {"type": "string"}},
    },
    "content": {
        "application/json": {
            "schema": {
                "type": "object",
                "properties": {
                    "success": {"type": "boolean", "example": False},
                    "error": {"type": "string", "example": "Rate limit exceeded"},
                    "code": {"type": "string", "example": "rate_limited"},
                    "retry_after": {"type": "string"},
                },
            }
        }
    },
}


def _path_item(rule, view) -> dict[str, Any]:
    """Build the OpenAPI path-item for a single Flask rule."""
    summary, description = _docstring_summary(view)
    tag = _tag_for_path(rule.rule)
    methods = sorted(rule.methods - {"HEAD", "OPTIONS"})
    item: dict[str, Any] = {}

    for method in methods:
        op: dict[str, Any] = {
            "tags": [tag],
            "summary": summary or rule.endpoint.replace("_", " ").title(),
            "operationId": rule.endpoint.replace(".", "_") + "_" + method.lower(),
            "responses": {
                "200": _OK_ENVELOPE,
                "400": _ERR_ENVELOPE,
                "429": _RATE_LIMIT_ENVELOPE,
                "500": _ERR_ENVELOPE,
            },
        }
        if description:
            op["description"] = description
        if method in {"POST", "PUT", "PATCH"}:
            op["requestBody"] = {
                "required": False,
                "content": {"application/json": {"schema": _GENERIC_BODY}},
            }
        item[method.lower()] = op

    return item


def build_spec() -> dict[str, Any]:
    """Build the full OpenAPI 3.0 document by reflecting on the live app."""
    paths: dict[str, dict] = {}
    for rule in current_app.url_map.iter_rules():
        # Only document the public API surface.
        if not (rule.rule.startswith("/api/v1/") or rule.rule.startswith("/api/status")):
            continue
        # Static endpoint isn't an API.
        if rule.endpoint == "static":
            continue
        view = current_app.view_functions.get(rule.endpoint)
        # Translate Flask rule syntax to OpenAPI's: /export/<format> → /export/{format}
        spec_path = re.sub(r"<(?:\w+:)?(\w+)>", r"{\1}", rule.rule)
        item = _path_item(rule, view)
        if spec_path in paths:
            # Multiple Flask rules can share a path (different methods).
            paths[spec_path].update(item)
        else:
            paths[spec_path] = item

    return {
        "openapi": "3.0.3",
        "info": _INFO,
        "servers": [{"url": "/", "description": "This server"}],
        "tags": _TAGS,
        "paths": paths,
        "components": {
            "schemas": {},
        },
    }


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #


@openapi_bp.route("/openapi.json")
def openapi_spec():
    return jsonify(build_spec())


# Swagger UI shipped from a pinned + SRI'd CDN. The HTML is small enough
# to live inline; keeping it here means we don't add another template.
_SWAGGER_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>NuGenRDKit API — Swagger UI</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5.17.14/swagger-ui.css"
          integrity="sha384-wxLW6kwyHktdDGr6Pv1zgm/VGJh99lfUbzSn6HNHBENZlCN7W602k9VkGdxuFvPn"
          crossorigin="anonymous">
    <link rel="icon" href="{{ url_for('static', filename='favicon.svg') }}">
    <style>
        html, body { margin: 0; padding: 0; }
        body { background: #fafafa; }
        .topbar { display: none; }
    </style>
</head>
<body>
    <div id="swagger-ui"></div>
    <script src="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5.17.14/swagger-ui-bundle.js"
            integrity="sha384-wmyclcVGX/WhUkdkATwhaK1X1JtiNrr2EoYJ+diV3vj4v6OC5yCeSu+yW13SYJep"
            crossorigin="anonymous"></script>
    <script>
        window.onload = function () {
            window.ui = SwaggerUIBundle({
                url: "{{ spec_url }}",
                dom_id: "#swagger-ui",
                deepLinking: true,
                docExpansion: "list",
                tagsSorter: "alpha",
                operationsSorter: "alpha",
                tryItOutEnabled: true,
            });
        };
    </script>
</body>
</html>"""


@openapi_bp.route("/docs")
def swagger_ui():
    return render_template_string(
        _SWAGGER_TEMPLATE,
        spec_url="/api/v1/openapi.json",
    )
