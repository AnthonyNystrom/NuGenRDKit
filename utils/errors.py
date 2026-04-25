"""
Standard response envelopes + exception classes.

Every /api/v1/** handler eventually returns one of:

    {"success": True, ...data...}
    {"success": False, "error": "human readable", "request_id": "..."}

`json_ok` / `json_err` keep that contract consistent across blueprints so
Phase-0 tests keep passing as Phase 2 migrates route handlers to
`@rdkit_route`.
"""

from __future__ import annotations

from typing import Any

from flask import g, jsonify


class UserInputError(ValueError):
    """
    The caller sent something the app cannot interpret (bad SMILES, bad
    JSON shape, missing required field, oversized input). Always maps to
    HTTP 400 via `json_err`.
    """


class RDKitRuntimeError(RuntimeError):
    """
    RDKit itself failed on what looked like valid input (embedding timed
    out, optimizer diverged, LazyPick blew up). Maps to 400 by default so
    the client sees a human-readable message instead of a 500.
    """


def json_ok(data: dict[str, Any] | None = None, *, status: int = 200, **extra: Any):
    """
    Serialize a success envelope. Keys from `data` and `extra` are merged
    at the top level so existing response shapes (e.g. `{"success": True,
    "mol_block": ...}`) are preserved byte-for-byte.
    """
    body: dict[str, Any] = {"success": True}
    if data:
        body.update(data)
    if extra:
        body.update(extra)
    return jsonify(body), status


def json_err(message: str, *, status: int = 400, code: str | None = None, **extra: Any):
    """
    Serialize an error envelope. Includes the request ID set by
    `utils.logging.install_request_logging()` when available so clients
    can quote it in bug reports.
    """
    body: dict[str, Any] = {"success": False, "error": message}
    request_id = getattr(g, "request_id", None)
    if request_id:
        body["request_id"] = request_id
    if code:
        body["code"] = code
    if extra:
        body.update(extra)
    return jsonify(body), status
