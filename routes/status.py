"""
Service-status probe blueprint — Phase 7 polish.

NuGenRDKit's only "service dependency" is RDKit itself (everything is
local — no DB, no upstream API). The probes here are mostly diagnostic;
they let an operator quickly answer "is the install healthy?" and
"how much disk has the upload area chewed through?" without needing
shell access.

Endpoints (all GET, no body):

    /api/status/rdkit    → RDKit version + roundtrip probe
    /api/status/disk     → free / total bytes for the static + tmp dirs
    /api/status          → both, plus uptime since worker boot

The Phase-1 `/health` endpoint stays (it's the contract for k8s
liveness / docker HEALTHCHECK and returns a stable JSON shape). These
new endpoints are richer and meant for human consumption.
"""

from __future__ import annotations

import logging
import os
import shutil
import time

from flask import Blueprint, jsonify

from utils.errors import json_err

logger = logging.getLogger(__name__)
status_bp = Blueprint("status", __name__)

# Process boot time — captured at import. Worker recycling
# (gunicorn.conf.py) resets this every `max_requests`.
_BOOT_TS = time.time()

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# --------------------------------------------------------------------------- #
# Probes
# --------------------------------------------------------------------------- #


def _rdkit_probe() -> dict:
    """Roundtrip a SMILES through RDKit and report the result."""
    try:
        import rdkit
        from rdkit import Chem

        start = time.monotonic()
        mol = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")  # aspirin
        canonical = Chem.MolToSmiles(mol) if mol is not None else None
        latency_ms = int((time.monotonic() - start) * 1000)

        ok = canonical == "CC(=O)Oc1ccccc1C(=O)O"
        return {
            "ok": ok,
            "version": rdkit.__version__,
            "roundtrip_ms": latency_ms,
            "canonical_aspirin": canonical,
        }
    except Exception as exc:
        logger.exception("RDKit probe failed")
        return {"ok": False, "error": str(exc)}


def _disk_probe() -> dict:
    """Report free + total bytes for the static + tmp dirs."""
    paths = {
        "project": _PROJECT_ROOT,
        "static": os.path.join(_PROJECT_ROOT, "static"),
        "tmp": os.environ.get("TMPDIR", "/tmp"),
    }
    out = {}
    for name, path in paths.items():
        try:
            usage = shutil.disk_usage(path)
            out[name] = {
                "path": path,
                "total_bytes": usage.total,
                "free_bytes": usage.free,
                "free_pct": round(100 * usage.free / usage.total, 1) if usage.total else 0.0,
            }
        except Exception as exc:
            out[name] = {"path": path, "error": str(exc)}
    return out


def _uptime_seconds() -> float:
    return time.time() - _BOOT_TS


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #


@status_bp.route("/rdkit")
def status_rdkit():
    probe = _rdkit_probe()
    status_code = 200 if probe.get("ok") else 503
    return jsonify({"success": probe.get("ok", False), "rdkit": probe}), status_code


@status_bp.route("/disk")
def status_disk():
    return jsonify({"success": True, "disk": _disk_probe()})


@status_bp.route("/")
def status_all():
    rdkit = _rdkit_probe()
    return (
        jsonify(
            {
                "success": rdkit.get("ok", False),
                "uptime_seconds": round(_uptime_seconds(), 1),
                "rdkit": rdkit,
                "disk": _disk_probe(),
            }
        ),
        200 if rdkit.get("ok") else 503,
    )


@status_bp.route("", methods=["GET"])
def status_all_no_slash():
    """Both /api/status and /api/status/ map to the combined probe."""
    return status_all()


# Keep the import line satisfied even if a future refactor changes the
# error envelope import; we currently don't return error envelopes from
# this blueprint but the helper is here for parity with the others.
_ = json_err
