"""
RDKit helper functions shared across route modules.

This module is the single source of truth for SMILES/SMARTS parsing,
fingerprint generation, and CSV escaping. It exists so individual route
handlers don't reimplement the same `Chem.MolFromSmiles(s); if mol is None:
return error` pattern dozens of times.

Phase-0 tests pin every JSON response shape — these helpers must keep their
return values byte-for-byte compatible with the inline code they replace.
"""

from __future__ import annotations

import functools
import logging
from collections.abc import Callable
from typing import Any

from rdkit import Chem
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect

from utils.errors import RDKitRuntimeError, UserInputError, json_err

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

# Hard caps; anything larger is almost certainly an attack/copy-paste mistake.
# Phase 1 added MAX_CONTENT_LENGTH at the Flask layer; this is the per-field cap.
MAX_SMILES_LEN = 4096
MAX_SMARTS_LEN = 4096


def parse_smiles(raw: Any, *, max_len: int = MAX_SMILES_LEN, field: str = "smiles") -> Chem.Mol:
    """
    Validate + parse a SMILES string. Raises UserInputError (→ HTTP 400) on
    any failure; never returns None.

    >>> mol = parse_smiles("CCO")
    >>> mol.GetNumAtoms()
    3
    """
    if raw is None or raw == "":
        raise UserInputError(f"{field} is required")
    if not isinstance(raw, str):
        raise UserInputError(f"{field} must be a string")
    if len(raw) > max_len:
        raise UserInputError(f"{field} exceeds maximum length ({max_len} chars)")
    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        raise UserInputError(f"Invalid {field} string")
    return mol


def parse_smarts(raw: Any, *, max_len: int = MAX_SMARTS_LEN, field: str = "smarts") -> Chem.Mol:
    if raw is None or raw == "":
        raise UserInputError(f"{field} is required")
    if not isinstance(raw, str):
        raise UserInputError(f"{field} must be a string")
    if len(raw) > max_len:
        raise UserInputError(f"{field} exceeds maximum length ({max_len} chars)")
    pattern = Chem.MolFromSmarts(raw)
    if pattern is None:
        raise UserInputError(f"Invalid {field} pattern")
    return pattern


def parse_optional_smiles(
    raw: Any, *, max_len: int = MAX_SMILES_LEN, field: str = "smiles"
) -> Chem.Mol | None:
    """Same as parse_smiles, but returns None when the input is empty."""
    if raw is None or raw == "":
        return None
    return parse_smiles(raw, max_len=max_len, field=field)


# --------------------------------------------------------------------------- #
# Similarity-input normalization (back-compat aliases)
# --------------------------------------------------------------------------- #


def get_query_target_smiles(data: dict) -> tuple[str, str]:
    """
    The similarity routes historically accepted multiple param names for the
    same value. Phase 2 keeps every alias working but reads them through one
    place so the `smiles1` double-check typo is fixed once.

    Accepted keys (in order of precedence):
        query  → query_smiles → smiles1
        target → target_smiles → smiles2
    """
    if not isinstance(data, dict):
        raise UserInputError("Request body must be a JSON object")
    query = data.get("query") or data.get("query_smiles") or data.get("smiles1")
    target = data.get("target") or data.get("target_smiles") or data.get("smiles2")
    if not query or not target:
        raise UserInputError("Both query and target SMILES required")
    return query, target


# --------------------------------------------------------------------------- #
# Fingerprints — single source of truth
# --------------------------------------------------------------------------- #

# Names accepted by the public API. Extending this list is a public-API change.
SUPPORTED_FP_TYPES = (
    "morgan",
    "rdkit",
    "maccs",
    "avalon",
    "atom_pairs",
    "topological_torsions",
    "pattern",
    "layered",
)


def make_fingerprint(
    mol: Chem.Mol,
    kind: str = "morgan",
    *,
    radius: int = 2,
    n_bits: int = 2048,
):
    """
    Build a fingerprint of the requested kind. Centralizes the dispatch so
    routes don't ship six identical if/elif ladders.

    Raises UserInputError on unknown `kind` so callers can let it propagate
    via @rdkit_route to a 400.
    """
    if mol is None:
        raise UserInputError("Cannot fingerprint a null molecule")
    kind = (kind or "morgan").lower()
    if kind == "morgan":
        return GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    if kind == "rdkit":
        return Chem.RDKFingerprint(mol, fpSize=n_bits)
    if kind == "maccs":
        return rdMolDescriptors.GetMACCSKeysFingerprint(mol)
    if kind == "avalon":
        return pyAvalonTools.GetAvalonFP(mol, nBits=n_bits)
    if kind == "atom_pairs":
        return rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits)
    if kind == "topological_torsions":
        return rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=n_bits)
    if kind == "pattern":
        return Chem.PatternFingerprint(mol, fpSize=n_bits)
    if kind == "layered":
        return Chem.LayeredFingerprint(mol, fpSize=n_bits)
    raise UserInputError(f"Unsupported fingerprint type: {kind}")


# --------------------------------------------------------------------------- #
# CSV-safe value formatting
# --------------------------------------------------------------------------- #

_CSV_INJECTION_PREFIXES = ("=", "+", "-", "@", "\t", "\r")


def csv_safe(value: Any) -> str:
    """
    Defang spreadsheet formula-injection. When a downloaded CSV is opened in
    Excel/Sheets, any cell starting with =+−@ is treated as a formula and
    can phone home / exfiltrate data. Prefixing with a single quote
    neutralizes the formula without changing the visible value.
    """
    if value is None:
        return ""
    s = str(value)
    if s and s[0] in _CSV_INJECTION_PREFIXES:
        return "'" + s
    return s


# --------------------------------------------------------------------------- #
# @rdkit_route decorator — uniform error handling
# --------------------------------------------------------------------------- #


def rdkit_route(view: Callable[..., Any]) -> Callable[..., Any]:
    """
    Wrap a Flask view so domain exceptions become consistent JSON envelopes:

        UserInputError       → 400 {success: False, error: "..."}
        RDKitRuntimeError    → 400 {success: False, error: "..."}
        Anything else        → 500 {success: False, error: "Internal error"}

    The full traceback is logged with the request_id so support can find it
    in logs, but is never returned to the client.

    Usage:
        @blueprint.route("/foo", methods=["POST"])
        @rdkit_route
        def foo():
            mol = parse_smiles(request.get_json().get("smiles"))
            return json_ok({"answer": Chem.MolToSmiles(mol)})
    """

    @functools.wraps(view)
    def wrapper(*args: Any, **kwargs: Any):
        try:
            return view(*args, **kwargs)
        except UserInputError as exc:
            # User input errors carry their own user-friendly text.
            return json_err(str(exc), status=400, code="user_input_error")
        except RDKitRuntimeError as exc:
            # Phase-7: route the underlying RDKit exception text through
            # the friendly_errors mapper so clients see actionable
            # messages instead of "Sanitization error: …".
            from utils.friendly_errors import friendly

            message, code = friendly(exc)
            return json_err(message, status=400, code=code)
        except Exception as exc:
            logger.exception("Unhandled error in %s", view.__name__)
            from utils.friendly_errors import friendly

            message, code = friendly(exc)
            # Surface the friendly version to the user but keep the
            # 500 status code — this path is for unexpected failures.
            return json_err(message, status=500, code=code)

    return wrapper
