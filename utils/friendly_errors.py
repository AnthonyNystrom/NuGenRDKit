"""
RDKit / file-parsing error → human-readable message mapping.

When `@rdkit_route` (utils/rdkit_helpers.py) catches a RDKit exception
or one of our own UserInputError / RDKitRuntimeError instances, the
default JSON envelope shows `error: str(exc)`. RDKit's native messages
are aimed at developers ("Sanitization error: Can't kekulize mol")
and are confusing to end users.

This module turns those into one-line, action-oriented sentences that
get returned in the same `error` field. A `code` discriminator is also
attached so frontend code can branch on the cause without parsing
prose.

Usage:

    from utils.friendly_errors import friendly

    try:
        mol = Chem.MolFromSmiles(raw)
    except Exception as exc:
        message, code = friendly(exc, raw=raw)
        return json_err(message, code=code)
"""

from __future__ import annotations

import re
from typing import Any

# --------------------------------------------------------------------------- #
# Cause classification
# --------------------------------------------------------------------------- #

# Each entry: (substring or regex, friendly_message_template, error_code).
# First match wins. Templates can reference {raw} (the original input)
# and {detail} (the original RDKit message). Use plain text — the
# Phase-1 envelope is JSON, no HTML.
_RULES: list[tuple[Any, str, str]] = [
    # SMILES / SMARTS parse errors
    (
        re.compile(r"can't kekulize|sanitization", re.I),
        "RDKit could not sanitize the molecule — check that aromatic rings "
        "have alternating single/double bonds and valences are valid.",
        "sanitization_failed",
    ),
    (
        re.compile(r"explicit valence.*greater than permitted", re.I),
        "An atom carries more bonds than its valence allows. Common causes: "
        "missing H specifier, incorrect charge, or wrong atom type.",
        "explicit_valence",
    ),
    (
        re.compile(r"unclosed ring|unmatched (\(|\))", re.I),
        "SMILES has unclosed ring or parenthesis. Re-check brackets and ring closure digits.",
        "smiles_syntax",
    ),
    (
        re.compile(r"smiles parse error|invalid smiles", re.I),
        "SMILES could not be parsed. Verify it is well-formed (e.g. CCO, c1ccccc1, CC(=O)O).",
        "smiles_invalid",
    ),
    (
        re.compile(r"invalid smarts|smarts parse error", re.I),
        "SMARTS pattern could not be parsed. Verify atom and bond tokens are "
        "well-formed (e.g. [#6]~[#7]).",
        "smarts_invalid",
    ),
    # 3D / conformer
    (
        re.compile(r"embed.*fail|could not generate 3d|no embedding", re.I),
        "RDKit could not embed the molecule in 3D. Try simplifying the "
        "structure or generating 2D coordinates first.",
        "embed_failed",
    ),
    (
        re.compile(r"force.?field.*could not|mmff.*fail|uff.*fail", re.I),
        "Force-field optimisation failed. The molecule may contain atoms "
        "MMFF94/UFF cannot parametrise — try a different force field.",
        "ff_failed",
    ),
    # Reaction
    (
        re.compile(r"no registered converter|reactionfromrxn", re.I),
        "RDKit could not parse the reaction. Verify the RXN file format "
        "or SMARTS reaction pattern.",
        "reaction_invalid",
    ),
    # Substructure / similarity
    (
        re.compile(r"no sub.?structure match", re.I),
        "The probe molecule does not share a substructure with the reference. "
        "Pick a probe that contains a fragment of the reference.",
        "no_substructure_match",
    ),
    # File parsing
    (
        re.compile(r"invalid (sdf|mol|pdb)", re.I),
        "The uploaded molecular file could not be parsed. Verify it is a "
        "valid SDF / MOL / PDB block (RDKit's MolFromMolBlock parser is "
        "strict about formatting).",
        "molblock_invalid",
    ),
    (
        re.compile(r"file size exceeds|too large", re.I),
        "Uploaded file exceeds the configured size limit. See MAX_UPLOAD_MB in your environment.",
        "payload_too_large",
    ),
    # Network / timeout
    (
        re.compile(r"timeout|timed out", re.I),
        "Operation timed out. Try a smaller input or fewer conformers.",
        "timeout",
    ),
]


_FALLBACK = (
    "RDKit reported an error processing your input. Original message: {detail}",
    "rdkit_error",
)


def friendly(exc: BaseException, *, raw: str | None = None) -> tuple[str, str]:
    """
    Return ``(message, code)`` for an exception.

    ``message`` is plain ASCII safe to embed in a JSON envelope.
    ``code`` is a short stable token a frontend can branch on.

    Never raises. Falls back to a generic wrapper if no rule matches.
    """
    detail = str(exc) or exc.__class__.__name__
    for matcher, template, code in _RULES:
        if isinstance(matcher, str):
            if matcher.lower() in detail.lower():
                return _format(template, detail, raw), code
        else:
            if matcher.search(detail):
                return _format(template, detail, raw), code
    return _format(_FALLBACK[0], detail, raw), _FALLBACK[1]


def _format(template: str, detail: str, raw: str | None) -> str:
    return template.format(detail=_truncate(detail, 200), raw=_truncate(raw or "", 80))


def _truncate(value: str, limit: int) -> str:
    if len(value) <= limit:
        return value
    return value[: limit - 1] + "…"
