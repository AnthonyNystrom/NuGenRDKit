"""
Central Flask configuration.

`configure_app(app)` is the single place where app-level concerns live:
- `.env` loading
- `MAX_CONTENT_LENGTH` + `SECRET_KEY`
- CORS allowlist (closed by default in production)
- Rate-limit tiers (default / expensive / bulk)
- Security headers + CSP
- Cache-busting static URL parameters (`?v=<mtime>`)
- Structured request logging

Env vars (see .env.example):
    HOST, PORT                     bind address (dev only; prod via gunicorn)
    LOG_LEVEL                      default INFO
    FLASK_ENV                      "development" relaxes CORS to "*"
    SECRET_KEY                     required in prod; random in dev
    CORS_ORIGINS                   comma-separated allowlist
    MAX_UPLOAD_MB                  default 50
    RATE_LIMIT_DEFAULT             default "100/minute;1000/hour"
    RATE_LIMIT_EXPENSIVE           default "10/minute;100/hour"
    RATE_LIMIT_BULK                default "5/minute;50/hour"
"""

from __future__ import annotations

import contextlib
import os
import secrets
import warnings

from flask import Flask
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

from utils.logging import configure_root_logger, install_request_logging

# Exported so blueprints can apply tiered decorators: @limiter.limit(EXPENSIVE)
limiter: Limiter = Limiter(key_func=get_remote_address)


def _load_dotenv_if_available() -> None:
    """Load `.env` into os.environ without requiring python-dotenv."""
    try:
        from dotenv import load_dotenv
    except ImportError:
        return
    load_dotenv()


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


def _env_csv(name: str) -> list[str]:
    raw = os.environ.get(name, "")
    return [item.strip() for item in raw.split(",") if item.strip()]


# --------------------------------------------------------------------------- #
# Rate-limit tiers — exported for blueprints
# --------------------------------------------------------------------------- #


def rate_limit_default() -> str:
    return os.environ.get("RATE_LIMIT_DEFAULT", "100/minute;1000/hour")


def rate_limit_expensive() -> str:
    """3D embed, conformer search, MMFF optimize, bulk similarity."""
    return os.environ.get("RATE_LIMIT_EXPENSIVE", "10/minute;100/hour")


def rate_limit_bulk() -> str:
    """File uploads that may contain thousands of molecules."""
    return os.environ.get("RATE_LIMIT_BULK", "5/minute;50/hour")


# --------------------------------------------------------------------------- #
# CSP
# --------------------------------------------------------------------------- #

# Phase 5 dropped 'unsafe-inline' from script-src — every page now loads
# its logic from external <script src=…> files; inline `onclick=` was
# migrated to NuGenUtils.registerAction.
#
# Phase 7 attempted to drop 'unsafe-inline' from style-src too, but
# RDKit's server-side SVG output emits hundreds of `<rect style='fill:
# #FFFFFF;...'>` and `<path style='fill:none;stroke:#000000;...'>`
# elements per molecule. CSP rejects those, leaving every depicted
# molecule with a black background. We control RDKit's output server-
# side so the security risk is minimal — keep style-src 'unsafe-inline'
# for the SVG fill case. Our own templates + JS are still inline-style
# free (Phase 7 audit removed all 18 occurrences).
#
# CDN allowances: jsdelivr (3Dmol, swagger-ui), unpkg (lucide), cdnjs.
_CSP = (
    "default-src 'self'; "
    "script-src 'self' "
    "https://cdn.jsdelivr.net "
    "https://cdnjs.cloudflare.com "
    "https://unpkg.com "
    "https://3Dmol.org; "
    "style-src 'self' 'unsafe-inline' "
    "https://cdn.jsdelivr.net "
    "https://fonts.googleapis.com; "
    "font-src 'self' data: https://fonts.gstatic.com https://cdnjs.cloudflare.com; "
    "img-src 'self' data: blob: https:; "
    # connect-src allows source-map fetches for the CDN-hosted scripts
    # we ship (lucide on unpkg, 3Dmol on jsdelivr). Source maps are
    # loaded via XHR by DevTools and fall under connect-src, not
    # script-src — leaving these off causes a CSP-violation message
    # in the console any time the user opens DevTools, even though
    # the app itself works fine.
    "connect-src 'self' https://unpkg.com https://cdn.jsdelivr.net; "
    # 3Dmol.js spawns Web Workers from blob: URLs.
    "worker-src 'self' blob:; "
    "object-src 'none'; "
    "base-uri 'self'; "
    "frame-ancestors 'none'"
)


# --------------------------------------------------------------------------- #
# Configuration entry point
# --------------------------------------------------------------------------- #


def configure_app(app: Flask) -> None:
    """Apply every cross-cutting concern to `app` in a deterministic order."""
    _load_dotenv_if_available()
    configure_root_logger()

    # --- Core Flask config -------------------------------------------------- #
    env = os.environ.get("FLASK_ENV", "production").lower()
    app.config.setdefault("ENV", env)
    app.config.setdefault(
        "MAX_CONTENT_LENGTH",
        _env_int("MAX_UPLOAD_MB", 50) * 1024 * 1024,
    )
    app.config.setdefault(
        "SECRET_KEY",
        os.environ.get("SECRET_KEY") or secrets.token_hex(32),
    )
    app.config.setdefault("JSON_SORT_KEYS", False)
    # Flask-Limiter reads these and also checks for the RATELIMIT_* variants.
    app.config.setdefault("RATELIMIT_DEFAULT", rate_limit_default())
    app.config.setdefault("RATELIMIT_STORAGE_URI", "memory://")
    app.config.setdefault("RATELIMIT_HEADERS_ENABLED", True)

    # --- Testing short-circuit --------------------------------------------- #
    # Phase-0 tests set app.config['TESTING'] = True. Rate limits + CSP are
    # kept in place so tests catch drift, but we skip bits that would make
    # tests non-deterministic (e.g. random SECRET_KEY reset is OK here).
    testing = app.config.get("TESTING", False)

    # --- CORS --------------------------------------------------------------- #
    origins = _env_csv("CORS_ORIGINS")
    if env == "development" and not origins:
        origins = ["*"]  # dev convenience — match legacy behavior
    if not origins:
        # Production with no allowlist → no cross-origin access permitted.
        origins = []
    CORS(app, resources={r"/api/*": {"origins": origins}}, supports_credentials=False)

    # --- Rate limiter ------------------------------------------------------- #
    # Suppress the in-memory storage warning we've already decided about
    # (see deployment docs: single-worker gunicorn is expected).
    warnings.filterwarnings("ignore", message=".*in-memory storage.*")
    limiter.init_app(app)
    # Apply the default to all routes. Blueprint-level overrides (expensive /
    # bulk) are added in Phase 2 via shared decorators on individual routes.

    # --- Structured request logging ---------------------------------------- #
    install_request_logging(app)

    # --- Security headers + CSP -------------------------------------------- #
    @app.after_request
    def _security_headers(response):
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("X-Frame-Options", "DENY")
        response.headers.setdefault("Referrer-Policy", "strict-origin-when-cross-origin")
        response.headers.setdefault(
            "Permissions-Policy", "geolocation=(), microphone=(), camera=()"
        )
        response.headers.setdefault("Content-Security-Policy", _CSP)
        # HSTS only when we know we're behind HTTPS (set by reverse proxy).
        if not testing and env != "development":
            response.headers.setdefault(
                "Strict-Transport-Security",
                "max-age=63072000; includeSubDomains",
            )
        return response

    # --- Cache-busting static URLs ----------------------------------------- #
    @app.url_defaults
    def _static_cache_bust(endpoint, values):
        if endpoint == "static" and "filename" in values and "v" not in values:
            static_folder = app.static_folder
            if not static_folder:
                return
            path = os.path.join(static_folder, values["filename"])
            with contextlib.suppress(OSError):
                values["v"] = int(os.path.getmtime(path))
