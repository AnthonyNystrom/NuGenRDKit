"""
Structured request logging + per-request correlation IDs.

`install_request_logging(app)` hooks before/after_request so every log line
carries the method, path, status, duration, and a request ID that is also
included in error envelopes (see utils/errors.py::json_err). Makes it
possible to correlate a user's bug report with a server-side stack trace.
"""

from __future__ import annotations

import logging
import os
import time
import uuid

from flask import Flask, g, request

_JSON_FMT = (
    '{"ts":"%(asctime)s","level":"%(levelname)s","logger":"%(name)s","message":%(message_json)s}'
)


class _JsonFormatter(logging.Formatter):
    """Lightweight JSON formatter with no extra deps."""

    def format(self, record: logging.LogRecord) -> str:
        # Escape the message for JSON embedding. `logging.Formatter.format`
        # handles `%(message)s` interpolation; we wrap the result in JSON.
        message = super().format(record)
        record.message_json = _as_json_string(message)
        return _JSON_FMT % record.__dict__


def _as_json_string(value: str) -> str:
    import json

    return json.dumps(value, ensure_ascii=False)


def configure_root_logger(level: str | None = None) -> None:
    """
    Configure the root logger. Called from utils.config.configure_app so
    the app never shadows the caller's logging setup.
    """
    resolved = (level or os.environ.get("LOG_LEVEL", "INFO")).upper()
    logging.basicConfig(
        level=getattr(logging, resolved, logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    # Silence Flask's werkzeug INFO chatter; keep warnings.
    logging.getLogger("werkzeug").setLevel(logging.WARNING)


def install_request_logging(app: Flask) -> None:
    """
    Add before_request + after_request hooks that populate `g.request_id`,
    time the handler, and log a single structured line per request at INFO.
    Error envelopes pick up `g.request_id` automatically.
    """
    logger = logging.getLogger("nugenrdkit.request")

    @app.before_request
    def _start_request():
        g.request_id = request.headers.get("X-Request-ID") or uuid.uuid4().hex[:16]
        g.request_start = time.monotonic()

    @app.after_request
    def _finish_request(response):
        duration_ms = int((time.monotonic() - getattr(g, "request_start", time.monotonic())) * 1000)
        # Propagate the id back to the client so they can quote it.
        response.headers.setdefault("X-Request-ID", getattr(g, "request_id", ""))
        # Skip noisy static-asset entries at INFO — log them at DEBUG.
        log = logger.debug if request.path.startswith("/static/") else logger.info
        log(
            "%s %s -> %d in %dms [req_id=%s]",
            request.method,
            request.path,
            response.status_code,
            duration_ms,
            getattr(g, "request_id", "-"),
        )
        return response
