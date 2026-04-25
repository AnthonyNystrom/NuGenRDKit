"""
Gunicorn configuration for NuGenRDKit.

Phase 6 deployment.

Single-worker constraint
------------------------
The Phase-1 rate limiter uses Flask-Limiter with `RATELIMIT_STORAGE_URI =
"memory://"` (decided in the upgrade plan: in-memory, document the
constraint, no Redis dep). With more than one worker each gunicorn
process maintains its own counter, so a 100/min limit becomes
N×100/min — which silently breaks the contract.

If/when you scale horizontally:
  1. Add Redis to docker-compose.yml
  2. Set RATELIMIT_STORAGE_URI=redis://redis:6379/0
  3. Bump `workers` here

Worker timeout
--------------
3D embed + conformer search routes can legitimately take 60+ seconds
(see Phase-2 cap of 200 conformers). Keep the gunicorn timeout above
the per-route budget so the worker isn't killed mid-RDKit call.

Logging
-------
We let gunicorn log to stdout/stderr; Phase-1 utils/logging.py adds
structured per-request lines on top. In Docker / k8s these stream to
the container collector.
"""

from __future__ import annotations

import os


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


# --- Bind ---------------------------------------------------------------- #

bind = os.environ.get("GUNICORN_BIND", f"0.0.0.0:{_env_int('PORT', 8000)}")

# --- Workers ------------------------------------------------------------- #

# Single sync worker by default — see "Single-worker constraint" above.
# Override via GUNICORN_WORKERS once an external rate-limit store is in
# place. With sync workers, threads handle concurrent slow-client cases.
workers = _env_int("GUNICORN_WORKERS", 1)
threads = _env_int("GUNICORN_THREADS", 4)
worker_class = os.environ.get("GUNICORN_WORKER_CLASS", "sync")

# 3 minutes — matches the longest legitimate RDKit operation (conformer
# search at MAX_CONFORMERS=200 with optimization).
timeout = _env_int("GUNICORN_TIMEOUT", 180)
graceful_timeout = _env_int("GUNICORN_GRACEFUL_TIMEOUT", 30)
keepalive = _env_int("GUNICORN_KEEPALIVE", 5)

# Worker recycling stops slow memory growth (RDKit holds C++ structures).
max_requests = _env_int("GUNICORN_MAX_REQUESTS", 2000)
max_requests_jitter = _env_int("GUNICORN_MAX_REQUESTS_JITTER", 200)

# Reject oversized request lines / headers early (defence in depth on top
# of the Flask MAX_CONTENT_LENGTH).
limit_request_line = 8190
limit_request_field_size = 16380

# --- Logging ------------------------------------------------------------- #

accesslog = os.environ.get("GUNICORN_ACCESS_LOG", "-")  # "-" → stdout
errorlog = os.environ.get("GUNICORN_ERROR_LOG", "-")
loglevel = os.environ.get("GUNICORN_LOG_LEVEL", "info")
access_log_format = '%(h)s %({X-Request-ID}o)s "%(r)s" %(s)s %(b)s "%(f)s" "%(a)s" %(L)ss'

# --- Process naming ------------------------------------------------------ #

proc_name = "nugenrdkit"

# --- Forwarded-protocol trust ------------------------------------------- #
# When running behind nginx / a load balancer the X-Forwarded-* headers
# carry the real scheme + remote addr. Set FORWARDED_ALLOW_IPS to the
# proxy's IP (or "*" if it's on a private network you trust).
forwarded_allow_ips = os.environ.get("FORWARDED_ALLOW_IPS", "127.0.0.1")


def _print_banner():
    """Sanity log on boot so it's obvious what mode we started in."""
    print(
        f"[gunicorn] bind={bind} workers={workers} threads={threads} "
        f"timeout={timeout}s class={worker_class} "
        f"max_requests={max_requests} (+/-{max_requests_jitter})"
    )


# Gunicorn calls these hooks if they exist.
def on_starting(_server):  # noqa: D401
    _print_banner()
