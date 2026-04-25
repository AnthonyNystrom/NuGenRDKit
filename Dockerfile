# syntax=docker/dockerfile:1.7
# ============================================================================
# NuGenRDKit production image — Phase 6.
#
# Multi-stage:
#   1. `builder`  — Debian slim with build tooling; installs Python deps
#                   into an isolated /opt/venv that we copy out.
#   2. `runtime`  — same base, no compilers; copies /opt/venv + app code,
#                   runs as a non-root user, exposes 8000, gunicorn entry.
#
# Why slim and not python:3.11-alpine? RDKit ships compiled wheels for
# manylinux but not musl, so alpine would force a source build (10+ min
# every CI run) without a meaningful size win.
# ============================================================================

ARG PYTHON_VERSION=3.11

# --- Stage 1: builder -------------------------------------------------------
FROM python:${PYTHON_VERSION}-slim AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_ROOT_USER_ACTION=ignore

# Build deps for any wheel that has to compile (some matplotlib backends,
# pillow image codecs).
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential \
        libfreetype6-dev \
        libpng-dev \
        libjpeg-dev \
        zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

WORKDIR /build
COPY requirements.txt ./
RUN pip install --upgrade pip wheel \
 && pip install -r requirements.txt

# --- Stage 2: runtime -------------------------------------------------------
FROM python:${PYTHON_VERSION}-slim AS runtime

# `.dockerignore` strips the `.git` directory from the build context, so
# the runtime container has no way to derive the commit SHA at request
# time. Pass `--build-arg GIT_SHA=$(git rev-parse HEAD)` (or the CI's
# ${{ github.sha }}) to bake the deploying commit into /health output.
ARG GIT_SHA=""

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    LOG_LEVEL=INFO \
    FLASK_ENV=production \
    HOST=0.0.0.0 \
    PORT=8000 \
    GIT_SHA=${GIT_SHA} \
    PATH="/opt/venv/bin:${PATH}"

# Runtime libraries that the wheels link against (no compilers).
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        libfreetype6 \
        libpng16-16 \
        libjpeg62-turbo \
        libgomp1 \
        curl \
        git \
 && rm -rf /var/lib/apt/lists/*

# Non-root app user. UID/GID 10001 avoids collisions with host accounts.
RUN groupadd --system --gid 10001 app \
 && useradd --system --uid 10001 --gid app --home /app --shell /usr/sbin/nologin app

# Copy the virtualenv + app code as the app user.
COPY --from=builder /opt/venv /opt/venv

WORKDIR /app
COPY --chown=app:app . /app

# Drop SELinux/AppArmor labels won't matter inside the container; what
# does matter is that gunicorn binds 8000 and runs as `app` (no root).
USER app

EXPOSE 8000

# Health probe that mirrors the Phase-1 /health endpoint contract.
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD curl -fsS http://127.0.0.1:8000/health || exit 1

# Phase-2 cap on conformers + Phase-1 rate-limit defaults all live inside
# utils/config.py. Override at deploy time via -e.
CMD ["gunicorn", "--config", "gunicorn.conf.py", "app:app"]
