#!/usr/bin/env python3
"""
Playwright visual baseline — Phase 0 safety net.

Captures a PNG of every page in `tests/fixtures/inventory.PAGES` against a
running instance of the app, and either writes it to
`tests/visual/baseline/` (first run / `--update`) or diffs against the
existing baseline (every subsequent run).

Usage:

    # 1. Start the app in a separate shell on port 8001 (not 8000, so we
    #    don't collide with dev):
    FLASK_APP=app.py PORT=8001 python app.py  # (Phase 1 will use $PORT)

    # 2a. First time: capture the baseline.
    APP_URL=http://127.0.0.1:8001 python tests/visual/take_baseline.py --update

    # 2b. Thereafter: diff against baseline and fail on regression.
    APP_URL=http://127.0.0.1:8001 python tests/visual/take_baseline.py

Install (once):

    pip install playwright pillow
    playwright install chromium

Phases 3 + 5 will INTENTIONALLY change the baseline. Rebase the baseline
with `--update` after a human has reviewed the diff PNGs.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

_THIS = Path(__file__).resolve()
_PROJECT_ROOT = _THIS.parent.parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from tests.fixtures.inventory import PAGES  # noqa: E402

BASELINE_DIR = _THIS.parent / "baseline"
DIFF_DIR = _THIS.parent / "diff"
VIEWPORT = {"width": 1440, "height": 900}


def _page_id(url: str) -> str:
    if url == "/":
        return "index"
    return url.strip("/").replace("/", "_")


def _capture(base_url: str, out_dir: Path) -> list[tuple[str, Path]]:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        sys.exit(
            "playwright not installed. Run:\n"
            "  pip install playwright pillow\n"
            "  playwright install chromium"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[tuple[str, Path]] = []

    with sync_playwright() as p:
        browser = p.chromium.launch()
        context = browser.new_context(viewport=VIEWPORT)
        page = context.new_page()
        for url in PAGES:
            page_id = _page_id(url)
            target = out_dir / f"{page_id}.png"
            full_url = base_url.rstrip("/") + url
            print(f"  {url}  ->  {target.name}")
            page.goto(full_url, wait_until="networkidle", timeout=30_000)
            page.screenshot(path=str(target), full_page=True)
            results.append((url, target))
        browser.close()
    return results


def _pixel_diff_ratio(a_path: Path, b_path: Path) -> float:
    from PIL import Image, ImageChops

    a = Image.open(a_path).convert("RGB")
    b = Image.open(b_path).convert("RGB")
    if a.size != b.size:
        return 1.0
    diff = ImageChops.difference(a, b)
    bbox = diff.getbbox()
    if bbox is None:
        return 0.0
    # Proportion of pixels changed by any amount
    area_changed = (bbox[2] - bbox[0]) * (bbox[3] - bbox[1])
    return area_changed / (a.size[0] * a.size[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--update",
        action="store_true",
        help="Rebase: overwrite the baseline PNGs with the current capture.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.02,
        help="Allowed fraction of pixels to differ (default 2%%).",
    )
    args = parser.parse_args()

    base_url = os.environ.get("APP_URL", "http://127.0.0.1:8001")
    print(f"Capturing pages from {base_url} ...")

    if args.update or not BASELINE_DIR.exists() or not any(BASELINE_DIR.iterdir()):
        _capture(base_url, BASELINE_DIR)
        print(f"Baseline written to {BASELINE_DIR}")
        return 0

    # Diff mode
    captures = _capture(base_url, DIFF_DIR)
    regressions: list[tuple[str, float]] = []
    for url, current in captures:
        page_id = _page_id(url)
        baseline = BASELINE_DIR / f"{page_id}.png"
        if not baseline.exists():
            print(f"  NEW: {url} has no baseline yet")
            continue
        ratio = _pixel_diff_ratio(baseline, current)
        if ratio > args.tolerance:
            regressions.append((url, ratio))
            print(f"  REGRESSED: {url}  {ratio:.2%}")
        else:
            print(f"  OK: {url}  {ratio:.2%}")

    if regressions:
        print("\nVisual regressions:")
        for url, ratio in regressions:
            print(f"  {url}: {ratio:.2%} (tolerance {args.tolerance:.2%})")
        print(f"\nDiff captures live in {DIFF_DIR}")
        print("If the change is intentional, rerun with --update.")
        return 1
    print("\nAll pages within tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
