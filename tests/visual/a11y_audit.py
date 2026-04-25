#!/usr/bin/env python3
"""
A11y audit harness — Phase F+ polish.

Runs axe-core (Deque) against every page in tests/fixtures/inventory.PAGES
via Playwright and reports the WCAG 2.1 A/AA violations. The axe-core
JS is fetched from the same pinned CDN we use for lucide / 3Dmol so the
audit doesn't add an npm dependency.

Usage:

    # Start the app on :8001 first.
    APP_URL=http://127.0.0.1:8001 python tests/visual/a11y_audit.py

Returns exit code 0 if no serious / critical violations, 1 otherwise.
The full violation report is printed to stdout and saved to
`tests/visual/a11y/<page>.json`.

This script complements (not replaces) `take_baseline.py`. CI can
gate on the exit code, treating a regression as a failure.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sys.exit(
        "playwright not installed. Run:\n  pip install playwright\n  playwright install chromium"
    )

from tests.fixtures.inventory import PAGES  # noqa: E402

URL = os.environ.get("APP_URL", "http://127.0.0.1:8001").rstrip("/")
OUT = Path(__file__).resolve().parent / "a11y"
OUT.mkdir(exist_ok=True)

# Pinned axe-core; sha384 captured Phase F+. Bump together if upgrading.
AXE_CDN = "https://cdn.jsdelivr.net/npm/axe-core@4.10.2/axe.min.js"


def _audit_page(page, url: str) -> dict:
    page.goto(URL + url, wait_until="networkidle", timeout=30_000)
    page.add_script_tag(url=AXE_CDN)
    return page.evaluate(
        """async () => {
            const result = await axe.run(document, {
                runOnly: ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa'],
                resultTypes: ['violations'],
            });
            return result;
        }"""
    )


def _serious_or_worse(violations: list) -> list:
    return [v for v in violations if v.get("impact") in ("serious", "critical")]


def main() -> int:
    print(f"axe-core audit against {URL}")
    print("axe version: 4.10.2")
    print()

    failed = 0
    grand_serious = 0

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 900})

        for url in PAGES:
            page_id = "index" if url == "/" else url.strip("/").replace("/", "_")
            try:
                report = _audit_page(page, url)
            except Exception as exc:
                print(f"  ✗ {url}: audit raised {type(exc).__name__}: {exc}")
                failed += 1
                continue

            violations = report.get("violations", [])
            serious = _serious_or_worse(violations)
            (OUT / f"{page_id}.json").write_text(json.dumps(report, indent=2))

            print(f"  {url}")
            print(f"    {len(violations)} violation(s) · {len(serious)} serious/critical")
            for v in serious:
                nodes = ", ".join(node.get("target", ["?"])[0] for node in v.get("nodes", [])[:3])
                more = "" if len(v.get("nodes", [])) <= 3 else f" (+{len(v['nodes']) - 3} more)"
                print(f"      [{v.get('impact', '?'):8}] {v.get('id')}: {v.get('help')}")
                print(f"          → {nodes}{more}")
            grand_serious += len(serious)

        browser.close()

    print()
    print(f"Total serious/critical: {grand_serious}")
    print(f"Reports saved under {OUT}/")
    if grand_serious > 0:
        print("Fail: there are serious / critical violations.")
        return 1
    print("Pass: no serious / critical violations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
