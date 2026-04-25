#!/usr/bin/env python3
"""
End-to-end Playwright flow tests — Phase F+ polish.

Complements `take_baseline.py` (visual snapshot) and `a11y_audit.py`
(static accessibility) by actually clicking through real workflows
and asserting the rendered output.

Each flow:
  1. Loads a tool page
  2. Verifies the auto-run results landed
  3. Triggers a state change (load an example, switch operation, etc.)
  4. Asserts the post-change UI is correct

Usage:

    APP_URL=http://127.0.0.1:8001 python tests/visual/e2e_flows.py

Returns 0 on all-green, 1 on any failure. Prints a per-flow line for
quick scanning.

CI can gate on this exit code in addition to pytest. Run after
`tests/visual/take_baseline.py --update` to confirm both contracts.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

try:
    from playwright.sync_api import expect, sync_playwright
except ImportError:
    sys.exit(
        "playwright not installed. Run:\n  pip install playwright\n  playwright install chromium"
    )


URL = os.environ.get("APP_URL", "http://127.0.0.1:8001").rstrip("/")
TIMEOUT = 30_000


# --------------------------------------------------------------------------- #
# Flow definitions — each returns (label, passed, detail) on completion.
# --------------------------------------------------------------------------- #


def flow_dashboard_quick_analysis(page):
    """Dashboard auto-loads aspirin and renders MW + LogP into the panel."""
    page.goto(URL + "/", wait_until="networkidle", timeout=TIMEOUT)
    # Quick analysis panel renders MW after the auto-run completes.
    page.wait_for_selector("#search-properties dt", timeout=TIMEOUT)
    mw_label = page.locator("#search-properties dt").first
    expect(mw_label).to_have_text("MW", timeout=TIMEOUT)
    # The structure SVG should be inside #search-structure.
    svg = page.locator("#search-structure svg").first
    expect(svg).to_be_visible(timeout=TIMEOUT)
    # Click an example chip; properties should refresh to caffeine's MW (~194).
    page.click('button[data-action="loadExample"][data-action-args*="Caffeine"]')
    page.wait_for_function(
        """() => {
            const dd = document.querySelector('#search-properties dd');
            if (!dd) return false;
            const v = parseFloat(dd.textContent);
            return Number.isFinite(v) && v > 190 && v < 200;
        }""",
        timeout=TIMEOUT,
    )


def flow_recent_populates_sidebar(page):
    """After dashboard auto-runs, the sidebar Recent list should be non-empty."""
    page.goto(URL + "/", wait_until="networkidle", timeout=TIMEOUT)
    # Wait for the auto-run to stamp something into Recent.
    page.wait_for_function(
        """() => {
            const list = document.getElementById('sidebar-recent');
            if (!list) return false;
            return !list.querySelector('.sidebar-recent-empty');
        }""",
        timeout=TIMEOUT,
    )
    # At least one button in the sidebar Recent list.
    btns = page.locator("#sidebar-recent button")
    expect(btns.first).to_be_visible(timeout=TIMEOUT)


def flow_structure_canonicalize(page):
    """/structure with default SMILES (CCO) should produce canonical_smiles output."""
    page.goto(URL + "/structure", wait_until="networkidle", timeout=TIMEOUT)
    # Structure Display panel should populate after auto-convert.
    page.wait_for_selector("#structure-display svg, #structure-output", timeout=TIMEOUT)
    # Switch to a different operation (InChI) and click Convert.
    page.select_option("#structure-format", "inchi")
    page.click("#convert-btn")
    page.wait_for_function(
        """() => {
            const out = document.body.innerText;
            return out.includes('InChI=') || out.includes('inchi');
        }""",
        timeout=TIMEOUT,
    )


def flow_descriptors_calculate(page):
    """/descriptors auto-runs CCO → Lipinski / 217+ values render in the grid."""
    page.goto(URL + "/descriptors", wait_until="networkidle", timeout=TIMEOUT)
    page.wait_for_selector("#descriptor-results", timeout=TIMEOUT)
    # The body should contain at least the molecular weight label.
    page.wait_for_function(
        """() => /Molecular Weight|MW|molecular_weight/i.test(document.body.innerText)""",
        timeout=TIMEOUT,
    )


def flow_send_to_navigates(page):
    """The Send-To submenu on the dashboard quick-action navigates with ?smiles="""
    page.goto(URL + "/", wait_until="networkidle", timeout=TIMEOUT)
    # Wait for default SMILES (aspirin) to be in the input.
    page.wait_for_function(
        "() => document.getElementById('global-search')?.value?.length > 5",
        timeout=TIMEOUT,
    )
    page.click('button[data-action="quickAction"][data-action-args*="descriptors"]')
    page.wait_for_url("**/descriptors?smiles=*", timeout=TIMEOUT)


def flow_sidebar_active_state(page):
    """Sidebar highlights the current page."""
    page.goto(URL + "/structure", wait_until="networkidle", timeout=TIMEOUT)
    active = page.locator(".sidebar-link.is-active").first
    expect(active).to_have_attribute("aria-current", "page", timeout=TIMEOUT)
    expect(active).to_contain_text("Structure", timeout=TIMEOUT)


def flow_density_toggle_persists(page):
    """Density toggle writes to localStorage + flips data-density."""
    page.goto(URL + "/", wait_until="networkidle", timeout=TIMEOUT)
    initial = page.evaluate("document.documentElement.getAttribute('data-density') || 'dense'")
    page.click('button[data-action="toggleDensity"]')
    page.wait_for_function(
        f"() => document.documentElement.getAttribute('data-density') !== '{initial}'",
        timeout=TIMEOUT,
    )
    stored = page.evaluate("localStorage.getItem('nug-density')")
    if stored not in ("dense", "comfortable"):
        raise AssertionError(f"density not stored: got {stored!r}")


FLOWS = [
    ("dashboard quick analysis (aspirin → MW)", flow_dashboard_quick_analysis),
    ("recent populates sidebar after auto-run", flow_recent_populates_sidebar),
    ("structure: canonicalize default + switch op", flow_structure_canonicalize),
    ("descriptors: auto-run renders MW", flow_descriptors_calculate),
    ("send-to: dashboard → /descriptors?smiles=", flow_send_to_navigates),
    ("sidebar marks current page as active", flow_sidebar_active_state),
    ("density toggle flips + persists", flow_density_toggle_persists),
]


def main() -> int:
    print(f"E2E flows against {URL}\n")
    failures: list[tuple[str, str]] = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        for label, flow in FLOWS:
            ctx = browser.new_context(viewport={"width": 1440, "height": 900})
            page = ctx.new_page()
            try:
                flow(page)
                print(f"  PASS  {label}")
            except Exception as exc:
                print(f"  FAIL  {label}")
                print(f"        {type(exc).__name__}: {exc}".replace("\n", " ")[:240])
                failures.append((label, str(exc)))
            finally:
                ctx.close()
        browser.close()

    print()
    if failures:
        print(f"{len(failures)} flow(s) failed.")
        return 1
    print(f"All {len(FLOWS)} flows passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
