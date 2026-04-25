/**
 * index.js — chemist dashboard module for "/".
 *
 * Drives:
 *   - Quick analysis (paste a SMILES → structure + properties)
 *   - "Send to tool" buttons (registered via data-action="quickAction")
 *   - "Load example" chips
 *   - System probe panel (RDKit / Python / Flask versions, status badge)
 *
 * Recent-molecules rendering lives in static/js/recent.js (Phase E).
 */
(function () {
  "use strict";
  if (window.location.pathname !== "/") return;

  const Utils = window.NuGenUtils;
  const Fmt = window.Formatters;
  if (!Utils || !Fmt) {
    console.warn("[index.js] NuGenUtils + Formatters required");
    return;
  }

  // ==========================================================================
  // data-action handlers
  // ==========================================================================

  Utils.registerAction("loadExample", function (smiles, name) {
    const input = document.getElementById("global-search");
    if (!input) return;
    input.value = smiles || "";
    input.dispatchEvent(new Event("input", { bubbles: true }));
    if (name) Utils.showAlert(`Loaded ${name}`, "info", 1500);
    runQuickAnalysis();
  });

  Utils.registerAction("runGlobalSearch", runQuickAnalysis);

  Utils.registerAction("quickAction", function (target) {
    const input = document.getElementById("global-search");
    const smiles = input ? input.value.trim() : "";
    if (!smiles) {
      Utils.showAlert("Enter a SMILES string first", "warning", 2000);
      return;
    }
    const url = new URL(`/${target}`, window.location.origin);
    url.searchParams.set("smiles", smiles);
    window.location.href = url.toString();
  });

  Utils.registerAction("clearRecent", function () {
    if (window.Recent && typeof window.Recent.clear === "function") {
      window.Recent.clear();
      Utils.showAlert("Recent cleared", "info", 1500);
    }
  });

  // ==========================================================================
  // Quick analysis
  // ==========================================================================

  let inFlight = null;

  async function runQuickAnalysis() {
    const input = document.getElementById("global-search");
    const structureDiv = document.getElementById("search-structure");
    const propertiesEl = document.getElementById("search-properties");
    if (!input || !structureDiv || !propertiesEl) return;

    const smiles = input.value.trim();
    if (!smiles) {
      Utils.showAlert("Please enter a molecule to search", "warning");
      return;
    }

    structureDiv.innerHTML = `<i data-lucide="loader-2" class="w-4 h-4 quick-results-spin"></i>`;
    refreshIcons();

    // Cancel any prior request via AbortController.
    if (inFlight) inFlight.abort();
    inFlight = new AbortController();

    let svg, props;
    try {
      [svg, props] = await Promise.all([
        Utils.fetchJSON("/api/v1/visualization/draw_svg", {
          method: "POST",
          body: { smiles, width: 220, height: 160 },
          signal: inFlight.signal,
        }).catch((err) => ({ success: false, error: err.message })),
        Utils.fetchJSON("/api/v1/properties/physicochemical", {
          method: "POST",
          body: { smiles },
          signal: inFlight.signal,
        }).catch((err) => ({ success: false, error: err.message })),
      ]);
    } catch (err) {
      if (err.name === "AbortError") return;
      console.error("[index] quick analysis failed:", err);
      Utils.showAlert("Quick analysis failed", "error", 2500);
      return;
    } finally {
      inFlight = null;
    }

    // Structure: server-rendered SVG, safe to inject.
    if (svg && svg.success && svg.svg) {
      structureDiv.innerHTML = svg.svg;
    } else {
      structureDiv.innerHTML = `<div class="alert-app-error">${Utils.escapeHtml(svg?.error || "Failed")}</div>`;
    }

    // Properties: format with Fmt.number for null safety.
    if (props && props.success && props.properties) {
      const p = props.properties;
      const pairs = {
        MW: Fmt.number(p.molecular_weight, { digits: 1 }),
        LogP: Fmt.number(p.logp, { digits: 2 }),
        HBD: Fmt.number(p.num_h_donors, { digits: 0 }),
        HBA: Fmt.number(p.num_h_acceptors, { digits: 0 }),
      };
      propertiesEl.innerHTML = Object.entries(pairs)
        .map(([k, v]) => `<div><dt>${Utils.escapeHtml(k)}</dt><dd>${Utils.escapeHtml(v)}</dd></div>`)
        .join("");
    } else {
      propertiesEl.innerHTML = `<div><dt>Error</dt><dd>${Utils.escapeHtml(props?.error || "Failed")}</dd></div>`;
    }

    // Stamp the recent-history ring buffer.
    if (window.Recent && typeof window.Recent.add === "function") {
      window.Recent.add({ smiles, page: "/" });
    }
  }

  // ==========================================================================
  // System probe panel
  // ==========================================================================

  async function loadSystemInfo() {
    const rdkit = document.getElementById("sys-rdkit");
    const python = document.getElementById("sys-python");
    const flask = document.getElementById("sys-flask");
    const status = document.getElementById("sys-status");
    if (!rdkit) return;
    try {
      const data = await Utils.fetchJSON("/health");
      rdkit.textContent = data?.versions?.rdkit || "?";
      python.textContent = data?.versions?.python || "?";
      flask.textContent = data?.versions?.flask || "?";
      if (data?.rdkit_working) {
        status.className = "badge-app-success";
        status.textContent = "healthy";
      } else {
        status.className = "badge-app-error";
        status.textContent = "degraded";
      }
    } catch (err) {
      status.className = "badge-app-error";
      status.textContent = "offline";
    }
  }

  // ==========================================================================
  // Boot
  // ==========================================================================

  function refreshIcons() {
    if (window.lucide) window.lucide.createIcons();
  }

  document.addEventListener("DOMContentLoaded", function () {
    const searchInput = document.getElementById("global-search");
    if (searchInput) {
      searchInput.addEventListener("keydown", (e) => {
        if (e.key === "Enter") {
          e.preventDefault();
          runQuickAnalysis();
        }
      });
      if (searchInput.value.trim()) runQuickAnalysis();
    }
    loadSystemInfo();
  });
})();
