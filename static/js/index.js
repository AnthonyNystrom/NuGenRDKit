/**
 * index.js — page module for "/" (the home dashboard).
 *
 * Phase-5 migration: replaces the inline <script> block that previously
 * lived at the bottom of templates/index.html. Uses
 *   NuGenUtils.fetchJSON     — for /api/v1/** calls (60 s timeout, JSON envelopes)
 *   NuGenUtils.escapeHtml    — for any user-supplied text inserted into DOM
 *   NuGenUtils.registerAction — wires `data-action` clicks (no inline onclick=)
 *   Formatters.number        — null-safe formatting of numeric properties
 *
 * Public actions registered (called from data-action="..." attributes):
 *   navigate(path)           — window.location = path
 *   loadExample(smiles)      — paste SMILES into #global-search and run
 *   runGlobalSearch()        — manually trigger the search
 *   quickAction(target)      — navigate to /<target>?smiles=...
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

  Utils.registerAction("navigate", function (path) {
    if (typeof path === "string" && path) window.location.href = path;
  });

  Utils.registerAction("loadExample", function (smiles, name) {
    const input = document.getElementById("global-search");
    if (!input) return;
    input.value = smiles || "";
    input.dispatchEvent(new Event("input", { bubbles: true }));
    if (name) Utils.showAlert(`Loaded ${name}`, "info", 1500);
    runGlobalSearch();
  });

  Utils.registerAction("runGlobalSearch", runGlobalSearch);

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

  // ==========================================================================
  // Global search (formerly performGlobalSearch in inline script)
  // ==========================================================================

  async function runGlobalSearch() {
    const input = document.getElementById("global-search");
    const resultsDiv = document.getElementById("search-results");
    const structureDiv = document.getElementById("search-structure");
    const propertiesDiv = document.getElementById("search-properties");
    if (!input || !resultsDiv || !structureDiv || !propertiesDiv) return;

    const smiles = input.value.trim();
    if (!smiles) {
      Utils.showAlert("Please enter a molecule to search", "warning");
      return;
    }

    resultsDiv.classList.remove("hidden");
    structureDiv.innerHTML = renderSpinner();
    propertiesDiv.innerHTML = renderPropertiesSkeleton();

    let structureData = null;
    let propertiesData = null;
    try {
      [structureData, propertiesData] = await Promise.all([
        Utils.fetchJSON("/api/v1/visualization/draw_svg", {
          method: "POST",
          body: { smiles, width: 200, height: 120 },
        }).catch((err) => ({ success: false, error: err.message })),
        Utils.fetchJSON("/api/v1/properties/physicochemical", {
          method: "POST",
          body: { smiles },
        }).catch((err) => ({ success: false, error: err.message })),
      ]);
    } catch (err) {
      console.error("Search error:", err);
      Utils.showAlert("Search failed. Please try again.", "error");
      structureDiv.innerHTML = renderError();
      propertiesDiv.innerHTML = renderPropertiesError();
      reInitIcons();
      return;
    }

    // Structure: server-rendered SVG. RDKit produces this on the
    // backend so the markup is trusted; no escapeHtml here would
    // break the SVG. Phase 7 may DOMPurify it for defense in depth.
    if (structureData.success && structureData.svg) {
      structureDiv.innerHTML = structureData.svg;
    } else {
      structureDiv.innerHTML = renderError();
    }

    // Properties: build via Formatters.number to defang any odd
    // numeric value the API might emit.
    if (propertiesData.success && propertiesData.properties) {
      const p = propertiesData.properties;
      propertiesDiv.innerHTML = renderProperties({
        MW: Fmt.number(p.molecular_weight, { digits: 1 }),
        LogP: Fmt.number(p.logp, { digits: 2 }),
        HBD: Fmt.number(p.num_h_donors, { digits: 0 }),
        HBA: Fmt.number(p.num_h_acceptors, { digits: 0 }),
      });
    } else {
      propertiesDiv.innerHTML = renderPropertiesError();
    }

    reInitIcons();
  }

  // ==========================================================================
  // DOM rendering helpers — strings are static or pre-formatted
  // ==========================================================================

  function renderSpinner() {
    return `
      <div class="flex items-center justify-center h-full">
        <div class="animate-spin rounded-full h-6 w-6 border-b-2 border-primary-600"></div>
      </div>`;
  }

  function renderError(msg = "Error") {
    return `
      <div class="flex items-center justify-center h-full text-red-500">
        <i data-lucide="alert-circle" class="w-6 h-6"></i>
        <span class="sr-only">${Utils.escapeHtml(msg)}</span>
      </div>`;
  }

  function renderPropertiesSkeleton() {
    return ["MW", "LogP", "HBD", "HBA"]
      .map(
        (label) => `
        <div class="flex justify-between">
          <span class="text-gray-600">${label}:</span>
          <span class="font-mono animate-pulse">…</span>
        </div>`
      )
      .join("");
  }

  function renderProperties(values) {
    return Object.entries(values)
      .map(
        ([label, value]) => `
        <div class="flex justify-between">
          <span class="text-gray-600">${Utils.escapeHtml(label)}:</span>
          <span class="font-mono">${Utils.escapeHtml(value)}</span>
        </div>`
      )
      .join("");
  }

  function renderPropertiesError() {
    return `
      <div class="text-center text-red-500 text-sm">
        <i data-lucide="alert-circle" class="w-4 h-4 mx-auto mb-1"></i>
        <p>Analysis failed</p>
      </div>`;
  }

  function reInitIcons() {
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
  }

  // ==========================================================================
  // Hero molecule (right side of the hero section).
  // Replaces the legacy app.js::loadHeroMolecule(). Renders a recognizable
  // molecule (caffeine) into #hero-molecule so the panel isn't empty on
  // first paint. The placeholder stays if the request fails.
  // ==========================================================================

  async function loadHeroMolecule() {
    const target = document.getElementById("hero-molecule");
    if (!target) return;
    try {
      const data = await Utils.fetchJSON("/api/v1/visualization/draw_svg", {
        method: "POST",
        body: {
          smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", // caffeine
          width: 320,
          height: 280,
        },
      });
      if (data && data.success && data.svg) {
        // SVG is server-rendered by RDKit, not user-supplied — safe to inject.
        // Solid white wrapper so the molecule renders crisp against the dark
        // hero gradient (bg-white/95 was undefined in app.css).
        target.innerHTML =
          `<div class="w-full h-full flex items-center justify-center bg-white rounded-xl p-4">` +
          data.svg +
          `</div>`;
      }
    } catch (err) {
      // Silent failure — placeholder remains visible. The /health probe
      // banner (Phase 1) will surface a global outage if the API is dead.
      console.warn("[index] hero molecule failed to load:", err.message);
    }
  }

  // ==========================================================================
  // Boot — wire keyboard + auto-search + hero molecule
  // ==========================================================================

  document.addEventListener("DOMContentLoaded", function () {
    loadHeroMolecule();

    const searchInput = document.getElementById("global-search");
    if (searchInput) {
      searchInput.addEventListener("keydown", (e) => {
        if (e.key === "Enter") {
          e.preventDefault();
          runGlobalSearch();
        }
      });
      // Auto-search if the field already has a value (page refresh/SSR).
      if (searchInput.value.trim()) runGlobalSearch();
    }
  });
})();
