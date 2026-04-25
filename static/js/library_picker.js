/**
 * library_picker.js — generic "Pick from library" modal.
 *
 * Renders a searchable, category-grouped browser over the
 * `NuGenLibrary` catalog and wires it to whatever the current page's
 * primary SMILES input is. Works on every page that ships with one of
 * the inputs in PAGE_INPUT below.
 *
 * What it does:
 *   1. Injects a small "Library" button next to the page's primary
 *      SMILES input (next to the Save → Workspace button shipped by
 *      workspace_save.js).
 *   2. Opens a modal with:
 *        - a search box (filters by name, synonyms, SMILES, category)
 *        - results grouped by category, virtualized via simple
 *          "show X · expand all" tabs
 *        - a special "Paste a MOL block" panel on /structure that
 *          calls /api/v1/structure/mol_to_smiles
 *
 * Replaces the original structure_picker.js (which only ran on
 * /structure and only knew about ~22 molecules).
 */
(function () {
  "use strict";

  const Utils = window.NuGenUtils;
  if (!Utils) return;

  // Same map workspace_save.js uses — keeps the two integrations in sync.
  const PAGE_INPUT = {
    "/":              "global-search",
    "/structure":     "structure-input",
    "/descriptors":   "descriptor-input",
    "/fingerprints":  "fingerprint-input",
    "/similarity":    "query-molecule",
    "/coordinates":   "coordinate-input",
    "/properties":    "property-input",
    "/visualization": "viz-input",
  };

  function targetInput() {
    const id = PAGE_INPUT[window.location.pathname];
    return id ? document.getElementById(id) : null;
  }

  // ---------------------------------------------------------------- //
  // Modal body
  // ---------------------------------------------------------------- //

  function buildBody() {
    const root = document.createElement("div");
    root.className = "nug-library-picker space-y-3";
    const total = window.NuGenLibrary?.count?.() || 0;

    root.innerHTML = `
      <div class="flex items-center gap-2 sticky top-0 bg-white pt-1 pb-2 z-10">
        <i data-lucide="search" class="w-4 h-4 text-gray-400"></i>
        <input id="lib-pick-search"
               type="search"
               placeholder="Search ${total} molecules · name, SMILES, synonym…"
               autocomplete="off"
               class="form-input flex-1"/>
        <span id="lib-pick-count" class="text-xs text-gray-500 whitespace-nowrap">${total}</span>
      </div>
      <div id="lib-pick-results" class="space-y-4 max-h-[60vh] overflow-y-auto pr-1"></div>
    `;

    if (window.location.pathname === "/structure") {
      const molSection = document.createElement("section");
      molSection.className = "border-t border-gray-200 pt-3 mt-3";
      molSection.innerHTML = `
        <h3 class="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">Or paste a MOL block</h3>
        <textarea id="lib-pick-molblock"
                  rows="4"
                  class="w-full px-3 py-2 border border-gray-300 rounded-lg font-mono text-xs focus:ring-2 focus:ring-primary-500"
                  placeholder="MOL block content here…"></textarea>
        <button type="button" id="lib-pick-convert-mol"
                class="mt-2 inline-flex items-center px-3 py-2 bg-primary-600 text-white rounded-lg text-sm hover:bg-primary-700">
          Convert to SMILES
        </button>
      `;
      root.appendChild(molSection);
    }
    return root;
  }

  function renderResults(query) {
    const out = document.getElementById("lib-pick-results");
    const countEl = document.getElementById("lib-pick-count");
    if (!out) return;
    const lib = window.NuGenLibrary;
    if (!lib) {
      out.innerHTML = `<p class="text-sm text-gray-500 px-2 py-4">Library not loaded.</p>`;
      return;
    }
    const matches = lib.search(query || "");
    countEl.textContent = String(matches.length);
    if (matches.length === 0) {
      out.innerHTML = `<p class="text-sm text-gray-500 px-2 py-6 text-center">No matches.</p>`;
      return;
    }
    // Group by category, preserving the catalog's section order.
    const order = lib.categories();
    const grouped = new Map();
    order.forEach((c) => grouped.set(c, []));
    matches.forEach((m) => {
      if (!grouped.has(m.category)) grouped.set(m.category, []);
      grouped.get(m.category).push(m);
    });
    let html = "";
    grouped.forEach((items, category) => {
      if (items.length === 0) return;
      html += `
        <section>
          <h3 class="text-[11px] font-semibold text-gray-500 uppercase tracking-wider mb-2 sticky -top-1 bg-white pb-1">
            ${Utils.escapeHtml(category)} · <span class="text-gray-400 font-normal">${items.length}</span>
          </h3>
          <div class="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-2">
            ${items.map(renderTile).join("")}
          </div>
        </section>
      `;
    });
    out.innerHTML = html;
    // Wire clicks
    out.querySelectorAll("[data-pick-smiles]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const smi = btn.getAttribute("data-pick-smiles");
        const name = btn.getAttribute("data-pick-name") || smi;
        pickSmiles(smi, name);
      });
    });
  }

  function renderTile(item) {
    return `
      <button type="button"
              data-pick-smiles="${Utils.escapeHtml(item.smiles)}"
              data-pick-name="${Utils.escapeHtml(item.name)}"
              class="text-left px-3 py-2 bg-gray-50 border border-gray-200 rounded-lg hover:bg-primary-50 hover:border-primary-300 transition-colors min-w-0">
        <div class="text-sm font-medium text-gray-900 truncate">${Utils.escapeHtml(item.name)}</div>
        <div class="text-[11px] text-gray-500 font-mono truncate">${Utils.escapeHtml(item.smiles)}</div>
      </button>
    `;
  }

  function pickSmiles(smiles, name) {
    const input = targetInput();
    if (input) {
      input.value = smiles;
      input.dispatchEvent(new Event("input", { bubbles: true }));
      input.dispatchEvent(new Event("change", { bubbles: true }));
    }
    Utils.showAlert(`Loaded ${name}`, "success", 1800);
    Utils.closeModal();
  }

  async function convertMolBlock() {
    const ta = document.getElementById("lib-pick-molblock");
    if (!ta || !ta.value.trim()) {
      Utils.showAlert("Paste a MOL block first", "warning", 2000);
      return;
    }
    try {
      const data = await Utils.fetchJSON("/api/v1/structure/mol_to_smiles", {
        method: "POST",
        body: { mol_block: ta.value },
      });
      if (data && data.success && data.canonical_smiles) {
        pickSmiles(data.canonical_smiles, "MOL conversion");
      } else {
        Utils.showAlert(data?.error || "Conversion failed", "error", 3000);
      }
    } catch (err) {
      Utils.showAlert(err.message || "Conversion failed", "error", 3000);
    }
  }

  function openPicker() {
    if (!targetInput()) {
      Utils.showAlert("This page has no primary SMILES input.", "warning", 2500);
      return;
    }
    Utils.openModal({
      title: "Pick a molecule from the library",
      body: buildBody(),
      size: "lg",
      buttons: [{ label: "Close", action: "closeModal", variant: "secondary" }],
    });
    setTimeout(() => {
      const search = document.getElementById("lib-pick-search");
      if (search) {
        search.addEventListener("input", () => renderResults(search.value));
        search.focus();
      }
      const mol = document.getElementById("lib-pick-convert-mol");
      if (mol) mol.addEventListener("click", convertMolBlock);
      renderResults("");
    }, 0);
  }

  // ---------------------------------------------------------------- //
  // Per-page integration: inject a "Library" button next to the input
  // and register both the legacy openStructurePicker action and the
  // new generic openLibraryPicker action so existing markup keeps working.
  // ---------------------------------------------------------------- //

  Utils.registerAction("openLibraryPicker", openPicker);
  Utils.registerAction("openStructurePicker", openPicker);

  function injectButton() {
    const input = targetInput();
    if (!input) return;
    if (input.dataset.libraryWired === "1") return;
    input.dataset.libraryWired = "1";

    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "nug-lib-pick-btn";
    btn.setAttribute("data-tooltip", "Pick from library (264 molecules)");
    btn.setAttribute("aria-label", "Pick a molecule from the library");
    btn.innerHTML = `<i data-lucide="library" class="w-4 h-4"></i><span>Library</span>`;
    btn.addEventListener("click", openPicker);

    // workspace_save.js may have already wrapped the input in
    // .nug-ws-save-wrap. If so, reuse that wrapper. Otherwise create one.
    const wrap = input.parentElement?.classList.contains("nug-ws-save-wrap")
      ? input.parentElement
      : (() => {
          const w = document.createElement("span");
          w.className = "nug-ws-save-wrap";
          input.parentNode.insertBefore(w, input);
          w.appendChild(input);
          return w;
        })();
    wrap.appendChild(btn);

    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", injectButton);
  } else {
    injectButton();
  }
})();
