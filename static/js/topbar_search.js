/**
 * topbar_search — ⌘K command palette.
 *
 * Phase 4 ships the search infrastructure (keyboard shortcut, fuzzy
 * matcher, palette UI). The topbar input element it focuses doesn't
 * exist yet — Phase 5's _topbar.html migration adds the
 * `#topbar-search-input` element. Until then, ⌘K opens the palette
 * directly.
 *
 * The catalog enumerates every tool page so the palette can navigate
 * to any of them without server round-trip.
 */
(function () {
  "use strict";

  // Catalog drives the palette listing. `desc` reads as a chemist-facing
  // tagline (what you'd actually do here) instead of a raw URL. `icon` is
  // a lucide name, matching the sidebar so the visual mapping is 1:1.
  const CATALOG = [
    { label: "Home",            href: "/",              icon: "home",       desc: "Dashboard · quick analysis · recent",  keywords: "dashboard index" },
    { label: "Structure",       href: "/structure",     icon: "box",        desc: "Convert SMILES ↔ MOL / SDF / InChI",   keywords: "smiles mol inchi convert validate canonical" },
    { label: "Descriptors",     href: "/descriptors",   icon: "bar-chart-2",desc: "Calculate 217+ molecular descriptors", keywords: "lipinski logp topological vsa" },
    { label: "Fingerprints",    href: "/fingerprints",  icon: "fingerprint",desc: "Morgan / RDKit / MACCS / Atom-pairs",  keywords: "morgan ecfp maccs avalon atom pairs torsion" },
    { label: "Similarity",      href: "/similarity",    icon: "search",     desc: "Tanimoto · MCS · diverse subsets",     keywords: "tanimoto dice cosine substructure mcs diverse" },
    { label: "3D Coordinates",  href: "/coordinates",   icon: "atom",       desc: "Embed conformers · UFF / MMFF",        keywords: "embed conformer optimize uff mmff geometry" },
    { label: "Properties",      href: "/properties",    icon: "beaker",     desc: "QED · ADMET · drug-likeness · BRICS",  keywords: "qed admet drug-likeness scaffold brics recap" },
    { label: "Reactions",       href: "/reactions",     icon: "flask-conical", desc: "Enumerate libraries from SMARTS",   keywords: "smarts enumeration library brics rxn" },
    { label: "Visualization",   href: "/visualization", icon: "image",      desc: "2D / 3D / pharmacophore / surface",    keywords: "draw svg png 3dmol pharmacophore similarity-map" },
    { label: "API reference",   href: "/api",           icon: "book-open",  desc: "OpenAPI / Swagger explorer",           keywords: "openapi swagger endpoints json" },
    { label: "Health",          href: "/health",        icon: "activity",   desc: "Status · version · diagnostics",       keywords: "status version git diagnostics" },
  ];

  // Quick actions surface chrome-level controls (workspace, theme, density)
  // alongside page navigation. Each runs an action via NuGenUtils so the
  // palette serves both navigation and command-bar use cases.
  const ACTIONS = [
    { label: "Add current SMILES to workspace", icon: "save",       desc: "Save the molecule from this page's input", action: "addCurrentToWorkspace", keywords: "save workspace add smiles current store" },
    { label: "Open workspace",                  icon: "briefcase",  desc: "Open / close the right-side drawer",       action: "toggleWorkspace",       keywords: "workspace clipboard drawer batch" },
    { label: "Copy all workspace SMILES",       icon: "copy",       desc: "All saved SMILES, one per line",           action: "copyAllWorkspace",      keywords: "workspace copy clipboard all smiles batch" },
    { label: "Toggle theme",                    icon: "sun-moon",   desc: "Switch between light and dark",            action: "toggleTheme",           keywords: "theme dark light mode" },
    { label: "Toggle density",                  icon: "rows-3",     desc: "Dense ↔ comfortable layout",                action: "toggleDensity",         keywords: "density spacing dense comfortable" },
  ];

  // ==========================================================================
  // Fuzzy matching — score each catalog entry for the query
  // ==========================================================================

  function score(query, entry) {
    if (!query) return 0.5;
    const q = query.toLowerCase();
    const haystack = `${entry.label} ${entry.keywords}`.toLowerCase();
    if (haystack.includes(q)) return 1 - haystack.indexOf(q) / haystack.length;
    // sub-sequence match: every char of q appears in order in haystack
    let i = 0;
    for (const ch of haystack) {
      if (ch === q[i]) i += 1;
      if (i === q.length) return 0.4;
    }
    return 0;
  }

  function rank(query, source) {
    return source
      .map((entry) => ({ entry, s: score(query, entry) }))
      .filter((r) => r.s > 0)
      .sort((a, b) => b.s - a.s)
      .slice(0, 8)
      .map((r) => r.entry);
  }

  // ==========================================================================
  // Palette UI (lazy)
  // ==========================================================================

  let palette = null;
  let resultsEl = null;
  // The topbar search input is the only text input — the dropdown is just
  // results + footer, so the user never sees two inputs at once.
  let topbarInputEl = null;
  let topbarInputBound = false;

  function ensurePalette() {
    if (palette) return palette;
    palette = document.createElement("div");
    palette.id = "nug-topbar-palette";
    palette.setAttribute("aria-hidden", "true");
    palette.setAttribute("role", "listbox");
    palette.setAttribute("aria-label", "Quick search results");
    // Anchored drop-down (no full-screen overlay). Position is set in
    // show() via getBoundingClientRect on the topbar input so the panel
    // sits flush under the search field rather than centering on screen.
    palette.className = "nug-palette hidden";
    palette.innerHTML = `
      <div class="bg-white rounded-xl shadow-xl border border-gray-200 overflow-hidden">
        <ul id="nug-topbar-palette-results" class="max-h-80 overflow-y-auto py-1" role="listbox"></ul>
        <div class="px-3 py-2 border-t border-gray-200 text-xs text-gray-500 flex items-center justify-between">
          <span><span class="kbd">↑</span> <span class="kbd">↓</span> navigate · <span class="kbd">↵</span> open · <span class="kbd">Esc</span> close</span>
        </div>
      </div>
    `;
    document.body.appendChild(palette);
    resultsEl = palette.querySelector("#nug-topbar-palette-results");
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
    return palette;
  }

  function bindTopbarInput() {
    if (topbarInputBound) return;
    topbarInputEl = document.getElementById("topbar-search-input");
    if (!topbarInputEl) return;
    topbarInputEl.addEventListener("input", () => {
      ensurePalette();
      render(topbarInputEl.value);
      if (palette.classList.contains("hidden")) show();
    });
    topbarInputEl.addEventListener("keydown", onKey);
    topbarInputEl.addEventListener("focus", () => {
      // Open on focus so click + ⌘K both surface the dropdown.
      ensurePalette();
      render(topbarInputEl.value || "");
      show();
    });
    topbarInputBound = true;
  }

  function positionUnderAnchor() {
    const anchor = document.getElementById("topbar-search-input")
                || document.querySelector(".topbar-search")
                || document.querySelector(".topbar");
    if (!anchor) return;
    const rect = anchor.getBoundingClientRect();
    // Cap the panel at the input's width up to a comfortable max so the
    // dropdown reads as belonging to that input.
    const left = Math.max(8, Math.round(rect.left));
    const top = Math.round(rect.bottom + 6);
    const minW = Math.max(rect.width, 360);
    palette.style.left = `${left}px`;
    palette.style.top = `${top}px`;
    palette.style.width = `${Math.min(minW, 520)}px`;
  }

  function onDocClick(e) {
    if (!palette || palette.classList.contains("hidden")) return;
    if (palette.contains(e.target)) return;
    // Clicking the topbar input itself shouldn't immediately close the
    // panel it just opened.
    const anchor = document.getElementById("topbar-search-input");
    if (anchor && (e.target === anchor || anchor.contains(e.target))) return;
    hide();
  }

  function show() {
    ensurePalette();
    bindTopbarInput();
    palette.setAttribute("aria-hidden", "false");
    palette.classList.remove("hidden");
    if (resultsEl.children.length === 0) render(topbarInputEl?.value || "");
    positionUnderAnchor();
    if (topbarInputEl && document.activeElement !== topbarInputEl) {
      setTimeout(() => topbarInputEl.focus(), 0);
    }
    document.addEventListener("click", onDocClick, true);
    window.addEventListener("resize", positionUnderAnchor);
    window.addEventListener("scroll", positionUnderAnchor, true);
  }

  function hide() {
    if (!palette) return;
    palette.setAttribute("aria-hidden", "true");
    palette.classList.add("hidden");
    document.removeEventListener("click", onDocClick, true);
    window.removeEventListener("resize", positionUnderAnchor);
    window.removeEventListener("scroll", positionUnderAnchor, true);
  }

  function renderRow(entry, isFirst) {
    const escape = window.NuGenUtils?.escapeHtml || ((s) => s);
    const isAction = !!entry.action;
    const dataAttrs = isAction
      ? `data-palette-action="${escape(entry.action)}"`
      : `data-palette-href="${escape(entry.href)}"`;
    return `
      <li role="option" ${dataAttrs} tabindex="0"
          class="palette-row flex items-center gap-3 px-3 py-2 text-sm cursor-pointer ${
            isFirst ? "is-active" : ""
          }">
        <span class="palette-row-icon">
          <i data-lucide="${escape(entry.icon || "circle")}" class="w-4 h-4"></i>
        </span>
        <span class="flex-1 min-w-0">
          <span class="block font-medium text-gray-900 truncate">${escape(entry.label)}</span>
          <span class="block text-xs text-gray-500 truncate">${escape(entry.desc || "")}</span>
        </span>
        ${isAction
          ? `<span class="text-[10px] uppercase tracking-wide text-gray-400">Action</span>`
          : `<i data-lucide="corner-down-left" class="w-3 h-3 text-gray-400"></i>`}
      </li>`;
  }

  function renderHeader(label) {
    const escape = window.NuGenUtils?.escapeHtml || ((s) => s);
    return `<li role="presentation" class="palette-header px-3 pt-2 pb-1 text-[10px] font-semibold uppercase tracking-wide text-gray-400">${escape(label)}</li>`;
  }

  function render(query) {
    const pages = rank(query, CATALOG);
    const actions = rank(query, ACTIONS);
    if (pages.length === 0 && actions.length === 0) {
      resultsEl.innerHTML = `<li class="px-3 py-3 text-sm text-gray-500">No matches.</li>`;
      return;
    }
    let html = "";
    let firstSeen = false;
    if (pages.length) {
      html += renderHeader("Pages");
      pages.forEach((entry) => {
        html += renderRow(entry, !firstSeen);
        firstSeen = true;
      });
    }
    if (actions.length) {
      html += renderHeader("Commands");
      actions.forEach((entry) => {
        html += renderRow(entry, !firstSeen);
        firstSeen = true;
      });
    }
    resultsEl.innerHTML = html;
    resultsEl.querySelectorAll("[data-palette-href]").forEach((li) => {
      li.addEventListener("click", () => navigate(li.dataset.paletteHref));
    });
    resultsEl.querySelectorAll("[data-palette-action]").forEach((li) => {
      li.addEventListener("click", () => runCommand(li.dataset.paletteAction));
    });
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
  }

  function navigate(href) {
    hide();
    if (href) window.location.href = href;
  }

  function runCommand(name) {
    hide();
    if (name && window.NuGenUtils && typeof window.NuGenUtils.runAction === "function") {
      window.NuGenUtils.runAction(name);
    }
  }

  function activatableRows() {
    return Array.from(resultsEl.querySelectorAll(".palette-row"));
  }

  function moveActive(delta) {
    const rows = activatableRows();
    if (rows.length === 0) return;
    const cur = rows.findIndex((r) => r.classList.contains("is-active"));
    const next = ((cur < 0 ? 0 : cur + delta) + rows.length) % rows.length;
    rows.forEach((r, i) => r.classList.toggle("is-active", i === next));
    rows[next].scrollIntoView({ block: "nearest" });
  }

  function activateCurrent() {
    const active = resultsEl.querySelector(".palette-row.is-active") ||
                   resultsEl.querySelector(".palette-row");
    if (!active) return;
    if (active.dataset.paletteHref) navigate(active.dataset.paletteHref);
    else if (active.dataset.paletteAction) runCommand(active.dataset.paletteAction);
  }

  function onKey(e) {
    if (e.key === "Escape") {
      e.preventDefault();
      hide();
    } else if (e.key === "Enter") {
      e.preventDefault();
      activateCurrent();
    } else if (e.key === "ArrowDown") {
      e.preventDefault();
      moveActive(1);
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      moveActive(-1);
    }
  }

  // ==========================================================================
  // Wire ⌘K shortcut
  // ==========================================================================

  if (window.NuGenUtils?.shortcuts) {
    window.NuGenUtils.shortcuts.add(
      "Mod+k",
      () => {
        // Focus the topbar input — its focus listener opens the dropdown.
        const topbarInput = document.getElementById("topbar-search-input");
        if (topbarInput) {
          topbarInput.focus();
          topbarInput.select();
          return;
        }
        show();
      },
      "Open quick search"
    );
  }

  if (window.NuGenUtils?.registerAction) {
    window.NuGenUtils.registerAction("openSearch", show);
  }

  // Bind the topbar input as soon as the DOM is ready so a focus on the
  // input opens the dropdown without first needing a ⌘K or click round-trip.
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", bindTopbarInput);
  } else {
    bindTopbarInput();
  }
})();
