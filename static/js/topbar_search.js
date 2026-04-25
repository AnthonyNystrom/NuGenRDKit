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

  const CATALOG = [
    { label: "Home", href: "/", keywords: "dashboard index" },
    { label: "Structure", href: "/structure", keywords: "smiles mol inchi convert validate canonical" },
    { label: "Descriptors", href: "/descriptors", keywords: "lipinski logp topological vsa" },
    { label: "Fingerprints", href: "/fingerprints", keywords: "morgan ecfp maccs avalon atom pairs torsion" },
    { label: "Similarity", href: "/similarity", keywords: "tanimoto dice cosine substructure mcs diverse" },
    { label: "3D Coordinates", href: "/coordinates", keywords: "embed conformer optimize uff mmff geometry" },
    { label: "Properties", href: "/properties", keywords: "qed admet drug-likeness scaffold brics recap" },
    { label: "Reactions", href: "/reactions", keywords: "smarts enumeration library brics rxn" },
    { label: "Visualization", href: "/visualization", keywords: "draw svg png 3dmol pharmacophore similarity-map" },
    { label: "API reference", href: "/api", keywords: "openapi swagger endpoints json" },
    { label: "Health", href: "/health", keywords: "status version git diagnostics" },
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

  function rank(query) {
    return CATALOG.map((entry) => ({ entry, s: score(query, entry) }))
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
  let inputEl = null;

  function ensurePalette() {
    if (palette) return palette;
    palette = document.createElement("div");
    palette.id = "nug-topbar-palette";
    palette.setAttribute("aria-hidden", "true");
    palette.setAttribute("role", "dialog");
    palette.setAttribute("aria-label", "Quick search");
    palette.className =
      "fixed inset-0 z-50 hidden items-start justify-center pt-24 px-4 bg-black bg-opacity-50";
    palette.innerHTML = `
      <div class="bg-white rounded-xl shadow-xl border border-gray-200 w-full max-w-lg overflow-hidden">
        <div class="flex items-center px-4 border-b border-gray-200">
          <i data-lucide="search" class="w-4 h-4 text-gray-400"></i>
          <input id="nug-topbar-palette-input"
                 type="text"
                 placeholder="Jump to a tool…"
                 autocomplete="off"
                 class="flex-1 py-3 px-3 text-sm bg-transparent focus:outline-none">
          <span class="text-xs text-gray-400 mr-1 hidden sm:inline">Esc</span>
        </div>
        <ul id="nug-topbar-palette-results" class="max-h-64 overflow-y-auto py-1" role="listbox"></ul>
      </div>
    `;
    palette.addEventListener("click", (e) => {
      if (e.target === palette) hide();
    });
    document.body.appendChild(palette);
    inputEl = palette.querySelector("#nug-topbar-palette-input");
    resultsEl = palette.querySelector("#nug-topbar-palette-results");
    inputEl.addEventListener("input", () => render(inputEl.value));
    inputEl.addEventListener("keydown", onKey);
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
    return palette;
  }

  function show() {
    ensurePalette();
    palette.setAttribute("aria-hidden", "false");
    palette.classList.remove("hidden");
    palette.classList.add("flex");
    inputEl.value = "";
    render("");
    setTimeout(() => inputEl.focus(), 0);
  }

  function hide() {
    if (!palette) return;
    palette.setAttribute("aria-hidden", "true");
    palette.classList.add("hidden");
    palette.classList.remove("flex");
  }

  function render(query) {
    const escape = window.NuGenUtils?.escapeHtml || ((s) => s);
    const items = rank(query);
    if (items.length === 0) {
      resultsEl.innerHTML = `<li class="px-4 py-3 text-sm text-gray-500">No matches.</li>`;
      return;
    }
    resultsEl.innerHTML = items
      .map(
        (entry, i) => `
          <li role="option" data-href="${escape(entry.href)}" tabindex="0"
              class="px-4 py-2 text-sm cursor-pointer ${
                i === 0 ? "bg-primary-50 text-primary-700" : "text-gray-700 hover:bg-gray-50"
              }">
            <div class="font-medium">${escape(entry.label)}</div>
            <div class="text-xs text-gray-500">${escape(entry.href)}</div>
          </li>`
      )
      .join("");
    resultsEl.querySelectorAll("[data-href]").forEach((li) => {
      li.addEventListener("click", () => navigate(li.dataset.href));
    });
  }

  function navigate(href) {
    hide();
    if (href) window.location.href = href;
  }

  function onKey(e) {
    if (e.key === "Escape") {
      e.preventDefault();
      hide();
    } else if (e.key === "Enter") {
      e.preventDefault();
      const first = resultsEl.querySelector("[data-href]");
      if (first) navigate(first.dataset.href);
    } else if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      // simple highlight cycling — Phase 7 polish can add proper cursor state
    }
  }

  // ==========================================================================
  // Wire ⌘K shortcut
  // ==========================================================================

  if (window.NuGenUtils?.shortcuts) {
    window.NuGenUtils.shortcuts.add(
      "Mod+k",
      () => {
        // Phase-5 hook: if a real input exists in the topbar, focus that
        // instead of opening the palette.
        const topbarInput = document.getElementById("topbar-search-input");
        if (topbarInput) {
          topbarInput.focus();
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
})();
