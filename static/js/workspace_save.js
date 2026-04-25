/**
 * workspace_save.js — make the workspace actually reachable.
 *
 * The Workspace API + drawer have shipped, but no per-page module ever
 * called Workspace.add(), so nothing populated the drawer. This module
 * closes the loop:
 *
 *   1. Injects a "+ Workspace" save button next to each page's primary
 *      SMILES input. Click → snapshot the current input value into
 *      Workspace + show a toast.
 *   2. Registers a palette command "Add current SMILES to workspace"
 *      so ⌘K → "save" → ↵ also works.
 *
 * Per-page input map mirrors what the topbar Recent / send-to flows
 * already use; reactions has no single primary input and is skipped.
 */
(function () {
  "use strict";

  const PAGE_INPUT = {
    "/":              { id: "global-search",  label: "molecule" },
    "/structure":     { id: "structure-input",label: "molecule" },
    "/descriptors":   { id: "descriptor-input",label: "molecule" },
    "/fingerprints":  { id: "fingerprint-input",label: "molecule" },
    "/similarity":    { id: "query-molecule", label: "query molecule" },
    "/coordinates":   { id: "coordinate-input",label: "molecule" },
    "/properties":    { id: "property-input", label: "molecule" },
    "/visualization": { id: "viz-input",      label: "molecule" },
  };

  function currentEntry() {
    return PAGE_INPUT[window.location.pathname] || null;
  }

  function getCurrentSmiles() {
    const entry = currentEntry();
    if (!entry) return null;
    const el = document.getElementById(entry.id);
    if (!el || typeof el.value !== "string") return null;
    const value = el.value.trim();
    return value ? { smiles: value, label: entry.label } : null;
  }

  function saveCurrent() {
    const cur = getCurrentSmiles();
    if (!cur) {
      if (window.NuGenUtils?.showAlert) {
        window.NuGenUtils.showAlert(
          "No SMILES in the current input to save.",
          "warning",
          2500
        );
      }
      return;
    }
    if (!window.Workspace || typeof window.Workspace.add !== "function") return;
    window.Workspace.add({
      type: "smiles",
      value: cur.smiles,
      name: cur.smiles,
      meta: { page: window.location.pathname },
    });
    if (window.NuGenUtils?.showAlert) {
      window.NuGenUtils.showAlert(
        `Added ${cur.label} to workspace`,
        "success",
        2500
      );
    }
  }

  function injectButton() {
    const entry = currentEntry();
    if (!entry) return;
    const input = document.getElementById(entry.id);
    if (!input) return;
    if (input.dataset.workspaceWired === "1") return;
    input.dataset.workspaceWired = "1";

    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "nug-ws-save-btn";
    btn.setAttribute("data-tooltip", "Save current SMILES to workspace");
    btn.setAttribute("aria-label", "Add current SMILES to workspace");
    btn.innerHTML = `<i data-lucide="briefcase" class="w-4 h-4"></i><span>Save</span>`;
    btn.addEventListener("click", saveCurrent);

    // Place the button to the immediate right of the input, in a wrapper
    // that keeps the input's flex / grow rules intact. Use a small inline
    // flex container so existing layouts that rely on the input's own
    // sizing continue to work.
    const wrap = document.createElement("span");
    wrap.className = "nug-ws-save-wrap";
    input.parentNode.insertBefore(wrap, input);
    wrap.appendChild(input);
    wrap.appendChild(btn);

    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
  }

  // Register palette command (topbar_search.js scans this catalog at first
  // ⌘K open, so we add to its ACTIONS list via the command-palette hook.)
  if (window.NuGenUtils?.registerAction) {
    window.NuGenUtils.registerAction("addCurrentToWorkspace", saveCurrent);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", injectButton);
  } else {
    injectButton();
  }
})();
