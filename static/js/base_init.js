/**
 * base_init.js — small bootstrap loaded on every page.
 *
 * Responsibilities:
 *   1. Run lucide.createIcons() once the SVGs land + every time
 *      a per-page module asks for an icon refresh.
 *   2. Wire the topbar sidebar-toggle (mobile hamburger).
 *   3. Register cross-cutting actions:
 *        - toggleSidebar  — shows/hides sidebar on mobile
 *        - toggleDensity  — swaps html[data-density] dense ↔ comfortable,
 *                          persisted in localStorage
 *      (toggleTheme + toggleWorkspace already live in utils.js +
 *       workspace.js respectively.)
 *   4. Reflect the current Workspace count into the topbar /
 *      sidebar badges by subscribing to Workspace.onChange.
 */
(function () {
  "use strict";

  // -------------------------------------------------------------- //
  // Lucide icon refresh
  // -------------------------------------------------------------- //
  function initLucide() {
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
  }

  // -------------------------------------------------------------- //
  // Sidebar toggle (mobile)
  // -------------------------------------------------------------- //
  function toggleSidebar() {
    const sidebar = document.getElementById("app-sidebar");
    const button = document.getElementById("sidebar-toggle");
    if (!sidebar) return;
    const open = sidebar.classList.toggle("is-open");
    if (button) button.setAttribute("aria-expanded", open ? "true" : "false");
  }

  // -------------------------------------------------------------- //
  // Density toggle (dense ↔ comfortable)
  // -------------------------------------------------------------- //
  const STORAGE_KEY = "nug-density";
  function applyStoredDensity() {
    const stored = localStorage.getItem(STORAGE_KEY);
    if (stored === "dense" || stored === "comfortable") {
      document.documentElement.setAttribute("data-density", stored);
    }
  }
  function toggleDensity() {
    const current = document.documentElement.getAttribute("data-density") || "dense";
    const next = current === "dense" ? "comfortable" : "dense";
    document.documentElement.setAttribute("data-density", next);
    try {
      localStorage.setItem(STORAGE_KEY, next);
    } catch {
      /* private browsing — ignore */
    }
    if (window.NuGenUtils?.showAlert) {
      window.NuGenUtils.showAlert(`Density: ${next}`, "info", 1500);
    }
  }

  // -------------------------------------------------------------- //
  // Workspace count badges (topbar + sidebar)
  // -------------------------------------------------------------- //
  function renderWorkspaceCount(items) {
    const count = Array.isArray(items) ? items.length : (window.Workspace?.count?.() || 0);
    const topbarBadge = document.getElementById("workspace-count-badge");
    const sidebarBadge = document.getElementById("sidebar-workspace-count");
    [topbarBadge, sidebarBadge].forEach((el) => {
      if (!el) return;
      el.textContent = String(count);
      el.hidden = count === 0;
    });
  }

  // -------------------------------------------------------------- //
  // Boot
  // -------------------------------------------------------------- //
  function boot() {
    initLucide();
    applyStoredDensity();

    // Register actions on NuGenUtils' delegator. utils.js already
    // wired toggleTheme + toggleWorkspace.
    if (window.NuGenUtils?.registerAction) {
      window.NuGenUtils.registerAction("toggleSidebar", toggleSidebar);
      window.NuGenUtils.registerAction("toggleDensity", toggleDensity);
    }

    // Live-update the workspace badge.
    if (window.Workspace?.onChange) {
      window.Workspace.onChange(renderWorkspaceCount);
      renderWorkspaceCount(window.Workspace.list?.() || []);
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
