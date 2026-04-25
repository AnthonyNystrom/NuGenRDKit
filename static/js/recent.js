/**
 * recent.js — recent-molecules ring buffer for the chemist UI shell.
 *
 * Persisted in sessionStorage (per-tab, gone on close — different
 * from Workspace which is meant for explicit "Save"). Designed to
 * answer "what was that thing I ran two tools ago?" without polluting
 * Workspace.
 *
 * Public API (window.Recent):
 *   Recent.add({ smiles, name?, page? })  add an entry; pushes oldest off
 *   Recent.list()                          → array, newest first
 *   Recent.clear()
 *   Recent.onChange(cb)                    subscribe; returns unsubscribe
 *
 * Renderers:
 *   - sidebar list (#sidebar-recent, every page)
 *   - dashboard grid (#recent-grid, only on /)
 *
 * Maximum 12 entries. Re-adding an existing SMILES bumps it to the
 * front (de-duplicates). Click a recent entry to load it back into
 * the current page's primary input.
 */
(function () {
  "use strict";

  const STORAGE_KEY = "nug-recent-v1";
  const MAX_ITEMS = 12;
  const subscribers = new Set();

  // ============================================================ //
  // Persistence
  // ============================================================ //

  function read() {
    try {
      const raw = sessionStorage.getItem(STORAGE_KEY);
      const parsed = raw ? JSON.parse(raw) : [];
      return Array.isArray(parsed) ? parsed : [];
    } catch {
      return [];
    }
  }

  function write(items) {
    try {
      sessionStorage.setItem(STORAGE_KEY, JSON.stringify(items));
    } catch {
      /* private browsing — ignore */
    }
    subscribers.forEach((cb) => {
      try {
        cb(items);
      } catch (err) {
        console.error("[Recent] subscriber threw:", err);
      }
    });
  }

  // ============================================================ //
  // Public API
  // ============================================================ //

  function add(item) {
    if (!item || typeof item !== "object") return null;
    const smiles = String(item.smiles || "").trim();
    if (!smiles) return null;

    const name = (item.name || "").trim();
    const page = (item.page || window.location.pathname || "").trim();
    const at = Date.now();

    const items = read();
    // De-dup by SMILES — bump existing entry to the front rather than duplicate.
    const filtered = items.filter((x) => x.smiles !== smiles);
    filtered.unshift({ smiles, name, page, at });
    if (filtered.length > MAX_ITEMS) filtered.length = MAX_ITEMS;
    write(filtered);
    return filtered[0];
  }

  function list() {
    return read();
  }

  function clear() {
    write([]);
  }

  function onChange(cb) {
    if (typeof cb !== "function") return () => {};
    subscribers.add(cb);
    return () => subscribers.delete(cb);
  }

  // ============================================================ //
  // Renderers
  // ============================================================ //

  const escape = (s) => (window.NuGenUtils?.escapeHtml || String)(s);

  function truncate(s, n) {
    if (!s) return "";
    return s.length > n ? s.slice(0, n - 1) + "…" : s;
  }

  function relativeTime(ts) {
    if (!ts) return "";
    const seconds = Math.floor((Date.now() - ts) / 1000);
    if (seconds < 60) return "just now";
    if (seconds < 3600) return `${Math.floor(seconds / 60)}m`;
    if (seconds < 86400) return `${Math.floor(seconds / 3600)}h`;
    return `${Math.floor(seconds / 86400)}d`;
  }

  function pageLabel(page) {
    if (!page || page === "/") return "Dashboard";
    return page.replace(/^\//, "");
  }

  function renderSidebar(items) {
    const ul = document.getElementById("sidebar-recent");
    if (!ul) return;
    if (!items || items.length === 0) {
      ul.innerHTML = `<li class="sidebar-recent-empty">No recent molecules.</li>`;
      return;
    }
    ul.innerHTML = items
      .slice(0, 5)
      .map((it) => {
        const display = it.name || it.smiles;
        return `
          <li>
            <button type="button"
                    data-action="loadRecent"
                    data-action-args='[${JSON.stringify(it.smiles)}]'
                    title="${escape(it.smiles)}">
              ${escape(truncate(display, 22))}
            </button>
          </li>`;
      })
      .join("");
    if (window.lucide?.createIcons) window.lucide.createIcons();
  }

  function renderDashboard(items) {
    const grid = document.getElementById("recent-grid");
    if (!grid) return;
    if (!items || items.length === 0) {
      grid.innerHTML = `
        <p class="rc-empty">
          <i data-lucide="inbox" class="rc-empty-icon"></i>
          <span>Run any tool to start a session history.</span>
        </p>`;
      if (window.lucide?.createIcons) window.lucide.createIcons();
      return;
    }
    grid.innerHTML = items
      .map(
        (it) => `
          <button type="button"
                  class="recent-card"
                  data-action="loadRecent"
                  data-action-args='[${JSON.stringify(it.smiles)}]'
                  title="${escape(it.smiles)}">
            <div class="recent-card-name">${escape(it.name || it.smiles)}</div>
            <div class="recent-card-smiles">${escape(it.smiles)}</div>
            <div class="recent-card-meta">
              <span>${escape(pageLabel(it.page))}</span>
              <span>${escape(relativeTime(it.at))}</span>
            </div>
          </button>`
      )
      .join("");
  }

  // ============================================================ //
  // Wiring
  // ============================================================ //

  function loadRecent(smiles) {
    if (!smiles) return;
    // Find the page's primary input. Most tool pages use one of these IDs.
    const candidates = [
      "structure-input",
      "descriptor-input",
      "fingerprint-input",
      "similarity-input",
      "coordinate-input",
      "properties-input",
      "reaction-input",
      "viz-input",
      "global-search",
    ];
    for (const id of candidates) {
      const el = document.getElementById(id);
      if (el) {
        el.value = smiles;
        el.dispatchEvent(new Event("input", { bubbles: true }));
        el.focus();
        if (window.NuGenUtils?.showAlert) {
          window.NuGenUtils.showAlert("Loaded into input", "info", 1500);
        }
        return;
      }
    }
    // Fallback: just put it on the clipboard / show alert
    if (window.NuGenUtils?.showAlert) {
      window.NuGenUtils.showAlert(`No input field found on this page for: ${smiles}`, "warning", 2500);
    }
  }

  function boot() {
    if (window.NuGenUtils?.registerAction) {
      window.NuGenUtils.registerAction("loadRecent", loadRecent);
    }
    onChange((items) => {
      renderSidebar(items);
      renderDashboard(items);
    });
    const initial = read();
    renderSidebar(initial);
    renderDashboard(initial);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }

  window.Recent = { add, list, clear, onChange };
})();
