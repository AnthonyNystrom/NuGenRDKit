/**
 * Workspace — cross-tool clipboard, persisted in sessionStorage.
 *
 * Public surface:
 *   Workspace.add(item)           // {id?, type, value, name, meta?}
 *   Workspace.list()              // → array
 *   Workspace.get(id)             // → item or null
 *   Workspace.remove(id)          // → bool
 *   Workspace.count()             // → integer
 *   Workspace.clear()
 *   Workspace.onChange(callback)  // returns unsubscribe fn
 *   Workspace.fillInto(elementId) // pastes the most recent matching value
 *
 * UI: a right-anchored drawer toggled by NuGenUtils.runAction("toggleWorkspace").
 * Drawer is created lazily on first toggle to keep first paint fast.
 *
 * Requires: window.NuGenUtils (utils.js loaded first).
 */
(function () {
  "use strict";

  const STORAGE_KEY = "nug-workspace-v1";
  const COUNTER_KEY = "nug-workspace-counter-v1";
  const subscribers = new Set();

  // ==========================================================================
  // Persistence
  // ==========================================================================

  function read() {
    try {
      const raw = sessionStorage.getItem(STORAGE_KEY);
      if (!raw) return [];
      const parsed = JSON.parse(raw);
      return Array.isArray(parsed) ? parsed : [];
    } catch {
      return [];
    }
  }

  function write(items) {
    try {
      sessionStorage.setItem(STORAGE_KEY, JSON.stringify(items));
    } catch {
      // QuotaExceeded or private mode — degrade silently.
    }
    subscribers.forEach((cb) => {
      try {
        cb(items);
      } catch (err) {
        console.error("[Workspace] subscriber threw:", err);
      }
    });
  }

  function nextId() {
    let n;
    try {
      n = parseInt(sessionStorage.getItem(COUNTER_KEY) || "0", 10) + 1;
      sessionStorage.setItem(COUNTER_KEY, String(n));
    } catch {
      n = Date.now();
    }
    return `ws-${n}`;
  }

  // ==========================================================================
  // Public API
  // ==========================================================================

  function add(item) {
    if (!item || typeof item !== "object") {
      throw new TypeError("Workspace.add: item must be an object");
    }
    const items = read();
    const next = {
      id: item.id || nextId(),
      type: item.type || "smiles",
      value: item.value ?? "",
      name: item.name || "",
      meta: item.meta || null,
      addedAt: Date.now(),
    };
    items.unshift(next);
    write(items);
    return next;
  }

  function list() {
    return read();
  }

  function get(id) {
    return read().find((i) => i.id === id) || null;
  }

  function remove(id) {
    const items = read();
    const filtered = items.filter((i) => i.id !== id);
    if (filtered.length === items.length) return false;
    write(filtered);
    return true;
  }

  function count() {
    return read().length;
  }

  function clear() {
    write([]);
  }

  function onChange(cb) {
    if (typeof cb !== "function") return () => {};
    subscribers.add(cb);
    return () => subscribers.delete(cb);
  }

  /**
   * Paste the most recent workspace item of compatible type into the named
   * input. Used by Phase-7 "Send to…" cross-tool flows.
   */
  function fillInto(elementId, { type } = {}) {
    const el = document.getElementById(elementId);
    if (!el) return false;
    const items = read();
    const candidate = type ? items.find((i) => i.type === type) : items[0];
    if (!candidate) return false;
    if ("value" in el) {
      el.value = candidate.value;
      el.dispatchEvent(new Event("input", { bubbles: true }));
      return true;
    }
    el.textContent = candidate.value;
    return true;
  }

  // ==========================================================================
  // Drawer UI (lazy)
  // ==========================================================================

  let drawer = null;

  function ensureDrawer() {
    if (drawer) return drawer;
    drawer = document.createElement("aside");
    drawer.id = "nug-workspace-drawer";
    drawer.setAttribute("aria-label", "Workspace");
    drawer.setAttribute("aria-hidden", "true");
    drawer.className =
      "fixed top-0 right-0 h-full w-80 bg-white border-l border-gray-200 shadow-xl transform translate-x-full transition-transform z-50";
    drawer.style.transitionDuration = "200ms";
    drawer.innerHTML = `
      <header class="flex items-center justify-between px-4 py-3 border-b border-gray-200">
        <h2 class="text-sm font-semibold text-gray-900">Workspace</h2>
        <div class="flex items-center space-x-2">
          <button type="button" data-action="clearWorkspace" class="text-xs text-gray-500 hover:text-red-600">Clear</button>
          <button type="button" data-action="toggleWorkspace" aria-label="Close workspace" class="text-gray-500 hover:text-gray-900">
            <i data-lucide="x" class="w-5 h-5"></i>
          </button>
        </div>
      </header>
      <div id="nug-workspace-list" class="overflow-y-auto h-full pb-12 divide-y divide-gray-200"></div>
      <p id="nug-workspace-empty" class="text-center text-sm text-gray-500 px-6 py-8">
        Workspace is empty. Use the "Send to workspace" button on a result card to add items.
      </p>
    `;
    document.body.appendChild(drawer);
    onChange(renderList);
    renderList(read());
    return drawer;
  }

  function renderList(items) {
    if (!drawer) return;
    const listEl = drawer.querySelector("#nug-workspace-list");
    const emptyEl = drawer.querySelector("#nug-workspace-empty");
    if (!items || items.length === 0) {
      listEl.innerHTML = "";
      emptyEl.classList.remove("hidden");
      return;
    }
    emptyEl.classList.add("hidden");
    listEl.innerHTML = items
      .map((item) => {
        const escape = window.NuGenUtils?.escapeHtml || ((s) => s);
        const label = item.name ? escape(item.name) : escape(item.value);
        const sub = item.name ? escape(item.value) : escape(item.type || "");
        return `
          <article class="px-4 py-3 hover:bg-gray-50">
            <div class="flex items-start justify-between gap-2">
              <div class="min-w-0 flex-1">
                <p class="text-sm font-medium text-gray-900 truncate">${label}</p>
                <p class="text-xs text-gray-500 font-mono truncate">${sub}</p>
              </div>
              <button type="button"
                      data-action="removeWorkspaceItem"
                      data-action-args='[${JSON.stringify(item.id)}]'
                      aria-label="Remove"
                      class="text-gray-400 hover:text-red-600">
                <i data-lucide="trash-2" class="w-4 h-4"></i>
              </button>
            </div>
          </article>`;
      })
      .join("");
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
  }

  function toggleDrawer() {
    const el = ensureDrawer();
    const open = el.getAttribute("aria-hidden") !== "false";
    el.setAttribute("aria-hidden", open ? "false" : "true");
    el.style.transform = open ? "translateX(0)" : "translateX(100%)";
  }

  // ==========================================================================
  // Wiring (data-action handlers + ⌘ shortcut)
  // ==========================================================================

  if (window.NuGenUtils?.registerAction) {
    window.NuGenUtils.registerAction("toggleWorkspace", toggleDrawer);
    window.NuGenUtils.registerAction("clearWorkspace", () => {
      clear();
      window.NuGenUtils.showAlert("Workspace cleared", "info", 2000);
    });
    window.NuGenUtils.registerAction("removeWorkspaceItem", (id) => {
      remove(id);
    });
    window.NuGenUtils.registerAction("addToWorkspace", (item) => {
      const stored = add(item || {});
      window.NuGenUtils.showAlert(
        `Added "${stored.name || stored.value}" to workspace`,
        "success",
        2500
      );
    });
    if (window.NuGenUtils.shortcuts) {
      window.NuGenUtils.shortcuts.add("g+w", toggleDrawer, "Toggle workspace");
    }
  }

  window.Workspace = {
    add,
    list,
    get,
    remove,
    count,
    clear,
    onChange,
    fillInto,
  };
})();
