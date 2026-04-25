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

  /**
   * Reorder the workspace by moving the item at `fromIndex` to `toIndex`.
   * Indices are clamped; out-of-range moves no-op.
   */
  function reorder(fromIndex, toIndex) {
    const items = read();
    if (fromIndex < 0 || fromIndex >= items.length) return false;
    const target = Math.max(0, Math.min(items.length - 1, toIndex));
    if (target === fromIndex) return false;
    const [moved] = items.splice(fromIndex, 1);
    items.splice(target, 0, moved);
    write(items);
    return true;
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
          <button type="button" data-action="copyAllWorkspace" class="text-xs text-gray-500 hover:text-primary-700" title="Copy all SMILES (one per line)">Copy all</button>
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

  // -------------------------------------------------------------------- //
  // Per-row utility: load SMILES into current page's input, copy, send-to.
  // The PAGE_INPUT map mirrors workspace_save.js / library_picker.js so
  // every page that surfaces a "Save" button also accepts a load-back.
  // -------------------------------------------------------------------- //

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

  const SEND_TO_DESTINATIONS = [
    { label: "Structure",     page: "/structure"    , icon: "box"          },
    { label: "Descriptors",   page: "/descriptors"  , icon: "bar-chart-3"  },
    { label: "Fingerprints",  page: "/fingerprints" , icon: "fingerprint"  },
    { label: "Similarity",    page: "/similarity"   , icon: "search"       },
    { label: "3D Coords",     page: "/coordinates"  , icon: "orbit"        },
    { label: "Properties",    page: "/properties"   , icon: "flask-conical"},
    { label: "Visualization", page: "/visualization", icon: "image"        },
  ];

  function loadIntoCurrentPage(item) {
    const inputId = PAGE_INPUT[window.location.pathname];
    const el = inputId ? document.getElementById(inputId) : null;
    if (!el) {
      // No primary input on this page (e.g. /reactions) — fall back to
      // sending the user to /structure with the SMILES pre-filled.
      window.location.href = `/structure?smiles=${encodeURIComponent(item.value)}`;
      return;
    }
    el.value = item.value;
    el.dispatchEvent(new Event("input", { bubbles: true }));
    el.dispatchEvent(new Event("change", { bubbles: true }));
    el.focus();
    if (window.NuGenUtils?.showAlert) {
      window.NuGenUtils.showAlert(`Loaded ${item.name || item.value}`, "success", 1800);
    }
  }

  async function copyItem(item) {
    if (!navigator.clipboard) {
      window.NuGenUtils?.showAlert("Clipboard not available", "error", 2500);
      return;
    }
    try {
      await navigator.clipboard.writeText(item.value);
      window.NuGenUtils?.showAlert("Copied SMILES", "success", 1500);
    } catch {
      window.NuGenUtils?.showAlert("Clipboard write blocked", "error", 2500);
    }
  }

  function buildSendToMenu(item) {
    const here = window.location.pathname;
    const escape = window.NuGenUtils?.escapeHtml || ((s) => s);
    const links = SEND_TO_DESTINATIONS
      .filter((d) => d.page !== here)
      .map(
        (d) => `
          <li role="menuitem">
            <a href="${d.page}?smiles=${encodeURIComponent(item.value)}"
               class="ws-sendto-link">
              <i data-lucide="${escape(d.icon)}" class="w-3 h-3"></i>
              <span>${escape(d.label)}</span>
            </a>
          </li>`
      )
      .join("");
    return `<ul class="ws-sendto-menu" role="menu" hidden>${links}</ul>`;
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
    const escape = window.NuGenUtils?.escapeHtml || ((s) => s);
    listEl.innerHTML = items
      .map((item, idx) => {
        const label = item.name ? escape(item.name) : escape(item.value);
        const sub = item.name ? escape(item.value) : escape(item.type || "");
        return `
          <article class="nug-workspace-item"
                   draggable="true"
                   data-index="${idx}"
                   data-id="${escape(item.id)}"
                   title="Click to load into current page">
            <span class="ws-drag-handle" aria-hidden="true" title="Drag to reorder">⋮⋮</span>
            <div class="ws-row-body">
              <p class="ws-row-name">${label}</p>
              <p class="ws-row-sub">${sub}</p>
            </div>
            <div class="ws-row-actions">
              <button type="button" class="ws-act ws-act-copy" aria-label="Copy SMILES" data-tooltip="Copy SMILES">
                <i data-lucide="clipboard" class="w-4 h-4"></i>
              </button>
              <span class="ws-sendto-wrap">
                <button type="button" class="ws-act ws-act-sendto" aria-haspopup="menu" aria-expanded="false" aria-label="Send to other tool" data-tooltip="Send to…">
                  <i data-lucide="send" class="w-4 h-4"></i>
                </button>
                ${buildSendToMenu(item)}
              </span>
              <button type="button" class="ws-act ws-act-remove" aria-label="Remove" data-tooltip="Remove">
                <i data-lucide="trash-2" class="w-4 h-4"></i>
              </button>
            </div>
          </article>`;
      })
      .join("");
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
    wireRowActions(listEl, items);
    wireDrag(listEl);
  }

  function wireRowActions(listEl, items) {
    listEl.querySelectorAll(".nug-workspace-item").forEach((row) => {
      const idx = parseInt(row.dataset.index, 10);
      const item = items[idx];
      if (!item) return;

      // Click anywhere on the body / name → load into current page.
      const body = row.querySelector(".ws-row-body");
      if (body) {
        body.addEventListener("click", (e) => {
          e.stopPropagation();
          loadIntoCurrentPage(item);
        });
      }

      // Copy
      const copyBtn = row.querySelector(".ws-act-copy");
      if (copyBtn) {
        copyBtn.addEventListener("click", (e) => {
          e.stopPropagation();
          copyItem(item);
        });
      }

      // Remove
      const removeBtn = row.querySelector(".ws-act-remove");
      if (removeBtn) {
        removeBtn.addEventListener("click", (e) => {
          e.stopPropagation();
          remove(item.id);
        });
      }

      // Send-to dropdown
      const sendBtn = row.querySelector(".ws-act-sendto");
      const menu = row.querySelector(".ws-sendto-menu");
      if (sendBtn && menu) {
        sendBtn.addEventListener("click", (e) => {
          e.stopPropagation();
          // Close any other open menus first.
          listEl.querySelectorAll(".ws-sendto-menu").forEach((m) => {
            if (m !== menu) m.hidden = true;
          });
          listEl.querySelectorAll(".ws-act-sendto").forEach((b) => {
            if (b !== sendBtn) b.setAttribute("aria-expanded", "false");
          });
          menu.hidden = !menu.hidden;
          sendBtn.setAttribute("aria-expanded", menu.hidden ? "false" : "true");
        });
      }
    });
  }

  // Document-level click closes any open send-to menu.
  document.addEventListener("click", () => {
    if (!drawer) return;
    drawer.querySelectorAll(".ws-sendto-menu").forEach((m) => (m.hidden = true));
    drawer.querySelectorAll(".ws-act-sendto").forEach((b) => b.setAttribute("aria-expanded", "false"));
  });

  function wireDrag(listEl) {
    let dragIdx = null;
    listEl.querySelectorAll(".nug-workspace-item").forEach((el) => {
      el.addEventListener("dragstart", (e) => {
        dragIdx = parseInt(el.dataset.index, 10);
        el.classList.add("opacity-50");
        if (e.dataTransfer) {
          e.dataTransfer.effectAllowed = "move";
          e.dataTransfer.setData("text/plain", String(dragIdx));
        }
      });
      el.addEventListener("dragend", () => {
        el.classList.remove("opacity-50");
        dragIdx = null;
      });
      el.addEventListener("dragover", (e) => {
        e.preventDefault();
        if (e.dataTransfer) e.dataTransfer.dropEffect = "move";
      });
      el.addEventListener("drop", (e) => {
        e.preventDefault();
        const targetIdx = parseInt(el.dataset.index, 10);
        const fromIdx = dragIdx ?? parseInt(e.dataTransfer?.getData("text/plain") || "-1", 10);
        if (Number.isFinite(fromIdx) && Number.isFinite(targetIdx)) {
          reorder(fromIdx, targetIdx);
        }
      });
    });
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
    window.NuGenUtils.registerAction("copyAllWorkspace", async () => {
      const items = list();
      if (items.length === 0) {
        window.NuGenUtils.showAlert("Workspace is empty", "info", 2000);
        return;
      }
      const text = items
        .filter((i) => (i.type || "smiles") === "smiles" && i.value)
        .map((i) => i.value)
        .join("\n");
      try {
        await navigator.clipboard.writeText(text);
        window.NuGenUtils.showAlert(`Copied ${items.length} SMILES to clipboard`, "success", 2500);
      } catch {
        window.NuGenUtils.showAlert("Clipboard write blocked by browser", "error", 3000);
      }
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
    reorder,
    onChange,
    fillInto,
  };
})();
