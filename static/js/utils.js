/**
 * NuGenUtils — Phase-4 foundation module.
 *
 * The other JS modules (formatters, workspace, results_card,
 * topbar_search) depend on this. Loaded with `defer` first in
 * base.html so subsequent modules can reference window.NuGenUtils
 * at parse time.
 *
 * Public surface (kept stable across phases):
 *
 *   NuGenUtils.escapeHtml(s)
 *   NuGenUtils.fetchJSON(url, options)
 *   NuGenUtils.showAlert(message, type)
 *   NuGenUtils.showLoading(elementId)
 *   NuGenUtils.hideLoading(elementId)
 *   NuGenUtils.openModal({ title, body, buttons })
 *   NuGenUtils.closeModal()
 *   NuGenUtils.registerAction(name, fn)
 *   NuGenUtils.trapFocus(element)
 *   NuGenUtils.shortcuts.add(key, handler, label)
 *
 * Inline `onclick=` handlers in templates migrate to
 *   <button data-action="myFn" data-action-args='[1, 2]'>
 * which `data-action` delegation resolves to a registered handler at
 * click time. The args attribute is JSON.
 */
(function () {
  "use strict";

  // ==========================================================================
  // String helpers
  // ==========================================================================

  function escapeHtml(value) {
    if (value === null || value === undefined) return "";
    return String(value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  // ==========================================================================
  // fetch wrapper
  //
  // - Default 60 s timeout (configurable per call). Routes that legitimately
  //   take longer (3D embed, conformer search) pass `timeoutMs: 180_000`.
  // - Always sends + accepts JSON on POST.
  // - Throws an Error whose `.status` is the HTTP code and whose `.body` is
  //   the parsed error envelope so callers can render `.error`.
  // ==========================================================================

  const DEFAULT_TIMEOUT_MS = 60_000;

  async function fetchJSON(url, options = {}) {
    const { method = "GET", body, headers = {}, timeoutMs = DEFAULT_TIMEOUT_MS, signal } =
      options;

    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), timeoutMs);
    const compositeSignal = signal
      ? anySignal([signal, controller.signal])
      : controller.signal;

    const init = {
      method,
      headers: {
        Accept: "application/json",
        ...(body !== undefined ? { "Content-Type": "application/json" } : {}),
        ...headers,
      },
      signal: compositeSignal,
    };
    if (body !== undefined) {
      init.body = typeof body === "string" ? body : JSON.stringify(body);
    }

    let response;
    try {
      response = await fetch(url, init);
    } catch (err) {
      clearTimeout(timer);
      if (err.name === "AbortError") {
        throw makeError(`Request to ${url} timed out after ${timeoutMs}ms`, 0, null);
      }
      throw makeError(`Network error: ${err.message}`, 0, null);
    }
    clearTimeout(timer);

    let payload = null;
    try {
      payload = await response.json();
    } catch {
      // Non-JSON body — leave payload as null.
    }

    if (!response.ok) {
      const message =
        (payload && (payload.error || payload.message)) ||
        `HTTP ${response.status}`;
      throw makeError(message, response.status, payload);
    }

    return payload;
  }

  function makeError(message, status, body) {
    const err = new Error(message);
    err.status = status;
    err.body = body;
    return err;
  }

  function anySignal(signals) {
    const ctl = new AbortController();
    for (const s of signals) {
      if (s.aborted) {
        ctl.abort();
        break;
      }
      s.addEventListener("abort", () => ctl.abort(), { once: true });
    }
    return ctl.signal;
  }

  // ==========================================================================
  // Toast / alert system
  //
  // Toasts stack in a fixed container at top-right. Each auto-dismisses
  // after `duration` ms. Click the × to close early. `type` controls color.
  // ==========================================================================

  const TOAST_CONTAINER_ID = "nug-toast-container";

  function ensureToastContainer() {
    let el = document.getElementById(TOAST_CONTAINER_ID);
    if (!el) {
      el = document.createElement("div");
      el.id = TOAST_CONTAINER_ID;
      el.className = "fixed top-4 right-4 z-50 space-y-2 pointer-events-none";
      el.setAttribute("role", "status");
      el.setAttribute("aria-live", "polite");
      document.body.appendChild(el);
    }
    return el;
  }

  function showAlert(message, type = "info", duration = 5000) {
    const container = ensureToastContainer();
    const palette = {
      success: "bg-green-100 text-green-800 border-green-200",
      error: "bg-red-50 text-red-700 border-red-200",
      warning: "bg-yellow-50 text-yellow-800 border-yellow-200",
      info: "bg-blue-50 text-blue-700 border-blue-200",
    }[type] || "bg-gray-100 text-gray-800 border-gray-200";

    const toast = document.createElement("div");
    toast.className =
      "pointer-events-auto inline-flex items-start max-w-md p-3 border rounded-lg shadow-sm " +
      palette;
    const text = document.createElement("span");
    text.className = "flex-1 text-sm pr-3";
    text.textContent = message; // textContent — safe by construction
    const close = document.createElement("button");
    close.type = "button";
    close.className = "ml-2 text-lg leading-none opacity-60 hover:opacity-100";
    close.setAttribute("aria-label", "Dismiss");
    close.textContent = "×";
    close.addEventListener("click", () => toast.remove());
    toast.appendChild(text);
    toast.appendChild(close);
    container.appendChild(toast);

    if (duration > 0) {
      setTimeout(() => toast.remove(), duration);
    }
    return toast;
  }

  // ==========================================================================
  // Loading state
  // ==========================================================================

  function showLoading(elementId) {
    const el = document.getElementById(elementId);
    if (el) el.classList.remove("hidden");
  }
  function hideLoading(elementId) {
    const el = document.getElementById(elementId);
    if (el) el.classList.add("hidden");
  }

  // ==========================================================================
  // Modal helpers (paired with templates/_partials/_modal.html)
  // ==========================================================================

  function openModal({ title = "", body = "", buttons = [], size = "md" } = {}) {
    const dialog = document.getElementById("app-modal");
    if (!dialog) {
      // Caller hasn't included _modal.html on this page; fall back to alert.
      showAlert(`${title ? title + ": " : ""}${typeof body === "string" ? body : ""}`, "info");
      return;
    }
    // Size dial — swap the dialog's max-width class. md is the default
    // 32rem from _modal.html; "lg" widens to 48rem for richer pickers
    // (library picker, multi-pane settings, etc.).
    const sizeClass = ({
      sm: "max-w-md",
      md: "max-w-lg",
      lg: "max-w-3xl",
      xl: "max-w-4xl",
    })[size] || "max-w-lg";
    ["max-w-md","max-w-lg","max-w-3xl","max-w-4xl"].forEach((c) => dialog.classList.remove(c));
    dialog.classList.add(sizeClass);
    document.getElementById("app-modal-title").textContent = title;
    const bodyEl = document.getElementById("app-modal-body");
    bodyEl.innerHTML = "";
    if (body instanceof Node) bodyEl.appendChild(body);
    else if (typeof body === "string") bodyEl.textContent = body;
    const footer = document.getElementById("app-modal-footer");
    footer.innerHTML = "";
    buttons.forEach(({ label, action, variant = "secondary", argv }) => {
      const btn = document.createElement("button");
      btn.type = "button";
      const cls =
        variant === "primary"
          ? "px-4 py-2 bg-primary-600 text-white rounded-lg hover:bg-primary-700"
          : variant === "danger"
          ? "px-4 py-2 bg-red-600 text-white rounded-lg hover:bg-red-700"
          : "px-4 py-2 bg-white text-gray-700 border border-gray-300 rounded-lg hover:bg-gray-50";
      btn.className = cls;
      btn.textContent = label;
      btn.addEventListener("click", () => {
        if (action) runAction(action, argv);
        if (variant !== "danger") closeModal();
      });
      footer.appendChild(btn);
    });
    if (typeof dialog.showModal === "function") dialog.showModal();
    else dialog.setAttribute("open", "");
  }

  function closeModal() {
    const dialog = document.getElementById("app-modal");
    if (!dialog) return;
    if (typeof dialog.close === "function") dialog.close();
    else dialog.removeAttribute("open");
  }

  // ==========================================================================
  // data-action delegator
  //
  // Single document-level click listener routes
  //     <button data-action="foo" data-action-args='[1, "x"]'>
  // to a function registered with `registerAction("foo", fn)`. Modules
  // populate the registry; templates stay free of inline `onclick=`.
  // ==========================================================================

  const actionRegistry = Object.create(null);

  function registerAction(name, fn) {
    if (typeof fn !== "function") throw new TypeError(`Action ${name} must be a function`);
    actionRegistry[name] = fn;
  }

  function runAction(name, argv) {
    const fn = actionRegistry[name];
    if (!fn) {
      console.warn(`[NuGenUtils] No registered action: ${name}`);
      return;
    }
    try {
      return fn.apply(null, Array.isArray(argv) ? argv : []);
    } catch (err) {
      console.error(`[NuGenUtils] Action "${name}" threw:`, err);
      showAlert(`Action error: ${err.message || err}`, "error");
    }
  }

  function onActionClick(event) {
    const trigger = event.target.closest("[data-action]");
    if (!trigger) return;
    event.preventDefault();
    const name = trigger.getAttribute("data-action");
    let argv = [];
    const raw = trigger.getAttribute("data-action-args");
    if (raw) {
      try {
        argv = JSON.parse(raw);
      } catch {
        console.warn(`[NuGenUtils] data-action-args on ${name} is not valid JSON: ${raw}`);
      }
    }
    runAction(name, argv);
  }

  document.addEventListener("click", onActionClick);

  // Built-in actions used by partials
  registerAction("closeModal", closeModal);

  // ==========================================================================
  // Focus trap (used by modal & any component that needs it)
  // ==========================================================================

  const FOCUSABLE =
    'a[href], button:not([disabled]), textarea:not([disabled]), input:not([disabled]):not([type="hidden"]), select:not([disabled]), [tabindex]:not([tabindex="-1"])';

  function trapFocus(container) {
    if (!container) return () => {};
    function onKey(e) {
      if (e.key !== "Tab") return;
      const items = container.querySelectorAll(FOCUSABLE);
      if (items.length === 0) return;
      const first = items[0];
      const last = items[items.length - 1];
      if (e.shiftKey && document.activeElement === first) {
        e.preventDefault();
        last.focus();
      } else if (!e.shiftKey && document.activeElement === last) {
        e.preventDefault();
        first.focus();
      }
    }
    container.addEventListener("keydown", onKey);
    return () => container.removeEventListener("keydown", onKey);
  }

  // ==========================================================================
  // Keyboard shortcuts catalog
  //
  // Catalog entries appear in the `?` help overlay (Phase 7 polish renders).
  //   ⌘K / Ctrl+K → focus topbar search    (registered by topbar_search.js)
  //   /           → focus first input on page
  //   ⌘Enter      → submit active form
  //   g s/d/f/y/c/p/r/v → navigate to tools
  //   ?           → toggle help overlay
  //   Esc         → close modal / overlays
  // ==========================================================================

  const shortcuts = {
    catalog: [],
    add(key, handler, label = "") {
      this.catalog.push({ key, label });
      const matcher = parseKey(key);
      document.addEventListener("keydown", (e) => {
        if (isInTextInput(e.target) && !matcher.allowsInputs) return;
        if (matchesShortcut(e, matcher)) {
          e.preventDefault();
          handler(e);
        }
      });
    },
  };

  function parseKey(key) {
    // Accepts e.g. "Mod+K", "/", "?", "Esc", "g s"
    const allowsInputs = key === "Esc" || key.startsWith("Mod+");
    const tokens = key.split("+").map((t) => t.trim());
    return { tokens, allowsInputs, raw: key };
  }

  function matchesShortcut(e, m) {
    const wantMod = m.tokens.includes("Mod");
    const hasMod = e.metaKey || e.ctrlKey;
    if (wantMod !== hasMod) return false;
    const main = m.tokens[m.tokens.length - 1];
    if (main === "Enter") return e.key === "Enter";
    if (main === "Esc") return e.key === "Escape";
    return e.key.toLowerCase() === main.toLowerCase();
  }

  function isInTextInput(target) {
    if (!target) return false;
    const tag = target.tagName;
    return tag === "INPUT" || tag === "TEXTAREA" || target.isContentEditable;
  }

  // Default global shortcuts. Per-page `g s/d/f/...` navigation registers
  // here too so Phase-7 help overlay can describe them all.
  shortcuts.add(
    "Esc",
    () => closeModal(),
    "Close modal"
  );
  shortcuts.add(
    "/",
    () => {
      const input = document.querySelector(
        "main input:not([type=hidden]):not([type=file]), main textarea"
      );
      if (input) input.focus();
    },
    "Focus first input on page"
  );
  shortcuts.add(
    "Mod+Enter",
    () => {
      const form = document.querySelector("main form");
      if (form && typeof form.requestSubmit === "function") form.requestSubmit();
      else {
        const btn = document.querySelector("main button[type=submit], main [data-action]");
        if (btn) btn.click();
      }
    },
    "Submit active form"
  );

  // g+letter navigation (Vim-ish). Pressing 'g' arms a 1.5 s window for
  // the next character to map to a tool URL.
  const navMap = {
    h: "/",
    s: "/structure",
    d: "/descriptors",
    f: "/fingerprints",
    y: "/similarity",
    c: "/coordinates",
    p: "/properties",
    r: "/reactions",
    v: "/visualization",
  };

  let gArmed = false;
  let gTimer = 0;
  document.addEventListener("keydown", (e) => {
    if (isInTextInput(e.target)) return;
    if (gArmed) {
      const dest = navMap[e.key.toLowerCase()];
      gArmed = false;
      clearTimeout(gTimer);
      if (dest) {
        e.preventDefault();
        window.location.href = dest;
      }
      return;
    }
    if (e.key === "g" && !e.metaKey && !e.ctrlKey && !e.altKey) {
      gArmed = true;
      clearTimeout(gTimer);
      gTimer = setTimeout(() => (gArmed = false), 1500);
      return;
    }
    // Phase E: Linear / Slack pattern — typing any printable letter
    // outside of an input opens the topbar search palette and
    // pre-fills the typed character. Only fires when no modifier
    // is held and the topbar search exists on the page.
    if (
      e.key.length === 1 &&
      /^[a-z0-9]$/i.test(e.key) &&
      !e.metaKey && !e.ctrlKey && !e.altKey
    ) {
      const topbarInput = document.getElementById("topbar-search-input");
      if (topbarInput) {
        e.preventDefault();
        topbarInput.value = e.key;
        topbarInput.focus();
        topbarInput.dispatchEvent(new Event("input", { bubbles: true }));
      }
    }
  });

  // ==========================================================================
  // Sortable tables (data-sort-by on <th> elements)
  //
  // Templates opt in by emitting:
  //   <th data-sort-by="numeric">MW</th>     numeric (font-variant tabular-nums)
  //   <th data-sort-by="text">Name</th>      lexicographic
  //   <th data-sort-by="date">Added</th>     ISO date / timestamp
  // Click toggles asc / desc; the sort-direction goes onto the th
  // so the CSS `::after` arrow renders.
  // ==========================================================================

  document.addEventListener("click", (e) => {
    const th = e.target.closest("th[data-sort-by]");
    if (!th) return;
    const table = th.closest("table");
    if (!table) return;

    const headerRow = th.parentElement;
    const headers = Array.from(headerRow.children);
    const idx = headers.indexOf(th);
    const kind = th.getAttribute("data-sort-by");
    const currentDir = th.getAttribute("data-sort-direction");
    const nextDir = currentDir === "asc" ? "desc" : "asc";

    // Clear other headers' indicators.
    headers.forEach((h) => h.removeAttribute("data-sort-direction"));
    th.setAttribute("data-sort-direction", nextDir);

    const tbody = table.tBodies[0];
    if (!tbody) return;
    const rows = Array.from(tbody.rows);

    function valueOf(row) {
      const cell = row.cells[idx];
      if (!cell) return null;
      const raw = (cell.dataset.sortValue ?? cell.textContent ?? "").trim();
      if (kind === "numeric") return parseFloat(raw);
      if (kind === "date") return Date.parse(raw) || 0;
      return raw.toLowerCase();
    }

    rows.sort((a, b) => {
      const va = valueOf(a);
      const vb = valueOf(b);
      if (va == null) return 1;
      if (vb == null) return -1;
      if (va < vb) return nextDir === "asc" ? -1 : 1;
      if (va > vb) return nextDir === "asc" ? 1 : -1;
      return 0;
    });

    const frag = document.createDocumentFragment();
    rows.forEach((r) => frag.appendChild(r));
    tbody.appendChild(frag);
  });

  // ==========================================================================
  // Theme toggle persistence
  // ==========================================================================

  function applyStoredTheme() {
    const stored = localStorage.getItem("nug-theme");
    if (stored === "dark" || stored === "light") {
      document.documentElement.setAttribute("data-theme", stored);
    }
  }
  function toggleTheme() {
    const current = document.documentElement.getAttribute("data-theme") || "light";
    const next = current === "dark" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", next);
    try {
      localStorage.setItem("nug-theme", next);
    } catch {
      /* private browsing, ignore */
    }
    return next;
  }
  registerAction("toggleTheme", toggleTheme);
  applyStoredTheme();

  // ==========================================================================
  // Public API
  // ==========================================================================

  window.NuGenUtils = {
    escapeHtml,
    fetchJSON,
    showAlert,
    showLoading,
    hideLoading,
    openModal,
    closeModal,
    registerAction,
    runAction,
    trapFocus,
    shortcuts,
    toggleTheme,
  };
})();
