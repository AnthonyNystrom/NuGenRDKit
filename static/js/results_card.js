/**
 * ResultsCard — replaces ad-hoc result rendering on every tool page.
 *
 * Public surface (kept stable across phases):
 *
 *   ResultsCard.mount(id, spec)
 *     spec keys: { title, meta, summary, details, raw, downloads,
 *                  workspaceItem, copyText, prefix, extraTabs }
 *
 *   ResultsCard.showLoading(id, { title })
 *   ResultsCard.showError(id, { title, message })
 *   ResultsCard.showEmpty(id, { icon, title, body })
 *
 * The mount target is whatever Jinja element has id=`id` (typically
 * inside _partials/_results_panel.html). This module owns the inner
 * markup; templates declare only the empty `<div id="…">` shell.
 *
 * `summary`, `details`, `raw`, and `extraTabs[].html` may contain HTML
 * but the caller is responsible for escaping. Use NuGenUtils.escapeHtml
 * for any user-supplied string.
 *
 * Requires: window.NuGenUtils (utils.js loaded first).
 */
(function () {
  "use strict";

  const escape = (s) => (window.NuGenUtils?.escapeHtml || String)(s);

  // ==========================================================================
  // Public methods
  // ==========================================================================

  function mount(id, spec = {}) {
    const target = document.getElementById(id);
    if (!target) {
      console.warn(`[ResultsCard] No element with id="${id}"`);
      return;
    }
    const {
      title = "Result",
      meta = "",
      summary = "",
      details = "",
      raw = "",
      downloads = [],
      workspaceItem = null,
      copyText = "",
      prefix = `${id}-rc`,
      extraTabs = [],
      // Phase F+: a list of cross-tool destinations the user can send
      // this molecule to. Either:
      //   sendTo: ["smiles string"]            ← shorthand: build default list
      //   sendTo: { smiles: "…", except: [] }  ← shorthand variant
      //   sendTo: [{label, page, smiles}]      ← explicit
      sendTo = null,
    } = spec;

    const tabs = buildTabs({ summary, details, raw, extraTabs }, prefix);

    target.innerHTML = `
      <article class="results-card border border-gray-200 rounded-lg bg-white">
        <header class="flex flex-wrap items-start justify-between gap-3 px-4 py-3 border-b border-gray-200">
          <div class="min-w-0">
            <h3 class="text-base font-semibold text-gray-900">${escape(title)}</h3>
            ${meta ? `<p class="text-xs text-gray-500 mt-1">${escape(meta)}</p>` : ""}
          </div>
          <div class="flex items-center gap-2 ml-auto" data-rc-actions></div>
        </header>
        ${
          tabs.headers
            ? `<nav class="flex gap-1 px-4 pt-2 border-b border-gray-200" role="tablist">${tabs.headers}</nav>`
            : ""
        }
        <div class="rc-body p-4">${tabs.panels}</div>
      </article>
    `;

    // Action buttons (right-aligned in header)
    const actions = target.querySelector("[data-rc-actions]");
    if (copyText) appendActionButton(actions, "Copy", "copy", () => copyToClipboard(copyText));
    if (workspaceItem)
      appendActionButton(actions, "Save", "save", () => {
        if (window.NuGenUtils?.runAction) {
          window.NuGenUtils.runAction("addToWorkspace", [workspaceItem]);
        }
      });
    if (sendTo) appendSendToMenu(actions, sendTo);
    downloads.forEach((d) => {
      const btn = appendActionButton(actions, d.label || "Download", "download");
      if (d.href) btn.setAttribute("href", d.href);
      if (d.onClick) btn.addEventListener("click", d.onClick);
      if (d.download) btn.setAttribute("download", d.download);
    });

    // Wire tab switching
    if (tabs.headers) {
      target.querySelectorAll("[data-rc-tab]").forEach((tab) => {
        tab.addEventListener("click", () => activateTab(target, tab.dataset.rcTab, prefix));
      });
    }

    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
    target.scrollIntoView({ behavior: "smooth", block: "nearest" });
  }

  function showLoading(id, { title = "Loading…" } = {}) {
    const target = document.getElementById(id);
    if (!target) return;
    target.innerHTML = `
      <article class="results-card border border-gray-200 rounded-lg bg-white p-8 text-center">
        <div class="inline-flex items-center text-gray-500">
          <span class="inline-block w-4 h-4 mr-2 border-2 border-gray-300 border-t-primary-500 rounded-full animate-spin"></span>
          <span class="text-sm">${escape(title)}</span>
        </div>
      </article>`;
  }

  function showError(id, { title = "Error", message = "" } = {}) {
    const target = document.getElementById(id);
    if (!target) return;
    target.innerHTML = `
      <article class="results-card border border-red-200 rounded-lg bg-red-50">
        <header class="px-4 py-3 border-b border-red-200">
          <h3 class="text-base font-semibold text-red-700">${escape(title)}</h3>
        </header>
        <div class="p-4 text-sm text-red-700">${escape(message)}</div>
      </article>`;
  }

  function showEmpty(id, { icon = "inbox", title = "No results", body = "" } = {}) {
    const target = document.getElementById(id);
    if (!target) return;
    target.innerHTML = `
      <div class="rc-empty text-center py-12">
        <i data-lucide="${escape(icon)}" class="w-12 h-12 mx-auto text-gray-300 mb-3"></i>
        <p class="text-sm font-medium text-gray-700">${escape(title)}</p>
        ${body ? `<p class="text-sm text-gray-500 mt-1">${escape(body)}</p>` : ""}
      </div>`;
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
  }

  // ==========================================================================
  // Internals
  // ==========================================================================

  function buildTabs({ summary, details, raw, extraTabs }, prefix) {
    const tabs = [];
    if (summary) tabs.push({ key: "summary", label: "Summary", html: summary });
    if (details) tabs.push({ key: "details", label: "Details", html: details });
    extraTabs.forEach((t) => tabs.push(t));
    if (raw) tabs.push({ key: "raw", label: "Raw", html: `<pre class="text-xs font-mono whitespace-pre-wrap break-all">${escape(raw)}</pre>` });

    if (tabs.length === 0) return { headers: "", panels: "" };
    if (tabs.length === 1) return { headers: "", panels: `<div>${tabs[0].html}</div>` };

    const headers = tabs
      .map(
        (t, i) => `
          <button type="button"
                  role="tab"
                  data-rc-tab="${t.key}"
                  id="${prefix}-tab-${t.key}"
                  aria-controls="${prefix}-panel-${t.key}"
                  aria-selected="${i === 0}"
                  class="rc-tab px-3 py-2 text-sm font-medium border-b-2 ${
                    i === 0
                      ? "border-primary-600 text-primary-700"
                      : "border-transparent text-gray-500 hover:text-gray-900"
                  }">
            ${escape(t.label)}
          </button>`
      )
      .join("");

    const panels = tabs
      .map(
        (t, i) => `
          <div role="tabpanel"
               id="${prefix}-panel-${t.key}"
               aria-labelledby="${prefix}-tab-${t.key}"
               class="rc-panel ${i === 0 ? "" : "hidden"}">
            ${t.html}
          </div>`
      )
      .join("");

    return { headers, panels };
  }

  function activateTab(target, key, prefix) {
    target.querySelectorAll(".rc-tab").forEach((tab) => {
      const active = tab.dataset.rcTab === key;
      tab.classList.toggle("border-primary-600", active);
      tab.classList.toggle("text-primary-700", active);
      tab.classList.toggle("border-transparent", !active);
      tab.classList.toggle("text-gray-500", !active);
      tab.setAttribute("aria-selected", String(active));
    });
    target.querySelectorAll(".rc-panel").forEach((panel) => {
      panel.classList.toggle("hidden", panel.id !== `${prefix}-panel-${key}`);
    });
  }

  function appendActionButton(container, label, icon, onClick) {
    const btn = document.createElement(onClick ? "button" : "a");
    btn.className = "rc-action-btn";
    if (btn.tagName === "BUTTON") btn.type = "button";
    btn.innerHTML = `<i data-lucide="${escape(icon)}" class="w-3 h-3"></i><span>${escape(label)}</span>`;
    if (onClick) btn.addEventListener("click", onClick);
    container.appendChild(btn);
    return btn;
  }

  // Build the "Send to <tool>" dropdown. `spec` may be a string (the
  // SMILES) → expand to default destinations, or {smiles, except: [pages]}
  // shorthand, or an explicit array of {label, page, smiles}.
  const DEFAULT_SEND_TO = [
    { label: "Structure",     page: "/structure"    , icon: "box"          },
    { label: "Descriptors",   page: "/descriptors"  , icon: "bar-chart-3"  },
    { label: "Fingerprints",  page: "/fingerprints" , icon: "fingerprint"  },
    { label: "Similarity",    page: "/similarity"   , icon: "search"       },
    { label: "3D Coords",     page: "/coordinates"  , icon: "orbit"        },
    { label: "Properties",    page: "/properties"   , icon: "flask-conical"},
    { label: "Visualization", page: "/visualization", icon: "image"        },
  ];

  function appendSendToMenu(container, spec) {
    let smiles, except = [], items = null;
    if (typeof spec === "string") {
      smiles = spec;
    } else if (Array.isArray(spec)) {
      items = spec;
    } else if (spec && typeof spec === "object") {
      smiles = spec.smiles;
      except = spec.except || [];
      items = spec.items;
    }
    if (!items) {
      const here = window.location.pathname;
      items = DEFAULT_SEND_TO
        .filter((d) => d.page !== here && !except.includes(d.page))
        .map((d) => ({ ...d, smiles }));
    }
    if (!items.length) return;

    const wrap = document.createElement("div");
    wrap.className = "rc-sendto-wrap";
    wrap.innerHTML = `
      <button type="button" class="rc-action-btn" aria-haspopup="menu" aria-expanded="false">
        <i data-lucide="send" class="w-3 h-3"></i><span>Send to</span>
        <i data-lucide="chevron-down" class="w-3 h-3"></i>
      </button>
      <ul class="rc-sendto-menu" role="menu" hidden></ul>
    `;
    const trigger = wrap.querySelector("button");
    const menu = wrap.querySelector(".rc-sendto-menu");

    items.forEach((it) => {
      const li = document.createElement("li");
      li.role = "menuitem";
      const a = document.createElement("a");
      a.href = it.page + (it.smiles ? `?smiles=${encodeURIComponent(it.smiles)}` : "");
      a.innerHTML = `<i data-lucide="${escape(it.icon || "arrow-right")}" class="w-3 h-3"></i><span>${escape(it.label)}</span>`;
      li.appendChild(a);
      menu.appendChild(li);
    });

    trigger.addEventListener("click", () => {
      const open = menu.hidden === false;
      menu.hidden = open;
      trigger.setAttribute("aria-expanded", open ? "false" : "true");
    });
    // Close on outside click.
    document.addEventListener("click", (e) => {
      if (!wrap.contains(e.target)) {
        menu.hidden = true;
        trigger.setAttribute("aria-expanded", "false");
      }
    });
    container.appendChild(wrap);
  }

  function copyToClipboard(text) {
    if (!navigator.clipboard) {
      window.NuGenUtils?.showAlert("Clipboard not available", "error", 2500);
      return;
    }
    navigator.clipboard
      .writeText(text)
      .then(() => window.NuGenUtils?.showAlert("Copied", "success", 1500))
      .catch(() => window.NuGenUtils?.showAlert("Copy failed", "error", 2500));
  }

  // ==========================================================================
  // Public API
  // ==========================================================================

  window.ResultsCard = { mount, showLoading, showError, showEmpty };
})();
