/**
 * pagination.js — shared paginated-list renderer.
 *
 * Replaces the per-page `slice(0, 20)` + "Showing X of Y results" footer
 * pattern with real prev/next + page-size controls. Each call mounts an
 * isolated paginator into a container and re-renders the visible window
 * without touching the rest of the page.
 *
 * Public API (window.NuGenPagination):
 *   render(opts) → { destroy(), refresh(items?) }
 *
 * opts:
 *   container  — HTMLElement to mount into
 *   items      — array of records
 *   pageSize   — initial page size (default 20)
 *   pageSizes  — selectable sizes (default [10, 20, 50, 100])
 *   renderItem — fn(item, index) → HTML string
 *   summary    — fn({start, end, total}) → optional header HTML
 */
(function () {
  "use strict";

  const escapeHtml = (s) =>
    (window.NuGenUtils && window.NuGenUtils.escapeHtml)
      ? window.NuGenUtils.escapeHtml(s)
      : String(s).replace(/[&<>"']/g, (c) => (
          { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
        ));

  function render(opts) {
    const {
      container,
      items: initialItems = [],
      pageSize: initialSize = 20,
      pageSizes = [10, 20, 50, 100],
      renderItem,
      summary,
      itemsClass = "space-y-3",
    } = opts;

    if (!container || typeof renderItem !== "function") {
      throw new Error("NuGenPagination.render: container + renderItem required");
    }

    let items = initialItems.slice();
    let size = initialSize;
    let page = 1;

    container.innerHTML = `
      <div class="nug-paginator-summary mb-2"></div>
      <div class="nug-paginator-items ${itemsClass}"></div>
      <div class="nug-paginator-controls mt-3 flex items-center justify-between gap-3 text-xs text-gray-600 border-t border-gray-200 pt-3"></div>
    `;
    const summaryEl = container.querySelector(".nug-paginator-summary");
    const listEl = container.querySelector(".nug-paginator-items");
    const ctrlEl = container.querySelector(".nug-paginator-controls");

    function pageCount() {
      return Math.max(1, Math.ceil(items.length / size));
    }

    function clampPage() {
      const max = pageCount();
      if (page > max) page = max;
      if (page < 1) page = 1;
    }

    function paint() {
      clampPage();
      const start = (page - 1) * size;
      const end = Math.min(items.length, start + size);
      const slice = items.slice(start, end);

      summaryEl.innerHTML = summary
        ? summary({ start: items.length ? start + 1 : 0, end, total: items.length })
        : "";

      listEl.innerHTML = slice.map((item, idx) => renderItem(item, start + idx)).join("");

      const max = pageCount();
      const sizeOpts = pageSizes
        .map((n) => `<option value="${n}"${n === size ? " selected" : ""}>${n}</option>`)
        .join("");

      ctrlEl.innerHTML = `
        <div>
          ${items.length === 0
            ? "No results"
            : `Showing ${escapeHtml(start + 1)}–${escapeHtml(end)} of ${escapeHtml(items.length)}`}
        </div>
        <div class="flex items-center gap-2">
          <label>
            <span class="sr-only">Rows per page</span>
            <select class="nug-paginator-size form-select" style="height:24px;padding:0 4px;font-size:11px;">
              ${sizeOpts}
            </select>
          </label>
          <button type="button" class="nug-paginator-prev btn btn-secondary"
                  style="height:24px;padding:0 8px;font-size:11px;"
                  ${page <= 1 ? "disabled" : ""}>Prev</button>
          <span>Page ${escapeHtml(page)} / ${escapeHtml(max)}</span>
          <button type="button" class="nug-paginator-next btn btn-secondary"
                  style="height:24px;padding:0 8px;font-size:11px;"
                  ${page >= max ? "disabled" : ""}>Next</button>
        </div>
      `;

      // Re-render any lucide icons inside freshly inserted markup.
      if (window.lucide && typeof window.lucide.createIcons === "function") {
        window.lucide.createIcons();
      }
    }

    function onClick(e) {
      const t = e.target;
      if (!t || !(t instanceof Element)) return;
      if (t.classList.contains("nug-paginator-prev")) {
        page -= 1;
        paint();
      } else if (t.classList.contains("nug-paginator-next")) {
        page += 1;
        paint();
      }
    }
    function onChange(e) {
      const t = e.target;
      if (t && t.classList && t.classList.contains("nug-paginator-size")) {
        const n = parseInt(t.value, 10);
        if (Number.isFinite(n) && n > 0) {
          size = n;
          page = 1;
          paint();
        }
      }
    }
    container.addEventListener("click", onClick);
    container.addEventListener("change", onChange);

    paint();

    return {
      refresh(next) {
        if (Array.isArray(next)) items = next.slice();
        page = 1;
        paint();
      },
      destroy() {
        container.removeEventListener("click", onClick);
        container.removeEventListener("change", onChange);
        container.innerHTML = "";
      },
    };
  }

  window.NuGenPagination = { render };
})();
