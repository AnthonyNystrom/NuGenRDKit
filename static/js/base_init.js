/**
 * base_init.js — small bootstrap loaded on every page.
 *
 * Two responsibilities:
 *   1. Tell lucide to scan the document and replace `<i data-lucide="…">`
 *      placeholder elements with the actual SVG icon. Lucide must already
 *      be on `window.lucide` (its CDN script tag is loaded just before us).
 *   2. Wire the topbar's mobile-menu toggle. Updates `aria-expanded`
 *      so the screen-reader announcement matches the visible state.
 *
 * Phase 5: extracted from the inline <script> at the bottom of
 * templates/base.html so CSP `script-src` no longer needs `'unsafe-inline'`.
 */
(function () {
  "use strict";

  function initLucide() {
    if (window.lucide && typeof window.lucide.createIcons === "function") {
      window.lucide.createIcons();
    }
  }

  function wireMobileMenu() {
    const button = document.getElementById("mobile-menu-button");
    const menu = document.getElementById("mobile-menu");
    if (!button || !menu) return;
    button.addEventListener("click", function () {
      const isHidden = menu.classList.contains("hidden");
      menu.classList.toggle("hidden");
      button.setAttribute("aria-expanded", isHidden ? "true" : "false");
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", function () {
      initLucide();
      wireMobileMenu();
    });
  } else {
    initLucide();
    wireMobileMenu();
  }
})();
