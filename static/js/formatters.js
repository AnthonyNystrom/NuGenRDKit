/**
 * Formatters — small pure functions used by ResultsCard and per-page
 * scripts. Keeping them in one module makes them trivially testable
 * (see tests/formatters/run.js).
 *
 * Public surface:
 *   Formatters.number(value, options)
 *   Formatters.percent(value, options)
 *   Formatters.molecularFormula(formula)
 *   Formatters.formatBoolean(value, { yes, no })
 *   Formatters.bytes(byteCount)
 *   Formatters.duration(milliseconds)
 *
 * Each function returns a string and never throws on null/undefined.
 */
(function () {
  "use strict";

  function isNum(v) {
    return typeof v === "number" && Number.isFinite(v);
  }

  function number(value, { digits = 4, fallback = "—" } = {}) {
    if (!isNum(value)) return fallback;
    if (Math.abs(value) >= 1_000_000) return value.toExponential(2);
    if (Number.isInteger(value)) return value.toString();
    return value.toFixed(digits).replace(/\.?0+$/, "");
  }

  function percent(value, { digits = 1, fallback = "—" } = {}) {
    if (!isNum(value)) return fallback;
    // Accepts 0..1 (proportion) or 0..100 (already a percent). > 1.5 ⇒
    // treat as already in percent.
    const v = value > 1.5 ? value : value * 100;
    return `${v.toFixed(digits)}%`;
  }

  /**
   * Add <sub> subscripting to digit runs in a Hill-order molecular formula:
   *   "C9H8O4" → "C<sub>9</sub>H<sub>8</sub>O<sub>4</sub>"
   * Returns HTML — caller is responsible for trusting the source.
   */
  function molecularFormula(formula) {
    if (!formula) return "";
    return String(formula).replace(/(\d+)/g, "<sub>$1</sub>");
  }

  function formatBoolean(value, { yes = "Yes", no = "No" } = {}) {
    return value ? yes : no;
  }

  function bytes(byteCount) {
    if (!isNum(byteCount) || byteCount < 0) return "—";
    const units = ["B", "KB", "MB", "GB", "TB"];
    let i = 0;
    let v = byteCount;
    while (v >= 1024 && i < units.length - 1) {
      v /= 1024;
      i += 1;
    }
    return `${v.toFixed(v >= 100 || v === Math.floor(v) ? 0 : 1)} ${units[i]}`;
  }

  function duration(milliseconds) {
    if (!isNum(milliseconds) || milliseconds < 0) return "—";
    if (milliseconds < 1000) return `${Math.round(milliseconds)} ms`;
    if (milliseconds < 60_000) return `${(milliseconds / 1000).toFixed(2)} s`;
    const minutes = Math.floor(milliseconds / 60_000);
    const seconds = Math.floor((milliseconds % 60_000) / 1000);
    return `${minutes}m ${seconds}s`;
  }

  window.Formatters = {
    number,
    percent,
    molecularFormula,
    formatBoolean,
    bytes,
    duration,
  };
})();
