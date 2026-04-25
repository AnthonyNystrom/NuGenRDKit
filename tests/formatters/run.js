#!/usr/bin/env node
/**
 * tests/formatters/run.js — Phase 4 JS unit harness.
 *
 * Runs in plain Node (no jsdom, no npm dependencies). The module
 * scripts in static/js/* are written as IIFEs that attach to a global
 * `window` object — this harness builds a minimal `window` shim, then
 * sources each module via Node's `vm` module so they can be exercised
 * in isolation.
 *
 * Usage:
 *
 *     node tests/formatters/run.js
 *
 * CI runs this in the same workflow step as pytest. A non-zero exit
 * code fails CI.
 *
 * Coverage:
 *   - Formatters: number / percent / molecularFormula / formatBoolean /
 *     bytes / duration
 *   - NuGenUtils: escapeHtml + action registry
 *   - Workspace: add / list / remove / count / fillInto
 *   - ResultsCard: mount / showLoading / showError / showEmpty (smoke)
 */

const fs = require("fs");
const path = require("path");
const vm = require("vm");

let passed = 0;
let failed = 0;
const failures = [];

function assert(cond, message) {
  if (cond) {
    passed += 1;
  } else {
    failed += 1;
    failures.push(message);
    console.error("  ✗", message);
  }
}

function eq(actual, expected, label) {
  const same =
    actual === expected ||
    (Number.isNaN(actual) && Number.isNaN(expected)) ||
    JSON.stringify(actual) === JSON.stringify(expected);
  assert(same, `${label}: expected ${JSON.stringify(expected)}, got ${JSON.stringify(actual)}`);
}

function header(name) {
  console.log(`\n— ${name} —`);
}

// ============================================================================
// jsdom-lite shim
// ============================================================================

function makeShim() {
  const listeners = new Map();
  const elementsById = new Map();

  function makeElement(tag = "div") {
    const el = {
      tagName: (tag || "div").toUpperCase(),
      classList: new Set(),
      attributes: new Map(),
      children: [],
      textContent: "",
      innerHTML: "",
      type: "",
      value: "",
      tabIndex: -1,
      addEventListener(name, fn) {
        const key = `${this._uid || "_"}-${name}`;
        listeners.set(key, fn);
      },
      setAttribute(k, v) {
        this.attributes.set(k, String(v));
        if (k === "id") elementsById.set(v, this);
      },
      getAttribute(k) {
        return this.attributes.has(k) ? this.attributes.get(k) : null;
      },
      removeAttribute(k) {
        this.attributes.delete(k);
      },
      hasAttribute(k) {
        return this.attributes.has(k);
      },
      appendChild(child) {
        this.children.push(child);
        if (child.attributes && child.attributes.has("id")) {
          elementsById.set(child.attributes.get("id"), child);
        }
        return child;
      },
      querySelector() {
        return null;
      },
      querySelectorAll() {
        return [];
      },
      closest() {
        return null;
      },
      remove() {
        this.children = [];
      },
      scrollIntoView() {},
      focus() {},
      click() {},
      dispatchEvent() {},
    };
    el._uid = Math.random().toString(36).slice(2);
    return el;
  }

  const documentEl = {
    body: makeElement("body"),
    documentElement: makeElement("html"),
    createElement: (tag) => makeElement(tag),
    getElementById(id) {
      return elementsById.get(id) || null;
    },
    querySelector() {
      return null;
    },
    querySelectorAll() {
      return [];
    },
    addEventListener() {},
  };
  documentEl.documentElement.setAttribute("data-theme", "light");

  // sessionStorage / localStorage in-memory shims
  function memoryStorage() {
    const store = new Map();
    return {
      getItem: (k) => (store.has(k) ? store.get(k) : null),
      setItem: (k, v) => store.set(k, String(v)),
      removeItem: (k) => store.delete(k),
      clear: () => store.clear(),
      get length() {
        return store.size;
      },
    };
  }

  return {
    document: documentEl,
    sessionStorage: memoryStorage(),
    localStorage: memoryStorage(),
    location: { pathname: "/", href: "http://localhost/" },
    setTimeout,
    clearTimeout,
    Event: class {
      constructor(type, opts) {
        this.type = type;
        this.bubbles = !!(opts && opts.bubbles);
      }
    },
    AbortController: class {
      constructor() {
        this.signal = { aborted: false, addEventListener() {} };
      }
      abort() {
        this.signal.aborted = true;
      }
    },
    fetch: () => Promise.reject(new Error("fetch not implemented in shim")),
    console,
  };
}

function loadModule(file, sandbox) {
  const code = fs.readFileSync(path.join(__dirname, "../..", "static/js", file), "utf8");
  vm.runInContext(code, sandbox, { filename: file });
}

// ============================================================================
// Build context, source modules
// ============================================================================

const win = makeShim();
const sandbox = vm.createContext({ window: win, ...win });
sandbox.window = win;
sandbox.document = win.document;
sandbox.sessionStorage = win.sessionStorage;
sandbox.localStorage = win.localStorage;
sandbox.location = win.location;
sandbox.console = console;
sandbox.setTimeout = setTimeout;
sandbox.clearTimeout = clearTimeout;
sandbox.fetch = win.fetch;
sandbox.AbortController = win.AbortController;
sandbox.Event = win.Event;

loadModule("utils.js", sandbox);
loadModule("formatters.js", sandbox);
loadModule("workspace.js", sandbox);
loadModule("results_card.js", sandbox);
loadModule("topbar_search.js", sandbox);

const { NuGenUtils, Formatters, Workspace, ResultsCard } = win;

// ============================================================================
// Tests
// ============================================================================

header("Formatters.number");
eq(Formatters.number(1.23456), "1.2346", "1.23456 → 4 digits trimmed");
eq(Formatters.number(1234), "1234", "integer stays integer");
eq(Formatters.number(0), "0", "zero");
eq(Formatters.number(null), "—", "null fallback");
eq(Formatters.number(NaN), "—", "NaN fallback");
eq(Formatters.number(undefined), "—", "undefined fallback");
eq(Formatters.number(-2.5), "-2.5", "negative float");
eq(Formatters.number(1e7).startsWith("1.00e"), true, "very large → scientific");

header("Formatters.percent");
eq(Formatters.percent(0.5), "50.0%", "0.5 → 50%");
eq(Formatters.percent(0.1234, { digits: 2 }), "12.34%", "0.1234 with 2 digits");
eq(Formatters.percent(75), "75.0%", "75 (already a percent)");
eq(Formatters.percent(null), "—", "null fallback");

header("Formatters.molecularFormula");
eq(Formatters.molecularFormula("C9H8O4"), "C<sub>9</sub>H<sub>8</sub>O<sub>4</sub>", "aspirin");
eq(Formatters.molecularFormula("CH4"), "CH<sub>4</sub>", "methane");
eq(Formatters.molecularFormula(""), "", "empty string");
eq(Formatters.molecularFormula(null), "", "null safe");

header("Formatters.formatBoolean");
eq(Formatters.formatBoolean(true), "Yes", "true default");
eq(Formatters.formatBoolean(false), "No", "false default");
eq(Formatters.formatBoolean(true, { yes: "On", no: "Off" }), "On", "custom yes label");
eq(Formatters.formatBoolean(0), "No", "falsy is no");

header("Formatters.bytes");
eq(Formatters.bytes(0), "0 B", "zero bytes");
eq(Formatters.bytes(512), "512 B", "small");
eq(Formatters.bytes(2048), "2 KB", "round KB");
eq(Formatters.bytes(1024 * 1024 * 3.5), "3.5 MB", "fractional MB");
eq(Formatters.bytes(-1), "—", "negative fallback");
eq(Formatters.bytes(null), "—", "null fallback");

header("Formatters.duration");
eq(Formatters.duration(0), "0 ms", "zero ms");
eq(Formatters.duration(450), "450 ms", "ms range");
eq(Formatters.duration(1500), "1.50 s", "seconds range");
eq(Formatters.duration(125_000), "2m 5s", "minutes:seconds");
eq(Formatters.duration(null), "—", "null fallback");

header("NuGenUtils.escapeHtml");
eq(NuGenUtils.escapeHtml("<script>alert(1)</script>"), "&lt;script&gt;alert(1)&lt;/script&gt;", "tags escaped");
eq(NuGenUtils.escapeHtml('hello "world"'), "hello &quot;world&quot;", "double quotes");
eq(NuGenUtils.escapeHtml("o'reilly"), "o&#39;reilly", "single quote");
eq(NuGenUtils.escapeHtml("a & b"), "a &amp; b", "ampersand");
eq(NuGenUtils.escapeHtml(null), "", "null safe");
eq(NuGenUtils.escapeHtml(undefined), "", "undefined safe");
eq(NuGenUtils.escapeHtml(42), "42", "non-string coerces");

header("NuGenUtils action registry");
let actionCalls = 0;
NuGenUtils.registerAction("nugTestAction", (a, b) => {
  actionCalls += a + b;
});
NuGenUtils.runAction("nugTestAction", [1, 2]);
eq(actionCalls, 3, "registered action runs with argv");
NuGenUtils.runAction("missingAction", []);
eq(actionCalls, 3, "missing action is no-op");

assert(typeof NuGenUtils.fetchJSON === "function", "fetchJSON exposed");
assert(typeof NuGenUtils.showAlert === "function", "showAlert exposed");
assert(typeof NuGenUtils.openModal === "function", "openModal exposed");
assert(typeof NuGenUtils.closeModal === "function", "closeModal exposed");
assert(typeof NuGenUtils.trapFocus === "function", "trapFocus exposed");
assert(NuGenUtils.shortcuts && Array.isArray(NuGenUtils.shortcuts.catalog), "shortcuts.catalog exposed");

header("Workspace persistence");
Workspace.clear();
eq(Workspace.count(), 0, "starts empty");
const a = Workspace.add({ type: "smiles", value: "CCO", name: "Ethanol" });
const b = Workspace.add({ type: "smiles", value: "CC(=O)OC1=CC=CC=C1C(=O)O", name: "Aspirin" });
eq(Workspace.count(), 2, "after two adds");
eq(Workspace.list()[0].name, "Aspirin", "most-recent first");
eq(Workspace.get(a.id).value, "CCO", "get by id");
assert(Workspace.remove(a.id), "remove returns true");
eq(Workspace.count(), 1, "after remove");
assert(!Workspace.remove("does-not-exist"), "remove of missing returns false");

let onChangeCalls = 0;
const unsubscribe = Workspace.onChange(() => {
  onChangeCalls += 1;
});
Workspace.add({ value: "CN" });
eq(onChangeCalls, 1, "onChange fired on add");
unsubscribe();
Workspace.add({ value: "C" });
eq(onChangeCalls, 1, "unsubscribed callback no longer fires");

Workspace.clear();
eq(Workspace.count(), 0, "clear resets");

header("ResultsCard surface");
assert(typeof ResultsCard.mount === "function", "mount exposed");
assert(typeof ResultsCard.showLoading === "function", "showLoading exposed");
assert(typeof ResultsCard.showError === "function", "showError exposed");
assert(typeof ResultsCard.showEmpty === "function", "showEmpty exposed");

// mount onto a stub element to ensure no thrown error
const stub = win.document.createElement("div");
stub.setAttribute("id", "rc-stub");
win.document.body.appendChild(stub);
let mountThrew = false;
try {
  ResultsCard.mount("rc-stub", { title: "Hello", summary: "<p>ok</p>", raw: "{}" });
} catch (err) {
  mountThrew = err;
}
assert(!mountThrew, `mount() runs without throwing (${mountThrew && mountThrew.message})`);
assert(typeof stub.innerHTML === "string", "mount wrote innerHTML");
assert(stub.innerHTML.includes("Hello"), "rendered title");

// ============================================================================
// Summary
// ============================================================================

console.log(`\n${passed} passed, ${failed} failed`);
if (failed > 0) {
  console.error(`\nFailures:\n${failures.map((f) => "  " + f).join("\n")}`);
  process.exit(1);
}
