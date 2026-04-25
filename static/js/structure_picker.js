/**
 * structure_picker.js — molecule input helper for /structure.
 *
 * Phase 7 polish. The original plan called for a JSME or Ketcher
 * sketcher; both either (a) ship a multi-megabyte bundle (Ketcher,
 * OpenChemLib) or (b) need a CSP-incompatible dynamic loader chain
 * (JSME). The 80% UX win is letting users pick from common molecules
 * + paste a MOL block, neither of which needs an editor canvas.
 *
 * Public action (registered with NuGenUtils):
 *   openStructurePicker()  — opens a modal letting the user
 *     • pick a curated drug/scaffold (button → pastes SMILES)
 *     • paste a MOL block (calls /api/v1/structure/mol_to_smiles
 *       and writes the canonical SMILES back into #structure-input)
 *
 * Wired into structure.html via:
 *   <button data-action="openStructurePicker">…</button>
 */
(function () {
  "use strict";
  if (window.location.pathname !== "/structure") return;
  const Utils = window.NuGenUtils;
  if (!Utils) {
    console.warn("[structure_picker] NuGenUtils unavailable");
    return;
  }

  // ==========================================================================
  // Curated catalog
  // ==========================================================================

  const CATALOG = {
    "Common drugs": [
      { name: "Aspirin",       smiles: "CC(=O)OC1=CC=CC=C1C(=O)O" },
      { name: "Caffeine",      smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
      { name: "Ibuprofen",     smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" },
      { name: "Paracetamol",   smiles: "CC(=O)NC1=CC=C(O)C=C1" },
      { name: "Tamoxifen",     smiles: "CC/C(=C(/c1ccc(OCCN(C)C)cc1)c1ccccc1)c1ccccc1" },
      { name: "Penicillin G",  smiles: "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O" },
      { name: "Sildenafil",    smiles: "CCCc1nn(C)c2c1NC(c1cc(S(=O)(=O)N3CCN(C)CC3)ccc1OCC)=NC2=O" },
    ],
    "Solvents & reagents": [
      { name: "Water",         smiles: "O" },
      { name: "Ethanol",       smiles: "CCO" },
      { name: "Acetone",       smiles: "CC(=O)C" },
      { name: "DMSO",          smiles: "CS(=O)C" },
      { name: "Acetic acid",   smiles: "CC(=O)O" },
      { name: "DCM",           smiles: "ClCCl" },
      { name: "Toluene",       smiles: "Cc1ccccc1" },
    ],
    "Scaffolds": [
      { name: "Benzene",       smiles: "c1ccccc1" },
      { name: "Cyclohexane",   smiles: "C1CCCCC1" },
      { name: "Naphthalene",   smiles: "c1ccc2ccccc2c1" },
      { name: "Indole",        smiles: "c1ccc2[nH]ccc2c1" },
      { name: "Pyridine",      smiles: "c1ccncc1" },
      { name: "Imidazole",     smiles: "c1nc[nH]c1" },
      { name: "Furan",         smiles: "c1ccoc1" },
      { name: "Quinoline",     smiles: "c1ccc2ncccc2c1" },
    ],
  };

  // ==========================================================================
  // Modal builder
  // ==========================================================================

  function buildBody() {
    const root = document.createElement("div");
    root.className = "space-y-5 max-h-96 overflow-y-auto";

    Object.entries(CATALOG).forEach(([category, items]) => {
      const section = document.createElement("section");
      section.innerHTML = `<h3 class="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">${Utils.escapeHtml(category)}</h3>`;
      const grid = document.createElement("div");
      grid.className = "grid grid-cols-2 sm:grid-cols-3 gap-2";
      items.forEach((item) => {
        const btn = document.createElement("button");
        btn.type = "button";
        btn.className =
          "px-3 py-2 text-left text-sm bg-gray-50 border border-gray-200 rounded-lg hover:bg-primary-50 hover:border-primary-300 transition-colors";
        btn.innerHTML =
          `<div class="font-medium text-gray-900">${Utils.escapeHtml(item.name)}</div>` +
          `<div class="text-xs text-gray-500 font-mono truncate">${Utils.escapeHtml(item.smiles)}</div>`;
        btn.addEventListener("click", () => pickSmiles(item.smiles, item.name));
        grid.appendChild(btn);
      });
      section.appendChild(grid);
      root.appendChild(section);
    });

    // MOL-block paste pane
    const molSection = document.createElement("section");
    molSection.innerHTML = `
      <h3 class="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">Paste a MOL block</h3>
      <textarea id="picker-molblock"
                rows="6"
                class="w-full px-3 py-2 border border-gray-300 rounded-lg font-mono text-xs focus:ring-2 focus:ring-primary-500"
                placeholder="MOL block content here…"></textarea>
      <button type="button" id="picker-convert-molblock"
              class="mt-2 inline-flex items-center px-3 py-2 bg-primary-600 text-white rounded-lg text-sm hover:bg-primary-700">
        Convert to SMILES
      </button>
    `;
    root.appendChild(molSection);
    return root;
  }

  function pickSmiles(smiles, name) {
    const input = document.getElementById("structure-input");
    if (input) {
      input.value = smiles;
      input.dispatchEvent(new Event("input", { bubbles: true }));
    }
    Utils.showAlert(`Loaded ${name}`, "success", 1800);
    Utils.closeModal();
  }

  async function convertMolBlock() {
    const ta = document.getElementById("picker-molblock");
    if (!ta || !ta.value.trim()) {
      Utils.showAlert("Paste a MOL block first", "warning", 2000);
      return;
    }
    try {
      const data = await Utils.fetchJSON("/api/v1/structure/mol_to_smiles", {
        method: "POST",
        body: { mol_block: ta.value },
      });
      if (data && data.success && data.canonical_smiles) {
        pickSmiles(data.canonical_smiles, "MOL conversion");
      } else {
        Utils.showAlert(data?.error || "Conversion failed", "error", 3000);
      }
    } catch (err) {
      Utils.showAlert(err.message || "Conversion failed", "error", 3000);
    }
  }

  // ==========================================================================
  // Action registration
  // ==========================================================================

  Utils.registerAction("openStructurePicker", function () {
    Utils.openModal({
      title: "Pick a molecule",
      body: buildBody(),
      buttons: [{ label: "Close", action: "closeModal", variant: "secondary" }],
    });
    // Defer so the modal is in the DOM when we hook the convert button.
    setTimeout(() => {
      const btn = document.getElementById("picker-convert-molblock");
      if (btn) btn.addEventListener("click", convertMolBlock);
    }, 0);
  });
})();
