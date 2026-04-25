/**
 * structure.js — page module for "/structure".
 *
 * Phase-5 mechanical extraction of the inline <script> that previously
 * lived in templates/structure.html. Logic preserved verbatim except that
 * the local `showAlert` function delegates to NuGenUtils.showAlert so
 * dynamic toast HTML no longer contains inline `onclick=` (CSP-safe).
 *
 * Phase 7 polish may refactor `fetch` calls to use NuGenUtils.fetchJSON.
 */
(function () {
  "use strict";
  if (window.location.pathname !== "/structure") return;

document.addEventListener('DOMContentLoaded', function() {
    // Check for SMILES parameter from URL
    const urlParams = new URLSearchParams(window.location.search);
    const smilesParam = urlParams.get('smiles');
    if (smilesParam) {
        document.getElementById('structure-input').value = smilesParam;
        // Auto-convert if SMILES was provided
        convertStructure();
    }
    
    // Convert button handler
    document.getElementById('convert-btn').addEventListener('click', function() {
        convertStructure();
    });
    
    // Example button handler
    document.getElementById('example-btn').addEventListener('click', function() {
        loadStructureExample();
    });
    
    // Example molecule handlers
    document.querySelectorAll('.example-molecule').forEach(button => {
        button.addEventListener('click', function() {
            const smiles = this.getAttribute('data-smiles');
            document.getElementById('structure-input').value = smiles;
            convertStructure();
            // Examples often sit deep in the input column. Scroll the
            // result panel into view so the user sees the loading state
            // (and then the result) instead of staring at the button.
            document.getElementById('structure-display')
                ?.scrollIntoView({ behavior: 'smooth', block: 'start' });
        });
    });

    // Auto-convert on page load
    convertStructure();
});

async function convertStructure() {
    const smiles = document.getElementById('structure-input').value.trim();
    const format = document.getElementById('structure-format').value;
    const fileInput = document.getElementById('batch-structure-file');

    // Check if file upload is being used
    if (fileInput && fileInput.files.length > 0) {
        await convertBatchStructures();
        return;
    }

    if (!smiles) {
        showAlert('Please enter a SMILES string', 'warning');
        return;
    }
    
    try {
        // Get structure visualization first
        const vizResult = await fetch('/api/v1/visualization/draw_svg', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles, width: 300, height: 300 })
        });
        
        if (vizResult.ok) {
            const vizData = await vizResult.json();
            if (vizData.success) {
                document.getElementById('structure-display').innerHTML = vizData.svg;
            }
        }
        
        // Get conversion result - map format to correct endpoint
        let endpoint;
        switch (format) {
            // Format conversions
            case 'inchi':
                endpoint = '/api/v1/structure/smiles_to_inchi';
                break;
            case 'mol':
                endpoint = '/api/v1/structure/smiles_to_mol';
                break;
            case 'sdf':
                endpoint = '/api/v1/structure/smiles_to_sdf';
                break;
            case 'pdb':
                endpoint = '/api/v1/structure/smiles_to_pdb';
                break;
            case 'canonicalize':
                endpoint = '/api/v1/structure/canonicalize';
                break;
            // Hydrogen operations
            case 'add_hydrogens':
                endpoint = '/api/v1/structure/add_hydrogens';
                break;
            case 'remove_hydrogens':
                endpoint = '/api/v1/structure/remove_hydrogens';
                break;
            // Stereochemistry
            case 'stereochemistry':
                endpoint = '/api/v1/structure/stereochemistry';
                break;
            case 'enumerate_stereoisomers':
                endpoint = '/api/v1/structure/enumerate_stereoisomers';
                break;
            // Aromaticity
            case 'kekulize':
                endpoint = '/api/v1/structure/kekulize';
                break;
            case 'aromaticity':
                endpoint = '/api/v1/structure/aromaticity';
                break;
            // Fragment operations
            case 'fragment':
                endpoint = '/api/v1/structure/fragment';
                break;
            case 'murcko_scaffold':
                endpoint = '/api/v1/structure/murcko_scaffold';
                break;
            // Ring analysis
            case 'ring_info':
                endpoint = '/api/v1/structure/ring_info';
                break;
            // Cleanup
            case 'neutralize':
                endpoint = '/api/v1/structure/neutralize';
                break;
            case 'remove_fragments':
                endpoint = '/api/v1/structure/remove_fragments';
                break;
            case 'cleanup':
                endpoint = '/api/v1/structure/cleanup';
                break;
            // Validation
            case 'validate':
                endpoint = '/api/v1/structure/validate';
                break;
            default:
                endpoint = '/api/v1/structure/canonicalize';
        }
        
        const response = await fetch(endpoint, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles })
        });
        
        const data = await response.json();
        
        if (data.success) {
            displayResult(data, format);
            // Store the result globally for quick links
            window.lastStructureData = data;
            window.Recent?.add({ smiles, page: '/structure' });
        } else {
            showAlert(data.error || 'Conversion failed', 'danger');
        }
        
    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
    }
}

function displayResult(data, format) {
    const resultDiv = document.getElementById('structure-result');
    let content = '';

    switch (format) {
        case 'canonicalize':
            content = `
                <div class="bg-gray-50 rounded-lg p-4">
                    <h4 class="font-medium text-gray-900 mb-2">Canonical SMILES</h4>
                    <div class="font-mono text-sm bg-white p-3 rounded border">${data.canonical_smiles}</div>
                </div>
            `;
            break;
        case 'inchi':
            content = `
                <div class="bg-gray-50 rounded-lg p-4">
                    <h4 class="font-medium text-gray-900 mb-2">InChI</h4>
                    <div class="font-mono text-sm bg-white p-3 rounded border break-all">${data.inchi}</div>
                    <h4 class="font-medium text-gray-900 mb-2 mt-4">InChI Key</h4>
                    <div class="font-mono text-sm bg-white p-3 rounded border">${data.inchi_key}</div>
                </div>
            `;
            break;
        case 'mol':
            content = `
                <div class="bg-gray-50 rounded-lg p-4">
                    <h4 class="font-medium text-gray-900 mb-2">MOL Block</h4>
                    <pre class="font-mono text-xs bg-white p-3 rounded border overflow-auto max-h-60">${data.mol_block}</pre>
                </div>
            `;
            break;
        case 'sdf':
            content = `
                <div class="bg-gray-50 rounded-lg p-4">
                    <h4 class="font-medium text-gray-900 mb-2">SDF Block</h4>
                    <pre class="font-mono text-xs bg-white p-3 rounded border overflow-auto max-h-60">${data.sdf_block}</pre>
                </div>
            `;
            break;
        case 'pdb':
            content = `
                <div class="bg-gray-50 rounded-lg p-4">
                    <h4 class="font-medium text-gray-900 mb-2">PDB Block</h4>
                    <pre class="font-mono text-xs bg-white p-3 rounded border overflow-auto max-h-60">${data.pdb_block}</pre>
                </div>
            `;
            break;
        case 'add_hydrogens':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Add Hydrogens Result</h4>
                    <div class="grid grid-cols-2 gap-3 text-sm">
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Original Atoms</div>
                            <div class="font-bold text-lg">${data.num_atoms_original}</div>
                        </div>
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">With Hydrogens</div>
                            <div class="font-bold text-lg text-green-600">${data.num_atoms_with_h}</div>
                        </div>
                    </div>
                    <div class="bg-green-50 border border-green-200 p-3 rounded">
                        <div class="text-sm text-green-800 mb-1">Added ${data.num_hydrogens_added} hydrogen(s)</div>
                        <div class="font-mono text-xs break-all">${data.smiles_with_hydrogens}</div>
                    </div>
                </div>
            `;
            break;
        case 'remove_hydrogens':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Remove Hydrogens Result</h4>
                    <div class="grid grid-cols-2 gap-3 text-sm">
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Original Atoms</div>
                            <div class="font-bold text-lg">${data.num_atoms_original}</div>
                        </div>
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Without Hydrogens</div>
                            <div class="font-bold text-lg text-blue-600">${data.num_atoms_without_h}</div>
                        </div>
                    </div>
                    <div class="bg-blue-50 border border-blue-200 p-3 rounded">
                        <div class="text-sm text-blue-800 mb-1">Removed ${data.num_hydrogens_removed} hydrogen(s)</div>
                        <div class="font-mono text-xs break-all">${data.smiles_without_hydrogens}</div>
                    </div>
                </div>
            `;
            break;
        case 'stereochemistry':
            const stereo = data.stereochemistry;
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Stereochemistry Analysis</h4>
                    <div class="bg-white p-3 rounded border">
                        <div class="text-sm text-gray-600 mb-1">Chiral Centers: <span class="font-bold text-gray-900">${stereo.num_chiral_centers}</span></div>
                        <div class="text-sm ${stereo.has_stereochemistry ? 'text-green-700' : 'text-gray-700'}">
                            ${stereo.has_stereochemistry ? '✓ Contains stereochemistry' : 'No stereochemistry detected'}
                        </div>
                    </div>
                    ${stereo.chiral_centers.length > 0 ? `
                        <div class="bg-white p-3 rounded border">
                            <div class="text-sm font-medium mb-2">Chiral Centers:</div>
                            <div class="space-y-1">
                                ${stereo.chiral_centers.map(cc => `
                                    <div class="text-xs font-mono">Atom ${cc.atom_idx}: ${cc.chirality}</div>
                                `).join('')}
                            </div>
                        </div>
                    ` : ''}
                    <div class="bg-blue-50 border border-blue-200 p-3 rounded">
                        <div class="text-xs text-blue-800 mb-1">With Stereo</div>
                        <div class="font-mono text-xs break-all">${stereo.smiles_with_stereo}</div>
                    </div>
                </div>
            `;
            break;
        case 'enumerate_stereoisomers':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Stereoisomer Enumeration</h4>
                    <div class="bg-white p-3 rounded border">
                        <div class="text-sm text-gray-600">Found <span class="font-bold text-gray-900">${data.num_stereoisomers}</span> stereoisomer(s)</div>
                    </div>
                    <div class="space-y-2 max-h-80 overflow-y-auto">
                        ${data.stereoisomers.map(iso => `
                            <div class="bg-white p-3 rounded border hover:border-blue-300">
                                <div class="text-xs text-gray-600 mb-1">Isomer ${iso.index + 1}</div>
                                <div class="font-mono text-xs break-all">${iso.smiles}</div>
                            </div>
                        `).join('')}
                    </div>
                </div>
            `;
            break;
        case 'kekulize':
            content = `
                <div class="bg-gray-50 rounded-lg p-4">
                    <h4 class="font-medium text-gray-900 mb-2">Kekulized Structure</h4>
                    <div class="font-mono text-sm bg-white p-3 rounded border">${data.kekule_smiles}</div>
                </div>
            `;
            break;
        case 'aromaticity':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Aromaticity Analysis</h4>
                    <div class="grid grid-cols-2 gap-3 text-sm">
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Aromatic Atoms</div>
                            <div class="font-bold text-lg">${data.num_aromatic_atoms}</div>
                        </div>
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Aromatic Bonds</div>
                            <div class="font-bold text-lg">${data.num_aromatic_bonds}</div>
                        </div>
                    </div>
                    <div class="bg-purple-50 border border-purple-200 p-3 rounded">
                        <div class="text-sm text-purple-800 mb-1">Aromatic SMILES</div>
                        <div class="font-mono text-xs break-all">${data.aromatic_smiles}</div>
                    </div>
                </div>
            `;
            break;
        case 'fragment':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Fragment Analysis</h4>
                    <div class="bg-white p-3 rounded border">
                        <div class="text-sm text-gray-600">Found <span class="font-bold text-gray-900">${data.num_fragments}</span> fragment(s)</div>
                        <div class="text-xs ${data.is_single_fragment ? 'text-green-700' : 'text-orange-700'}">
                            ${data.is_single_fragment ? 'Single connected molecule' : 'Multiple disconnected fragments'}
                        </div>
                    </div>
                    <div class="space-y-2 max-h-80 overflow-y-auto">
                        ${data.fragments.map(frag => `
                            <div class="bg-white p-3 rounded border hover:border-blue-300">
                                <div class="flex justify-between mb-1">
                                    <div class="text-xs text-gray-600">Fragment ${frag.index + 1}</div>
                                    <div class="text-xs text-gray-600">${frag.num_atoms} atoms, MW: ${frag.molecular_weight.toFixed(2)}</div>
                                </div>
                                <div class="font-mono text-xs break-all">${frag.smiles}</div>
                            </div>
                        `).join('')}
                    </div>
                </div>
            `;
            break;
        case 'murcko_scaffold':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Murcko Scaffold</h4>
                    <div class="bg-white p-3 rounded border">
                        <div class="text-sm text-gray-600 mb-1">Scaffold (${data.scaffold_num_atoms} atoms)</div>
                        <div class="font-mono text-xs break-all">${data.scaffold_smiles}</div>
                    </div>
                    ${data.generic_scaffold_smiles ? `
                        <div class="bg-blue-50 border border-blue-200 p-3 rounded">
                            <div class="text-sm text-blue-800 mb-1">Generic Scaffold</div>
                            <div class="font-mono text-xs break-all">${data.generic_scaffold_smiles}</div>
                        </div>
                    ` : ''}
                </div>
            `;
            break;
        case 'ring_info':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Ring Analysis</h4>
                    <div class="grid grid-cols-2 gap-3 text-sm">
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Total Rings</div>
                            <div class="font-bold text-lg">${data.num_rings}</div>
                        </div>
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Aromatic Rings</div>
                            <div class="font-bold text-lg text-purple-600">${data.num_aromatic_rings}</div>
                        </div>
                    </div>
                    <div class="grid grid-cols-2 gap-3 text-sm">
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Aliphatic Rings</div>
                            <div class="font-bold text-lg">${data.num_aliphatic_rings}</div>
                        </div>
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Saturated Rings</div>
                            <div class="font-bold text-lg">${data.num_saturated_rings}</div>
                        </div>
                    </div>
                    ${data.ring_sizes.length > 0 ? `
                        <div class="bg-white p-3 rounded border">
                            <div class="text-sm font-medium mb-2">Ring Sizes</div>
                            <div class="text-xs">${data.ring_sizes.join(', ')}-membered rings</div>
                        </div>
                    ` : ''}
                </div>
            `;
            break;
        case 'neutralize':
            content = `
                <div class="bg-gray-50 rounded-lg p-4">
                    <h4 class="font-medium text-gray-900 mb-2">Neutralized Structure</h4>
                    <div class="font-mono text-sm bg-white p-3 rounded border">${data.neutralized_smiles}</div>
                </div>
            `;
            break;
        case 'remove_fragments':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Fragment Removal</h4>
                    <div class="grid grid-cols-2 gap-3 text-sm">
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Original Atoms</div>
                            <div class="font-bold text-lg">${data.original_num_atoms}</div>
                        </div>
                        <div class="bg-white p-3 rounded border">
                            <div class="text-gray-600">Cleaned Atoms</div>
                            <div class="font-bold text-lg text-green-600">${data.cleaned_num_atoms}</div>
                        </div>
                    </div>
                    <div class="bg-green-50 border border-green-200 p-3 rounded">
                        <div class="text-sm text-green-800 mb-1">Cleaned SMILES</div>
                        <div class="font-mono text-xs break-all">${data.cleaned_smiles}</div>
                    </div>
                </div>
            `;
            break;
        case 'cleanup':
            content = `
                <div class="bg-gray-50 rounded-lg p-4 space-y-3">
                    <h4 class="font-medium text-gray-900 mb-2">Full Molecular Cleanup</h4>
                    <div class="bg-green-50 border border-green-200 p-3 rounded">
                        <div class="text-sm text-green-800 mb-1">Cleaned SMILES</div>
                        <div class="font-mono text-xs break-all">${data.cleaned_smiles}</div>
                    </div>
                    <div class="bg-white p-3 rounded border">
                        <div class="text-sm text-gray-600">Formula: <span class="font-mono font-medium">${data.molecular_formula}</span></div>
                    </div>
                </div>
            `;
            break;
        case 'validate':
            const isValid = data.is_valid;
            content = `
                <div class="bg-gray-50 rounded-lg p-4">
                    <div class="flex items-center space-x-2 mb-2">
                        <i data-lucide="${isValid ? 'check-circle' : 'x-circle'}" class="w-5 h-5 ${isValid ? 'text-green-600' : 'text-red-600'}"></i>
                        <h4 class="font-medium text-gray-900">Validation Result</h4>
                    </div>
                    <div class="text-sm ${isValid ? 'text-green-700' : 'text-red-700'}">
                        ${isValid ? 'Valid molecular structure' : data.error || 'Invalid structure'}
                    </div>
                </div>
            `;
            break;
    }

    resultDiv.innerHTML = content;
    // Re-initialize icons
    lucide.createIcons();
}

function loadStructureExample() {
    const examples = ['CCO', 'CC(=O)OC1=CC=CC=C1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'];
    const randomExample = examples[Math.floor(Math.random() * examples.length)];
    document.getElementById('structure-input').value = randomExample;
    convertStructure();
}

async function convertBatchStructures() {
    const fileInput = document.getElementById('batch-structure-file');
    const format = document.getElementById('structure-format').value;
    const file = fileInput.files[0];

    if (!file) {
        showAlert('Please select a file', 'warning');
        return;
    }

    try {
        document.getElementById('structure-result').innerHTML = `
            <div class="flex flex-col items-center justify-center space-y-3">
                <div class="animate-spin rounded-full h-12 w-12 border-b-2 border-primary-600"></div>
                <div class="text-center">
                    <p class="text-gray-700 font-medium">Processing ${file.name}...</p>
                    <p class="text-sm text-gray-500">Converting structures</p>
                </div>
            </div>
        `;

        const formData = new FormData();
        formData.append('file', file);
        formData.append('output_format', format);

        const response = await fetch('/api/v1/structure/batch_convert', {
            method: 'POST',
            body: formData
        });

        const data = await response.json();

        if (data.success) {
            displayBatchStructureResults(data);
        } else {
            showAlert(data.error || 'Batch conversion failed', 'danger');
            document.getElementById('structure-result').innerHTML = '';
        }

    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
        document.getElementById('structure-result').innerHTML = '';
    }
}

function displayBatchStructureResults(data) {
    const results = data.results || [];
    const successful = data.num_successful || 0;
    const total = data.num_molecules || 0;
    const format = data.output_format || 'unknown';

    let html = `
        <div class="mb-4 bg-gradient-to-r from-blue-50 to-blue-100 p-4 rounded-lg border border-blue-200">
            <h3 class="text-lg font-semibold text-gray-900 mb-2">Batch Conversion Results</h3>
            <div class="grid grid-cols-3 gap-4 text-sm">
                <div>
                    <span class="text-gray-600">File:</span>
                    <span class="ml-2 font-medium">${data.file_name}</span>
                </div>
                <div>
                    <span class="text-gray-600">Format:</span>
                    <span class="ml-2 font-medium">${format.toUpperCase()}</span>
                </div>
                <div>
                    <span class="text-gray-600">Success:</span>
                    <span class="ml-2 font-bold text-green-700">${successful} / ${total}</span>
                </div>
            </div>
        </div>

        <div class="space-y-3 max-h-96 overflow-y-auto">
            ${results.slice(0, 20).map((result, idx) => `
                <div class="p-4 border border-gray-200 rounded-lg hover:border-blue-300 hover:shadow-sm transition-all">
                    <div class="flex justify-between items-start mb-2">
                        <div class="flex-1">
                            <h4 class="font-medium text-gray-900">${result.name}</h4>
                            <code class="text-xs text-gray-600 break-all">${result.input_smiles || result.smiles}</code>
                        </div>
                        ${result.error ? `
                            <span class="text-xs text-red-600">Error</span>
                        ` : `
                            <span class="text-xs text-green-600">✓ Converted</span>
                        `}
                    </div>
                    ${!result.error ? `
                        <div class="mt-2 bg-gray-50 rounded p-2">
                            <div class="text-xs font-medium text-gray-700 mb-1">${format.toUpperCase()} Output:</div>
                            <code class="text-xs text-gray-800 break-all block max-h-20 overflow-y-auto">
                                ${result.canonical_smiles || result.inchi || result.mol_block?.substring(0, 200) || result.output || 'N/A'}
                            </code>
                        </div>
                    ` : `
                        <div class="text-xs text-red-600 mt-1">Error: ${result.error}</div>
                    `}
                </div>
            `).join('')}
            ${results.length > 20 ? `
                <div class="text-center text-sm text-gray-500 py-3">
                    Showing 20 of ${results.length} results
                </div>
            ` : ''}
        </div>
    `;

    document.getElementById('structure-result').innerHTML = html;

    // Clear structure display when showing batch results
    document.getElementById('structure-display').innerHTML = `
        <div class="text-center text-gray-500">
            <i data-lucide="layers" class="w-16 h-16 mx-auto mb-4 text-gray-300"></i>
            <p>Batch conversion - see results below</p>
        </div>
    `;

    lucide.createIcons();
}

function showLoading(show) {
    const overlay = document.getElementById('loading-overlay');
    if (show) {
        overlay.classList.remove('hidden');
    } else {
        overlay.classList.add('hidden');
    }
}

function showAlert(message, type) {
  // Delegated to NuGenUtils.showAlert (Phase 5). The legacy `type`
  // values include 'danger' which we translate to 'error' for the
  // shared API. Falls back to alert() if NuGenUtils is unavailable
  // (e.g. utils.js failed to load).
  const mapped = (type === "danger") ? "error" : (type || "info");
  if (window.NuGenUtils && typeof window.NuGenUtils.showAlert === "function") {
    window.NuGenUtils.showAlert(message, mapped);
  } else {
    console.warn("[showAlert] NuGenUtils unavailable:", message);
  }
}

function navigateWithMolecule(targetPage) {
    const inputSmiles = document.getElementById('structure-input')?.value?.trim();
    
    // Try to get converted SMILES from the result display if available
    let convertedSmiles = null;
    const resultDiv = document.getElementById('structure-result');
    if (resultDiv && window.lastStructureData) {
        // Use the last conversion result if available
        if (window.lastStructureData.canonical_smiles) {
            convertedSmiles = window.lastStructureData.canonical_smiles;
        } else if (window.lastStructureData.smiles) {
            convertedSmiles = window.lastStructureData.smiles;
        }
    }
    
    // Use converted SMILES if available, otherwise use input SMILES
    const smiles = convertedSmiles || inputSmiles;
    
    if (smiles) {
        const url = new URL(`/${targetPage}`, window.location.origin);
        url.searchParams.set('smiles', smiles);
        window.location.href = url.toString();
    } else {
        showAlert('Please enter a molecule first', 'warning');
    }
}

// Phase 5: register data-action so the template buttons no longer
// need inline onclick="navigateWithMolecule(...)" handlers.
if (window.NuGenUtils && typeof window.NuGenUtils.registerAction === "function") {
  window.NuGenUtils.registerAction("navigateWithMolecule", navigateWithMolecule);
}
})();
