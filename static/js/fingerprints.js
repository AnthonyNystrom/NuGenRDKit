/**
 * fingerprints.js — page module for "/fingerprints".
 *
 * Phase-5 mechanical extraction of the inline <script> that previously
 * lived in templates/fingerprints.html. Logic preserved verbatim except that
 * the local `showAlert` function delegates to NuGenUtils.showAlert so
 * dynamic toast HTML no longer contains inline `onclick=` (CSP-safe).
 *
 * Phase 7 polish may refactor `fetch` calls to use NuGenUtils.fetchJSON.
 */
(function () {
  "use strict";
  if (window.location.pathname !== "/fingerprints") return;

document.addEventListener('DOMContentLoaded', function() {
    // Check for SMILES parameter from URL
    const urlParams = new URLSearchParams(window.location.search);
    const smilesParam = urlParams.get('smiles');
    if (smilesParam) {
        document.getElementById('fingerprint-input').value = smilesParam;
        // Auto-generate if SMILES was provided
        generateFingerprint();
    }
    
    // Generate button handler
    document.getElementById('generate-btn').addEventListener('click', function() {
        generateFingerprint();
    });
    
    // Compare button handler
    document.getElementById('compare-btn').addEventListener('click', function() {
        compareFingerprints();
    });
    
    // Export button handler
    document.getElementById('export-btn').addEventListener('click', function() {
        exportFingerprint();
    });
    
    // Example molecule handlers
    document.querySelectorAll('.example-molecule').forEach(button => {
        button.addEventListener('click', function() {
            const smiles = this.getAttribute('data-smiles');
            document.getElementById('fingerprint-input').value = smiles;
            generateFingerprint();
            document.getElementById('fingerprint-display')
                ?.scrollIntoView({ behavior: 'smooth', block: 'start' });
        });
    });
    
    // Fingerprint type change handler
    document.getElementById('fingerprint-type').addEventListener('change', function() {
        updateFingerprintParams();
    });
    
    updateFingerprintParams();
    generateFingerprint(); // Auto-generate on page load
});

function updateFingerprintParams() {
    const fpType = document.getElementById('fingerprint-type').value;
    const paramsDiv = document.getElementById('fingerprint-params');
    
    // Show/hide parameters based on fingerprint type
    if (fpType === 'maccs') {
        paramsDiv.classList.add('hidden');
    } else {
        paramsDiv.classList.remove('hidden');
        
        if (fpType === 'avalon') {
            document.getElementById('fp-radius').parentElement.classList.add('hidden');
        } else {
            document.getElementById('fp-radius').parentElement.classList.remove('hidden');
        }
    }
}

async function generateFingerprint() {
    const smiles = document.getElementById('fingerprint-input').value.trim();
    const fpType = document.getElementById('fingerprint-type').value;
    const fileInput = document.getElementById('batch-fingerprints-file');

    // Check if file upload is being used
    if (fileInput && fileInput.files.length > 0) {
        await calculateBatchFingerprints();
        return;
    }

    if (!smiles) {
        showAlert('Please enter a SMILES string', 'warning');
        return;
    }
    
    try {
        const params = {
            smiles: smiles
        };
        
        // Add parameters based on fingerprint type
        if (fpType !== 'maccs') {
            if (fpType !== 'avalon') {
                params.radius = parseInt(document.getElementById('fp-radius').value);
            }
            params.n_bits = parseInt(document.getElementById('fp-bits').value);
        }
        
        const response = await fetch(`/api/v1/fingerprints/${fpType}`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(params)
        });
        
        const data = await response.json();
        
        if (data.success) {
            displayFingerprint(data, smiles, fpType);
            document.getElementById('export-btn').disabled = false;
            window.lastFingerprintData = data; // Store for export and comparison
            window.Recent?.add({ smiles, page: '/fingerprints' });
        } else {
            showAlert(data.error || 'Generation failed', 'danger');
        }
        
    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
    }
}

function displayFingerprint(data, smiles, fpType) {
    const displayDiv = document.getElementById('fingerprint-display');
    const statsDiv = document.getElementById('fingerprint-stats');
    
    // Create fingerprint visualization
    const fingerprint = data.fingerprint;
    const onBits = data.on_bits || fingerprint.filter(bit => bit === 1).length;
    const totalBits = fingerprint.length;
    
    // Create bit pattern visualization
    const bitSize = Math.max(2, Math.floor(400 / Math.sqrt(totalBits)));
    const bitsPerRow = Math.floor(400 / bitSize);
    
    let visualization = `
        <div class="mb-4">
            <h4 class="text-sm font-medium text-gray-900 mb-2 break-words">
                <span class="block sm:inline">${fpType.replace('_', ' ').toUpperCase()} Fingerprint for:</span>
                <span class="font-mono text-primary-600 break-all">${smiles}</span>
            </h4>
            <div class="bg-white border rounded-lg p-4 overflow-hidden">
                <div id="fp-bit-grid" class="grid gap-1" data-cols="${bitsPerRow}" data-cell-px="${bitSize}">
    `;

    fingerprint.forEach((bit, index) => {
        visualization += `<div class="w-${bitSize === 2 ? '2' : '3'} h-${bitSize === 2 ? '2' : '3'} ${bit ? 'bg-primary-600' : 'bg-gray-200'} rounded-sm" title="Bit ${index}: ${bit}"></div>`;
    });

    visualization += `
                </div>
                <div class="mt-3 text-center">
                    <span class="text-xs text-gray-600">
                        Dark squares represent set bits (${onBits}/${totalBits} bits set, ${(onBits/totalBits*100).toFixed(1)}%)
                    </span>
                </div>
            </div>
        </div>
    `;

    displayDiv.innerHTML = visualization;
    // Phase-7 CSP: grid-template-columns is set after innerHTML so the
    // value isn't an inline `style=` attribute (blocked by style-src).
    // DOM property writes are still allowed.
    const bitGrid = document.getElementById("fp-bit-grid");
    if (bitGrid) {
        bitGrid.style.gridTemplateColumns = `repeat(${bitsPerRow}, ${bitSize}px)`;
        bitGrid.style.justifyContent = "center";
    }
    
    // Show statistics
    const stats = `
        <div class="grid grid-cols-2 sm:grid-cols-4 gap-4">
            <div class="bg-gray-50 rounded-lg p-3 text-center">
                <div class="text-2xl font-bold text-primary-600">${totalBits}</div>
                <div class="text-xs text-gray-600">Total Bits</div>
            </div>
            <div class="bg-gray-50 rounded-lg p-3 text-center">
                <div class="text-2xl font-bold text-green-600">${onBits}</div>
                <div class="text-xs text-gray-600">Set Bits</div>
            </div>
            <div class="bg-gray-50 rounded-lg p-3 text-center">
                <div class="text-2xl font-bold text-blue-600">${(onBits/totalBits*100).toFixed(1)}%</div>
                <div class="text-xs text-gray-600">Density</div>
            </div>
            <div class="bg-gray-50 rounded-lg p-3 text-center">
                <div class="text-sm font-bold text-purple-600 break-words leading-tight">${fpType.replace('_', ' ').toUpperCase()}</div>
                <div class="text-xs text-gray-600 mt-1">Type</div>
            </div>
        </div>
    `;
    
    statsDiv.innerHTML = stats;
    statsDiv.classList.remove('hidden');
}

async function compareFingerprints() {
    const smiles1 = document.getElementById('fingerprint-input').value.trim();
    const smiles2 = document.getElementById('compare-input').value.trim();
    const fpType = document.getElementById('fingerprint-type').value;
    
    if (!smiles1 || !smiles2) {
        showAlert('Please enter both SMILES strings', 'warning');
        return;
    }
    
    try {
        const response = await fetch('/api/v1/similarity/tanimoto', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                smiles1: smiles1,
                smiles2: smiles2,
                fingerprint_type: fpType
            })
        });
        
        const data = await response.json();
        
        if (data.success) {
            displaySimilarityResults(data, smiles1, smiles2);
        } else {
            showAlert(data.error || 'Comparison failed', 'danger');
        }
        
    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
    }
}

function displaySimilarityResults(data, smiles1, smiles2) {
    const resultsDiv = document.getElementById('similarity-results');
    const similarity = data.tanimoto_similarity;
    
    const getSimularityColor = (sim) => {
        if (sim >= 0.8) return 'text-green-600 bg-green-50';
        if (sim >= 0.6) return 'text-yellow-600 bg-yellow-50';
        if (sim >= 0.4) return 'text-orange-600 bg-orange-50';
        return 'text-red-600 bg-red-50';
    };
    
    const getSimilarityDescription = (sim) => {
        if (sim >= 0.8) return 'Very Similar';
        if (sim >= 0.6) return 'Similar';
        if (sim >= 0.4) return 'Somewhat Similar';
        return 'Different';
    };
    
    resultsDiv.innerHTML = `
        <div class="border border-gray-200 rounded-lg p-4">
            <h4 class="font-medium text-gray-900 mb-3">Tanimoto Similarity</h4>
            <div class="flex items-center justify-between mb-3">
                <div class="flex-1">
                    <div class="text-xs text-gray-600 mb-1">Molecule 1</div>
                    <div class="font-mono text-sm truncate">${smiles1}</div>
                </div>
                <div class="px-4 py-2 mx-4 rounded-lg ${getSimularityColor(similarity)}">
                    <div class="text-2xl font-bold">${(similarity * 100).toFixed(1)}%</div>
                </div>
                <div class="flex-1 text-right">
                    <div class="text-xs text-gray-600 mb-1">Molecule 2</div>
                    <div class="font-mono text-sm truncate">${smiles2}</div>
                </div>
            </div>
            <div class="text-center">
                <span class="inline-flex items-center px-3 py-1 rounded-full text-sm font-medium ${getSimularityColor(similarity)}">
                    ${getSimilarityDescription(similarity)}
                </span>
            </div>
        </div>
    `;
    
    resultsDiv.classList.remove('hidden');
}

function exportFingerprint() {
    if (!window.lastFingerprintData) return;
    
    const data = window.lastFingerprintData;
    
    // Create text content with fingerprint data
    let content = `Molecule: ${data.smiles}\\n`;
    content += `Fingerprint Type: ${document.getElementById('fingerprint-type').value}\\n`;
    content += `Total Bits: ${data.fingerprint.length}\\n`;
    content += `Set Bits: ${data.on_bits || data.fingerprint.filter(bit => bit === 1).length}\\n`;
    content += `Fingerprint Vector:\\n`;
    content += data.fingerprint.join(',');
    
    // Download file
    const blob = new Blob([content], { type: 'text/plain' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `fingerprint_${data.smiles.replace(/[^a-zA-Z0-9]/g, '_')}.txt`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
}

async function calculateBatchFingerprints() {
    const fileInput = document.getElementById('batch-fingerprints-file');
    const fpType = document.getElementById('fingerprint-type').value;
    const radius = document.getElementById('fp-radius').value;
    const nBits = document.getElementById('fp-bits').value;
    const exportFormat = document.getElementById('fp-export-format').value;
    const file = fileInput.files[0];

    if (!file) {
        showAlert('Please select a file', 'warning');
        return;
    }

    try {
        const resultsDiv = document.getElementById('fingerprint-results');
        resultsDiv.innerHTML = `
            <div class="flex flex-col items-center justify-center h-64 space-y-3">
                <div class="animate-spin rounded-full h-12 w-12 border-b-2 border-primary-600"></div>
                <div class="text-center">
                    <p class="text-gray-700 font-medium">Processing ${file.name}...</p>
                    <p class="text-sm text-gray-500">Generating fingerprints</p>
                    <p class="text-xs text-gray-400 mt-2">This may take a moment for large files</p>
                </div>
            </div>
        `;

        const formData = new FormData();
        formData.append('file', file);
        formData.append('fp_type', fpType);
        formData.append('radius', radius);
        formData.append('n_bits', nBits);
        formData.append('output_format', exportFormat);

        const response = await fetch('/api/v1/fingerprints/batch_file', {
            method: 'POST',
            body: formData
        });

        if (exportFormat === 'csv' || exportFormat === 'numpy') {
            // Download file
            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            const ext = exportFormat === 'numpy' ? 'npy' : 'csv';
            a.download = `fingerprints_${file.name.replace(/\\.[^.]+$/, '')}.${ext}`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);

            resultsDiv.innerHTML = `
                <div class="text-center py-12">
                    <div class="inline-flex items-center justify-center w-16 h-16 bg-green-100 rounded-full mb-4">
                        <i data-lucide="download" class="w-8 h-8 text-green-600"></i>
                    </div>
                    <h3 class="text-lg font-semibold text-gray-900 mb-2">${exportFormat.toUpperCase()} Downloaded!</h3>
                    <p class="text-gray-600">Check your downloads folder for the results</p>
                    <p class="text-sm text-gray-500 mt-2">File: fingerprints_${file.name.replace(/\\.[^.]+$/, '')}.${ext}</p>
                </div>
            `;
            lucide.createIcons();
        } else {
            // JSON display
            const data = await response.json();
            if (data.success) {
                displayBatchFingerprintResults(data);
            } else {
                showAlert(data.error || 'Batch generation failed', 'danger');
                resultsDiv.innerHTML = '';
            }
        }

    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
    }
}

function displayBatchFingerprintResults(data) {
    const results = data.results || [];
    const successful = data.num_successful || 0;
    const total = data.num_molecules || 0;

    let html = `
        <div class="mb-6 bg-gradient-to-r from-indigo-50 to-indigo-100 p-4 rounded-lg border border-indigo-200">
            <h3 class="text-lg font-semibold text-gray-900 mb-2">Batch Results</h3>
            <div class="grid grid-cols-3 gap-4 text-sm">
                <div>
                    <span class="text-gray-600">File:</span>
                    <span class="ml-2 font-medium">${data.file_name}</span>
                </div>
                <div>
                    <span class="text-gray-600">Type:</span>
                    <span class="ml-2 font-medium">${data.fingerprint_type}</span>
                </div>
                <div>
                    <span class="text-gray-600">Success:</span>
                    <span class="ml-2 font-bold text-green-700">${successful} / ${total}</span>
                </div>
            </div>
        </div>

        <div id="fingerprint-batch-paginator"></div>
    `;

    const container = document.getElementById('fingerprint-results');
    container.innerHTML = html;

    const renderRow = (result) => `
        <div class="p-4 border border-gray-200 rounded-lg hover:border-indigo-300 hover:shadow-sm transition-all">
            <div class="flex justify-between items-start mb-2">
                <div class="flex-1">
                    <h4 class="font-medium text-gray-900">${result.name}</h4>
                    <code class="text-xs text-gray-600">${result.smiles}</code>
                </div>
                <div class="ml-4">
                    ${result.error ? `
                        <span class="text-xs text-red-600">Error</span>
                    ` : `
                        <span class="text-xs text-green-600">${result.on_bits} bits set</span>
                    `}
                </div>
            </div>
            ${result.error ? `
                <div class="text-xs text-red-600 mt-1">Error: ${result.error}</div>
            ` : ''}
        </div>
    `;

    if (window.NuGenPagination) {
        window.NuGenPagination.render({
            container: container.querySelector('#fingerprint-batch-paginator'),
            items: results,
            renderItem: renderRow,
            pageSize: 20,
        });
    } else {
        container.querySelector('#fingerprint-batch-paginator').innerHTML =
            `<div class="space-y-3">${results.slice(0, 20).map(renderRow).join('')}</div>`;
    }
    if (window.lucide) lucide.createIcons();
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
})();
