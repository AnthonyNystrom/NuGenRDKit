/**
 * descriptors.js — page module for "/descriptors".
 *
 * Phase-5 mechanical extraction of the inline <script> that previously
 * lived in templates/descriptors.html. Logic preserved verbatim except that
 * the local `showAlert` function delegates to NuGenUtils.showAlert so
 * dynamic toast HTML no longer contains inline `onclick=` (CSP-safe).
 *
 * Phase 7 polish may refactor `fetch` calls to use NuGenUtils.fetchJSON.
 */
(function () {
  "use strict";
  if (window.location.pathname !== "/descriptors") return;

document.addEventListener('DOMContentLoaded', function() {
    // Check for SMILES parameter from URL
    const urlParams = new URLSearchParams(window.location.search);
    const smilesParam = urlParams.get('smiles');
    if (smilesParam) {
        document.getElementById('descriptor-input').value = smilesParam;
        // Auto-calculate if SMILES was provided
        calculateDescriptors();
        return; // Skip the normal auto-calculate at the bottom
    }
    
    // Calculate button handler
    document.getElementById('calculate-btn').addEventListener('click', function() {
        calculateDescriptors();
    });
    
    // Export button handler
    document.getElementById('export-btn').addEventListener('click', function() {
        exportDescriptors();
    });
    
    // Example molecule handlers
    document.querySelectorAll('.example-molecule').forEach(button => {
        button.addEventListener('click', function() {
            const smiles = this.getAttribute('data-smiles');
            document.getElementById('descriptor-input').value = smiles;
            calculateDescriptors();
            // Scroll the result panel into view so the user sees
            // loading-then-result rather than staring at the button.
            document.getElementById('descriptor-results')
                ?.scrollIntoView({ behavior: 'smooth', block: 'start' });
        });
    });
    
    // Auto-calculate on page load
    calculateDescriptors();
});

async function calculateDescriptors() {
    const smiles = document.getElementById('descriptor-input').value.trim();
    const descriptorType = document.getElementById('descriptor-type').value;
    const fileInput = document.getElementById('batch-descriptors-file');

    // Check if file upload is being used
    if (fileInput && fileInput.files.length > 0) {
        await calculateBatchDescriptors();
        return;
    }

    if (!smiles) {
        showAlert('Please enter a SMILES string', 'warning');
        return;
    }
    
    try {
        document.getElementById('descriptor-loading').classList.remove('hidden');
        document.getElementById('descriptor-results').innerHTML = '';
        document.getElementById('export-btn').disabled = true;
        
        const endpoint = `/api/v1/descriptors/${descriptorType}`;
        const response = await fetch(endpoint, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles })
        });
        
        const data = await response.json();
        
        if (data.success) {
            displayDescriptors(data.descriptors, smiles, descriptorType);
            document.getElementById('export-btn').disabled = false;
            window.lastDescriptorData = data; // Store for export
            window.Recent?.add({ smiles, page: '/descriptors' });
        } else {
            showAlert(data.error || 'Calculation failed', 'danger');
            document.getElementById('descriptor-results').innerHTML = `
                <div class="flex items-center justify-center h-96 text-red-500">
                    <div class="text-center">
                        <i data-lucide="alert-circle" class="w-16 h-16 mx-auto mb-4"></i>
                        <p class="text-lg font-medium">Calculation Failed</p>
                        <p class="text-sm">${data.error || 'Unknown error occurred'}</p>
                    </div>
                </div>
            `;
        }
        
    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
    } finally {
        document.getElementById('descriptor-loading').classList.add('hidden');
        lucide.createIcons();
    }
}

function displayDescriptors(descriptors, smiles, type) {
    const resultsDiv = document.getElementById('descriptor-results');
    
    // Create molecule header
    let html = `
        <div class="mb-6 p-4 bg-gray-50 rounded-lg">
            <div class="flex items-center justify-between">
                <div>
                    <h3 class="font-medium text-gray-900">Molecule: <span class="font-mono text-primary-600">${smiles}</span></h3>
                    <p class="text-sm text-gray-600">${getDescriptorTypeDescription(type)}</p>
                </div>
                <div class="text-right">
                    <div class="text-2xl font-bold text-primary-600">${Object.keys(descriptors).length}</div>
                    <div class="text-xs text-gray-600">Descriptors</div>
                </div>
            </div>
        </div>
    `;
    
    // Special formatting for different descriptor types
    if (type === 'lipinski') {
        html += formatLipinskiResults(descriptors);
    } else if (type === 'basic') {
        html += formatBasicResults(descriptors);
    } else {
        html += formatGeneralResults(descriptors);
    }
    
    resultsDiv.innerHTML = html;
}

function formatLipinskiResults(descriptors) {
    const violations = descriptors.lipinski_violations || 0;
    const passes = descriptors.passes_lipinski || violations <= 1;
    
    return `
        <div class="mb-6 p-4 rounded-lg ${passes ? 'bg-green-50 border border-green-200' : 'bg-red-50 border border-red-200'}">
            <div class="flex items-center space-x-2 mb-3">
                <i data-lucide="${passes ? 'check-circle' : 'x-circle'}" class="w-5 h-5 ${passes ? 'text-green-600' : 'text-red-600'}"></i>
                <h3 class="font-medium ${passes ? 'text-green-800' : 'text-red-800'}">
                    Lipinski Rule of 5: ${passes ? 'PASSED' : 'FAILED'} (${violations} violation${violations !== 1 ? 's' : ''})
                </h3>
            </div>
            <p class="text-sm ${passes ? 'text-green-700' : 'text-red-700'} mb-4">
                ${passes ? 'This compound shows good drug-like properties.' : 'This compound may have poor oral bioavailability.'}
            </p>
        </div>
        
        <div class="overflow-x-auto">
            <table class="min-w-full divide-y divide-gray-200">
                <thead class="bg-gray-50">
                    <tr>
                        <th class="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Property</th>
                        <th class="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Value</th>
                        <th class="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Limit</th>
                        <th class="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Status</th>
                    </tr>
                </thead>
                <tbody class="bg-white divide-y divide-gray-200">
                    ${formatLipinskiRow('Molecular Weight', descriptors.molecular_weight, '≤ 500 Da', descriptors.molecular_weight <= 500)}
                    ${formatLipinskiRow('LogP', descriptors.logp, '≤ 5', descriptors.logp <= 5)}
                    ${formatLipinskiRow('H-bond Donors', descriptors.num_h_donors, '≤ 5', descriptors.num_h_donors <= 5)}
                    ${formatLipinskiRow('H-bond Acceptors', descriptors.num_h_acceptors, '≤ 10', descriptors.num_h_acceptors <= 10)}
                </tbody>
            </table>
        </div>
    `;
}

function formatLipinskiRow(property, value, limit, passes) {
    return `
        <tr>
            <td class="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900">${property}</td>
            <td class="px-6 py-4 whitespace-nowrap text-sm text-gray-900">${typeof value === 'number' ? value.toFixed(2) : value}</td>
            <td class="px-6 py-4 whitespace-nowrap text-sm text-gray-500">${limit}</td>
            <td class="px-6 py-4 whitespace-nowrap">
                <span class="inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium ${
                    passes ? 'bg-green-100 text-green-800' : 'bg-red-100 text-red-800'
                }">
                    <i data-lucide="${passes ? 'check' : 'x'}" class="w-3 h-3 mr-1"></i>
                    ${passes ? 'Pass' : 'Fail'}
                </span>
            </td>
        </tr>
    `;
}

function formatBasicResults(descriptors) {
    return `
        <div class="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 mb-6">
            ${Object.entries(descriptors).map(([key, value]) => `
                <div class="bg-gray-50 rounded-lg p-4">
                    <div class="text-2xl font-bold text-primary-600">${typeof value === 'number' ? value.toFixed(value < 1 ? 3 : 1) : (value || 'N/A')}</div>
                    <div class="text-sm font-medium text-gray-900">${formatDescriptorName(key)}</div>
                </div>
            `).join('')}
        </div>
    `;
}

function formatGeneralResults(descriptors) {
    return `
        <div class="overflow-x-auto">
            <table class="min-w-full divide-y divide-gray-200">
                <thead class="bg-gray-50">
                    <tr>
                        <th class="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Descriptor</th>
                        <th class="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Value</th>
                    </tr>
                </thead>
                <tbody class="bg-white divide-y divide-gray-200">
                    ${Object.entries(descriptors).map(([key, value], index) => `
                        <tr class="${index % 2 === 0 ? 'bg-white' : 'bg-gray-50'}">
                            <td class="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900">${formatDescriptorName(key)}</td>
                            <td class="px-6 py-4 whitespace-nowrap text-sm text-gray-900 font-mono">${
                                typeof value === 'number' ? value.toFixed(4) : (value || 'N/A')
                            }</td>
                        </tr>
                    `).join('')}
                </tbody>
            </table>
        </div>
    `;
}

function formatDescriptorName(name) {
    return name.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
}

function getDescriptorTypeDescription(type) {
    const descriptions = {
        'basic': 'Essential molecular properties and counts',
        'lipinski': 'Drug-likeness assessment using Lipinski Rule of 5',
        'logp': 'Lipophilicity and partition coefficient descriptors',
        'topological': 'Graph-based molecular topology indices',
        'all': 'Complete set of available molecular descriptors'
    };
    return descriptions[type] || 'Molecular descriptors';
}

function exportDescriptors() {
    if (!window.lastDescriptorData) return;
    
    const data = window.lastDescriptorData;
    const descriptors = data.descriptors;
    
    // Create CSV content
    let csv = 'Descriptor,Value\\n';
    Object.entries(descriptors).forEach(([key, value]) => {
        csv += `"${formatDescriptorName(key)}","${value}"\\n`;
    });
    
    // Download CSV
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `descriptors_${data.smiles.replace(/[^a-zA-Z0-9]/g, '_')}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
}

async function calculateBatchDescriptors() {
    const fileInput = document.getElementById('batch-descriptors-file');
    const descriptorType = document.getElementById('descriptor-type').value;
    const exportFormat = document.getElementById('export-format').value;
    const file = fileInput.files[0];

    if (!file) {
        showAlert('Please select a file', 'warning');
        return;
    }

    try {
        document.getElementById('descriptor-loading').classList.remove('hidden');
        document.getElementById('descriptor-results').innerHTML = `
            <div class="flex flex-col items-center justify-center h-64 space-y-3">
                <div class="animate-spin rounded-full h-12 w-12 border-b-2 border-primary-600"></div>
                <div class="text-center">
                    <p class="text-gray-700 font-medium">Processing ${file.name}...</p>
                    <p class="text-sm text-gray-500">Calculating descriptors for molecules</p>
                    <p class="text-xs text-gray-400 mt-2">This may take a moment for large files</p>
                </div>
            </div>
        `;

        const formData = new FormData();
        formData.append('file', file);
        formData.append('descriptor_type', descriptorType);
        formData.append('output_format', exportFormat);

        const response = await fetch('/api/v1/descriptors/batch_file', {
            method: 'POST',
            body: formData
        });

        if (exportFormat === 'csv') {
            // CSV download
            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `descriptors_${file.name.replace(/\\.[^.]+$/, '')}.csv`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);

            document.getElementById('descriptor-results').innerHTML = `
                <div class="text-center py-12">
                    <div class="inline-flex items-center justify-center w-16 h-16 bg-green-100 rounded-full mb-4">
                        <i data-lucide="download" class="w-8 h-8 text-green-600"></i>
                    </div>
                    <h3 class="text-lg font-semibold text-gray-900 mb-2">CSV Downloaded!</h3>
                    <p class="text-gray-600">Check your downloads folder for the results</p>
                    <p class="text-sm text-gray-500 mt-2">File: descriptors_${file.name.replace(/\\.[^.]+$/, '')}.csv</p>
                </div>
            `;
            lucide.createIcons();
        } else {
            // JSON display
            const data = await response.json();

            if (data.success) {
                displayBatchResults(data);
            } else {
                showAlert(data.error || 'Batch calculation failed', 'danger');
                document.getElementById('descriptor-results').innerHTML = '';
            }
        }

    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
        document.getElementById('descriptor-results').innerHTML = '';
    } finally {
        document.getElementById('descriptor-loading').classList.add('hidden');
    }
}

function displayBatchResults(data) {
    const results = data.results || [];
    const successful = data.num_successful || 0;
    const total = data.num_molecules || 0;

    const container = document.getElementById('descriptor-results');
    container.innerHTML = `
        <div class="mb-6 bg-gradient-to-r from-purple-50 to-purple-100 p-4 rounded-lg border border-purple-200">
            <h3 class="text-lg font-semibold text-gray-900 mb-2">Batch Results</h3>
            <div class="grid grid-cols-3 gap-4 text-sm">
                <div>
                    <span class="text-gray-600">File:</span>
                    <span class="ml-2 font-medium">${data.file_name}</span>
                </div>
                <div>
                    <span class="text-gray-600">Type:</span>
                    <span class="ml-2 font-medium">${data.descriptor_type}</span>
                </div>
                <div>
                    <span class="text-gray-600">Success:</span>
                    <span class="ml-2 font-bold text-green-700">${successful} / ${total}</span>
                </div>
            </div>
        </div>
        <div id="descriptor-batch-paginator"></div>
    `;

    const renderRow = (result) => `
        <div class="p-4 border border-gray-200 rounded-lg hover:border-purple-300 hover:shadow-sm transition-all">
            <div class="flex justify-between items-start mb-2">
                <div class="flex-1">
                    <h4 class="font-medium text-gray-900">${result.name}</h4>
                    <code class="text-xs text-gray-600">${result.smiles}</code>
                </div>
                ${result.error ? `
                    <span class="ml-4 text-xs text-red-600">Error: ${result.error}</span>
                ` : `
                    <span class="ml-4 text-xs text-green-600">✓ Success</span>
                `}
            </div>
            ${result.descriptors ? `
                <div class="mt-2 grid grid-cols-2 gap-2 text-xs">
                    ${Object.entries(result.descriptors).slice(0, 6).map(([key, val]) => `
                        <div>
                            <span class="text-gray-600">${key}:</span>
                            <span class="font-medium">${typeof val === 'number' ? val.toFixed(2) : val}</span>
                        </div>
                    `).join('')}
                </div>
            ` : ''}
        </div>
    `;

    if (window.NuGenPagination) {
        window.NuGenPagination.render({
            container: container.querySelector('#descriptor-batch-paginator'),
            items: results,
            renderItem: renderRow,
            pageSize: 20,
        });
    } else {
        container.querySelector('#descriptor-batch-paginator').innerHTML =
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
