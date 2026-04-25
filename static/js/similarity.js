/**
 * similarity.js — page module for "/similarity".
 *
 * Phase-5 mechanical extraction of the inline <script> that previously
 * lived in templates/similarity.html. Logic preserved verbatim except that
 * the local `showAlert` function delegates to NuGenUtils.showAlert so
 * dynamic toast HTML no longer contains inline `onclick=` (CSP-safe).
 *
 * Phase 7 polish may refactor `fetch` calls to use NuGenUtils.fetchJSON.
 */
(function () {
  "use strict";
  if (window.location.pathname !== "/similarity") return;

document.addEventListener('DOMContentLoaded', function() {
    // Calculate similarity button handler
    document.getElementById('calculate-similarity-btn').addEventListener('click', function() {
        calculateSimilarity();
    });
    
    // Search substructure button handler
    document.getElementById('search-substructure-btn').addEventListener('click', function() {
        searchSubstructure();
    });

    // Select diverse subset button handler
    document.getElementById('select-diverse-subset-btn').addEventListener('click', function() {
        selectDiverseSubset();
    });

    // Example pair handlers
    document.querySelectorAll('.example-pair').forEach(button => {
        button.addEventListener('click', function() {
            const query = this.getAttribute('data-query');
            const target = this.getAttribute('data-target');
            document.getElementById('query-molecule').value = query;
            document.getElementById('target-molecule').value = target;
            calculateSimilarity();
        });
    });
    
    // Example substructure handlers
    document.querySelectorAll('.example-substructure').forEach(button => {
        button.addEventListener('click', function() {
            const pattern = this.getAttribute('data-pattern');
            document.getElementById('substructure-pattern').value = pattern;
            searchSubstructure();
        });
    });
    
    // Auto-calculate on page load
    calculateSimilarity();
});

async function calculateSimilarity() {
    const querySmiles = document.getElementById('query-molecule').value.trim();
    const targetSmiles = document.getElementById('target-molecule').value.trim();
    const fpType = document.getElementById('similarity-fp-type').value;
    const metric = document.getElementById('similarity-metric').value;
    const fileInput = document.getElementById('bulk-similarity-file');

    // Check if file upload is being used
    if (fileInput && fileInput.files.length > 0) {
        await calculateBulkSimilarity();
        return;
    }

    if (!querySmiles || !targetSmiles) {
        showAlert('Please enter both query and target molecules', 'warning');
        return;
    }
    
    try {
        // Show loading state
        document.getElementById('similarity-results').innerHTML = `
            <div class="flex items-center justify-center h-32">
                <div class="animate-spin rounded-full h-8 w-8 border-b-2 border-primary-600"></div>
                <span class="ml-3 text-gray-600">Calculating similarity...</span>
            </div>
        `;
        
        // Calculate similarity
        const response = await fetch(`/api/v1/similarity/${metric}`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                smiles1: querySmiles,
                smiles2: targetSmiles,
                fingerprint_type: fpType
            })
        });
        
        const data = await response.json();
        
        if (data.success) {
            displaySimilarityResults(data, querySmiles, targetSmiles, metric);
            loadMolecularStructures(querySmiles, targetSmiles);
        } else {
            showAlert(data.error || 'Similarity calculation failed', 'danger');
        }
        
    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
    }
}

function displaySimilarityResults(data, query, target, metric) {
    const resultsDiv = document.getElementById('similarity-results');
    const similarity = data[`${metric}_similarity`] || data.similarity || 0;
    
    const getColor = (sim) => {
        if (sim >= 0.8) return { bg: 'bg-green-100', text: 'text-green-800', bar: 'bg-green-500' };
        if (sim >= 0.6) return { bg: 'bg-yellow-100', text: 'text-yellow-800', bar: 'bg-yellow-500' };
        if (sim >= 0.4) return { bg: 'bg-orange-100', text: 'text-orange-800', bar: 'bg-orange-500' };
        return { bg: 'bg-red-100', text: 'text-red-800', bar: 'bg-red-500' };
    };
    
    const getDescription = (sim) => {
        if (sim >= 0.8) return 'Very Similar';
        if (sim >= 0.6) return 'Similar';
        if (sim >= 0.4) return 'Somewhat Similar';
        return 'Different';
    };
    
    const colors = getColor(similarity);
    
    resultsDiv.innerHTML = `
        <div class="space-y-6">
            <!-- Main Similarity Score -->
            <div class="text-center">
                <div class="inline-flex items-center justify-center w-32 h-32 ${colors.bg} rounded-full mb-4">
                    <div class="text-center">
                        <div class="text-3xl font-bold ${colors.text}">${(similarity * 100).toFixed(1)}%</div>
                        <div class="text-sm ${colors.text} font-medium">${metric.toUpperCase()}</div>
                    </div>
                </div>
                <h3 class="text-lg font-semibold text-gray-900">${getDescription(similarity)}</h3>
                <p class="text-gray-600">Based on ${data.fingerprint_type.toUpperCase()} fingerprints</p>
            </div>
            
            <!-- Similarity Bar -->
            <div class="space-y-2">
                <div class="flex justify-between text-sm text-gray-600">
                    <span>Similarity Score</span>
                    <span>${similarity.toFixed(4)}</span>
                </div>
                <div class="w-full bg-gray-200 rounded-full h-3">
                    <div id="similarity-bar-fill"
                         class="${colors.bar} h-3 rounded-full transition-all duration-1000 ease-out"
                         data-fill-pct="${similarity * 100}"></div>
                </div>
                <div class="flex justify-between text-xs text-gray-500">
                    <span>0% (Different)</span>
                    <span>100% (Identical)</span>
                </div>
            </div>
            
            <!-- Molecule Comparison -->
            <div class="bg-gray-50 rounded-lg p-4">
                <div class="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div class="text-center">
                        <h4 class="font-medium text-gray-900 mb-2">Query Molecule</h4>
                        <code class="text-xs bg-white px-2 py-1 rounded border break-all">${query}</code>
                    </div>
                    <div class="text-center">
                        <h4 class="font-medium text-gray-900 mb-2">Target Molecule</h4>
                        <code class="text-xs bg-white px-2 py-1 rounded border break-all">${target}</code>
                    </div>
                </div>
            </div>
            
            <!-- Additional Metrics -->
            ${data.additional_metrics ? `
                <div class="grid grid-cols-2 md:grid-cols-4 gap-4">
                    ${Object.entries(data.additional_metrics).map(([key, value]) => `
                        <div class="bg-gray-50 rounded-lg p-3 text-center">
                            <div class="text-lg font-semibold text-gray-900">${typeof value === 'number' ? value.toFixed(3) : value}</div>
                            <div class="text-xs text-gray-600">${key.replace(/_/g, ' ').toUpperCase()}</div>
                        </div>
                    `).join('')}
                </div>
            ` : ''}
        </div>
    `;
}

async function loadMolecularStructures(query, target) {
    const queryDiv = document.getElementById('query-structure');
    const targetDiv = document.getElementById('target-structure');
    
    // Load query structure
    try {
        const queryResponse = await fetch('/api/v1/visualization/draw_svg', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles: query, width: 250, height: 200 })
        });
        
        const queryData = await queryResponse.json();
        if (queryData.success) {
            queryDiv.innerHTML = queryData.svg;
        }
    } catch (error) {
        console.warn('Failed to load query structure:', error);
    }
    
    // Load target structure
    try {
        const targetResponse = await fetch('/api/v1/visualization/draw_svg', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles: target, width: 250, height: 200 })
        });
        
        const targetData = await targetResponse.json();
        if (targetData.success) {
            targetDiv.innerHTML = targetData.svg;
        }
    } catch (error) {
        console.warn('Failed to load target structure:', error);
    }
}

async function searchSubstructure() {
    const pattern = document.getElementById('substructure-pattern').value.trim();
    const moleculesText = document.getElementById('target-molecules').value.trim();
    
    if (!pattern || !moleculesText) {
        showAlert('Please enter both pattern and target molecules', 'warning');
        return;
    }
    
    const molecules = moleculesText.split('\\n').map(s => s.trim()).filter(s => s);
    
    try {
        document.getElementById('substructure-results').innerHTML = `
            <div class="flex items-center justify-center h-32">
                <div class="animate-spin rounded-full h-8 w-8 border-b-2 border-primary-600"></div>
                <span class="ml-3 text-gray-600">Searching substructures...</span>
            </div>
        `;
        
        const response = await fetch('/api/v1/similarity/substructure', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                pattern: pattern,
                molecules: molecules
            })
        });
        
        const data = await response.json();
        
        if (data.success) {
            displaySubstructureResults(data, pattern);
        } else {
            showAlert(data.error || 'Substructure search failed', 'danger');
        }
        
    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
    }
}

function displaySubstructureResults(data, pattern) {
    const resultsDiv = document.getElementById('substructure-results');
    const matches = data.matches || [];
    const total = data.num_targets || 0;
    
    let html = `
        <div class="mb-6">
            <div class="flex items-center justify-between">
                <h3 class="text-lg font-semibold text-gray-900">
                    Pattern: <code class="text-sm bg-gray-100 px-2 py-1 rounded">${pattern}</code>
                </h3>
                <span class="inline-flex items-center px-3 py-1 rounded-full text-sm font-medium ${
                    matches.length > 0 ? 'bg-green-100 text-green-800' : 'bg-gray-100 text-gray-800'
                }">
                    ${matches.length}/${total} matches
                </span>
            </div>
        </div>
    `;
    
    if (matches.length === 0) {
        html += `
            <div class="text-center py-8">
                <i data-lucide="x-circle" class="w-12 h-12 mx-auto mb-4 text-gray-400"></i>
                <p class="text-gray-600">No molecules contain the specified substructure</p>
            </div>
        `;
    } else {
        html += `
            <div class="space-y-4">
                ${matches.map((match, index) => `
                    <div class="border border-gray-200 rounded-lg p-4 hover:bg-gray-50 transition-colors">
                        <div class="flex items-center justify-between">
                            <div class="flex-1">
                                <h4 class="font-medium text-gray-900">Match ${index + 1}</h4>
                                <code class="text-sm text-gray-600 break-all">${match.target_smiles}</code>
                                ${match.name ? `<p class="text-sm text-gray-500 mt-1">${match.name}</p>` : ''}
                            </div>
                            <div class="ml-4 flex-shrink-0">
                                <span class="inline-flex items-center px-2 py-1 rounded-full text-xs font-medium bg-green-100 text-green-800">
                                    <i data-lucide="check" class="w-3 h-3 mr-1"></i>
                                    Match
                                </span>
                            </div>
                        </div>
                        ${match.match_atoms && match.match_atoms.length > 0 ? `
                            <div class="mt-3 text-xs text-gray-500">
                                Matching atoms: ${match.match_atoms.join(', ')}
                            </div>
                        ` : ''}
                    </div>
                `).join('')}
            </div>
        `;
    }
    
    resultsDiv.innerHTML = html;

    // Phase-7 CSP: width was an inline `style=` attribute; set it via
    // DOM property after innerHTML so style-src no longer needs
    // 'unsafe-inline'. The value comes from the API response (a number),
    // so coerce + clamp to defend against unexpected payload shapes.
    const fill = document.getElementById("similarity-bar-fill");
    if (fill) {
        const pct = Math.max(0, Math.min(100, parseFloat(fill.dataset.fillPct) || 0));
        fill.style.width = `${pct}%`;
    }

    // Re-initialize icons
    if (window.lucide) {
        lucide.createIcons();
    }
}

async function calculateBulkSimilarity() {
    const querySmiles = document.getElementById('query-molecule').value.trim();
    const fileInput = document.getElementById('bulk-similarity-file');
    const fpType = document.getElementById('similarity-fp-type').value;
    const metric = document.getElementById('similarity-metric').value;

    if (!querySmiles) {
        showAlert('Please enter a query molecule', 'warning');
        return;
    }

    if (!fileInput.files.length) {
        showAlert('Please select a file to upload', 'warning');
        return;
    }

    const file = fileInput.files[0];

    try {
        document.getElementById('similarity-results').innerHTML = `
            <div class="flex flex-col items-center justify-center h-32 space-y-2">
                <div class="animate-spin rounded-full h-8 w-8 border-b-2 border-primary-600"></div>
                <span class="text-gray-600">Uploading and analyzing ${file.name}...</span>
                <span class="text-sm text-gray-500">This may take a moment for large files</span>
            </div>
        `;

        const formData = new FormData();
        formData.append('file', file);
        formData.append('query_smiles', querySmiles);
        formData.append('fp_type', fpType);
        formData.append('similarity_metric', metric);
        formData.append('threshold', '0.0');

        const response = await fetch('/api/v1/similarity/bulk_similarity_file', {
            method: 'POST',
            body: formData
        });

        const data = await response.json();

        if (data.success) {
            displayBulkSimilarityResults(data, querySmiles);
        } else {
            showAlert(data.error || 'Bulk similarity calculation failed', 'danger');
        }

    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
    }
}

function displayBulkSimilarityResults(data, querySmiles) {
    const resultsDiv = document.getElementById('similarity-results');
    const results = data.results || [];
    const hits = data.num_hits || 0;
    const total = data.num_targets || 0;

    let html = `
        <div class="mb-6">
            <div class="bg-gradient-to-r from-green-50 to-green-100 p-4 rounded-lg border border-green-200">
                <h3 class="text-lg font-semibold text-gray-900 mb-2">Bulk Similarity Results</h3>
                <div class="grid grid-cols-2 gap-4 text-sm">
                    <div>
                        <span class="text-gray-600">Query:</span>
                        <code class="ml-2 bg-white px-2 py-1 rounded text-xs">${querySmiles}</code>
                    </div>
                    <div>
                        <span class="text-gray-600">File:</span>
                        <span class="ml-2 font-medium">${data.file_name}</span>
                    </div>
                    <div>
                        <span class="text-gray-600">Matches Found:</span>
                        <span class="ml-2 font-bold text-green-700">${hits} / ${total}</span>
                    </div>
                    <div>
                        <span class="text-gray-600">Fingerprint:</span>
                        <span class="ml-2 font-medium">${data.fingerprint_type}</span>
                    </div>
                </div>
            </div>
        </div>

        <div class="space-y-3 max-h-96 overflow-y-auto">
            ${results.slice(0, 50).map(result => `
                <div class="p-4 border border-gray-200 rounded-lg hover:border-primary-300 hover:shadow-sm transition-all">
                    <div class="flex justify-between items-start mb-2">
                        <div class="flex-1">
                            <h4 class="font-medium text-gray-900">${result.name}</h4>
                            <code class="text-xs text-gray-600 break-all">${result.target_smiles}</code>
                        </div>
                        <div class="ml-4">
                            ${result.similarity !== null ? `
                                <span class="inline-flex items-center px-3 py-1 rounded-full text-sm font-medium ${
                                    result.similarity >= 0.8 ? 'bg-green-100 text-green-800' :
                                    result.similarity >= 0.5 ? 'bg-yellow-100 text-yellow-800' :
                                    'bg-gray-100 text-gray-800'
                                }">
                                    ${(result.similarity * 100).toFixed(1)}%
                                </span>
                            ` : '<span class="text-red-500 text-sm">Error</span>'}
                        </div>
                    </div>
                    ${result.error ? `
                        <div class="text-xs text-red-600 mt-1">Error: ${result.error}</div>
                    ` : ''}
                </div>
            `).join('')}
            ${results.length > 50 ? `
                <div class="text-center text-sm text-gray-500 py-3">
                    Showing top 50 of ${results.length} results
                </div>
            ` : ''}
        </div>
    `;

    resultsDiv.innerHTML = html;

    if (window.lucide) {
        lucide.createIcons();
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

async function selectDiverseSubset() {
    const fileInput = document.getElementById('diverse-subset-file-upload');
    const subsetSize = parseInt(document.getElementById('diverse-subset-size').value);
    const fpType = document.getElementById('diverse-fp-type').value;

    if (!fileInput || fileInput.files.length === 0) {
        showAlert('Please upload a molecule library file', 'warning');
        return;
    }

    const file = fileInput.files[0];

    try {
        // Show loading state
        document.getElementById('diverse-subset-results').innerHTML = `
            <div class="flex items-center justify-center h-32">
                <div class="animate-spin rounded-full h-8 w-8 border-b-2 border-purple-600"></div>
                <span class="ml-3 text-gray-600">Selecting diverse subset...</span>
            </div>
        `;

        const formData = new FormData();
        formData.append('file', file);
        formData.append('num_molecules', subsetSize);
        formData.append('fp_type', fpType);

        const response = await fetch('/api/v1/similarity/diverse_subset_file', {
            method: 'POST',
            body: formData
        });

        const data = await response.json();

        if (data.success) {
            displayDiverseSubsetResults(data);
        } else {
            showAlert(data.error || 'Diverse subset selection failed', 'danger');
            document.getElementById('diverse-subset-results').innerHTML = `
                <div class="flex items-center justify-center h-32 text-red-500">
                    <div class="text-center">
                        <i data-lucide="alert-circle" class="w-12 h-12 mx-auto mb-2"></i>
                        <p class="text-sm">${data.error || 'Error occurred'}</p>
                    </div>
                </div>
            `;
        }

    } catch (error) {
        console.error('Error:', error);
        showAlert('Network error occurred', 'danger');
        document.getElementById('diverse-subset-results').innerHTML = `
            <div class="flex items-center justify-center h-32 text-red-500">
                <div class="text-center">
                    <i data-lucide="alert-circle" class="w-12 h-12 mx-auto mb-2"></i>
                    <p class="text-sm">Network error</p>
                </div>
            </div>
        `;
    }
}

function displayDiverseSubsetResults(data) {
    const resultsDiv = document.getElementById('diverse-subset-results');
    const selected = data.selected_molecules || [];
    const totalInput = data.total_input_molecules || 0;

    let html = `
        <div class="mb-6">
            <div class="bg-gradient-to-r from-purple-50 to-purple-100 p-4 rounded-lg border border-purple-200">
                <h3 class="text-lg font-semibold text-gray-900 mb-2">Diverse Subset Selection</h3>
                <div class="grid grid-cols-2 gap-4 text-sm">
                    <div>
                        <span class="text-gray-600">Input Molecules:</span>
                        <span class="ml-2 font-bold text-purple-700">${totalInput}</span>
                    </div>
                    <div>
                        <span class="text-gray-600">Selected:</span>
                        <span class="ml-2 font-bold text-purple-700">${selected.length}</span>
                    </div>
                    <div>
                        <span class="text-gray-600">Fingerprint:</span>
                        <span class="ml-2 font-medium">${data.fingerprint_type || 'morgan'}</span>
                    </div>
                    <div>
                        <span class="text-gray-600">Method:</span>
                        <span class="ml-2 font-medium">MaxMin Picker</span>
                    </div>
                </div>
            </div>
        </div>

        <div class="space-y-3 max-h-96 overflow-y-auto">
            ${selected.map((mol, idx) => `
                <div class="p-4 border border-gray-200 rounded-lg hover:border-purple-300 hover:shadow-sm transition-all">
                    <div class="flex justify-between items-start mb-2">
                        <div class="flex-1">
                            <div class="flex items-center space-x-2 mb-1">
                                <span class="inline-flex items-center justify-center w-6 h-6 rounded-full bg-purple-100 text-purple-700 text-xs font-bold">
                                    ${idx + 1}
                                </span>
                                <h4 class="font-medium text-gray-900">${mol.name}</h4>
                            </div>
                            <code class="text-xs text-gray-600 break-all">${mol.smiles}</code>
                        </div>
                        ${mol.diversity_score !== undefined ? `
                            <div class="ml-4">
                                <span class="inline-flex items-center px-3 py-1 rounded-full text-sm font-medium bg-purple-100 text-purple-800">
                                    Diversity: ${(mol.diversity_score * 100).toFixed(1)}%
                                </span>
                            </div>
                        ` : ''}
                    </div>
                    ${mol.error ? `
                        <div class="text-xs text-red-600 mt-1">Error: ${mol.error}</div>
                    ` : ''}
                </div>
            `).join('')}
        </div>

        ${selected.length > 0 ? `
            <div class="mt-4 p-3 bg-gray-50 rounded-lg border border-gray-200">
                <p class="text-sm text-gray-600">
                    <i data-lucide="info" class="w-4 h-4 inline mr-1"></i>
                    These ${selected.length} molecules were selected using the MaxMin picker algorithm to maximize structural diversity.
                </p>
            </div>
        ` : ''}
    `;

    resultsDiv.innerHTML = html;

    if (window.lucide) {
        lucide.createIcons();
    }
}
})();
