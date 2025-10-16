// Modern NuGenRDKit Web Application JavaScript

class NuGenRDKitApp {
    constructor() {
        this.baseURL = '';
        this.init();
    }

    init() {
        this.checkServerHealth();
        this.setupGlobalEventListeners();
        this.initializeCurrentPage();
    }

    // Global setup methods
    setupGlobalEventListeners() {
        // Smooth scrolling for anchor links
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                const target = document.querySelector(this.getAttribute('href'));
                if (target) {
                    target.scrollIntoView({ behavior: 'smooth' });
                }
            });
        });

        // Mobile menu toggle
        const mobileMenuButton = document.getElementById('mobile-menu-button');
        const mobileMenu = document.getElementById('mobile-menu');
        
        if (mobileMenuButton && mobileMenu) {
            mobileMenuButton.addEventListener('click', () => {
                mobileMenu.classList.toggle('hidden');
            });
        }

        // Close mobile menu when clicking outside
        document.addEventListener('click', (e) => {
            if (mobileMenu && !mobileMenu.contains(e.target) && !mobileMenuButton.contains(e.target)) {
                mobileMenu.classList.add('hidden');
            }
        });

        // Keyboard navigation
        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                this.closeAllModals();
                if (mobileMenu) mobileMenu.classList.add('hidden');
            }
        });
    }

    initializeCurrentPage() {
        const path = window.location.pathname;
        
        switch (path) {
            case '/':
                this.initHomePage();
                break;
            case '/structure':
                this.initStructurePage();
                break;
            case '/descriptors':
                this.initDescriptorsPage();
                break;
            case '/fingerprints':
                this.initFingerprintsPage();
                break;
            default:
                // Initialize common functionality for other pages
                break;
        }
    }

    // Page-specific initialization
    initHomePage() {
        this.loadHeroMolecule();
        this.animateStats();
        this.setupFeatureCards();
    }

    initStructurePage() {
        // Structure page is handled by inline scripts
    }

    initDescriptorsPage() {
        // Descriptors page is handled by inline scripts
    }

    initFingerprintsPage() {
        // Fingerprints page is handled by inline scripts
    }

    // Utility Methods
    async checkServerHealth() {
        try {
            const response = await fetch('/health');
            const data = await response.json();
            
            const indicator = document.getElementById('status-indicator');
            if (indicator) {
                if (data.success && data.rdkit_working) {
                    indicator.className = 'inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium bg-green-100 text-green-800';
                    indicator.innerHTML = '<div class="w-2 h-2 bg-green-400 rounded-full mr-1.5 animate-pulse"></div>Connected';
                } else {
                    indicator.className = 'inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium bg-red-100 text-red-800';
                    indicator.innerHTML = '<div class="w-2 h-2 bg-red-400 rounded-full mr-1.5"></div>Disconnected';
                }
            }
        } catch (error) {
            console.warn('Health check failed:', error);
            const indicator = document.getElementById('status-indicator');
            if (indicator) {
                indicator.className = 'inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium bg-yellow-100 text-yellow-800';
                indicator.innerHTML = '<div class="w-2 h-2 bg-yellow-400 rounded-full mr-1.5"></div>Checking...';
            }
        }
    }

    async loadHeroMolecule() {
        const heroContainer = document.getElementById('hero-molecule');
        if (!heroContainer) return;

        const exampleMolecules = [
            'CC(=O)OC1=CC=CC=C1C(=O)O', // Aspirin
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', // Caffeine  
            'CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C', // Testosterone
            'CCO' // Ethanol (fallback)
        ];

        const randomSmiles = exampleMolecules[Math.floor(Math.random() * exampleMolecules.length)];

        try {
            const response = await fetch('/api/v1/visualization/draw_svg', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ 
                    smiles: randomSmiles, 
                    width: 320, 
                    height: 320 
                })
            });

            const data = await response.json();
            if (data.success) {
                heroContainer.innerHTML = data.svg;
                heroContainer.classList.add('animate-fade-in');
            }
        } catch (error) {
            console.warn('Failed to load hero molecule:', error);
        }
    }

    animateStats() {
        const stats = document.querySelectorAll('[class*="text-4xl font-bold"]');
        
        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    const target = entry.target;
                    const finalValue = target.textContent;
                    
                    if (finalValue.includes('+') || finalValue.includes('∞')) {
                        // Handle special cases
                        target.style.opacity = '0';
                        setTimeout(() => {
                            target.style.opacity = '1';
                            target.classList.add('animate-fade-in');
                        }, Math.random() * 500);
                    } else if (!isNaN(parseInt(finalValue))) {
                        // Animate numbers
                        const numValue = parseInt(finalValue);
                        this.animateNumber(target, numValue);
                    }
                }
            });
        }, { threshold: 0.5 });

        stats.forEach(stat => observer.observe(stat));
    }

    animateNumber(element, target) {
        const duration = 2000;
        const start = 0;
        const startTime = performance.now();

        const animate = (currentTime) => {
            const elapsed = currentTime - startTime;
            const progress = Math.min(elapsed / duration, 1);
            
            // Easing function
            const easeOutCubic = 1 - Math.pow(1 - progress, 3);
            const current = Math.floor(start + (target - start) * easeOutCubic);
            
            element.textContent = current;
            
            if (progress < 1) {
                requestAnimationFrame(animate);
            } else {
                element.textContent = target;
            }
        };

        requestAnimationFrame(animate);
    }

    setupFeatureCards() {
        const cards = document.querySelectorAll('.feature-card');
        
        cards.forEach((card, index) => {
            // Stagger the animation
            setTimeout(() => {
                card.classList.add('animate-fade-in');
            }, index * 100);
            
            // Add enhanced hover effects
            card.addEventListener('mouseenter', function() {
                this.style.transform = 'translateY(-4px) scale(1.02)';
                this.style.boxShadow = '0 10px 25px rgba(0, 0, 0, 0.15)';
            });
            
            card.addEventListener('mouseleave', function() {
                this.style.transform = 'translateY(0) scale(1)';
                this.style.boxShadow = '';
            });
        });
    }

    // API utility methods
    async apiCall(endpoint, method = 'POST', data = null) {
        try {
            const options = {
                method: method,
                headers: {
                    'Content-Type': 'application/json',
                }
            };

            if (data) {
                options.body = JSON.stringify(data);
            }

            const response = await fetch(`${this.baseURL}${endpoint}`, options);
            const result = await response.json();

            if (!result.success) {
                throw new Error(result.error || 'Unknown error occurred');
            }

            return result;
        } catch (error) {
            console.error('API Error:', error);
            throw error;
        }
    }

    // UI utility methods
    showToast(message, type = 'info', duration = 5000) {
        const toast = document.createElement('div');
        const toastId = 'toast-' + Date.now();
        
        const typeClasses = {
            'success': 'bg-green-50 border-green-200 text-green-800',
            'error': 'bg-red-50 border-red-200 text-red-800',
            'warning': 'bg-yellow-50 border-yellow-200 text-yellow-800',
            'info': 'bg-blue-50 border-blue-200 text-blue-800'
        };

        const iconNames = {
            'success': 'check-circle',
            'error': 'alert-circle',
            'warning': 'alert-triangle',  
            'info': 'info'
        };

        toast.id = toastId;
        toast.className = `fixed top-4 right-4 z-50 max-w-sm p-4 rounded-lg shadow-lg border ${typeClasses[type] || typeClasses.info} animate-slide-up`;
        
        toast.innerHTML = `
            <div class="flex items-center space-x-2">
                <i data-lucide="${iconNames[type] || iconNames.info}" class="w-4 h-4 flex-shrink-0"></i>
                <span class="text-sm font-medium flex-1">${message}</span>
                <button onclick="this.closest('.fixed').remove()" class="flex-shrink-0 ml-2 p-1 rounded hover:bg-black hover:bg-opacity-10 transition-colors">
                    <i data-lucide="x" class="w-4 h-4"></i>
                </button>
            </div>
        `;

        document.body.appendChild(toast);
        
        // Initialize Lucide icons
        if (window.lucide) {
            lucide.createIcons();
        }

        // Auto-remove after duration
        setTimeout(() => {
            const toastElement = document.getElementById(toastId);
            if (toastElement) {
                toastElement.style.animation = 'slideUp 0.3s ease-out reverse';
                setTimeout(() => toastElement.remove(), 300);
            }
        }, duration);

        return toast;
    }

    showModal(title, content, actions = []) {
        const modal = document.createElement('div');
        modal.className = 'fixed inset-0 z-50 overflow-y-auto';
        modal.innerHTML = `
            <div class="flex items-center justify-center min-h-screen pt-4 px-4 pb-20 text-center sm:block sm:p-0">
                <div class="fixed inset-0 bg-gray-500 bg-opacity-75 transition-opacity" onclick="this.closest('.fixed').remove()"></div>
                <div class="inline-block align-bottom bg-white rounded-lg text-left overflow-hidden shadow-xl transform transition-all sm:my-8 sm:align-middle sm:max-w-lg sm:w-full">
                    <div class="bg-white px-4 pt-5 pb-4 sm:p-6 sm:pb-4">
                        <div class="sm:flex sm:items-start">
                            <div class="mt-3 text-center sm:mt-0 sm:text-left w-full">
                                <h3 class="text-lg leading-6 font-medium text-gray-900 mb-4">${title}</h3>
                                <div class="text-gray-600">${content}</div>
                            </div>
                        </div>
                    </div>
                    ${actions.length > 0 ? `
                        <div class="bg-gray-50 px-4 py-3 sm:px-6 sm:flex sm:flex-row-reverse">
                            ${actions.map(action => `
                                <button class="${action.class || 'bg-gray-300 text-gray-700'} ml-3 inline-flex justify-center rounded-md border border-transparent px-4 py-2 text-base font-medium shadow-sm hover:opacity-90 focus:outline-none focus:ring-2 focus:ring-offset-2 sm:ml-3 sm:w-auto sm:text-sm" 
                                        onclick="${action.onclick || 'this.closest(\'.fixed\').remove()'}">
                                    ${action.text}
                                </button>
                            `).join('')}
                        </div>
                    ` : ''}
                </div>
            </div>
        `;

        document.body.appendChild(modal);
        return modal;
    }

    closeAllModals() {
        document.querySelectorAll('.fixed.inset-0.z-50').forEach(modal => {
            modal.remove();
        });
    }

    // Data export utilities
    exportToCSV(data, filename) {
        let csv = '';
        
        if (Array.isArray(data) && data.length > 0) {
            // Array of objects
            const headers = Object.keys(data[0]);
            csv = headers.join(',') + '\\n';
            
            data.forEach(row => {
                csv += headers.map(header => {
                    const value = row[header];
                    return typeof value === 'string' ? `"${value.replace(/"/g, '""')}"` : value;
                }).join(',') + '\\n';
            });
        } else if (typeof data === 'object') {
            // Single object
            csv = 'Property,Value\\n';
            Object.entries(data).forEach(([key, value]) => {
                csv += `"${key}","${value}"\\n`;
            });
        }

        this.downloadFile(csv, filename, 'text/csv');
    }

    downloadFile(content, filename, mimeType = 'text/plain') {
        const blob = new Blob([content], { type: mimeType });
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);
    }

    // Format utilities
    formatMolecularName(smiles) {
        // Simple molecule name mapping
        const names = {
            'CCO': 'Ethanol',
            'CC(=O)OC1=CC=CC=C1C(=O)O': 'Aspirin',
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': 'Caffeine',
            'CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C': 'Testosterone'
        };
        return names[smiles] || 'Unknown Compound';
    }

    formatNumber(num, decimals = 2) {
        if (typeof num !== 'number') return 'N/A';
        if (Math.abs(num) < 0.001) return num.toExponential(2);
        return num.toFixed(decimals);
    }

    // Performance monitoring
    measurePerformance(name, fn) {
        const start = performance.now();
        const result = fn();
        const end = performance.now();
        return result;
    }

    // Enhanced error handling
    handleError(error, context = 'Unknown') {
        console.error(`Error in ${context}:`, error);
        
        let message = 'An unexpected error occurred';
        let severity = 'error';
        
        // Classify error types
        if (error.name === 'TypeError') {
            message = 'Invalid data format or structure';
        } else if (error.name === 'NetworkError' || error.message.includes('fetch')) {
            message = 'Network connection failed - please check your connection';
            severity = 'warning';
        } else if (error.name === 'SyntaxError' || error.message.includes('JSON')) {
            message = 'Server response format error - please try again';
        } else if (error.message.includes('timeout')) {
            message = 'Request timed out - please try again';
            severity = 'warning';
        } else if (error.message.includes('404')) {
            message = 'API endpoint not found - please contact support';
        } else if (error.message.includes('500')) {
            message = 'Server error - please try again later';
        } else if (error.message.includes('403')) {
            message = 'Access denied - please check permissions';
        } else if (error.message) {
            message = error.message;
        }

        // Log to persistent error storage
        this.logError(error, context);
        
        this.showToast(message, severity);
        return { error: true, message, context };
    }

    // Error logging for debugging
    logError(error, context) {
        try {
            const errorLog = {
                timestamp: new Date().toISOString(),
                context,
                message: error.message,
                stack: error.stack,
                userAgent: navigator.userAgent,
                url: window.location.href
            };

            // Store in local storage (keep last 50 errors)
            let errors = JSON.parse(localStorage.getItem('rdkit_errors') || '[]');
            errors.unshift(errorLog);
            errors = errors.slice(0, 50);
            localStorage.setItem('rdkit_errors', JSON.stringify(errors));
        } catch (e) {
            console.warn('Failed to log error:', e);
        }
    }

    // Get error logs for debugging
    getErrorLogs() {
        try {
            return JSON.parse(localStorage.getItem('rdkit_errors') || '[]');
        } catch (e) {
            return [];
        }
    }

    // Clear error logs
    clearErrorLogs() {
        localStorage.removeItem('rdkit_errors');
    }

    // Data persistence utilities
    saveToHistory(type, data) {
        try {
            const historyKey = `rdkit_history_${type}`;
            let history = JSON.parse(localStorage.getItem(historyKey) || '[]');
            
            const entry = {
                id: Date.now(),
                timestamp: new Date().toISOString(),
                data,
                type
            };
            
            // Add to beginning and limit to 100 entries
            history.unshift(entry);
            history = history.slice(0, 100);
            
            localStorage.setItem(historyKey, JSON.stringify(history));
            return entry.id;
        } catch (e) {
            console.warn('Failed to save to history:', e);
            return null;
        }
    }

    getHistory(type, limit = 50) {
        try {
            const historyKey = `rdkit_history_${type}`;
            const history = JSON.parse(localStorage.getItem(historyKey) || '[]');
            return history.slice(0, limit);
        } catch (e) {
            console.warn('Failed to get history:', e);
            return [];
        }
    }

    clearHistory(type = null) {
        try {
            if (type) {
                localStorage.removeItem(`rdkit_history_${type}`);
            } else {
                // Clear all history
                const keys = Object.keys(localStorage).filter(key => key.startsWith('rdkit_history_'));
                keys.forEach(key => localStorage.removeItem(key));
            }
        } catch (e) {
            console.warn('Failed to clear history:', e);
        }
    }

    // Save user preferences
    savePreference(key, value) {
        try {
            const prefs = JSON.parse(localStorage.getItem('rdkit_preferences') || '{}');
            prefs[key] = value;
            localStorage.setItem('rdkit_preferences', JSON.stringify(prefs));
        } catch (e) {
            console.warn('Failed to save preference:', e);
        }
    }

    getPreference(key, defaultValue = null) {
        try {
            const prefs = JSON.parse(localStorage.getItem('rdkit_preferences') || '{}');
            return prefs.hasOwnProperty(key) ? prefs[key] : defaultValue;
        } catch (e) {
            console.warn('Failed to get preference:', e);
            return defaultValue;
        }
    }

    // Session data management
    saveSessionData(key, data, expireMinutes = 60) {
        try {
            const sessionData = {
                data,
                timestamp: Date.now(),
                expires: Date.now() + (expireMinutes * 60 * 1000)
            };
            sessionStorage.setItem(`rdkit_session_${key}`, JSON.stringify(sessionData));
        } catch (e) {
            console.warn('Failed to save session data:', e);
        }
    }

    getSessionData(key) {
        try {
            const stored = sessionStorage.getItem(`rdkit_session_${key}`);
            if (!stored) return null;
            
            const sessionData = JSON.parse(stored);
            if (Date.now() > sessionData.expires) {
                sessionStorage.removeItem(`rdkit_session_${key}`);
                return null;
            }
            
            return sessionData.data;
        } catch (e) {
            console.warn('Failed to get session data:', e);
            return null;
        }
    }

    // Export all user data
    exportUserData() {
        const userData = {
            preferences: JSON.parse(localStorage.getItem('rdkit_preferences') || '{}'),
            history: {},
            errors: this.getErrorLogs(),
            exported: new Date().toISOString(),
            version: '1.0'
        };

        // Get all history types
        const historyTypes = ['search', 'descriptors', 'fingerprints', 'similarity', 'coordinates', 'properties', 'reactions', 'visualization'];
        historyTypes.forEach(type => {
            userData.history[type] = this.getHistory(type);
        });

        return userData;
    }

    // Import user data
    importUserData(userData) {
        try {
            if (userData.preferences) {
                localStorage.setItem('rdkit_preferences', JSON.stringify(userData.preferences));
            }

            if (userData.history) {
                Object.entries(userData.history).forEach(([type, history]) => {
                    localStorage.setItem(`rdkit_history_${type}`, JSON.stringify(history));
                });
            }

            this.showToast('Data imported successfully', 'success');
            return true;
        } catch (e) {
            this.handleError(e, 'Data Import');
            return false;
        }
    }

    // Get storage usage statistics
    getStorageStats() {
        const stats = {
            localStorage: {
                used: 0,
                keys: 0,
                rdkitKeys: 0
            },
            sessionStorage: {
                used: 0,
                keys: 0,
                rdkitKeys: 0
            }
        };

        // Calculate localStorage usage
        for (let key in localStorage) {
            if (localStorage.hasOwnProperty(key)) {
                const value = localStorage.getItem(key);
                stats.localStorage.used += key.length + (value ? value.length : 0);
                stats.localStorage.keys++;
                
                if (key.startsWith('rdkit_')) {
                    stats.localStorage.rdkitKeys++;
                }
            }
        }

        // Calculate sessionStorage usage
        for (let key in sessionStorage) {
            if (sessionStorage.hasOwnProperty(key)) {
                const value = sessionStorage.getItem(key);
                stats.sessionStorage.used += key.length + (value ? value.length : 0);
                stats.sessionStorage.keys++;
                
                if (key.startsWith('rdkit_')) {
                    stats.sessionStorage.rdkitKeys++;
                }
            }
        }

        // Convert bytes to human readable
        stats.localStorage.usedMB = (stats.localStorage.used / 1024 / 1024).toFixed(2);
        stats.sessionStorage.usedMB = (stats.sessionStorage.used / 1024 / 1024).toFixed(2);

        return stats;
    }
}

// Initialize the application
document.addEventListener('DOMContentLoaded', function() {
    window.app = new NuGenRDKitApp();
    
    // Initialize Lucide icons
    if (window.lucide) {
        lucide.createIcons();
    }
});

// Global utility functions for backward compatibility
function showAlert(message, type) {
    if (window.app) {
        window.app.showToast(message, type);
    }
}

function showLoading(show) {
    const overlay = document.getElementById('loading-overlay');
    if (overlay) {
        if (show) {
            overlay.classList.remove('hidden');
        } else {
            overlay.classList.add('hidden');
        }
    }
}

// Export for use in other scripts
if (typeof module !== 'undefined' && module.exports) {
    module.exports = NuGenRDKitApp;
}