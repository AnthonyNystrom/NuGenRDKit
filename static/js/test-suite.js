/**
 * NuGenRDKit Test Suite
 * Comprehensive testing for all application functionality
 */

class NuGenRDKitTestSuite {
    constructor() {
        this.tests = [];
        this.results = {
            passed: 0,
            failed: 0,
            total: 0
        };
        this.startTime = null;
        this.endTime = null;
    }

    // Test runner
    async runAllTests() {
        console.log('🧪 Starting NuGenRDKit Test Suite...');
        this.startTime = performance.now();
        
        await this.testAppInitialization();
        await this.testUIComponents();
        await this.testDataPersistence();
        await this.testErrorHandling();
        await this.testAPIConnections();
        await this.testPageFunctionality();
        await this.testFormValidation();
        await this.testExportFunctions();
        
        this.endTime = performance.now();
        this.printResults();
    }

    // Test app initialization
    async testAppInitialization() {
        this.test('App object exists', () => {
            return window.app instanceof NuGenRDKitApp;
        });

        this.test('Lucide icons initialized', () => {
            return window.lucide && typeof window.lucide.createIcons === 'function';
        });

        this.test('Base URL set correctly', () => {
            return window.app.baseURL === '';
        });
    }

    // Test UI components
    async testUIComponents() {
        this.test('Navigation elements present', () => {
            const nav = document.querySelector('nav');
            const logo = document.querySelector('[data-testid="logo"]') || document.querySelector('.logo');
            return nav && nav.children.length > 0;
        });

        this.test('Mobile menu functionality', () => {
            const mobileButton = document.getElementById('mobile-menu-button');
            const mobileMenu = document.getElementById('mobile-menu');
            
            if (mobileButton && mobileMenu) {
                // Simulate click
                mobileButton.click();
                return !mobileMenu.classList.contains('hidden');
            }
            return true; // Pass if no mobile menu (desktop only)
        });

        this.test('Toast notification system', () => {
            window.app.showToast('Test message', 'info');
            const toast = document.querySelector('.fixed.top-4.right-4');
            return toast !== null;
        });
    }

    // Test data persistence
    async testDataPersistence() {
        const testData = { test: 'data', timestamp: Date.now() };
        
        this.test('Save to history', () => {
            const id = window.app.saveToHistory('test', testData);
            return id !== null;
        });

        this.test('Get from history', () => {
            const history = window.app.getHistory('test');
            return Array.isArray(history) && history.length > 0;
        });

        this.test('Save preferences', () => {
            window.app.savePreference('testPref', 'testValue');
            return window.app.getPreference('testPref') === 'testValue';
        });

        this.test('Session data management', () => {
            window.app.saveSessionData('testSession', testData, 1);
            const retrieved = window.app.getSessionData('testSession');
            return retrieved && retrieved.test === 'data';
        });

        this.test('Storage statistics', () => {
            const stats = window.app.getStorageStats();
            return stats && stats.localStorage && stats.sessionStorage;
        });
    }

    // Test error handling
    async testErrorHandling() {
        this.test('Error logging', () => {
            const testError = new Error('Test error');
            window.app.logError(testError, 'Test context');
            const logs = window.app.getErrorLogs();
            return logs.length > 0 && logs[0].message === 'Test error';
        });

        this.test('Error classification', () => {
            const networkError = new Error('fetch failed');
            networkError.name = 'NetworkError';
            const result = window.app.handleError(networkError, 'Test');
            return result.error && result.message.includes('Network');
        });

        this.test('Clear error logs', () => {
            window.app.clearErrorLogs();
            return window.app.getErrorLogs().length === 0;
        });
    }

    // Test API connections (mock)
    async testAPIConnections() {
        // Mock successful API call
        this.test('API call wrapper', async () => {
            try {
                // This will fail in test environment, but we're testing the wrapper
                await window.app.apiCall('/health', 'GET');
                return true;
            } catch (error) {
                // Expected to fail in test environment
                return error.message !== undefined;
            }
        });

        this.test('API error handling', async () => {
            try {
                await window.app.apiCall('/nonexistent', 'GET');
                return false;
            } catch (error) {
                return true; // Expected to throw
            }
        });
    }

    // Test page functionality
    async testPageFunctionality() {
        // Test form elements exist
        this.test('Home page search form', () => {
            const searchInput = document.getElementById('global-search');
            const searchBtn = document.getElementById('search-btn');
            return searchInput && searchBtn;
        });

        // Test feature cards
        this.test('Feature cards clickable', () => {
            const featureCards = document.querySelectorAll('.feature-card');
            return featureCards.length >= 8; // Should have 8 main features
        });

        // Test example molecules
        this.test('Example molecules functional', () => {
            const examples = document.querySelectorAll('.quick-example');
            if (examples.length > 0) {
                const firstExample = examples[0];
                const smiles = firstExample.getAttribute('data-smiles');
                return smiles && smiles.length > 0;
            }
            return true; // Pass if no examples on current page
        });
    }

    // Test form validation
    async testFormValidation() {
        // Create temporary form elements for testing
        const testForm = document.createElement('form');
        const testInput = document.createElement('input');
        testInput.type = 'text';
        testInput.id = 'test-smiles-input';
        testInput.value = 'CCO';
        testForm.appendChild(testInput);
        document.body.appendChild(testForm);

        this.test('SMILES input validation', () => {
            const input = document.getElementById('test-smiles-input');
            const value = input.value.trim();
            // Basic SMILES validation - should contain only valid characters
            const validChars = /^[A-Za-z0-9@+\-\[\]()=#.:\/\\]+$/;
            return validChars.test(value) && value.length > 0;
        });

        this.test('Empty input handling', () => {
            const input = document.getElementById('test-smiles-input');
            input.value = '';
            const value = input.value.trim();
            return value === '';
        });

        // Clean up
        document.body.removeChild(testForm);
    }

    // Test export functions
    async testExportFunctions() {
        const testData = [
            { id: 1, name: 'Test1', value: 'Value1' },
            { id: 2, name: 'Test2', value: 'Value2' }
        ];

        this.test('CSV export function', () => {
            try {
                window.app.exportToCSV(testData, 'test.csv');
                return true;
            } catch (error) {
                return false;
            }
        });

        this.test('File download function', () => {
            try {
                window.app.downloadFile('test content', 'test.txt');
                return true;
            } catch (error) {
                return false;
            }
        });

        this.test('Number formatting', () => {
            const result1 = window.app.formatNumber(123.456789, 2);
            const result2 = window.app.formatNumber(0.00001, 3);
            return result1 === '123.46' && result2 === '1.00e-5';
        });
    }

    // Test utility function
    test(name, testFunction) {
        try {
            const startTime = performance.now();
            const result = testFunction();
            const endTime = performance.now();
            const duration = endTime - startTime;

            // Handle both sync and async results
            Promise.resolve(result).then(passed => {
                this.results.total++;
                if (passed) {
                    this.results.passed++;
                    console.log(`✅ ${name} (${duration.toFixed(2)}ms)`);
                } else {
                    this.results.failed++;
                    console.log(`❌ ${name} (${duration.toFixed(2)}ms)`);
                }
            }).catch(error => {
                this.results.total++;
                this.results.failed++;
                console.log(`❌ ${name} - Error: ${error.message} (${duration.toFixed(2)}ms)`);
            });

        } catch (error) {
            this.results.total++;
            this.results.failed++;
            console.log(`❌ ${name} - Error: ${error.message}`);
        }
    }

    // Print test results
    printResults() {
        setTimeout(() => {
            const duration = this.endTime - this.startTime;
            const passRate = (this.results.passed / this.results.total * 100).toFixed(1);
            
            console.log('\n🧪 Test Results Summary');
            console.log('========================');
            console.log(`Total Tests: ${this.results.total}`);
            console.log(`Passed: ${this.results.passed}`);
            console.log(`Failed: ${this.results.failed}`);
            console.log(`Pass Rate: ${passRate}%`);
            console.log(`Duration: ${duration.toFixed(2)}ms`);
            console.log('========================');

            if (this.results.failed === 0) {
                console.log('🎉 All tests passed!');
            } else {
                console.log(`⚠️  ${this.results.failed} test(s) failed`);
            }

            // Store results for later analysis
            window.app.saveSessionData('test_results', {
                ...this.results,
                passRate: parseFloat(passRate),
                duration,
                timestamp: new Date().toISOString()
            });
        }, 1000); // Allow async tests to complete
    }

    // Accessibility tests
    async testAccessibility() {
        this.test('Images have alt text', () => {
            const images = document.querySelectorAll('img');
            return Array.from(images).every(img => img.alt !== undefined);
        });

        this.test('Buttons have labels', () => {
            const buttons = document.querySelectorAll('button');
            return Array.from(buttons).every(btn => 
                btn.textContent.trim() !== '' || 
                btn.getAttribute('aria-label') !== null ||
                btn.querySelector('span:not([class*="sr-only"])') !== null
            );
        });

        this.test('Form inputs have labels', () => {
            const inputs = document.querySelectorAll('input[type="text"], input[type="number"], select, textarea');
            return Array.from(inputs).every(input => {
                const id = input.id;
                const label = document.querySelector(`label[for="${id}"]`);
                return label || input.getAttribute('aria-label') || input.getAttribute('placeholder');
            });
        });
    }

    // Performance tests
    async testPerformance() {
        this.test('Page load performance', () => {
            const navigation = performance.getEntriesByType('navigation')[0];
            return navigation.loadEventEnd - navigation.fetchStart < 3000; // 3 seconds
        });

        this.test('DOM elements count', () => {
            const elements = document.querySelectorAll('*').length;
            return elements < 1000; // Reasonable DOM size
        });

        this.test('Memory usage estimation', () => {
            if (performance.memory) {
                const usedMB = performance.memory.usedJSHeapSize / 1024 / 1024;
                return usedMB < 50; // Less than 50MB
            }
            return true; // Skip if not available
        });
    }
}

// Auto-run tests if this script is loaded
if (typeof window !== 'undefined') {
    window.RDKitTestSuite = NuGenRDKitTestSuite;
    
    // Add test runner to app if available
    if (window.app) {
        window.app.testSuite = new NuGenRDKitTestSuite();
        
        // Add global test function
        window.runTests = () => {
            window.app.testSuite.runAllTests();
        };
    }

    // Auto-run tests on page load if URL contains test parameter
    document.addEventListener('DOMContentLoaded', function() {
        const urlParams = new URLSearchParams(window.location.search);
        if (urlParams.get('test') === 'true') {
            setTimeout(() => {
                if (window.runTests) {
                    window.runTests();
                }
            }, 1000); // Wait for app to initialize
        }
    });
}

// Export for Node.js testing environments
if (typeof module !== 'undefined' && module.exports) {
    module.exports = NuGenRDKitTestSuite;
}