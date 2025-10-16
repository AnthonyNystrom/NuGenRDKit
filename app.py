from flask import Flask, request, jsonify, send_file, render_template
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
import logging
import os
import io
import base64
from functools import wraps
import traceback
import warnings

# Suppress RDKit warnings and deprecation messages
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', message='.*MorganGenerator.*')
warnings.filterwarnings('ignore', message='.*DEPRECATION.*')

# Suppress RDKit logger warnings
from rdkit import RDLogger
rdkit_logger = RDLogger.logger()
rdkit_logger.setLevel(RDLogger.ERROR)

from routes.molecular_structure import molecular_structure_bp
from routes.descriptors import descriptors_bp
from routes.fingerprints import fingerprints_bp
from routes.similarity import similarity_bp
from routes.coordinates import coordinates_bp
from routes.properties import properties_bp
from routes.reactions import reactions_bp
from routes.visualization import visualization_bp

app = Flask(__name__)
CORS(app)

# Suppress Flask-Limiter warning about in-memory storage
warnings.filterwarnings('ignore', message='.*in-memory storage.*')

# Rate limiting
limiter = Limiter(
    key_func=get_remote_address,
    app=app,
    default_limits=["1000 per hour", "100 per minute"],
    storage_uri="memory://"
)

# Configure logging - reduce verbosity
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress werkzeug (Flask) INFO logs
logging.getLogger('werkzeug').setLevel(logging.WARNING)

app.register_blueprint(molecular_structure_bp, url_prefix='/api/v1/structure')
app.register_blueprint(descriptors_bp, url_prefix='/api/v1/descriptors')
app.register_blueprint(fingerprints_bp, url_prefix='/api/v1/fingerprints')
app.register_blueprint(similarity_bp, url_prefix='/api/v1/similarity')
app.register_blueprint(coordinates_bp, url_prefix='/api/v1/coordinates')
app.register_blueprint(properties_bp, url_prefix='/api/v1/properties')
app.register_blueprint(reactions_bp, url_prefix='/api/v1/reactions')
app.register_blueprint(visualization_bp, url_prefix='/api/v1/visualization')

def handle_rdkit_errors(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as e:
            logger.error(f"RDKit error in {f.__name__}: {str(e)}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            return jsonify({'error': str(e), 'success': False}), 400
    return decorated_function

# Global error handlers
@app.errorhandler(400)
def bad_request(error):
    logger.warning(f"Bad request: {error}")
    return jsonify({'error': 'Bad request', 'success': False}), 400

@app.errorhandler(404)
def not_found(error):
    logger.warning(f"Not found: {error}")
    return jsonify({'error': 'Endpoint not found', 'success': False}), 404

@app.errorhandler(405)
def method_not_allowed(error):
    logger.warning(f"Method not allowed: {error}")
    return jsonify({'error': 'Method not allowed', 'success': False}), 405

@app.errorhandler(429)
def ratelimit_handler(e):
    logger.warning(f"Rate limit exceeded: {e}")
    return jsonify({'error': 'Rate limit exceeded', 'retry_after': str(e.retry_after), 'success': False}), 429

@app.errorhandler(500)
def internal_error(error):
    logger.error(f"Internal server error: {error}")
    logger.error(f"Traceback: {traceback.format_exc()}")
    return jsonify({'error': 'Internal server error', 'success': False}), 500

@app.before_request
def validate_json():
    """Validate JSON for POST requests"""
    if request.method == 'POST' and request.content_type == 'application/json':
        try:
            if not request.get_json():
                return jsonify({'error': 'Invalid or empty JSON', 'success': False}), 400
        except Exception as e:
            logger.warning(f"JSON parsing error: {str(e)}")
            return jsonify({'error': 'Invalid JSON format', 'success': False}), 400

@app.after_request
def add_security_headers(response):
    """Add security headers"""
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['X-Frame-Options'] = 'DENY'
    response.headers['X-XSS-Protection'] = '1; mode=block'
    return response

@app.route('/')
def index():
    try:
        return render_template('index.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500

@app.route('/structure')
def structure():
    try:
        return render_template('structure.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500

@app.route('/descriptors')
def descriptors():
    try:
        return render_template('descriptors.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500

@app.route('/fingerprints')
def fingerprints():
    try:
        return render_template('fingerprints.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500

@app.route('/similarity')
def similarity():
    try:
        return render_template('similarity.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500

@app.route('/coordinates')
def coordinates():
    try:
        return render_template('coordinates.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500

@app.route('/properties')
def properties():
    try:
        return render_template('properties.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500

@app.route('/reactions')
def reactions():
    try:
        return render_template('reactions.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500

@app.route('/visualization')
def visualization():
    try:
        return render_template('visualization.html')
    except Exception as e:
        logger.error(f"Error rendering template: {str(e)}")
        return jsonify({'error': f'Template error: {str(e)}'}), 500


@app.route('/api')
def api_info():
    return jsonify({
        'message': 'RDKit Flask API',
        'version': '1.0.0',
        'success': True,
        'endpoints': {
            'structure': '/api/v1/structure',
            'descriptors': '/api/v1/descriptors', 
            'fingerprints': '/api/v1/fingerprints',
            'similarity': '/api/v1/similarity',
            'coordinates': '/api/v1/coordinates',
            'properties': '/api/v1/properties',
            'reactions': '/api/v1/reactions',
            'visualization': '/api/v1/visualization'
        }
    })

@app.route('/health')
def health():
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles('C')
        return jsonify({
            'status': 'healthy', 
            'rdkit_working': mol is not None,
            'success': True
        })
    except Exception as e:
        return jsonify({
            'status': 'unhealthy', 
            'error': str(e),
            'success': False
        }), 500

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=8000)