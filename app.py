# -*- coding: utf-8 -*-
"""
BioToolkit – Flask API Server
==============================
Serves the static frontend and acts as a secure proxy for external
bioinformatics databases (NCBI, Ensembl, UniProt, KEGG) to bypass
browser CORS restrictions and protect API keys.

Deploy target: Render.com (Web Service, Python runtime)
"""

import os
import time
import json
import logging
import requests
import re
from functools import wraps
from flask import Flask, jsonify, request, send_from_directory, Response, g, url_for
from flask_cors import CORS
from flask_caching import Cache
from flask_login import LoginManager, UserMixin, login_required, current_user
from models import db, User, Workspace, AnalysisResult
import io
try:
    from Bio import SeqIO
    _BIOPYTHON_AVAILABLE = True
except ImportError:
    _BIOPYTHON_AVAILABLE = False

try:
    import redis
    _REDIS_AVAILABLE = True
except ImportError:
    _REDIS_AVAILABLE = False

# ─────────────────────────────────────────────────────────────
# App Setup
# ─────────────────────────────────────────────────────────────
app = Flask(__name__, static_folder=".", static_url_path="")
app.config['SECRET_KEY'] = os.environ.get("SECRET_KEY", "dev-secret-key-12345")
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get("DATABASE_URL", "sqlite:///biotoolkit.db")
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# Initialize Extensions
db.init_app(app)
login_manager = LoginManager()
login_manager.init_app(app)

REDIS_CACHE_URL = os.environ.get("REDIS_URL", "redis://localhost:6379/0")
cache_config = {'CACHE_TYPE': 'SimpleCache', 'CACHE_DEFAULT_TIMEOUT': 86400}
if _REDIS_AVAILABLE:
    cache_config = {
        'CACHE_TYPE': 'RedisCache',
        'CACHE_REDIS_URL': REDIS_CACHE_URL,
        'CACHE_DEFAULT_TIMEOUT': 86400
    }

cache = Cache(app, config=cache_config)

# Database table creation (Preparation Phase)
with app.app_context():
    try:
        db.create_all()
    except Exception as e:
        app.logger.warning("Database initialization deferred: {}".format(e))

# 1. CORS Enforcement (Restrict to exact frontend domain)
ALLOWED_ORIGIN = os.environ.get("ALLOWED_ORIGIN", "http://localhost:5173") # Default to local dev
CORS(app, resources={r"/api/*": {"origins": ALLOWED_ORIGIN}})

# 2. Payload size Limit (5MB) - Prevents OOM and large sequence abuse
app.config['MAX_CONTENT_LENGTH'] = 5 * 1024 * 1024
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("BioToolkit")

# ─────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_BASE = "https://rest.ensembl.org"
UNIPROT_BASE = "https://rest.uniprot.org"
KEGG_BASE = "https://rest.kegg.jp"

REDIS_URL = os.environ.get("REDIS_URL", "")
redis_client = None
if _REDIS_AVAILABLE:
    try:
        redis_client = redis.from_url(REDIS_URL) if REDIS_URL else redis.Redis(host='localhost', port=6379, db=0)
        redis_client.ping()
    except Exception as e:
        logger.warning("Redis unavailable for rate limiting, falling back to in-memory: %s", e)
        redis_client = None

# Server-side rate limit tracking (per database)
_last_request_time = {}
RATE_LIMITS = {
    "ncbi": 0.34,      # ~3 req/sec without key, ~10 req/sec with key
    "ensembl": 1.0,    # 1 req/sec per Ensembl guidelines
    "uniprot": 0.5,
    "kegg": 0.5
}

def enforce_rate_limit(db_name):
    """Scientific standard rate limiting to prevent API blocking."""
    if db_name not in RATE_LIMITS:
        return
    
    now = time.time()
    last = _last_request_time.get(db_name, 0)
    wait = RATE_LIMITS[db_name] - (now - last)
    
    if wait > 0:
        time.sleep(wait)
    _last_request_time[db_name] = time.time()

# ─────────────────────────────────────────────────────────────
# Auth / Session
# ─────────────────────────────────────────────────────────────
@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))

# ─────────────────────────────────────────────────────────────
# Static Content Routes
# ─────────────────────────────────────────────────────────────
@app.route("/")
def index():
    return send_from_directory(app.static_folder, "index.html")

@app.route("/js/<path:path>")
def send_js(path):
    return send_from_directory(os.path.join(app.static_folder, "js"), path)

@app.route("/css/<path:path>")
def send_css(path):
    return send_from_directory(os.path.join(app.static_folder, "css"), path)

# ─────────────────────────────────────────────────────────────
# API Endpoints
# ─────────────────────────────────────────────────────────────

@app.route("/api/fetch", methods=["GET"])
@cache.cached(timeout=86400, query_string=True)
def fetch_external():
    """Proxy for external bioinformatics databases."""
    db_type = request.args.get("db")
    acc_id = request.args.get("id")
    
    if not db_type or not acc_id:
        return jsonify({"error": "Missing db or id parameter"}), 400

    try:
        if db_type == "ncbi":
            enforce_rate_limit("ncbi")
            is_nt = db_type == "ncbi" and request.args.get("db", "nucleotide") == "nucleotide"
            params = {
                "db": "nucleotide" if is_nt else "protein",
                "id": acc_id,
                "rettype": "gb" if is_nt else "fasta",
                "retmode": "text",
                "api_key": NCBI_API_KEY
            }
            resp = requests.get("{}/efetch.fcgi".format(NCBI_BASE), params=params, timeout=15)
            resp.raise_for_status()
            raw_data = resp.text

            if is_nt and _BIOPYTHON_AVAILABLE:
                try:
                    record = SeqIO.read(io.StringIO(raw_data), "genbank")
                    transcripts = {}
                    
                    # 1. Primary pass: Discover all mRNA/transcript models
                    for feat in record.features:
                        if feat.type in ["mRNA", "transcript", "mRNA"]:
                            tid = feat.qualifiers.get("transcript_id", [None])[0] or \
                                  feat.qualifiers.get("locus_tag", [None])[0] or \
                                  str(feat.location)
                            transcripts[tid] = {
                                "id": tid,
                                "strand": "+" if feat.location.strand >= 0 else "-",
                                "exon_list": [],
                                "cds_regions": [],
                                "gene_range": [int(feat.location.start) + 1, int(feat.location.end)],
                                "description": feat.qualifiers.get("product", [record.description])[0]
                            }

                    # 2. Association pass: Link exons and CDS to transcripts
                    for feat in record.features:
                        if feat.type in ["exon", "CDS"]:
                            tid = feat.qualifiers.get("transcript_id", [None])[0]
                            if not tid:
                                # Fallback: Positional overlap with known transcript
                                for mtid, mdata in transcripts.items():
                                    if int(feat.location.start) >= mdata["gene_range"][0]-1 and \
                                       int(feat.location.end) <= mdata["gene_range"][1]:
                                        tid = mtid
                                        break
                            
                            if tid in transcripts:
                                if feat.type == "exon":
                                    transcripts[tid]["exon_list"].append([int(feat.location.start) + 1, int(feat.location.end)])
                                else:
                                    parts = feat.location.parts if hasattr(feat.location, 'parts') else [feat.location]
                                    transcripts[tid]["cds_regions"] = [[int(p.start) + 1, int(p.end)] for p in parts]

                    # 3. Finalize transcript list with TSS and sorting
                    tx_list = list(transcripts.values())
                    for tx in tx_list:
                        if tx["exon_list"]:
                            tx["exon_list"].sort(key=lambda x: x[0])
                            tx["tss"] = tx["exon_list"][0][0] if tx["strand"] == "+" else tx["exon_list"][-1][1]
                        else:
                            tx["tss"] = tx["gene_range"][0] if tx["strand"] == "+" else tx["gene_range"][1]
                    
                    # If no transcripts found (prokaryotic or simple FASTA-like GB), fallback to gene-level
                    if not tx_list:
                        metadata = {
                            "strand": "+", "tss": 1, "exon_list": [], "cds_regions": [],
                            "gene_range": [1, len(record.seq)], "description": record.description
                        }
                        for feat in record.features:
                            if feat.type == "gene":
                                metadata["gene_range"] = [int(feat.location.start) + 1, int(feat.location.end)]
                                metadata["strand"] = "+" if feat.location.strand >= 0 else "-"
                            elif feat.type == "CDS":
                                parts = feat.location.parts if hasattr(feat.location, 'parts') else [feat.location]
                                metadata["cds_regions"] = [[int(p.start) + 1, int(p.end)] for p in parts]
                        tx_list = [metadata]

                    return jsonify({
                        "sequence": str(record.seq),
                        "metadata": tx_list[0], # Primary transcript for backward compat
                        "transcripts": tx_list
                    })
                except Exception as e:
                    logger.warning("Isoform parsing failed: %s", e)

            return Response(raw_data, mimetype="text/plain")
        
        elif db_type == "ncbiprotein":
            enforce_rate_limit("ncbi")
            params = {
                "db": "protein",
                "id": acc_id,
                "rettype": "fasta",
                "retmode": "text",
                "api_key": NCBI_API_KEY
            }
            resp = requests.get("{}/efetch.fcgi".format(NCBI_BASE), params=params, timeout=15)

        elif db_type == "uniprot":
            enforce_rate_limit("uniprot")
            resp = requests.get("{}/{}.fasta".format(UNIPROT_BASE, acc_id), timeout=15)

        elif db_type == "ensembl":
            enforce_rate_limit("ensembl")
            # Ensembl metadata is easier via JSON, but the user expects sequence.
            # We'll fetch sequence and try to get metadata separately or from headers.
            headers = {"Content-Type": "application/json"}
            resp = requests.get("{}/sequence/id/{}?expand_5prime=2000;expand_3prime=2000".format(ENSEMBL_BASE, acc_id), headers=headers, timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                # Return structured if Ensembl gives it
                return jsonify({
                    "sequence": data.get("seq"),
                    "metadata": {
                        "strand": "+" if data.get("strand") >= 0 else "-",
                        "tss": 2001, # We expanded by 2000
                        "exon_list": [],
                        "description": data.get("desc", "Ensembl Sequence")
                    }
                })
            # Fallback to plain FASTA if JSON fails or not found
            headers = {"Content-Type": "text/x-fasta"}
            resp = requests.get("{}/sequence/id/{}".format(ENSEMBL_BASE, acc_id), headers=headers, timeout=15)

        else:
            return jsonify({"error": "Unsupported database: {}".format(db_type)}), 400

        resp.raise_for_status()
        return Response(resp.text, mimetype="text/plain")

    except requests.exceptions.RequestException as e:
        logger.error("External API error: %s", str(e))
        return jsonify({"error": "Failed to fetch from biological database: {}".format(str(e))}), 502

@app.route("/api/search", methods=["POST"])
def search_ncbi():
    """Proxy for NCBI ESearch."""
    enforce_rate_limit("ncbi")
    data = request.get_json()
    query = data.get("query")
    db_name = data.get("db", "nucleotide")
    
    if not query:
        return jsonify({"error": "Query required"}), 400

    try:
        params = {
            "db": db_name,
            "term": query,
            "retmode": "json",
            "api_key": NCBI_API_KEY
        }
        resp = requests.get("{}/esearch.fcgi".format(NCBI_BASE), params=params, timeout=10)
        resp.raise_for_status()
        return jsonify(resp.json())
    except Exception as e:
        return jsonify({"error": str(e)}), 502

# ─────────────────────────────────────────────────────────────
# BioKit Local (Non-Worker) Processing
# ─────────────────────────────────────────────────────────────
@app.route("/api/analyze/basic", methods=["POST"])
def analyze_basic():
    """Handle basic sequence metrics on server side if Biopython available."""
    if not _BIOPYTHON_AVAILABLE:
        return jsonify({"error": "Biopython not configured on host"}), 503
    
    data = request.get_json()
    seq_str = data.get("sequence", "").upper()
    if not seq_str:
        return jsonify({"error": "No sequence"}), 400

    try:
        bio_seq = BioSeq(seq_str)
        # Handle nucleotide vs protein
        is_protein = any(aa in seq_str for aa in "EFILPQZ")
        
        results = {"length": len(seq_str)}
        if not is_protein:
            results["gc_content"] = (seq_str.count("G") + seq_str.count("C")) / len(seq_str) * 100
            results["complement"] = str(bio_seq.complement())
            results["rev_complement"] = str(bio_seq.reverse_complement())
        
        return jsonify(results)
    except Exception as e:
        return jsonify({"error": str(e)}), 400

# ─────────────────────────────────────────────────────────────
# Error Handlers
# ─────────────────────────────────────────────────────────────
@app.errorhandler(413)
def request_entity_too_large(error):
    return jsonify({"error": "Sequence too large (Limit 5MB)"}), 413

@app.errorhandler(404)
def not_found(error):
    return jsonify({"error": "Resource not found"}), 404

# ─────────────────────────────────────────────────────────────
# Entry Point
# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Create DB on startup if not existing
    with app.app_context():
        db.create_all()
    
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)
