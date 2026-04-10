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
    from Bio.Seq import Seq as BioSeq
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

REDIS_CACHE_URL = os.environ.get("REDIS_URL")
cache_config = {'CACHE_TYPE': 'SimpleCache', 'CACHE_DEFAULT_TIMEOUT': 86400}
if _REDIS_AVAILABLE and REDIS_CACHE_URL:
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
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
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

@app.route("/api/fetch", methods=["GET", "POST"])
def fetch_external():
    """Proxy for external bioinformatics databases."""
    # Cache stores serializable dicts: {"json": {...}} or {"text": "..."}
    cache_key = "proxy_cache_{}".format(request.full_path)
    cached_val = cache.get(cache_key)
    if cached_val:
        if isinstance(cached_val, dict) and "json" in cached_val:
            return jsonify(cached_val["json"])
        if isinstance(cached_val, dict) and "text" in cached_val:
            return Response(cached_val["text"], mimetype="text/plain")

    if request.method == "POST":
        data = request.get_json() or {}
        db_type = data.get("db", "ncbi")
        acc_id = data.get("id") or data.get("accession")
    else:
        db_type = request.args.get("db")
        acc_id = request.args.get("id")
    
    if not db_type or not acc_id:
        return jsonify({"error": "Missing db or id/accession parameter"}), 400

    try:
        if db_type == "ncbi":
            enforce_rate_limit("ncbi")
            # For primer design, we usually want nucleotide
            params = {
                "db": "nucleotide",
                "id": acc_id,
                "rettype": "gb",
                "retmode": "text"
            }
            # Determine which API Key to use (Header > ENV)
            user_key = request.headers.get("X-NCBI-API-Key", "").strip()
            effective_key = user_key or NCBI_API_KEY
            
            if effective_key and str(effective_key).strip():
                params["api_key"] = effective_key

            resp = requests.get("{}/efetch.fcgi".format(NCBI_BASE), params=params, timeout=15)
            resp.raise_for_status()
            raw_data = resp.text
            
            final_resp = None

            if _BIOPYTHON_AVAILABLE:
                try:
                    record = SeqIO.read(io.StringIO(raw_data), "genbank")
                    features = []
                    transcripts = {}
                    
                    for feat in record.features:
                        # Prompt 3 flat features
                        if feat.type in ["exon", "CDS", "variation", "SNP"]:
                            f_type = "SNP" if feat.type in ["variation", "SNP"] else feat.type
                            features.append({
                                "type": f_type,
                                "start": int(feat.location.start) + 1,
                                "end": int(feat.location.end)
                            })
                        
                        # Legacy transcript parsing
                        if feat.type in ["mRNA", "transcript"]:
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

                    # Association for legacy transcripts
                    for feat in record.features:
                        if feat.type in ["exon", "CDS"]:
                            tid = feat.qualifiers.get("transcript_id", [None])[0]
                            if not tid:
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

                    tx_list = list(transcripts.values())
                    if not tx_list:
                        tx_list = [{"exon_list": [], "cds_regions": [], "gene_range": [1, len(record.seq)], "strand": "+"}]

                    result = {
                        "sequence": str(record.seq),
                        "accession": acc_id,
                        "organism": record.annotations.get("organism", "Unknown"),
                        "length": len(record.seq),
                        "features": features,
                        "metadata": tx_list[0],
                        "transcripts": tx_list
                    }
                    cache.set(cache_key, {"json": result}, timeout=86400)
                    return jsonify(result)
                except Exception as e:
                    logger.warning("Biopython parsing failed, falling back to raw GenBank text: %s", e)
            
            # Fallback: return raw GenBank text
            cache.set(cache_key, {"text": raw_data}, timeout=86400)
            return Response(raw_data, mimetype="text/plain")
        
        elif db_type == "ncbiprotein":
            enforce_rate_limit("ncbi")
            params = {
                "db": "protein",
                "id": acc_id,
                "rettype": "fasta",
                "retmode": "text"
            }
            user_key = request.headers.get("X-NCBI-API-Key", "").strip()
            effective_key = user_key or NCBI_API_KEY
            if effective_key and str(effective_key).strip():
                params["api_key"] = effective_key
            resp = requests.get("{}/efetch.fcgi".format(NCBI_BASE), params=params, timeout=15)

        elif db_type == "uniprot":
            enforce_rate_limit("uniprot")
            # New UniProt REST API endpoint (uniprotkb)
            resp = requests.get("{}/uniprotkb/{}.fasta".format(UNIPROT_BASE, acc_id), timeout=15)

        elif db_type == "ensembl":
            enforce_rate_limit("ensembl")
            # Ensembl REST: use Accept header (not Content-Type) for GET requests
            json_headers = {"Accept": "application/json"}
            resp = requests.get(
                "{}/sequence/id/{}".format(ENSEMBL_BASE, acc_id),
                headers=json_headers,
                params={"expand_5prime": 2000, "expand_3prime": 2000},
                timeout=15
            )
            if resp.status_code == 200:
                data = resp.json()
                result = {
                    "sequence": data.get("seq"),
                    "metadata": {
                        "strand": "+" if data.get("strand", 1) >= 0 else "-",
                        "tss": 2001,
                        "exon_list": [],
                        "description": data.get("desc", "Ensembl Sequence")
                    }
                }
                cache.set(cache_key, {"json": result}, timeout=86400)
                return jsonify(result)
            # Fallback: plain FASTA
            fasta_headers = {"Accept": "text/x-fasta"}
            resp = requests.get("{}/sequence/id/{}".format(ENSEMBL_BASE, acc_id), headers=fasta_headers, timeout=15)

        elif db_type == "kegg":
            enforce_rate_limit("kegg")
            seq_type = request.args.get("type", "ntseq")
            if seq_type not in ("ntseq", "aaseq"):
                return jsonify({"error": "Invalid KEGG type. Use 'ntseq' or 'aaseq'."}), 400
            resp = requests.get("{}/get/{}/{}".format(KEGG_BASE, acc_id, seq_type), timeout=15)
            if resp.status_code == 200:
                text = resp.text.strip()
                if not text or text.startswith("<!"):
                    return jsonify({"error": "KEGG ID '{}' returned no sequence.".format(acc_id)}), 404
                return Response(text, mimetype="text/plain")
            elif resp.status_code == 404:
                return jsonify({"error": "KEGG ID '{}' not found.".format(acc_id)}), 404
            else:
                return jsonify({"error": "KEGG returned HTTP {}.".format(resp.status_code)}), resp.status_code

        elif db_type == "ncbisymbol":
            enforce_rate_limit("ncbi")
            species = request.args.get("species", "Homo sapiens")
            user_key = request.headers.get("X-NCBI-API-Key", "").strip()
            effective_key = user_key or NCBI_API_KEY

            # Step 1: ESearch to find gene ID
            search_params = {
                "db": "gene",
                "term": "{}[Gene] AND {}[Organism]".format(acc_id, species),
                "retmax": 1,
                "retmode": "json",
                "tool": "BioToolkit",
                "email": "admin@biotoolkit.dev"
            }
            if effective_key and str(effective_key).strip():
                search_params["api_key"] = effective_key

            s_resp = requests.get("{}/esearch.fcgi".format(NCBI_BASE), params=search_params, timeout=15)
            s_resp.raise_for_status()
            s_data = s_resp.json()
            ids = s_data.get("esearchresult", {}).get("idlist", [])

            if not ids:
                return jsonify({"error": "No results for symbol '{}' in '{}'.".format(acc_id, species)}), 404

            gene_id = str(ids[0])

            # Step 2: ELink to find associated nucleotide accession
            enforce_rate_limit("ncbi")
            link_params = {
                "dbfrom": "gene",
                "db": "nuccore",
                "id": gene_id,
                "retmode": "json",
                "tool": "BioToolkit",
                "email": "admin@biotoolkit.dev"
            }
            if effective_key and str(effective_key).strip():
                link_params["api_key"] = effective_key

            l_resp = requests.get("{}/elink.fcgi".format(NCBI_BASE), params=link_params, timeout=15)
            l_resp.raise_for_status()
            l_data = l_resp.json()

            nuc_id = None
            try:
                for linkset in l_data.get("linksets", []):
                    for ldb in linkset.get("linksetdbs", []):
                        if ldb.get("dbto") == "nuccore" and ldb.get("links"):
                            nuc_id = str(ldb["links"][0])
                            break
                    if nuc_id:
                        break
            except (KeyError, IndexError):
                pass

            if not nuc_id:
                return jsonify({"error": "No nucleotide record linked to gene '{}'. Try using an accession ID directly.".format(acc_id)}), 404

            # Step 3: EFetch the nucleotide record
            enforce_rate_limit("ncbi")
            fetch_params = {
                "db": "nucleotide",
                "id": nuc_id,
                "rettype": "fasta",
                "retmode": "text"
            }
            if effective_key and str(effective_key).strip():
                fetch_params["api_key"] = effective_key

            resp = requests.get("{}/efetch.fcgi".format(NCBI_BASE), params=fetch_params, timeout=15)

        else:
            return jsonify({"error": "Unsupported database: {}".format(db_type)}), 400

        resp.raise_for_status()
        cache.set(cache_key, {"text": resp.text}, timeout=86400)
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
            "retmode": "json"
        }
        if NCBI_API_KEY and str(NCBI_API_KEY).strip():
            params["api_key"] = NCBI_API_KEY
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
