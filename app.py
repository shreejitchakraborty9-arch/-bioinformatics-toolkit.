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
try:
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
        app.logger.warning(f"Database initialization deferred: {e}")

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
    "kegg": 0.5,
}

REQUEST_TIMEOUT = 15  # seconds

# ─────────────────────────────────────────────────────────────
# Resilience State (Circuit Breaker)
# ─────────────────────────────────────────────────────────────
_cb_failure_count = 0
_cb_open_until = 0

# ─────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────
def error_response(message: str, status_code: int = 400, details: dict = None):
    """Standardized scientific error responder."""
    payload = {"error": message, "status": "failure"}
    if details:
        payload["details"] = details
    return jsonify(payload), status_code

def enforce_rate_limit(db_name: str):
    """
    Redis-backed server-side per-database rate limiter for distributed workers.
    Gracefully degrades to in-memory if Redis is unavailable.
    """
    min_interval = RATE_LIMITS.get(db_name, 0.5)
    
    # Specific adjustment for NCBI with API KEY
    if db_name == "ncbi" and NCBI_API_KEY:
        min_interval = 0.1  # 10 req/sec
    if redis_client:
        try:
            key = f"rate_limit:{db_name}"
            # Atomic Lua script for rate limit calculation
            lua_script = """
            local last = redis.call('GET', KEYS[1])
            local current_time = tonumber(ARGV[1])
            local min_interval = tonumber(ARGV[2])
            
            if not last then
                redis.call('SET', KEYS[1], current_time)
                return 0
            end
            local elapsed = current_time - tonumber(last)
            if elapsed < min_interval then
                local wait_time = min_interval - elapsed
                redis.call('SET', KEYS[1], current_time + wait_time)
                return wait_time
            else
                redis.call('SET', KEYS[1], current_time)
                return 0
            end
            """
            
            wait_time = redis_client.eval(lua_script, 1, key, time.time(), min_interval)
            
            if wait_time and float(wait_time) > 0:
                time.sleep(float(wait_time))
            return
        except Exception as e:
            logger.error("Redis rate limit failed for %s, falling back to memory: %s", db_name, e)
            
    # In-memory fallback
    last = _last_request_time.get(db_name, 0)
    elapsed = time.time() - last
    if elapsed < min_interval:
        time.sleep(min_interval - elapsed)
        _last_request_time[db_name] = time.time() + (min_interval - elapsed)
    else:
        _last_request_time[db_name] = time.time()


# ─────────────────────────────────────────────────────────────
# Auth & Session Preparation (Placeholders)
# ─────────────────────────────────────────────────────────────
@login_manager.user_loader
def load_user(user_id):
    # Prepared for future DB-backed auth
    return User.query.get(int(user_id)) if user_id else None

def optional_login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        # Currently allows all users (Guest mode enabled)
        return f(*args, **kwargs)
    return decorated_function

def proxy_get(url: str, headers: dict = None, params: dict = None) -> requests.Response:
    """
    Make a proxied GET request with a consistent User-Agent,
    retry logic (exponential backoff), and a circuit breaker.
    """
    global _cb_failure_count, _cb_open_until
    
    service = "NCBI" if "ncbi" in url else "Ensembl" if "ensembl" in url else "UniProt" if "uniprot" in url else "Biological service"
    
    # 2. Check Circuit Breaker
    if time.time() < _cb_open_until:
        wait_left = int(_cb_open_until - time.time())
        logger.warning("[Circuit Breaker] Blocking request to %s (Open for %ds more)", service, wait_left)
        raise requests.RequestException(f"{service} is temporarily unavailable due to consecutive failures. Please try again in {wait_left} seconds.")
    
    h = {
        "User-Agent": "BioToolkit/1.0-MVP (contact: admin@biotoolkit.dev)",
        "Accept": "text/plain, application/json, text/x-fasta",
    }
    if headers:
        h.update(headers)
    
    backoffs = [0.5, 1.0, 2.0]
    total_attempts = len(backoffs) + 1
    for i in range(total_attempts):
        try:
            resp = requests.get(url, headers=h, params=params, timeout=15)
            
            if resp.status_code == 200 or (resp.status_code < 500 and resp.status_code not in (429,)):
                _cb_failure_count = 0 # Reset on server-reachable response
                return resp
            
            if i < len(backoffs):
                logger.warning("[Retry %d] %s returned HTTP %d. Backing off %ss.", i + 1, service, resp.status_code, backoffs[i])
                time.sleep(backoffs[i])
                continue
            else:
                break
        except (requests.Timeout, requests.RequestException) as e:
            if i < len(backoffs):
                logger.warning("[Retry %d] Connection to %s failed: %s. Backing off %ss.", i + 1, service, e, backoffs[i])
                time.sleep(backoffs[i])
                continue
            break
            
    _cb_failure_count += 1
    if _cb_failure_count >= 5:
        _cb_open_until = time.time() + 30
        logger.error("[Circuit Breaker] OPENING for 30s after 5 consecutive failures.")
        
    msg = f"{service} temporarily unavailable. Please try again soon."
    if _cb_failure_count >= 5:
        msg = f"{service} service is under heavy load (circuit breaker active). Please try again in 30 seconds."
        
    raise requests.RequestException(msg)


def sanitize_input(val: str) -> str:
    """Sanitizes alphanumeric strings and accessions, blocking script/SQL injection."""
    if not val: return ""
    return re.sub(r"[^a-zA-Z0-9_\-\.\:\*\s\>\+]", "", val).strip()


def validate_sequence(raw_text: str, seq_type: str = "dna") -> tuple[str, list[str], tuple]:
    """
    Validates and normalizes raw sequence text (plain or FASTA).
    Returns (cleaned_seq, warnings_list, error_response_tuple).
    """
    if len(raw_text) > 5 * 1024 * 1024:
        return "", [], error_response("Sequence payload exceeds 5MB limit.", 413)

    raw_text = raw_text.strip()
    if not raw_text:
        return "", [], error_response("Sequence is empty.", 400)

    lines = raw_text.splitlines()
    if lines and lines[0].startswith(">"):
        seq_lines = []
        for line in lines[1:]:
            line = line.strip()
            if line.startswith(">"):
                return "", [], error_response("Malformed FASTA or multiple sequences not supported.", 400)
            seq_lines.append(line)
        cleaned_seq = "".join(seq_lines).upper()
    else:
        if ">" in raw_text:
            return "", [], error_response("Malformed FASTA format.", 400)
        cleaned_seq = "".join(raw_text.split()).upper()

    warnings = []
    length = len(cleaned_seq)

    if length == 0:
        return "", [], error_response("Sequence is empty after removing whitespace.", 400)

    if length < 15:
        warnings.append(f"Short sequence ({length} bases/aa). Results may lack statistical significance.")

    if seq_type in ("dna", "rna"):
        valid_strict = set("ATGCNU")
        ambiguous = set("RYSWKMBDHV")
        allowed = valid_strict | ambiguous
        found = set(cleaned_seq)
        
        invalid = found - allowed
        if invalid:
            return "", [], error_response(
                f"Sequence contains invalid {seq_type.upper()} characters: {', '.join(sorted(invalid))}.", 422
            )
            
        found_ambig = found & set("RYSWKMBDHV")
        if found_ambig:
            warnings.append(f"Sequence contains ambiguous bases: {', '.join(sorted(found_ambig))}.")

    elif seq_type == "protein":
        allowed = set("ACDEFGHIKLMNPQRSTVWY*")
        found = set(cleaned_seq)
        invalid = found - allowed
        if invalid:
            return "", [], error_response(
                f"Sequence contains invalid protein characters: {', '.join(sorted(invalid))}.", 422
            )

    return cleaned_seq, warnings, None


def _cacheable_response(response) -> bool:
    try:
        if isinstance(response, tuple):
            status_code = response[1] if len(response) > 1 else 200
            data = response[0] if response else b""
        else:
            status_code = getattr(response, "status_code", 200)
            data = response.get_data() if hasattr(response, "get_data") else b""
        
        if status_code != 200:
            return False
        if not data:
            return False
        return True
    except Exception:
        return False


@app.before_request
def _security_gate():
    ip = request.remote_addr
    if not ip: return
    
    if _REDIS_AVAILABLE and redis_client:
        try:
            if redis_client.get(f"block_ip:{ip}"):
                return error_response("Your IP is temporarily blocked.", 403)
            
            key = f"ip_rate:{ip}:{int(time.time()/60)}"
            count = redis_client.incr(key)
            if count == 1:
                redis_client.expire(key, 60)
            if count > 30:
                return error_response("IP Rate Limit Exceeded.", 429)
        except Exception:
            pass

    if not request.path.startswith("/api/"):
        return
    g._btk_request_start = time.monotonic()


@app.after_request
def _security_monitor(response):
    ip = request.remote_addr
    if not ip or not _REDIS_AVAILABLE or not redis_client:
        return response

    if response.status_code >= 400 and response.status_code != 429:
        try:
            fail_key = f"fail_count:{ip}"
            fails = redis_client.incr(fail_key)
            if fails == 1:
                redis_client.expire(fail_key, 600)
            if fails >= 15:
                redis_client.setex(f"block_ip:{ip}", 600, "1")
        except Exception:
            pass
    return response


@app.route("/")
def serve_index():
    return send_from_directory(".", "index.html")


@app.route("/<path:path>")
def serve_static(path):
    return send_from_directory(".", path)

# NCBI Accession Map
ACCESSION_PREFIX_MAP = {
    "NM_": ("nucleotide", "mRNA"), "NR_": ("nucleotide", "non-coding RNA"),
    "XM_": ("nucleotide", "predicted mRNA"), "XR_": ("nucleotide", "predicted non-coding RNA"),
    "NG_": ("nucleotide", "genomic region"), "NC_": ("nucleotide", "chromosome"),
    "NT_": ("nucleotide", "genomic contig"), "NW_": ("nucleotide", "genomic scaffold"),
    "AC_": ("nucleotide", "alternate chromosome"),
    "NP_": ("protein", "protein"), "XP_": ("protein", "predicted protein"),
    "WP_": ("protein", "non-redundant protein"), "YP_": ("protein", "protein"),
}

@app.route("/api/fetch/ncbi")
@cache.cached(timeout=21600, query_string=True, response_filter=_cacheable_response)
def fetch_ncbi():
    acc_id = sanitize_input(request.args.get("id", ""))
    fmt    = sanitize_input(request.args.get("format", "fasta")).lower()
    requested_db = sanitize_input(request.args.get("db", "")).lower() or None

    if not acc_id:
        return error_response("Missing 'id' parameter.", 400)

    db_name = requested_db or "nucleotide"
    rettype = "fasta" if fmt == "fasta" else "gb"

    enforce_rate_limit("ncbi")
    params = {
        "db": db_name, "id": acc_id, "rettype": rettype, "retmode": "text",
        "tool": "BioToolkit", "email": "admin@biotoolkit.dev",
    }
    user_key = request.headers.get("X-NCBI-API-Key", "").strip()
    if user_key or NCBI_API_KEY:
        params["api_key"] = user_key or NCBI_API_KEY

    try:
        resp = proxy_get(f"{NCBI_BASE}/efetch.fcgi", params=params)
        if resp.status_code == 200:
            return Response(resp.text.strip(), content_type="text/plain")
        return error_response("NCBI retrieval failed.", resp.status_code)
    except Exception as e:
        return error_response(str(e), 502)

@app.route("/api/fetch/ncbi/search")
@cache.cached(timeout=86400, query_string=True, response_filter=_cacheable_response)
def search_ncbi():
    symbol = sanitize_input(request.args.get("symbol", ""))
    species = sanitize_input(request.args.get("species", "Homo sapiens"))
    if not symbol: return error_response("Missing 'symbol'.", 400)

    enforce_rate_limit("ncbi")
    try:
        s_resp = proxy_get(f"{NCBI_BASE}/esearch.fcgi", params={
            "db": "gene", "term": f"{symbol}[Gene] AND {species}[Organism]",
            "retmax": 1, "retmode": "json"
        })
        ids = s_resp.json().get("esearchresult", {}).get("idlist", [])
        if not ids: return error_response("Gene not found.", 404)
        
        gene_id = ids[0]
        sum_resp = proxy_get(f"{NCBI_BASE}/esummary.fcgi", params={"db": "gene", "id": gene_id, "retmode": "json"})
        gene_info = sum_resp.json().get("result", {}).get(gene_id, {})
        
        return jsonify({
            "gene_symbol": gene_info.get("name", symbol),
            "description": gene_info.get("description", ""),
            "gene_id": gene_id,
            "organism": species
        })
    except Exception as e:
        return error_response(str(e), 502)

@app.route("/api/fetch/ensembl")
@cache.cached(timeout=86400, query_string=True, response_filter=_cacheable_response)
def fetch_ensembl():
    raw_id = sanitize_input(request.args.get("id", ""))
    species = sanitize_input(request.args.get("species", "homo_sapiens")).lower()
    if not raw_id: return error_response("Missing 'id'.", 400)

    enforce_rate_limit("ensembl")
    try:
        url = f"{ENSEMBL_BASE}/lookup/id/{raw_id}?expand=1;content-type=application/json"
        resp = proxy_get(url)
        if resp.status_code == 200:
            return jsonify(resp.json())
        return error_response("Ensembl lookup failed.", resp.status_code)
    except Exception as e:
        return error_response(str(e), 502)

@app.route("/api/fetch/uniprot")
def fetch_uniprot():
    acc = request.args.get("id", "").strip()
    if not acc: return error_response("Missing 'id'.", 400)
    enforce_rate_limit("uniprot")
    try:
        resp = proxy_get(f"{UNIPROT_BASE}/uniprotkb/{acc}.fasta")
        if resp.status_code == 200:
            return Response(resp.text, content_type="text/plain")
        return error_response("UniProt fetch failed.", resp.status_code)
    except Exception as e:
        return error_response(str(e), 502)

@app.route("/api/fetch/kegg")
def fetch_kegg():
    kegg_id = request.args.get("id", "").strip()
    if not kegg_id: return error_response("Missing 'id'.", 400)
    enforce_rate_limit("kegg")
    try:
        resp = proxy_get(f"{KEGG_BASE}/get/{kegg_id}/ntseq")
        if resp.status_code == 200:
            return Response(resp.text, content_type="text/plain")
        return error_response("KEGG fetch failed.", resp.status_code)
    except Exception as e:
        return error_response(str(e), 502)

@app.route("/api/health")
def health():
    return jsonify({"status": "ok", "version": "1.1.0"})

@app.route("/api/jobs/create", methods=["POST"])
def create_job():
    from worker import process_sequence_task
    data = request.get_json()
    if not data or "tool_type" not in data:
        return error_response("Missing 'tool_type'.", 400)
    
    task = process_sequence_task.delay(data["tool_type"], data.get("params", {}))
    return jsonify({"status": "accepted", "task_id": task.id}), 202

@app.route("/api/jobs/status/<task_id>")
def get_task_status(task_id):
    from worker import celery
    task = celery.AsyncResult(task_id)
    response = {"state": task.state, "task_id": task_id}
    if task.state == 'SUCCESS':
        response["result"] = task.result
    return jsonify(response)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5000)))
