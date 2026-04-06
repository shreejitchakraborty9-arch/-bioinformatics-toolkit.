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
    _BIOPYTHON_AVAILABLE = Falsetry:
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
# Pre-flight check
with app.app_context():
    try:
        db.create_all()
    except Exception as e:
        logger.warning("DB deferred.")

# 1. CORS Enforcement (Restrict to exact frontend domain)
ALLOWED_ORIGIN = os.environ.get("ALLOWED_ORIGIN", "http://localhost:5173") # Default to local dev
CORS(app, resources={r"/api/*": {"origins": ALLOWED_ORIGIN}})

# 2. Payload size Limit (5MB) - Prevents OOM and large sequence abuse
app.config['MAX_CONTENT_LENGTH'] = 5 * 1024 * 1024logging.basicConfig(level=logging.INFO)
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
        payload["details"] = detailsreturn jsonify(payload), status_codedef enforce_rate_limit(db: str):
    """
    Redis-backed server-side per-database rate limiter for distributed workers.
    Gracefully degrades to in-memory if Redis is unavailable.
    """
    min_interval = RATE_LIMITS.get(db, 0.5)
    
    # Specific adjustment for NCBI with API KEYif db == "ncbi" and NCBI_API_KEY:
        min_interval = 0.1  # 10 req/secif redis_client:
        try:
            key = f"rate_limit:{db}"
            # Atomic Lua script for rate limit calculationlua_script = """
            local last = redis.call('GET', KEYS[1])
            local current_time = tonumber(ARGV[1])
            local min_interval = tonumber(ARGV[2])
            
            if not last thenredis.call('SET', KEYS[1], current_time)
                return 0endlocal elapsed = current_time - tonumber(last)
            if elapsed < min_interval thenlocal wait_time = min_interval - elapsedredis.call('SET', KEYS[1], current_time + wait_time)
                return wait_timeelseredis.call('SET', KEYS[1], current_time)
                return 0end
            """
            
            wait_time = redis_client.eval(lua_script, 1, key, time.time(), min_interval)
            
            if wait_time and float(wait_time) > 0:
                time.sleep(float(wait_time))
            returnexcept Exception as e:
            logger.error("Redis rate limit failed for %s, falling back to memory: %s", db, e)
            
    # In-memory fallbacklast = _last_request_time.get(db, 0)
    elapsed = time.time() - lastif elapsed < min_interval:
        time.sleep(min_interval - elapsed)
        _last_request_time[db] = time.time() + (min_interval - elapsed)
    else:
        _last_request_time[db] = time.time()


# ─────────────────────────────────────────────────────────────
# Auth & Session Preparation (Placeholders)
# ─────────────────────────────────────────────────────────────
@login_manager.user_loader
def load_user(user_id):
    # Prepared for future DB-backed authreturn User.query.get(int(user_id)) if user_id else Nonedef optional_login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        # Currently allows all users (Guest mode enabled)
        # Flip this to @login_required later to enforce authenticationreturn f(*args, **kwargs)
    return decorated_functiondef proxy_get(url: str, headers: dict = None, params: dict = None) -> requests.Response:
    """
    Make a proxied GET request with a consistent User-Agent,
    retry logic (exponential backoff), and a circuit breaker.
    """
    global _cb_failure_count, _cb_open_until
    
    # 1. Determine service for reportingservice = "NCBI" if "ncbi" in url else "Ensembl" if "ensembl" in url else "UniProt" if "uniprot" in url else "Biological service"
    
    # 2. Check Circuit Breakerif time.time() < _cb_open_until:
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
    total_attempts = len(backoffs) + 1for i in range(total_attempts):
        try:
            resp = requests.get(url, headers=h, params=params, timeout=REQUEST_TIMEOUT)
            
            # success (200) or client errors (400/404) that aren't service-level failuresif resp.status_code == 200 or (resp.status_code < 500 and resp.status_code not in (429,)):
                _cb_failure_count = 0 # Reset on server-reachable responsereturn resp
            
            # retry on 429 (rate limit) or 500+ (server-side issues)
            if i < len(backoffs):
                logger.warning("[Retry %d] %s returned HTTP %d. Backing off %ss.", i + 1, service, resp.status_code, backoffs[i])
                time.sleep(backoffs[i])
                continueelse:
                # Max retries reachedbreakexcept (requests.Timeout, requests.RequestException) as e:
            if i < len(backoffs):
                logger.warning("[Retry %d] Connection to %s failed: %s. Backing off %ss.", i + 1, service, e, backoffs[i])
                time.sleep(backoffs[i])
                continuebreak
            
    # If we reached here, it's a persistent failure
    _cb_failure_count += 1if _cb_failure_count >= 5:
        _cb_open_until = time.time() + 30logger.error("[Circuit Breaker] OPENING for 30s after 5 consecutive failures.")
        
    # User-friendly message constructionmsg = f"{service} temporarily unavailable. Please try again soon."
    if _cb_failure_count >= 5:
        msg = f"{service} service is under heavy load (circuit breaker active). Please try again in 30 seconds."
        
    raise requests.RequestException(msg)


def sanitize_input(val: str) -> str:
    """Sanitizes alphanumeric strings and accessions, blocking script/SQL injection."""
    if not val: return ""
    # Only allow safe biological/standard charactersreturn re.sub(r"[^a-zA-Z0-9_\-\.\:\*\s\>\+]", "", val).strip()


def validate_sequence(raw_text: str, seq_type: str = "dna") -> tuple[str, list[str], tuple]:
    """
    Validates and normalizes raw sequence text (plain or FASTA).
    seq_type: "dna", "rna", or "protein".
    Returns (cleaned_seq, warnings_list, error_response_tuple).
    If error_response_tuple is not None, validation completely failed.
    """
    if len(raw_text) > 5 * 1024 * 1024:
        return "", [], error_response("Sequence payload exceeds 5MB limit.", 413)

    raw_text = raw_text.strip()
    if not raw_text:
        return "", [], error_response("Sequence is empty.", 400)

    # Simple FASTA handlinglines = raw_text.splitlines()
    if lines[0].startswith(">"):
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
        # Match real bioinformatics tool standardsvalid_strict = set("ATGCNU")
        ambiguous = set("RYSWKMBDHV")
        allowed = valid_strict | ambiguousfound = set(cleaned_seq)
        
        invalid = found - allowedif invalid:
            return "", [], error_response(
                f"Sequence contains invalid {seq_type.upper()} characters: {', '.join(sorted(invalid))}. "
                "Only standard and IUPAC nucleotide bases are permitted.",
                422
            )
            
        found_ambig = found & set("RYSWKMBDHV")
        if found_ambig:
            warnings.append(f"Sequence contains ambiguous bases: {', '.join(sorted(found_ambig))}.")

    elif seq_type == "protein":
        allowed = set("ACDEFGHIKLMNPQRSTVWY*")
        found = set(cleaned_seq)
        invalid = found - allowedif invalid:
            return "", [], error_response(
                f"Sequence contains invalid protein characters: {', '.join(sorted(invalid))}. "
                "Only standard 20 amino acids and '*' are permitted.",
                422
            )

    return cleaned_seq, warnings, None


# ─────────────────────────────────────────────────────────────
# Cache Utilities
# ─────────────────────────────────────────────────────────────
def _cacheable_response(response) -> bool:
    """
    Flask-Caching response_filter — only cache successful, non-empty,
    non-error responses. Prevents stale errors from being served.
    Handles both Flask Response objects and raw (body, status) tuples.
    """
    try:
        if isinstance(response, tuple):
            status_code = response[1] if len(response) > 1 else 200data = response[0] if response else b""
            if isinstance(data, str):
                data = data.encode()
            content_type = ""
            if len(response) > 2 and isinstance(response[2], dict):
                content_type = response[2].get("Content-Type", "")
        else:
            status_code = getattr(response, "status_code", 200)
            data = response.get_data() if hasattr(response, "get_data") else b""
            content_type = getattr(response, "content_type", "") or ""

        # Must be HTTP 200if status_code != 200:
            return False
        # Must have a bodyif not data:
            return False
        # Must not be a JSON error payload  {"error": "..."}
        if "json" in content_type:
            try:
                parsed = json.loads(data)
                if isinstance(parsed, dict) and "error" in parsed:
                    return Falseexcept (ValueError, TypeError):
                passreturn Trueexcept Exception:
        return False


@app.before_request
def _security_gate():
    """Centralized security hardening: IP rate limiting and block/ban checking."""
    # 3. IP Rate Limiting (30 requests/min per IP)
    ip = request.remote_addrif not ip: return
    
    # Check if IP is explicitly blocked (Abuse/Multiple failures)
    if _REDIS_AVAILABLE and redis_client:
        try:
            if redis_client.get(f"block_ip:{ip}"):
                logger.warning("[Security] Blocked request from banned IP: %s", ip)
                return error_response("Your IP is temporarily blocked due to repeated failed requests (10m ban). Please check your queries and try again later.", 403)
            
            # Minute-bucket rate limit (30 per IP per minute)
            key = f"ip_rate:{ip}:{int(time.time()/60)}"
            count = redis_client.incr(key)
            if count == 1:
                redis_client.expire(key, 60)
            if count > 30:
                logger.warning("[Security] IP Rate limit reached: %s", ip)
                return error_response("IP Rate Limit Exceeded (30 req/min). Public API exposure is protected.", 429)
        except Exception as e:
            logger.debug("Redis security gate failure: %s", e)

    # ─────────────────────────────────────────────────────────────
    # Cache Logs & Tracking 
    # ─────────────────────────────────────────────────────────────
    if not request.path.startswith("/api/"):
        returng._btk_request_start = time.monotonic()
    try:
        g._btk_cache_had_entry = cache.get(request.url) is not Noneexcept Exception:
        g._btk_cache_had_entry = False


@app.after_request
def _security_monitor(response):
    """Monitors failures to automatically block abusive IPs."""
    ip = request.remote_addrif not ip or not _REDIS_AVAILABLE or not redis_client:
        return response

    # 4. Prevent Abuse: Block if > 15 failures in 10 minutesif response.status_code >= 400 and response.status_code != 429:
        try:
            fail_key = f"fail_count:{ip}"
            fails = redis_client.incr(fail_key)
            if fails == 1:
                redis_client.expire(fail_key, 600)  # 10m windowif fails >= 15:
                logger.error("[Security] BANNING IP %s for 10 minutes (Too many failures).", ip)
                redis_client.setex(f"block_ip:{ip}", 600, "1")
        except Exception:
            passif not request.path.startswith("/api/"):
        return responsetry:
        elapsed_ms = (time.monotonic() - getattr(g, "_btk_request_start", time.monotonic())) * 1000was_cached = getattr(g, "_btk_cache_had_entry", False)
        status = "HIT " if was_cached else "MISS"
        logger.info(
            "[Cache %s] %s %s \u2192 HTTP %d (%.1f ms)",
            status, request.method, request.full_path.rstrip("?"), response.status_code, elapsed_ms,
        )
    except Exception:
        passreturn response


# ─────────────────────────────────────────────────────────────
# Static Frontend Serving
# ─────────────────────────────────────────────────────────────
@app.route("/")
def serve_index():
    return send_from_directory(".", "index.html")


@app.route("/<path:path>")
def serve_static(path):
    """Serve css/, js/, lib/ directories and any top-level static files."""
    return send_from_directory(".", path)


# ─────────────────────────────────────────────────────────────
# API – NCBI Proxy
# ─────────────────────────────────────────────────────────────

# Accession prefix → (database, molecule_type) mapping
# Covers RefSeq curated + INSDC accession namespaces
ACCESSION_PREFIX_MAP = {
    # mRNA / transcript
    "NM_": ("nucleotide", "mRNA"),
    "NR_": ("nucleotide", "non-coding RNA"),
    "XM_": ("nucleotide", "predicted mRNA"),
    "XR_": ("nucleotide", "predicted non-coding RNA"),
    # Genomic
    "NG_": ("nucleotide", "genomic region"),
    "NC_": ("nucleotide", "chromosome"),
    "NT_": ("nucleotide", "genomic contig"),
    "NW_": ("nucleotide", "genomic scaffold"),
    "AC_": ("nucleotide", "alternate chromosome"),
    # Protein
    "NP_": ("protein", "protein"),
    "XP_": ("protein", "predicted protein"),
    "WP_": ("protein", "non-redundant protein"),
    "YP_": ("protein", "protein"),
    "AP_": ("protein", "protein"),
    # INSDC (legacy GenBank/EMBL/DDBJ) — alphanumeric prefixes
    "U":   ("nucleotide", "sequence"),
    "AF":  ("nucleotide", "sequence"),
    "AY":  ("nucleotide", "sequence"),
    "DQ":  ("nucleotide", "sequence"),
    "EF":  ("nucleotide", "sequence"),
    "EU":  ("nucleotide", "sequence"),
    "FJ":  ("nucleotide", "sequence"),
    "GU":  ("nucleotide", "sequence"),
    "HM":  ("nucleotide", "sequence"),
    "JN":  ("nucleotide", "sequence"),
    "JQ":  ("nucleotide", "sequence"),
    "KC":  ("nucleotide", "sequence"),
    "KF":  ("nucleotide", "sequence"),
}

# GenBank format: allowed rettype tokens per db
_GENBANK_RETTYPE = {
    "nucleotide": "gb",
    "protein": "gp",
}

def _infer_db_from_accession(acc: str):
    """
    Infer the NCBI database and molecule type from an accession prefix.
    Returns (db, mol_type) or (None, None) if prefix is unrecognised.
    """
    upper = acc.upper()
    # Try longest match first (e.g. NM_ before N)
    for prefix, (db, mol) in sorted(ACCESSION_PREFIX_MAP.items(), key=lambda x: -len(x[0])):
        if upper.startswith(prefix.upper()):
            return db, molreturn None, Nonedef _validate_accession(acc: str, requested_db: str):
    """
    Validate the accession against the requested database.
    Returns (is_valid, inferred_db, mol_type, error_message).
    """
    inferred_db, mol_type = _infer_db_from_accession(acc)
    if inferred_db is None:
        # Unknown prefix — allow it through with a warning; NCBI will reject if truly invalidreturn True, requested_db, "sequence", Noneif requested_db != inferred_db:
        prefix = acc.split("_")[0] + ("_" if "_" in acc else "")
        return (
            False,
            inferred_db,
            mol_type,
            (
                f"Accession '{acc}' ({prefix}* = {mol_type}) belongs to the "
                f"'{inferred_db}' database, but '{requested_db}' was requested. "
                f"Switch to db={inferred_db} or use a matching accession."
            ),
        )
    return True, inferred_db, mol_type, None


@app.route("/api/fetch/ncbi")
@cache.cached(timeout=21600, query_string=True, response_filter=_cacheable_response)
def fetch_ncbi():
    """
    Proxy for NCBI E-utilities eFetch.
    Supports dual scientific output formats.

    Query params:
      - id      : Accession (e.g. NM_007294, NP_009225, NC_000017)
      - db      : 'nucleotide' | 'protein'  (auto-inferred from accession if omitted)
      - format  : 'fasta' | 'genbank'  (default: fasta)
                  'genbank' returns the full annotated GenBank / GenPept recordincluding FEATURES table (CDS, exon, gene, etc.)

    Returns:
      text/plain — raw FASTA or GenBank/GenPept recordapplication/json — on error, with scientific context
    """
    acc_id = sanitize_input(request.args.get("id", ""))
    fmt    = sanitize_input(request.args.get("format", "fasta")).lower()
    # Legacy 'rettype' param still honoured for backwards compatif "rettype" in request.args and "format" not in request.args:
        fmt = "genbank" if sanitize_input(request.args["rettype"]) == "gb" else "fasta"

    requested_db = sanitize_input(request.args.get("db", "")).lower() or None

    # ── 1. Presence check ─────────────────────────────────────
    if not acc_id:
        return error_response("Missing 'id' parameter. Provide an NCBI accession, e.g. NM_007294.", 400)

    # ── 2. Format check ───────────────────────────────────────
    if fmt not in ("fasta", "genbank"):
        return error_response(
            f"Invalid format '{fmt}'. Use 'fasta' (raw sequence) or "
            "'genbank' (annotated record with CDS/exon features).",
            400,
        )

    # ── 3. Accession validation & db inference ─────────────────
    inferred_db, mol_type = _infer_db_from_accession(acc_id)
    db = requested_db or inferred_db or "nucleotide"

    if requested_db and inferred_db:
        valid, _, _, val_err = _validate_accession(acc_id, requested_db)
        if not valid:
            return error_response(val_err, 400)

    if db not in ("nucleotide", "protein"):
        return error_response("Invalid 'db'. Use 'nucleotide' or 'protein'.", 400)

    # ── 4. Map format → NCBI rettype ──────────────────────────
    if fmt == "genbank":
        rettype = _GENBANK_RETTYPE.get(db, "gb")
    else:
        rettype = "fasta"

    # ── 5. Fetch from NCBI ────────────────────────────────────
    enforce_rate_limit("ncbi")

    params = {
        "db":      db,
        "id":      acc_id,
        "rettype": rettype,
        "retmode": "text",
        "tool":    "BioToolkit",
        "email":   "admin@biotoolkit.dev",
    }
    user_key = request.headers.get("X-NCBI-API-Key", "").strip()
    effective_key = user_key or NCBI_API_KEYif effective_key:
        params["api_key"] = effective_keytry:
        resp = proxy_get(f"{NCBI_BASE}/efetch.fcgi", params=params)

        if resp.status_code == 200:
            text = resp.text.strip()

            # Catch NCBI HTML error pages masquerading as 200if text.startswith("<!DOCTYPE") or text.startswith("<html"):
                return error_response(
                    f"NCBI could not retrieve '{acc_id}'. The accession may be "
                    "invalid, suppressed, or not available in the requested database.",
                    404,
                )

            # Catch NCBI "not found" text responsesif "Error" in text[:30] or "Nothing" in text[:40]:
                return error_response(
                    f"NCBI returned no data for accession '{acc_id}' in database '{db}'.",
                    404,
                )

            content_type = "text/plain"
            return Response(
                text,
                content_type=content_type,
                headers={
                    "X-BioToolkit-Format":    fmt,
                    "X-BioToolkit-DB":        db,
                    "X-BioToolkit-MolType":   mol_type or "unknown",
                    "X-BioToolkit-Accession": acc_id,
                },
            )

        elif resp.status_code == 400:
            return error_response(
                f"NCBI rejected the request for '{acc_id}'. "
                "Check that the accession is correct and the database matches the accession prefix.",
                400,
            )
        elif resp.status_code == 404:
            return error_response(
                f"Accession '{acc_id}' not found in NCBI '{db}'. "
                "Verify the accession exists and has not been suppressed or retracted.",
                404,
            )
        elif resp.status_code == 429:
            return error_response(
                "NCBI rate limit exceeded. The server is throttling requests. "
                "Please wait a few seconds and try again.",
                429,
            )
        else:
            return error_response(f"NCBI returned an unexpected HTTP {resp.status_code}.", resp.status_code)

    except requests.RequestException as e:
        logger.error("NCBI fetch error: %s", e)
        return error_response(str(e), 502)


@app.route("/api/fetch/ncbi/search")
@cache.cached(timeout=86400, query_string=True, response_filter=_cacheable_response)
def search_ncbi():
    """
    Scientific-grade gene resolver.
    Query params:
      - symbol: Gene symbol (e.g. BRCA1, TP53)
      - species: Species name (default: Homo sapiens)
    """
    symbol = sanitize_input(request.args.get("symbol", ""))
    species = sanitize_input(request.args.get("species", "Homo sapiens"))

    if not symbol:
        return error_response("Missing 'symbol' parameter.", 400)

    # --- 1. ESEARCH to get gene ID ---
    enforce_rate_limit("ncbi")
    search_term = f"{symbol}[Gene] AND {species}[Organism]"
    search_params = {
        "db": "gene",
        "term": search_term,
        "retmax": 1,
        "retmode": "json",
        "tool": "BioToolkit",
        "email": "admin@biotoolkit.dev"
    }
    user_key = request.headers.get("X-NCBI-API-Key", "").strip()
    effective_key = user_key or NCBI_API_KEYif effective_key:
        search_params["api_key"] = effective_keytry:
        s_resp = proxy_get(f"{NCBI_BASE}/esearch.fcgi", params=search_params)
        if s_resp.status_code != 200:
            return error_response(f"NCBI esearch returned HTTP {s_resp.status_code}.", s_resp.status_code)
            
        s_data = s_resp.json()
        ids = s_data.get("esearchresult", {}).get("idlist", [])
        
        if not ids:
            return error_response(f"No results found for gene symbol '{symbol}' in '{species}'.", 404)
            
        gene_id = str(ids[0])

        # --- 2. ESUMMARY to get metadata ---
        enforce_rate_limit("ncbi")
        sum_params = {
            "db": "gene",
            "id": gene_id,
            "retmode": "json",
            "tool": "BioToolkit",
            "email": "admin@biotoolkit.dev"
        }
        user_key = request.headers.get("X-NCBI-API-Key", "").strip()
        effective_key = user_key or NCBI_API_KEYif effective_key:
            sum_params["api_key"] = effective_keysum_resp = proxy_get(f"{NCBI_BASE}/esummary.fcgi", params=sum_params)
        sum_data = sum_resp.json()
        
        gene_info = sum_data.get("result", {}).get(gene_id, {})

        # --- 3. ELINK to get linked nucleotide & protein accessions ---
        enforce_rate_limit("ncbi")
        link_params = {
            "dbfrom": "gene",
            "db": "nuccore,protein",
            "id": gene_id,
            "retmode": "json",
            "tool": "BioToolkit",
            "email": "admin@biotoolkit.dev"
        }
        user_key = request.headers.get("X-NCBI-API-Key", "").strip()
        effective_key = user_key or NCBI_API_KEYif effective_key:
            link_params["api_key"] = effective_keylink_resp = proxy_get(f"{NCBI_BASE}/elink.fcgi", params=link_params)
        link_data = link_resp.json()
        
        nucleotide_accessions = []
        protein_accessions = []
        
        linksets = link_data.get("linksets", [])
        if linksets:
            linksetdbs = linksets[0].get("linksetdbs", [])
            for db_links in linksetdbs:
                db_to = db_links.get("dbto", "").lower()
                links_data = db_links.get("links", [])
                if isinstance(links_data, dict):
                    links_data = [links_data]
                
                links = [str(link.get("id")) for link in links_data if isinstance(link, dict) and link.get("id")]
                if db_to == "nuccore":
                    nucleotide_accessions.extend(links)
                elif db_to == "protein":
                    protein_accessions.extend(links)

        # Structure the final scientific-grade JSONreturn jsonify({
            "gene_symbol": gene_info.get("name", symbol),
            "organism": gene_info.get("organism", {}).get("scientificname", species),
            "description": gene_info.get("description", ""),
            "chromosome": gene_info.get("chromosome", ""),
            "map_location": gene_info.get("maplocation", ""),
            "gene_id": gene_id,
            "summary": gene_info.get("summary", ""),
            "nucleotide_accessions": nucleotide_accessions,
            "protein_accessions": protein_accessions
        })

    except requests.RequestException as e:
        logger.error("NCBI search error: %s", e)
        return error_response(str(e), 502)


# ─────────────────────────────────────────────────────────────
# API – Ensembl Proxy
# ─────────────────────────────────────────────────────────────
@app.route("/api/fetch/ensembl")
@cache.cached(timeout=86400, query_string=True, response_filter=_cacheable_response)
def fetch_ensembl():
    """
    Full gene information endpoint backed by Ensembl REST API.
    Supports gene symbol or Ensembl stable ID as input.

    Query params:
      - id               : Gene symbol (e.g. BRCA2) OR Ensembl stable ID (e.g. ENSG00000139618)
      - species          : Species slug for symbol lookup (default: homo_sapiens)
      - include_sequence : 'true' to additionally fetch the genomic FASTA sequence (default: false)

    Returns structured JSON:
      {
        "ensembl_id", "gene_name", "biotype", "strand",
        "chromosome", "start", "end",
        "transcripts": [{ "id", "biotype", "exons": [...] }],
        "sequence": "ATGC..."   ← only when include_sequence=true
      }
    """
    raw_id           = sanitize_input(request.args.get("id", ""))
    species          = sanitize_input(request.args.get("species", "homo_sapiens")).lower()
    include_seq      = request.args.get("include_sequence", "false").lower() == "true"

    if not raw_id:
        return error_response(
            "Missing 'id' parameter. Provide a gene symbol (e.g. BRCA2) "
            "or an Ensembl stable ID (e.g. ENSG00000139618).",
            400,
        )

    json_headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # ── 1. Symbol → Ensembl ID resolution ─────────────────────
    ens_id = raw_idif not raw_id.upper().startswith("ENS"):
        enforce_rate_limit("ensembl")
        try:
            xref_resp = proxy_get(
                f"{ENSEMBL_BASE}/xrefs/symbol/{species}/{raw_id}",
                headers=json_headers,
                params={"content-type": "application/json"},
            )
            if xref_resp.status_code != 200 or not xref_resp.json():
                return error_response(
                    f"Gene symbol '{raw_id}' not found in Ensembl for species '{species}'. "
                    "Check the spelling or try a direct Ensembl stable ID (ENSG...).",
                    404,
                )
            # Pick the first gene hit (object_type == Gene preferred)
            hits = xref_resp.json()
            gene_hits = [h for h in hits if h.get("type", "").lower() == "gene"]
            ens_id = (gene_hits or hits)[0]["id"]
        except requests.RequestException as e:
            logger.error("Ensembl xrefs lookup error: %s", e)
            return error_response(str(e), 502)

    # ── 2. /lookup/id?expand=1 — full gene + transcript metadata ──
    enforce_rate_limit("ensembl")
    try:
        lookup_resp = proxy_get(
            f"{ENSEMBL_BASE}/lookup/id/{ens_id}",
            headers=json_headers,
            params={"expand": 1, "content-type": "application/json"},
        )

        if lookup_resp.status_code == 404:
            return error_response(
                f"Ensembl ID '{ens_id}' not found. "
                "The ID may be retired or belong to a different assembly.",
                404,
            )
        if lookup_resp.status_code == 429:
            return error_response(
                "Ensembl rate limit exceeded. Please wait a few seconds and retry.",
                429,
            )
        if lookup_resp.status_code != 200:
            return error_response(
                f"Ensembl /lookup returned HTTP {lookup_resp.status_code}.",
                lookup_resp.status_code,
            )

        g = lookup_resp.json()

        # ── 3. Build transcript + exon structure ──────────────
        transcripts_out = []
        for tx in g.get("Transcript", []):
            exons_out = [
                {
                    "id":     ex.get("id"),
                    "start":  ex.get("start"),
                    "end":    ex.get("end"),
                    "strand": ex.get("strand"),
                    "rank":   ex.get("rank"),
                }
                for ex in tx.get("Exon", [])
            ]
            transcripts_out.append({
                "id":           tx.get("id"),
                "display_name": tx.get("display_name"),
                "biotype":      tx.get("biotype"),
                "start":        tx.get("start"),
                "end":          tx.get("end"),
                "strand":       tx.get("strand"),
                "is_canonical": bool(tx.get("is_canonical")),
                "exon_count":   len(exons_out),
                "exons":        exons_out,
            })

        # Sort: canonical first, then by starttranscripts_out.sort(key=lambda t: (not t["is_canonical"], t.get("start") or 0))

        result = {
            "ensembl_id":       ens_id,
            "gene_name":        g.get("display_name") or g.get("external_name") or raw_id,
            "description":      g.get("description", ""),
            "biotype":          g.get("biotype", ""),
            "strand":           g.get("strand"),
            "chromosome":       g.get("seq_region_name", ""),
            "start":            g.get("start"),
            "end":              g.get("end"),
            "assembly_name":    g.get("assembly_name", ""),
            "species":          g.get("species", species),
            "transcript_count": len(transcripts_out),
            "transcripts":      transcripts_out,
        }

        # ── 4. Optional: fetch genomic sequence ───────────────
        if include_seq:
            enforce_rate_limit("ensembl")
            try:
                seq_resp = proxy_get(
                    f"{ENSEMBL_BASE}/sequence/id/{ens_id}",
                    params={"content-type": "text/plain", "type": "genomic"},
                )
                if seq_resp.status_code == 200:
                    result["sequence"] = seq_resp.text.strip()
                else:
                    result["sequence_error"] = (
                        f"Sequence fetch returned HTTP {seq_resp.status_code}."
                    )
            except requests.RequestException as e:
                logger.warning("Ensembl sequence fetch failed: %s", e)
                result["sequence_error"] = str(e)

        return jsonify(result)

    except requests.Timeout:
        return error_response(
            "Ensembl request timed out. The REST API may be under load — try again.",
            504,
        )
    except requests.RequestException as e:
        logger.error("Ensembl fetch error: %s", e)
        return error_response("Failed to reach Ensembl REST API.", 502)



# ─────────────────────────────────────────────────────────────
# API – UniProt Proxy
# ─────────────────────────────────────────────────────────────
@app.route("/api/fetch/uniprot")
def fetch_uniprot():
    """
    Proxy for UniProt REST API.
    Query params:
      - id: UniProt accession (e.g. P04637 for TP53)
    """
    acc = request.args.get("id", "").strip()
    if not acc:
        return error_response("Missing 'id' parameter.", 400)

    enforce_rate_limit("uniprot")

    try:
        resp = proxy_get(f"{UNIPROT_BASE}/uniprotkb/{acc}.fasta")
        if resp.status_code == 200:
            return Response(resp.text, content_type="text/plain")
        elif resp.status_code == 404:
            return error_response(f"UniProt accession '{acc}' not found.", 404)
        else:
            return error_response(f"UniProt returned HTTP {resp.status_code}.", resp.status_code)
    except requests.Timeout:
        return error_response("UniProt request timed out.", 504)
    except requests.RequestException as e:
        logger.error("UniProt fetch error: %s", e)
        return error_response("Failed to reach UniProt.", 502)


# ─────────────────────────────────────────────────────────────
# API – KEGG Proxy
# ─────────────────────────────────────────────────────────────
@app.route("/api/fetch/kegg")
def fetch_kegg():
    """
    Proxy for KEGG REST API.
    Query params:
      - id: KEGG entry ID (e.g. hsa:7157 for human TP53)
      - type: 'ntseq' or 'aaseq' (default: ntseq)
    """
    kegg_id = request.args.get("id", "").strip()
    seq_type = request.args.get("type", "ntseq")

    if not kegg_id:
        return error_response("Missing 'id' parameter.", 400)
    if seq_type not in ("ntseq", "aaseq"):
        return error_response("Invalid 'type'. Use 'ntseq' or 'aaseq'.", 400)

    enforce_rate_limit("kegg")

    try:
        resp = proxy_get(f"{KEGG_BASE}/get/{kegg_id}/{seq_type}")
        if resp.status_code == 200:
            text = resp.text.strip()
            if not text or text.startswith("<!"):
                return error_response(f"KEGG ID '{kegg_id}' returned no sequence.", 404)
            return Response(text, content_type="text/plain")
        elif resp.status_code == 404:
            return error_response(f"KEGG ID '{kegg_id}' not found.", 404)
        else:
            return error_response(f"KEGG returned HTTP {resp.status_code}.", resp.status_code)
    except requests.Timeout:
        return error_response("KEGG request timed out.", 504)
    except requests.RequestException as e:
        logger.error("KEGG fetch error: %s", e)
        return error_response("Failed to reach KEGG.", 502)


# ─────────────────────────────────────────────────────────────
# API – Ensembl CDS → Protein Translation
# ─────────────────────────────────────────────────────────────
@app.route("/api/ensembl/translate")
@cache.cached(timeout=86400, query_string=True, response_filter=_cacheable_response)
def ensembl_translate():
    """
    Fetches the CDS for an Ensembl gene or transcript ID and performsserver-side translation using Biopython (standard codon table 1).

    Query params:
      - id : Ensembl gene ID (ENSG...) or transcript ID (ENST...)

    Returns:
      {
        "ensembl_id"      : "ENST...",
        "cds_sequence"    : "ATGAAA...",
        "protein_sequence": "MK...",
        "cds_length"      : 1234,
        "protein_length"  : 411,
        "stop_codons"     : [{ "codon": "TGA", "position": 411 }],
        "has_start_codon" : true,
        "has_stop_codon"  : true
      }
    """
    if not _BIOPYTHON_AVAILABLE:
        return error_response(
            "Biopython is not installed on this server. "
            "Add 'biopython>=1.83' to requirements.txt and redeploy.",
            501,
        )

    ens_id = request.args.get("id", "").strip()
    if not ens_id:
        return error_response(
            "Missing 'id' parameter. Provide an Ensembl gene (ENSG...) "
            "or transcript (ENST...) ID.",
            400,
        )

    upper_id = ens_id.upper()
    if not (upper_id.startswith("ENSG") or upper_id.startswith("ENST")):
        return error_response(
            f"Invalid Ensembl ID '{ens_id}'. "
            "Expected a gene (ENSG...) or transcript (ENST...) stable ID.",
            400,
        )

    # ── 1. Fetch CDS from Ensembl ─────────────────────────────
    enforce_rate_limit("ensembl")
    try:
        resp = proxy_get(
            f"{ENSEMBL_BASE}/sequence/id/{ens_id}",
            params={"content-type": "text/plain", "type": "cds"},
        )

        if resp.status_code == 404:
            return error_response(
                f"No CDS sequence found for '{ens_id}'. "
                "The ID may be a non-coding gene, retired, or lack annotated CDS. "
                "Try using a protein-coding transcript ID (ENST...).",
                404,
            )
        if resp.status_code == 400:
            return error_response(
                f"Ensembl rejected the CDS request for '{ens_id}'. "
                "Ensure the ID is a valid, protein-coding Ensembl gene or transcript.",
                400,
            )
        if resp.status_code == 429:
            return error_response(
                "Ensembl rate limit exceeded. Please wait a few seconds and retry.",
                429,
            )
        if resp.status_code != 200:
            return error_response(
                f"Ensembl CDS fetch returned HTTP {resp.status_code}.",
                resp.status_code,
            )

        cds_raw = resp.text.strip()
        if not cds_raw:
            return error_response(
                f"Ensembl returned an empty CDS for '{ens_id}'.", 404
            )

    except requests.RequestException as e:
        logger.error("Ensembl translation fetch error: %s", e)
        return error_response(str(e), 502)

    # ── 2. Translate using Biopython ──────────────────────────
    cds_seq, warnings, err = validate_sequence(cds_raw, seq_type="dna")
    if err:
        return errremainder = len(cds_seq) % 3if remainder != 0:
        cds_seq = cds_seq[:-remainder]
        warnings.append(f"CDS length was not a multiple of 3 (remainder: {remainder}). Extra bases ignored.")

    bio_seq = BioSeq(cds_seq)
    has_start = cds_seq[:3] == "ATG"

    # Translate — stop at first stop codon (*), to_stop=False so we can report allprotein_with_stops = str(bio_seq.translate(table=1, to_stop=False))

    # Locate all stop codons in the CDSstop_codons = []
    for i in range(0, len(cds_seq) - 2, 3):
        codon = cds_seq[i:i + 3]
        if len(codon) == 3 and str(BioSeq(codon).translate(table=1)) == "*":
            stop_codons.append({
                "codon":            codon,
                "position":         (i // 3) + 1,   # 1-indexed amino acid position
                "nucleotide_offset": i,
            })

    # Clean protein: everything up to (but not including) the first * for the primary ORFprotein_seq = protein_with_stops.split("*")[0]
    has_stop = "*" in protein_with_stopsreturn jsonify({
        "ensembl_id":       ens_id,
        "cds_sequence":     cds_seq,
        "protein_sequence": protein_seq,
        "cds_length":       len(cds_seq),
        "protein_length":   len(protein_seq),
        "has_start_codon":  has_start,
        "has_stop_codon":   has_stop,
        "stop_codons":      stop_codons,
        "warnings":         warnings,
    })

# ─────────────────────────────────────────────────────────────
# API – Client Sequence Analysis
# ─────────────────────────────────────────────────────────────
@app.route("/api/analyze/dna", methods=["POST"])
def analyze_dna():
    """
    Performs heavy server-side bioinformatic computation on a provided DNA sequence.
    Offloading this to Python prevents browser crashes on large contiguous sequences.
    
    Expected JSON payload:
      {
        "sequence": "ATGC..."
      }
    
    Returns:
      {
        "length"            : 1234,
        "gc_content"        : 45.2,
        "reverse_complement": "...",
        "rna"               : "...",
        "protein_frames"    : ["...", "...", "..."]
      }
    """
    if not _BIOPYTHON_AVAILABLE:
        return error_response(
            "Biopython is not installed on this server.",
            501,
        )

    data = request.get_json()
    if not data or "sequence" not in data:
        return error_response("Missing 'sequence' in JSON body.", 400)

    raw_seq = str(data["sequence"])
    
    cln_seq, warnings, err = validate_sequence(raw_seq, seq_type="dna")
    if err:
        return errlength = len(cln_seq)
    
    # Manual GC computation avoids version-specific Bio.SeqUtils importsgc_count = cln_seq.count('G') + cln_seq.count('C') + cln_seq.count('S')
    gc_content = round((gc_count / length) * 100, 2) if length > 0 else 0.0b_seq = BioSeq(cln_seq)
    
    rev_comp = str(b_seq.reverse_complement())
    
    # Transcription (T -> U). Bio.Seq.transcribe() handles DNA->RNA strictly.
    # If the user passed RNA natively, transcribe() will raise an error if T is missing and U is present
    # We gracefully manually handle it if it contains U alreadyif "U" in cln_seq and "T" not in cln_seq:
        rna = cln_seqelse:
        rna = str(b_seq.transcribe())

    # 3-Frame conceptual translationframes = []
    for frame in range(3):
        # Slice sequence to start at the frame offsetframe_seq = cln_seq[frame:]
        # Biopython translate requires sequence length multiple of 3remainder = len(frame_seq) % 3if remainder != 0:
            frame_seq = frame_seq[:-remainder]
            
        if len(frame_seq) >= 3:
            frames.append(str(BioSeq(frame_seq).translate(table=1, to_stop=False)))
        else:
            frames.append("")

    return jsonify({
        "length"            : length,
        "gc_content"        : gc_content,
        "reverse_complement": rev_comp,
        "rna"               : rna,
        "protein_frames"    : frames,
        "warnings"          : warnings
    })


# ─────────────────────────────────────────────────────────────
# API – Unified Identity Mapping
# ─────────────────────────────────────────────────────────────
@app.route("/api/gene/map")
@cache.cached(timeout=86400, query_string=True, response_filter=_cacheable_response)
def gene_map():
    """
    Creates a unified identity layer mapping a gene query (symbol or accession)
    across NCBI, Ensembl, UniProt, and KEGG.
    """
    query = sanitize_input(request.args.get("query", ""))
    species = sanitize_input(request.args.get("species", "human"))
    
    if not query:
        return error_response("Missing 'query' parameter (symbol or accession required).", 400)
        
    out = {
        "ncbi_gene_id": None,
        "ensembl_id": None,
        "uniprot_id": None,
        "kegg_id": None,
        "organism": None,
        "confidence": 0.0
    }
    
    user_key = request.headers.get("X-NCBI-API-Key", "").strip()
    effective_key = user_key or NCBI_API_KEY

    # ── 1. Resolve NCBI Gene ID (via eSearch/eLink) ─────
    is_acc = "_" in queryif is_acc:
        db = "protein" if query.startswith("NP_") or query.startswith("XP_") or query.startswith("WP_") else "nuccore"
        req_url = f"{NCBI_BASE}/esearch.fcgi?db={db}&term={query}&retmode=json"
        if NCBI_API_KEY: req_url += f"&api_key={NCBI_API_KEY}"
        
        enforce_rate_limit("ncbi")
        try:
            r1 = proxy_get(req_url).json()
            uids = r1.get("esearchresult", {}).get("idlist", [])
            if uids:
                link_url = f"{NCBI_BASE}/elink.fcgi?dbfrom={db}&db=gene&id={uids[0]}&retmode=json"
                if effective_key: link_url += f"&api_key={effective_key}"
                enforce_rate_limit("ncbi")
                r2 = proxy_get(link_url).json()
                try:
                    out["ncbi_gene_id"] = r2["linksets"][0]["linksetdbs"][0]["links"][0]
                except (KeyError, IndexError):
                    passexcept Exception as e:
            logger.error("NCBI accession to gene mapping failed: %s", e)
    else:
        term = f"{query}[Gene Name] AND {species}[Organism]"
        req_url = f"{NCBI_BASE}/esearch.fcgi?db=gene&term={term}&retmode=json"
        if NCBI_API_KEY: req_url += f"&api_key={NCBI_API_KEY}"
        enforce_rate_limit("ncbi")
        try:
            r1 = proxy_get(req_url).json()
            uids = r1.get("esearchresult", {}).get("idlist", [])
            if uids:
                out["ncbi_gene_id"] = uids[0]
        except Exception as e:
            logger.error("NCBI symbol to gene mapping failed: %s", e)
            
    # ── 2. Resolve Ensembl ID (via xrefs) ───────────────
    try:
        # Avoid empty queries if species parsing is weird, but we validated queryif is_acc and out["ncbi_gene_id"]:
            url = f"{ENSEMBL_BASE}/xrefs/name/{species}/{out['ncbi_gene_id']}"
        else:
            url = f"{ENSEMBL_BASE}/xrefs/symbol/{species}/{query}"
            
        enforce_rate_limit("ensembl")
        r_ens = proxy_get(url, headers={"Content-type": "application/json"})
        if r_ens.ok:
            data = r_ens.json()
            if data and isinstance(data, list):
                out["ensembl_id"] = data[0].get("id")
    except Exception as e:
        logger.error("Ensembl xref mapping failed: %s", e)
            
    # ── 3. Resolve UniProt & KEGG (via UniProt API) ─────
    try:
        search_query = Noneif out["ensembl_id"]:
            search_query = out["ensembl_id"]
        elif out["ncbi_gene_id"]:
            search_query = out["ncbi_gene_id"]
        elif not is_acc:
            search_query = f"(gene:{query}) AND (reviewed:true)"
            
        if search_query:
            url = f"{UNIPROT_BASE}/uniprotkb/search?query={search_query}&format=json"
            enforce_rate_limit("uniprot")
            r_up = proxy_get(url)
            if r_up.ok:
                data = r_up.json().get("results", [])
                if data:
                    # Prefer reviewed entries if multiple existreviewed = [x for x in data if x.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"]
                    entry = reviewed[0] if reviewed else data[0]
                    
                    out["uniprot_id"] = entry.get("primaryAccession")
                    out["organism"]   = entry.get("organism", {}).get("scientificName", species)
                    # Parse KEGGfor d in entry.get("uniProtKBCrossReferences", []):
                        if d.get("database") == "KEGG":
                            out["kegg_id"] = d.get("id")
                            breakexcept Exception as e:
        logger.error("UniProt API mapping failed: %s", e)

    # ── 4. Confidence Score ─────────────────────────────
    score = 0.0if out["ncbi_gene_id"]: score += 0.25if out["ensembl_id"]:   score += 0.25if out["uniprot_id"]:   score += 0.25if out["kegg_id"]:      score += 0.25out["confidence"] = round(score, 2)
    
    if out["confidence"] == 0.0:
        return error_response(f"Could not map '{query}' to any known gene identities in specified databases.", 404)
        
    return jsonify(out)


# ─────────────────────────────────────────────────────────────
# Health Check (for Render)
# ─────────────────────────────────────────────────────────────
@app.route("/api/health")
def health():
    return jsonify({"status": "ok", "version": "1.1.0-scaling-prep", "db": "ready"})


# ─────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────
# API – Async Job Processing
# ─────────────────────────────────────────────────────────────
@app.route("/api/jobs/create", methods=["POST"])
def create_job():
    """
    Submits a bioinformatics task for asynchronous processing with intelligent caching.
    """
    from worker import process_sequence_taskdata = request.get_json()
    if not data or "tool_type" not in data:
        return error_response("Missing 'tool_type' in JSON body.", 400)
    
    tool_type = sanitize_input(data["tool_type"]).lower()
    params = data.get("params", {})
    
    # 1. Server-side Sequence Validation (Scientific grade)
    if "sequence" in params:
        stype = "protein" if tool_type in ("protein_analysis", "hydrophobicity") else "dna"
        cleaned, warnings, error = validate_sequence(params["sequence"], seq_type=stype)
        if error: return errorparams["sequence"] = cleaned
    
    # 2. Intelligent Task Caching (Deduplication)
    # Reduces server load for identical scientific queriescache_key = Noneif redis_client:
        try:
            # Generate a stable hash of the task requestimport hashlibpayload_str = json.dumps({"t": tool_type, "p": params}, sort_keys=True)
            task_hash = hashlib.sha256(payload_str.encode()).hexdigest()
            cache_key = f"task_cache:{task_hash}"
            
            existing_id = redis_client.get(cache_key)
            if existing_id:
                logger.info("[Cache HIT] Reusing task %s for matching request", existing_id.decode())
                return jsonify({
                    "status": "accepted",
                    "task_id": existing_id.decode(),
                    "message": "Found cached task for this request.",
                    "cached": True
                }), 202except Exception as e:
            logger.debug("Task caching failed: %s", e)
    
    # 3. Queue the Tasktask = process_sequence_task.delay(tool_type, params)
    
    # Store caching entry if possibleif cache_key and redis_client:
        try:
            redis_client.setex(cache_key, 86400, task.id) # Cache task mapping for 24hexcept Exception:
            passreturn jsonify({
        "status": "accepted",
        "task_id": task.id,
        "message": f"Task {tool_type} has been queued for processing.",
        "cached": False
    }), 202


@app.route("/api/jobs/status/<task_id>")
def get_task_status(task_id):
    """
    Endpoint to check the status of a background bioinformatics task.
    """
    from worker import celerytask = celery.AsyncResult(task_id)
    
    response = {"state": task.state, "task_id": task_id}
    
    if task.state == 'PENDING':
        response["status"] = "Queued..."
    elif task.state == 'STARTED':
        response["status"] = "Processing..."
    elif task.state == 'PROGRESS':
        response["status"] = "Processing..."
        response["progress"] = task.info.get('progress', 0) if task.info else 0elif task.state == 'SUCCESS':
        # task.result contains the dictionary returned by the workerresult = task.resultresponse["status"] = "Completed"
        response["result"] = resultelif task.state == 'FAILURE':
        response["status"] = "Failed"
        # Extract scientific error message if possibleinfo = task.inforesponse["error"] = str(info) if info else "Internal worker failure."
    else:
        response["status"] = task.statereturn jsonify(response)



# ─────────────────────────────────────────────────────────────
# Entry Point
# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    debug = os.environ.get("FLASK_ENV", "production") != "production"
    app.run(host="0.0.0.0", port=port, debug=debug)
