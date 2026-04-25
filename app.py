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
import json
import logging
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
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
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.SeqUtils import molecular_weight as bio_molecular_weight
    from Bio import motifs as Bio_motifs
    from Bio import Restriction as BioRestriction
    _BIOPYTHON_AVAILABLE = True
except ImportError:
    _BIOPYTHON_AVAILABLE = False

try:
    import redis
    _REDIS_AVAILABLE = True
except ImportError:
    _REDIS_AVAILABLE = False

# Celery / batch task queue
try:
    from worker import celery as _celery_app, process_batch_task
    from celery.result import AsyncResult
    _CELERY_AVAILABLE = True
except Exception as _celery_import_err:
    _CELERY_AVAILABLE = False
    logger = logging.getLogger("BioToolkit")
    logging.basicConfig(level=logging.INFO)
    logging.getLogger("BioToolkit").warning(
        "Celery worker not importable — batch endpoints will return 503: %s",
        _celery_import_err
    )

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
        logger.warning("Redis connection failed: %s", e)
        redis_client = None

# ─────────────────────────────────────────────────────────────
# Resilient HTTP Session — exponential backoff on 429/503
# ─────────────────────────────────────────────────────────────
_retry = Retry(
    total=3,
    backoff_factor=1.5,
    status_forcelist=[429, 503],
    allowed_methods=["GET"],
    raise_on_status=False,
)
_adapter = HTTPAdapter(max_retries=_retry)
http_session = requests.Session()
http_session.mount("https://", _adapter)
http_session.mount("http://", _adapter)

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
            # For primer design, we usually want nucleotide
            params = {
                "db": "nucleotide",
                "id": acc_id,
                "rettype": "gb",
                "retmode": "text",
                "tool": "BioToolkit",
                "email": "admin@biotoolkit.dev"
            }
            # Determine which API Key to use (Header > ENV)
            user_key = request.headers.get("X-NCBI-API-Key", "").strip()
            effective_key = user_key or NCBI_API_KEY
            
            if effective_key and str(effective_key).strip():
                params["api_key"] = effective_key

            resp = http_session.get("{}/efetch.fcgi".format(NCBI_BASE), params=params, timeout=15)
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
            params = {
                "db": "protein",
                "id": acc_id,
                "rettype": "fasta",
                "retmode": "text",
                "tool": "BioToolkit",
                "email": "admin@biotoolkit.dev"
            }
            user_key = request.headers.get("X-NCBI-API-Key", "").strip()
            effective_key = user_key or NCBI_API_KEY
            if effective_key and str(effective_key).strip():
                params["api_key"] = effective_key
            resp = http_session.get("{}/efetch.fcgi".format(NCBI_BASE), params=params, timeout=15)

        elif db_type == "uniprot":
            resp = http_session.get("{}/uniprotkb/{}.fasta".format(UNIPROT_BASE, acc_id), timeout=15)

        elif db_type == "ensembl":
            # Ensembl REST: use Accept header (not Content-Type) for GET requests
            json_headers = {"Accept": "application/json"}
            resp = http_session.get(
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
            resp = http_session.get("{}/sequence/id/{}".format(ENSEMBL_BASE, acc_id), headers=fasta_headers, timeout=15)

        elif db_type == "kegg":
            seq_type = request.args.get("type", "ntseq")
            if seq_type not in ("ntseq", "aaseq"):
                return jsonify({"error": "Invalid KEGG type. Use 'ntseq' or 'aaseq'."}), 400
            resp = http_session.get("{}/get/{}/{}".format(KEGG_BASE, acc_id, seq_type), timeout=15)
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

            s_resp = http_session.get("{}/esearch.fcgi".format(NCBI_BASE), params=search_params, timeout=15)
            s_resp.raise_for_status()
            s_data = s_resp.json()
            ids = s_data.get("esearchresult", {}).get("idlist", [])

            if not ids:
                return jsonify({"error": "No results for symbol '{}' in '{}'.".format(acc_id, species)}), 404

            gene_id = str(ids[0])

            # Step 2: ELink to find associated nucleotide accession
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

            l_resp = http_session.get("{}/elink.fcgi".format(NCBI_BASE), params=link_params, timeout=15)
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
            fetch_params = {
                "db": "nucleotide",
                "id": nuc_id,
                "rettype": "fasta",
                "retmode": "text"
            }
            if effective_key and str(effective_key).strip():
                fetch_params["api_key"] = effective_key

            resp = http_session.get("{}/efetch.fcgi".format(NCBI_BASE), params=fetch_params, timeout=15)

        else:
            return jsonify({"error": "Unsupported database: {}".format(db_type)}), 400

        resp.raise_for_status()
        cache.set(cache_key, {"text": resp.text}, timeout=86400)
        return Response(resp.text, mimetype="text/plain")

    except requests.exceptions.RequestException as e:
        logger.error("External API error: %s", str(e))
        upstream_status = getattr(getattr(e, "response", None), "status_code", None)
        if upstream_status == 429:
            return jsonify({"error": "Rate limit exceeded by biological database. Please try again in a few seconds."}), 429
        if upstream_status in (400, 404):
            return jsonify({"error": "Sequence or ID not found in the database."}), 404
        return jsonify({"error": "Failed to fetch from biological database: {}".format(str(e))}), 502

@app.route("/api/search", methods=["POST"])
def search_ncbi():
    """Proxy for NCBI ESearch."""
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
        resp = http_session.get("{}/esearch.fcgi".format(NCBI_BASE), params=params, timeout=10)
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
# Protein Analysis — BioPython ProtParam Backend
# ─────────────────────────────────────────────────────────────

# Allowed protein residues (IUPAC + ambiguous + selenocysteine + stop)
_VALID_AA = frozenset("ACDEFGHIKLMNPQRSTVWYXBZU*")
# Residues that confuse ProtParam — strip before analysis
_AMBIGUOUS_AA = frozenset("XBZ")

@app.route("/api/analyze/protein", methods=["POST"])
def analyze_protein():
    """Server-side protein physicochemical analysis via BioPython ProtParam."""
    if not _BIOPYTHON_AVAILABLE:
        return jsonify({"error": "Biopython not available on this server"}), 503

    data = request.get_json(silent=True)
    if not data:
        return jsonify({"error": "Invalid or missing JSON payload"}), 400

    raw = data.get("sequence", "")
    mass_type = data.get("mass_type", "average")
    redox_state = data.get("redox_state", "reduced")

    if mass_type not in ("average", "monoisotopic"):
        return jsonify({"error": "mass_type must be 'average' or 'monoisotopic'"}), 422
    if redox_state not in ("reduced", "oxidized"):
        return jsonify({"error": "redox_state must be 'reduced' or 'oxidized'"}), 422

    # Strip FASTA header lines
    lines = raw.strip().splitlines()
    seq_lines = [l.strip() for l in lines if not l.strip().startswith(">")]
    seq = "".join(seq_lines).upper()

    # Reject empty or suspiciously short sequences
    if not seq:
        return jsonify({"error": "No sequence provided"}), 400
    if len(seq) < 4:
        return jsonify({"error": "Sequence too short (minimum 4 residues)"}), 400
    if len(seq) > 50000:
        return jsonify({"error": "Sequence too long (maximum 50,000 residues)"}), 400

    # Reject characters that are not valid protein residues
    invalid = set(seq) - _VALID_AA
    if invalid:
        return jsonify({"error": "Invalid characters for a protein sequence: {}".format("".join(sorted(invalid)))}), 422

    # Truncate at first internal stop codon
    seq_clean = seq.split("*")[0] if "*" in seq else seq
    truncated = len(seq_clean) < len(seq)

    # Ambiguous residues (X/B/Z) cannot be assigned a definite molecular formula.
    # Silent removal would produce a scientifically incorrect result — abort instead.
    if set(seq_clean) & _AMBIGUOUS_AA:
        return jsonify({
            "error": "Ambiguous residue detected. Precise calculation aborted to maintain scientific integrity.",
            "code": "AMBIGUOUS_RESIDUE"
        }), 422

    # Selenocysteine → Cysteine is a standard, documented approximation.
    analysis_seq = seq_clean.replace("U", "C")

    try:
        pa = ProteinAnalysis(analysis_seq)

        # Molecular weight
        if mass_type == "monoisotopic":
            mw_da = bio_molecular_weight(analysis_seq, seq_type="protein", monoisotopic=True)
        else:
            mw_da = pa.molecular_weight()
        mw_kda = round(mw_da / 1000.0, 3)

        # Isoelectric point
        pi = round(pa.isoelectric_point(), 2)

        # Extinction coefficient (reduced, oxidized)
        ec_reduced, ec_oxidized = pa.molar_extinction_coefficient()
        ext_coeff = ec_oxidized if redox_state == "oxidized" else ec_reduced

        # GRAVY hydrophobicity index
        gravy = round(pa.gravy(), 3)

        # Aliphatic index (Ikai 1980)
        aliphatic = round(pa.aliphatic_index(), 2)

        # Instability index (Guruprasad 1990)
        instability = round(pa.instability_index(), 2)

        # Net charge at physiological pH 7.4
        net_charge = round(pa.charge_at_pH(7.4), 2)

        return jsonify({
            "length": len(analysis_seq),
            "molecular_weight_kda": mw_kda,
            "mass_type": mass_type,
            "isoelectric_point": pi,
            "extinction_coefficient": ext_coeff,
            "redox_state": redox_state,
            "gravy": gravy,
            "aliphatic_index": aliphatic,
            "instability_index": instability,
            "instability_label": "Stable" if instability < 40 else "Unstable",
            "net_charge_ph74": net_charge,
            "truncated": truncated
        })

    except Exception as e:
        logger.error("ProtParam analysis error: %s", str(e))
        return jsonify({"error": "Analysis failed: {}".format(str(e))}), 500


# ─────────────────────────────────────────────────────────────
# UniProt Biological Context — Real Metadata Endpoint
# ─────────────────────────────────────────────────────────────

# UniProt accession format: [A-Z][0-9][A-Z0-9]{3}[0-9] or [O,P,Q][0-9][A-Z0-9]{3}[0-9]
_UNIPROT_ACC_RE = re.compile(r'^[A-Za-z][A-Za-z0-9_\-]{1,19}$')

@app.route("/api/uniprot/<accession_id>", methods=["GET"])
def get_uniprot_entry(accession_id):
    """Fetch rich biological metadata for a UniProt accession from UniProt REST API."""
    if not _UNIPROT_ACC_RE.match(accession_id):
        return jsonify({"error": "Invalid UniProt accession format"}), 400

    cache_key = "uniprot_meta_{}".format(accession_id.upper())
    cached = cache.get(cache_key)
    if cached:
        return jsonify(cached)

    try:
        resp = http_session.get(
            "{}/uniprotkb/{}.json".format(UNIPROT_BASE, accession_id),
            headers={"Accept": "application/json"},
            timeout=15
        )
        if resp.status_code == 404:
            return jsonify({"error": "UniProt entry '{}' not found".format(accession_id)}), 404
        resp.raise_for_status()
        raw = resp.json()

        # Protein name
        protein_name = ""
        rec_name = raw.get("proteinDescription", {}).get("recommendedName", {})
        if rec_name:
            protein_name = rec_name.get("fullName", {}).get("value", "")
        if not protein_name:
            sub_names = raw.get("proteinDescription", {}).get("submissionNames", [])
            if sub_names:
                protein_name = sub_names[0].get("fullName", {}).get("value", "")

        # Gene symbol
        gene_symbol = ""
        genes = raw.get("genes", [])
        if genes:
            gene_symbol = genes[0].get("geneName", {}).get("value", "")

        # Organism
        organism = raw.get("organism", {}).get("scientificName", "")

        # Sequence
        sequence = raw.get("sequence", {}).get("value", "")
        sequence_length = raw.get("sequence", {}).get("length", 0)

        # PTM text descriptions (commentType == "PTM")
        ptms = []
        for comment in raw.get("comments", []):
            if comment.get("commentType") == "PTM":
                for t in comment.get("texts", []):
                    val = t.get("value", "").strip()
                    if val:
                        ptms.append(val)

        # Subcellular locations (commentType == "SUBCELLULAR LOCATION")
        locations = []
        for comment in raw.get("comments", []):
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for sub in comment.get("subcellularLocations", []):
                    loc = sub.get("location", {}).get("value", "")
                    if loc and loc not in locations:
                        locations.append(loc)

        # Structural features: TM helices, signal peptide, modified residues
        _FEATURE_TYPES = {
            "Transmembrane", "Signal", "Propeptide",
            "Modified residue", "Glycosylation", "Lipidation", "Disulfide bond"
        }
        features = []
        for feat in raw.get("features", []):
            ftype = feat.get("type", "")
            if ftype in _FEATURE_TYPES:
                loc_obj = feat.get("location", {})
                entry = {
                    "type": ftype,
                    "start": loc_obj.get("start", {}).get("value"),
                    "end": loc_obj.get("end", {}).get("value"),
                }
                desc = feat.get("description", "").strip()
                if desc:
                    entry["description"] = desc
                features.append(entry)

        result = {
            "accession": accession_id.upper(),
            "protein_name": protein_name,
            "gene_symbol": gene_symbol,
            "organism": organism,
            "sequence": sequence,
            "sequence_length": sequence_length,
            "ptms": ptms[:5],
            "subcellular_locations": locations,
            "features": features
        }

        cache.set(cache_key, result, timeout=86400)
        return jsonify(result)

    except requests.exceptions.RequestException as e:
        logger.error("UniProt metadata fetch error: %s", str(e))
        return jsonify({"error": "Failed to fetch from UniProt: {}".format(str(e))}), 502


# ─────────────────────────────────────────────────────────────
# High-Throughput Batch Processing (Step 8)
# ─────────────────────────────────────────────────────────────

@app.route("/api/analyze/batch", methods=["POST"])
def analyze_batch():
    """
    Accept a multi-FASTA file upload (multipart/form-data, field name 'file').
    Dispatches an async Celery task and immediately returns a job_id.
    The client should poll GET /api/status/<job_id> for progress.
    """
    if not _CELERY_AVAILABLE:
        return jsonify({
            "error": "Batch processing requires the background worker (Celery + Redis). "
                     "Start it with: celery -A worker worker --loglevel=info"
        }), 503

    if "file" not in request.files:
        return jsonify({"error": "No file field in request. POST with multipart/form-data, field name 'file'."}), 400

    uploaded = request.files["file"]
    if not uploaded.filename:
        return jsonify({"error": "No file selected."}), 400

    try:
        content = uploaded.read().decode("utf-8", errors="replace").strip()
    except Exception as e:
        return jsonify({"error": "Could not read uploaded file: {}".format(str(e))}), 400

    if not content:
        return jsonify({"error": "Uploaded file is empty."}), 400

    # Lightweight FASTA sanity check — must contain at least one header line
    if ">" not in content:
        return jsonify({
            "error": "File does not appear to be in FASTA format (no '>' header lines found)."
        }), 422

    task = process_batch_task.apply_async(args=[content])
    logger.info("Batch job queued: %s", task.id)

    return jsonify({"job_id": task.id, "status": "queued"}), 202


@app.route("/api/status/<job_id>", methods=["GET"])
def batch_job_status(job_id):
    """
    Poll for async batch job progress.
    Returns JSON: { status, processed, total, progress, [failed], [error] }
    """
    if not _CELERY_AVAILABLE:
        return jsonify({"error": "Celery not available"}), 503

    # Sanitize job_id (Celery UUIDs are hex + hyphens)
    if not re.match(r'^[a-f0-9\-]{8,64}$', job_id, re.IGNORECASE):
        return jsonify({"error": "Invalid job ID format"}), 400

    result = AsyncResult(job_id, app=_celery_app)
    state  = result.state

    if state == "PENDING":
        return jsonify({"status": "pending", "processed": 0, "total": 0, "progress": 0})

    if state == "STARTED":
        return jsonify({"status": "running", "processed": 0, "total": 0, "progress": 0})

    if state == "PROGRESS":
        meta = result.info or {}
        return jsonify({
            "status":    "running",
            "processed": meta.get("processed", 0),
            "total":     meta.get("total", 0),
            "progress":  meta.get("progress", 0),
        })

    if state == "SUCCESS":
        res = result.result or {}
        # If the worker itself returned an error dict
        if res.get("status") == "error":
            return jsonify({"status": "failed", "error": res.get("message", "Worker error")}), 200
        return jsonify({
            "status":   "complete",
            "total":    res.get("total", 0),
            "failed":   res.get("failed", 0),
            "progress": 100,
        })

    if state in ("FAILURE", "REVOKED"):
        err_msg = str(result.info) if result.info else state
        return jsonify({"status": "failed", "error": err_msg})

    return jsonify({"status": state.lower(), "progress": 0})


@app.route("/api/batch/<job_id>/download", methods=["GET"])
def download_batch_results(job_id):
    """
    Return completed batch results as a CSV file attachment.
    The CSV is read from the Celery result stored in Redis.
    """
    if not _CELERY_AVAILABLE:
        return jsonify({"error": "Celery not available"}), 503

    if not re.match(r'^[a-f0-9\-]{8,64}$', job_id, re.IGNORECASE):
        return jsonify({"error": "Invalid job ID format"}), 400

    result = AsyncResult(job_id, app=_celery_app)
    if result.state != "SUCCESS":
        return jsonify({
            "error": "Job is not complete (state: {}).".format(result.state)
        }), 400

    res = result.result or {}
    csv_data = res.get("csv", "")
    if not csv_data:
        return jsonify({"error": "No CSV data found for this job."}), 404

    short_id = job_id.replace("-", "")[:8].upper()
    filename = "batch_results_{}.csv".format(short_id)

    return Response(
        csv_data,
        mimetype="text/csv",
        headers={
            "Content-Disposition": 'attachment; filename="{}"'.format(filename),
            "Content-Type":        "text/csv; charset=utf-8",
        },
    )


# ─────────────────────────────────────────────────────────────
# Promoter & Motif Analysis — PSSM-based, Python backend
# ─────────────────────────────────────────────────────────────

# Eukaryotic core promoter element PSSMs.
# Sources: Bucher 1990 (NAR 18:6299); Breathnach & Chambon 1981.
_EUK_MOTIF_DEFS = [
    {
        "name": "TATA-box",
        "consensus": "TATAAATA",
        "typical_tss_offset": -30,
        # 8-position matrix (Bucher 1990, Table 2)
        "counts": {
            "A": [14, 99, 100, 81, 99, 16, 56, 63],
            "C": [ 2,  0,   0,  1,  0,  2, 11,  8],
            "G": [ 2,  0,   0, 15,  0,  1, 29, 25],
            "T": [82,  1,   0,  3,  1, 81,  4,  4],
        },
        "threshold_fraction": 0.80,
    },
    {
        "name": "CAAT-box",
        "consensus": "CCAAT",
        "typical_tss_offset": -80,
        # 5-position core matrix (Breathnach & Chambon 1981)
        "counts": {
            "A": [ 0,  3, 100, 100,  4],
            "C": [71, 79,   0,   0, 14],
            "G": [ 8,  5,   0,   0, 24],
            "T": [21, 13,   0,   0, 58],
        },
        "threshold_fraction": 0.78,
    },
    {
        "name": "GC-box",
        "consensus": "GGGCGG",
        "typical_tss_offset": -90,
        # 6-position SP1-site matrix (Bucher 1990, Table 4)
        "counts": {
            "A": [ 2,  1,  1,  3,  1,  1],
            "C": [ 3,  4,  4, 82,  5,  5],
            "G": [91, 91, 91,  8, 90, 90],
            "T": [ 4,  4,  4,  7,  4,  4],
        },
        "threshold_fraction": 0.82,
    },
]

# Prokaryotic σ70 promoter element PSSMs.
# Count matrices compiled from 168 E. coli σ70 promoters.
# Source: Harley & Reynolds 1987 (NAR 15:2343).
_PROK_MOTIF_DEFS = [
    {
        "name": "-35 box",
        "consensus": "TTGACA",
        "typical_tss_offset": -35,
        "counts": {
            "A": [ 9,  9,  8, 62, 13, 54],
            "C": [ 7, 10,  8, 11, 60,  5],
            "G": [ 7,  6, 69, 13, 11, 10],
            "T": [71, 68, 10,  9, 16, 28],
        },
        "threshold_fraction": 0.75,
    },
    {
        "name": "-10 box",
        "consensus": "TATAAT",
        "typical_tss_offset": -10,
        "counts": {
            "A": [10, 85, 10, 63, 48, 11],
            "C": [ 7,  8,  6, 10,  8,  6],
            "G": [ 3,  6,  2, 12, 16,  8],
            "T": [80,  1, 82, 15, 28, 75],
        },
        "threshold_fraction": 0.75,
    },
]

# Valid spacer (bp) between 3' end of -35 hexamer and 5' end of -10 hexamer.
# Source: Harley & Reynolds 1987 (NAR 15:2343).
_PROK_SPACER_MIN = 15
_PROK_SPACER_MAX = 19

# Distance (bp) from TATA box 5'-most base to TSS proxy (ORF ATG).
# Source: Bucher 1990 (NAR 18:6299).
_EUK_TATA_TSS_MIN = 25
_EUK_TATA_TSS_MAX = 35

_MIN_ORF_NT = 90   # 30 amino acids minimum
_MAX_PROMOTER_SEQ = 100_000


# ── Promoter helpers ─────────────────────────────────────────

def _pssm_scan_motif(mdef, scan_str, bio_seq):
    """Run PSSM scan for one motif definition. Returns a list of raw hit dicts."""
    hits = []
    motif_obj = Bio_motifs.Motif(counts=mdef["counts"])
    motif_obj.pseudocounts = 0.5
    pssm = motif_obj.pssm
    if pssm.max <= 0:
        return hits
    threshold = mdef["threshold_fraction"] * pssm.max
    motif_len = motif_obj.length
    seq_len = len(scan_str)
    for strand_label, search_seq in (("+", bio_seq), ("-", bio_seq.reverse_complement())):
        scores = pssm.calculate(search_seq)
        if scores is None:
            continue
        for i, score in enumerate(scores):
            # NaN values (N bases) compare False — skipped automatically
            if score >= threshold:
                if strand_label == "+":
                    hit_start = i + 1
                    hit_end   = i + motif_len
                    hit_seq   = scan_str[i: i + motif_len]
                else:
                    hit_start = seq_len - (i + motif_len) + 1
                    hit_end   = seq_len - i
                    hit_seq   = str(search_seq[i: i + motif_len])
                hits.append({
                    "motif":              mdef["name"],
                    "strand":             strand_label,
                    "start":              hit_start,
                    "end":                hit_end,
                    "score":              round(float(score), 3),
                    "max_score":          round(float(pssm.max), 3),
                    "score_pct":          round(float(score) / float(pssm.max) * 100, 1),
                    "sequence":           hit_seq,
                    "consensus":          mdef["consensus"],
                    "typical_tss_offset": mdef.get("typical_tss_offset", 0),
                })
    return hits


def _scan_prokaryotic(scan_str, bio_seq):
    """
    Find σ70 -35/-10 box pairs whose inter-element spacer is 15–19 bp.

    Spacer is measured as the gap (in bp) between the 3' end of the -35
    hexamer and the 5' end of the -10 hexamer, following the convention
    of Harley & Reynolds 1987 (NAR 15:2343).

    Only motifs that participate in at least one valid pair are returned.
    Each returned hit carries a 'spacer_to_partner_bp' field.
    """
    hits_35 = _pssm_scan_motif(_PROK_MOTIF_DEFS[0], scan_str, bio_seq)
    hits_10 = _pssm_scan_motif(_PROK_MOTIF_DEFS[1], scan_str, bio_seq)

    validated = []
    for strand in ("+", "-"):
        s35 = [h for h in hits_35 if h["strand"] == strand]
        s10 = [h for h in hits_10 if h["strand"] == strand]
        for h35 in s35:
            for h10 in s10:
                if strand == "+":
                    # spacer = gap between 3' end of -35 and 5' end of -10
                    spacer = h10["start"] - h35["end"] - 1
                else:
                    # Minus-strand hits are reported in fwd coords.
                    # Biologically, the -35 box lies at higher fwd coordinates
                    # than the -10 box (gene reads right→left).
                    # hit["start"] is the downstream (3') tip in minus-strand
                    # orientation; hit["end"] is the upstream (5') tip.
                    spacer = h35["start"] - h10["end"] - 1
                if _PROK_SPACER_MIN <= spacer <= _PROK_SPACER_MAX:
                    h35c = dict(h35)
                    h10c = dict(h10)
                    h35c["spacer_to_partner_bp"] = spacer
                    h10c["spacer_to_partner_bp"] = spacer
                    validated.extend([h35c, h10c])

    # Deduplicate: a motif may satisfy multiple valid pairings
    seen = set()
    result = []
    for m in validated:
        key = (m["motif"], m["strand"], m["start"])
        if key not in seen:
            seen.add(key)
            result.append(m)
    return result


def _tata_tss_anchored(tata_hit, orfs):
    """
    Return True iff an ORF ATG codon lies 25–35 bp downstream of the
    TATA box 5'-most base (Bucher 1990).

    For + strand: TSS window = [tata_start + 25, tata_start + 35].
    For - strand: TSS window mirrors symmetrically using fwd coordinates
    (TATA 5' tip = tata_hit["end"]; downstream = lower fwd coords).
    """
    if tata_hit["strand"] == "+":
        lo = tata_hit["start"] + _EUK_TATA_TSS_MIN
        hi = tata_hit["start"] + _EUK_TATA_TSS_MAX
        return any(lo <= o["start"] <= hi for o in orfs if o["strand"] == "+")
    else:
        # orf["end"] holds the fwd coordinate of the ATG for minus-strand ORFs
        lo = tata_hit["end"] - _EUK_TATA_TSS_MAX
        hi = tata_hit["end"] - _EUK_TATA_TSS_MIN
        return any(lo <= o["end"] <= hi for o in orfs if o["strand"] == "-")


def _scan_eukaryotic(scan_str, bio_seq, orfs):
    """
    Scan for TATA-box, CAAT-box, and GC-box using eukaryotic PSSMs.
    TATA-box hits are filtered to those whose TSS proxy (nearest downstream
    ORF ATG) falls 25–35 bp downstream (Bucher 1990).
    """
    found = []
    for mdef in _EUK_MOTIF_DEFS:
        hits = _pssm_scan_motif(mdef, scan_str, bio_seq)
        if mdef["name"] == "TATA-box":
            hits = [h for h in hits if _tata_tss_anchored(h, orfs)]
        found.extend(hits)
    return found


def _dist_to_nearest_orf_start(motif, orfs):
    """
    Distance in bp from the motif's 3' end (in the direction of transcription)
    to the nearest downstream ORF start codon.

    + strand: 3' end of motif = motif["end"]; ORF ATG = orf["start"].
    - strand: 3' end of motif = motif["start"] (lower fwd coord);
              ORF ATG (fwd coord) = orf["end"].

    Returns None when no downstream ORF exists.
    """
    if motif["strand"] == "+":
        dists = [o["start"] - motif["end"]
                 for o in orfs if o["strand"] == "+" and o["start"] > motif["end"]]
    else:
        dists = [motif["start"] - o["end"]
                 for o in orfs if o["strand"] == "-" and o["end"] < motif["start"]]
    return min(dists) if dists else None


@app.route("/api/analyze/promoter", methods=["POST"])
def analyze_promoter():
    """
    Organism-aware PSSM promoter motif scanning and 6-frame ORF detection.

    Accepts:
      {
        "sequence":        "ATCG...",
        "organism_type":   "prokaryote" | "eukaryote",
        "upstream_length": 500          (optional)
      }

    Returns: { "orfs": [...], "motifs": [...], "error": null }

    Prokaryote mode: scans for σ70 -35/-10 box pairs with a valid 15–19 bp
    spacer (Harley & Reynolds 1987). Eukaryotic elements are not searched.

    Eukaryote mode: scans for TATA, CAAT, GC-box elements (Bucher 1990).
    TATA-box hits are further filtered to those where an ORF ATG lies
    25–35 bp downstream (TSS anchoring). Prokaryotic elements are not searched.

    Every returned motif includes 'dist_to_nearest_orf_start' (bp to the
    nearest downstream ORF ATG, or null if none found).
    """
    if not _BIOPYTHON_AVAILABLE:
        return jsonify({"orfs": [], "motifs": [], "error": "Biopython not available on this server"}), 503

    data = request.get_json(silent=True)
    if not data:
        return jsonify({"orfs": [], "motifs": [], "error": "Invalid or missing JSON payload"}), 400

    raw_seq = data.get("sequence", "")
    upstream_length = data.get("upstream_length", None)
    organism_type = data.get("organism_type", "").lower().strip()

    if organism_type not in ("prokaryote", "eukaryote"):
        return jsonify({
            "orfs": [], "motifs": [],
            "error": "organism_type must be 'prokaryote' or 'eukaryote'"
        }), 400

    lines = raw_seq.strip().splitlines()
    seq_str = "".join(l.strip() for l in lines if not l.strip().startswith(">")).upper()

    if not seq_str:
        return jsonify({"orfs": [], "motifs": [], "error": "No sequence provided"}), 400
    if len(seq_str) > _MAX_PROMOTER_SEQ:
        return jsonify({"orfs": [], "motifs": [], "error": "Sequence too long (max {:,} bp)".format(_MAX_PROMOTER_SEQ)}), 400

    invalid_chars = set(seq_str) - set("ACGTN")
    if invalid_chars:
        return jsonify({"orfs": [], "motifs": [], "error": "Invalid DNA characters: {}".format("".join(sorted(invalid_chars)))}), 422

    scan_str = seq_str
    if upstream_length and isinstance(upstream_length, int) and 0 < upstream_length < len(seq_str):
        scan_str = seq_str[:upstream_length]

    try:
        bio_seq = BioSeq(scan_str)

        # ── 6-Frame ORF Detection (full sequence, both strands) ───────────
        orfs = []
        stop_codons = {"TAA", "TAG", "TGA"}

        for strand_label, strand_seq in (("+", seq_str), ("-", str(BioSeq(seq_str).reverse_complement()))):
            seq_len = len(seq_str)
            for frame in range(3):
                start_pos = None
                i = frame
                while i + 2 < len(strand_seq):
                    codon = strand_seq[i: i + 3]
                    if codon == "ATG" and start_pos is None:
                        start_pos = i
                    elif codon in stop_codons and start_pos is not None:
                        orf_nt_len = (i + 3) - start_pos
                        if orf_nt_len >= _MIN_ORF_NT:
                            orf_dna = strand_seq[start_pos: i + 3]
                            try:
                                protein = str(BioSeq(orf_dna).translate(to_stop=True))
                            except Exception:
                                protein = ""
                            if strand_label == "+":
                                fwd_start = start_pos + 1
                                fwd_end   = i + 3
                            else:
                                fwd_start = seq_len - (i + 3) + 1
                                fwd_end   = seq_len - start_pos
                            orfs.append({
                                "frame":     (frame + 1) if strand_label == "+" else -(frame + 1),
                                "strand":    strand_label,
                                "start":     fwd_start,
                                "end":       fwd_end,
                                "length_nt": orf_nt_len,
                                "length_aa": len(protein),
                                "protein":   protein,
                            })
                        start_pos = None
                    i += 3

        orfs.sort(key=lambda x: x["start"])

        # ── Organism-Specific Motif Scanning ─────────────────────────────
        if organism_type == "prokaryote":
            found_motifs = _scan_prokaryotic(scan_str, bio_seq)
        else:
            found_motifs = _scan_eukaryotic(scan_str, bio_seq, orfs)

        # ── Distance to Nearest Downstream ORF Start Codon ───────────────
        for motif in found_motifs:
            motif["dist_to_nearest_orf_start"] = _dist_to_nearest_orf_start(motif, orfs)

        found_motifs.sort(key=lambda x: x["score"], reverse=True)

        return jsonify({"orfs": orfs, "motifs": found_motifs, "error": None})

    except Exception as e:
        logger.error("Promoter analysis error: %s", str(e))
        return jsonify({"orfs": [], "motifs": [], "error": "Analysis failed: {}".format(str(e))}), 500


# ─────────────────────────────────────────────────────────────
# Restriction Enzyme Analysis — Biopython Bio.Restriction backend
# ─────────────────────────────────────────────────────────────

_MAX_RESTRICTION_SEQ = 100_000

@app.route("/api/analyze/restriction", methods=["POST"])
def analyze_restriction():
    """
    Restriction enzyme digest analysis via Biopython Bio.Restriction.

    Accepts:
      {
        "sequence":  "ATCG...",
        "topology":  "linear" | "circular",
        "enzymes":   ["EcoRI", "BamHI"] | "ALL"
      }

    Returns:
      {
        "results": [{"enzyme": "EcoRI", "cuts": [105, 400], "fragments": [295, 4705]}],
        "error": null
      }
    """
    if not _BIOPYTHON_AVAILABLE:
        return jsonify({"results": [], "error": "Biopython not available on this server"}), 503

    data = request.get_json(silent=True)
    if not data:
        return jsonify({"results": [], "error": "Invalid or missing JSON payload"}), 400

    # Parse and clean sequence (strip FASTA headers)
    raw_seq = data.get("sequence", "")
    lines = raw_seq.strip().splitlines()
    seq_str = "".join(l.strip() for l in lines if not l.strip().startswith(">"))
    # Strip whitespace, digits, and any non-alphabet character, then normalise case
    seq_str = re.sub(r'[^A-Za-z]', '', seq_str).upper()
    # Reject sequences containing letters outside the IUPAC nucleotide alphabet
    if re.search(r'[^ACGTURYSWKMBDHVN]', seq_str):
        return jsonify({"results": [], "error": "Sequence contains invalid characters. Only IUPAC nucleotide codes are accepted."}), 400

    if not seq_str:
        return jsonify({"results": [], "error": "No sequence provided"}), 400
    if len(seq_str) > _MAX_RESTRICTION_SEQ:
        return jsonify({"error": "Sequence exceeds the 100,000 bp limit for the free research tier."}), 413

    # Validate topology
    topology = data.get("topology", "").lower().strip()
    if topology not in ("linear", "circular"):
        return jsonify({"results": [], "error": "topology must be 'linear' or 'circular'"}), 400

    # Parse enzymes
    enzymes_input = data.get("enzymes")
    if enzymes_input is None:
        return jsonify({"results": [], "error": "enzymes field is required"}), 400

    try:
        invalid_enzymes = []
        if isinstance(enzymes_input, str):
            if enzymes_input == "ALL":
                batch = BioRestriction.CommOnly
            else:
                return jsonify({"results": [], "error": "enzymes must be a list of enzyme names or the string 'ALL'"}), 400
        elif isinstance(enzymes_input, list):
            if not enzymes_input:
                return jsonify({"results": [], "error": "enzymes list must not be empty"}), 400
            valid_enzymes = []
            invalid_enzymes = []
            for name in enzymes_input:
                enz_name = str(name)
                enz = getattr(BioRestriction, enz_name, None)
                if enz is None or not hasattr(BioRestriction, enz_name):
                    invalid_enzymes.append(enz_name)
                else:
                    valid_enzymes.append(enz)
            if not valid_enzymes:
                return jsonify({"results": [], "warnings": {"unknown_enzymes": invalid_enzymes}, "error": "None of the requested enzymes are recognized by the Biopython database."}), 400
            batch = BioRestriction.RestrictionBatch(valid_enzymes)
        else:
            return jsonify({"results": [], "error": "enzymes must be a list of enzyme names or the string 'ALL'"}), 400

        seq_obj = BioSeq(seq_str)
        search_results = batch.search(seq_obj, linear=(topology == "linear"))

        seq_len = len(seq_str)
        results = []

        for enzyme, cuts in search_results.items():
            if not cuts:
                continue

            cuts = sorted(cuts)
            enzyme_name = str(enzyme)

            if topology == "linear":
                positions = [0] + cuts + [seq_len]
                fragments = [positions[i + 1] - positions[i] for i in range(len(positions) - 1)]
            else:
                if len(cuts) == 1:
                    fragments = [seq_len]
                else:
                    fragments = [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]
                    wrap_fragment = (seq_len - cuts[-1]) + cuts[0]
                    fragments.append(wrap_fragment)

            results.append({
                "enzyme": enzyme_name,
                "cuts": cuts,
                "fragments": sorted(fragments),
            })

        results.sort(key=lambda x: x["enzyme"])
        warnings = {"unknown_enzymes": invalid_enzymes} if isinstance(enzymes_input, list) and invalid_enzymes else {"unknown_enzymes": []}
        return jsonify({"results": results, "warnings": warnings, "error": None})

    except Exception as e:
        logger.error("Restriction analysis error: %s", str(e))
        return jsonify({"results": [], "error": "Analysis failed: {}".format(str(e))}), 500


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
