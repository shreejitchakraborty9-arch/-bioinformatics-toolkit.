import os
import io as _io
import csv
from celery import Celery

# Celery Configuration
CELERY_BROKER_URL    = os.environ.get("REDIS_URL", "redis://localhost:6379/0")
CELERY_RESULT_BACKEND = os.environ.get("REDIS_URL", "redis://localhost:6379/0")

celery = Celery(
    "biotoolkit-worker",
    broker=CELERY_BROKER_URL,
    backend=CELERY_RESULT_BACKEND,
)

# Standard Celery settings for reliability
celery.conf.update(
    task_serializer='json',
    accept_content=['json'],
    result_serializer='json',
    timezone='UTC',
    enable_utc=True,
    task_track_started=True,
    task_time_limit=600,       # 10-minute ceiling for large FASTA batches
    result_expires=86400,      # keep results in Redis for 24 h
)

try:
    from Bio import Align, motifs, SeqIO, Restriction
    from Bio.Seq import Seq
    from Bio.SeqUtils import GC_skew, MeltingTemp as mt
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as _ProteinAnalysis
    from Bio.SeqUtils import molecular_weight as _bio_mw
    _BIOPYTHON_AVAILABLE = True
except ImportError:
    _BIOPYTHON_AVAILABLE = False

# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────

# Valid protein residues (standard 20 + selenocysteine mapped to C)
_PROTEIN_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")

def _detect_seq_type(seq_str):
    """Return 'dna', 'rna', or 'protein' from sequence composition."""
    s   = seq_str.upper()
    tot = len(s)
    if tot == 0:
        return "unknown"

    dna_chars = frozenset("ACGTN")
    rna_chars = frozenset("ACGUN")

    dna_hits = sum(1 for c in s if c in dna_chars)
    rna_hits = sum(1 for c in s if c in rna_chars)

    if "U" in s and rna_hits / tot > 0.85:
        return "rna"
    if dna_hits / tot > 0.85:
        return "dna"
    return "protein"

@celery.task(bind=True, max_retries=3)
def process_sequence_task(self, tool_type, params):
    """
    Scientific-grade bioinformatics worker task.
    """
    if not _BIOPYTHON_AVAILABLE:
        return {"status": "error", "message": "Biopython is not available on this worker."}

    # Use .format() for Python 2.7 compatibility if needed, though target is Python 3
    print("Executing {} task with params: {}".format(tool_type, params.keys()))
    
    try:
        t_type = tool_type.lower().replace('-', '_')
        
        # 1. Alignment (Global/Local SMITH-WATERMAN / NEEDLEMAN-WUNSCH)
        if t_type in ("alignment", "align_global", "align_local"):
            seq1 = params.get("seq1", params.get("a", ""))
            seq2 = params.get("seq2", params.get("b", ""))
            mode = "local" if t_type == "align_local" else "global"
            
            aligner = Align.PairwiseAligner()
            aligner.mode = mode
            if "match" in params: aligner.match_score = float(params["match"])
            if "mismatch" in params: aligner.mismatch_score = float(params["mismatch"])
            if "gap" in params: aligner.open_gap_score = float(params["gap"])
            
            alignments = aligner.align(seq1, seq2)
            if not alignments:
                return {"score": 0, "alnA": "", "alnB": "", "mid": ""}
            
            best = alignments[0]
            alnA = str(best[0])
            alnB = str(best[1])
            mid = "".join("|" if a == b else (" " if a == "-" or b == "-" else "·") for a, b in zip(alnA, alnB))
            
            return {
                "score": float(best.score),
                "alnA": alnA,
                "alnB": alnB,
                "mid": mid
            }
            
        # 2. GC-Skew Analysis
        elif t_type in ("gc_skew", "gc-skew"):
            sequence = params.get("sequence", "")
            window = int(params.get("window", 100))
            if not sequence: return {"error": "No sequence"}
            
            skew = GC_skew(sequence, window)
            return {"skew": skew}

        # 3. Restriction Enzyme Mapping
        elif t_type in ("restriction_search", "restriction-search", "restriction_map"):
            sequence = params.get("sequence", "")
            enzymes = params.get("enzymes", []) # List of enzyme names
            if not sequence: return {"error": "No sequence"}
            
            seq_obj = Seq(sequence)
            # Filter enzymes to only those that exist in Biopython
            valid_enzymes = [e for e in enzymes if hasattr(Restriction, e)]
            rb = Restriction.RestrictionBatch(valid_enzymes) if valid_enzymes else Restriction.CommOnly
            results = rb.search(seq_obj)
            
            flattened = []
            for enz, sites in results.items():
                for s in sites:
                    pos = int(s)
                    flattened.append({
                        "name": str(enz),
                        "pos": pos,
                        "site": str(enz.site)
                    })
            
            flattened.sort(key=lambda x: x["pos"])
            return {"sites": flattened}

        else:
            return {"error": "Unknown tool type: {}".format(tool_type)}

    except Exception as e:
        self.retry(exc=e, countdown=5)
        return {"status": "error", "message": str(e)}


# ─────────────────────────────────────────────────────────────
# Batch FASTA Processor
# ─────────────────────────────────────────────────────────────

# CSV column order
_BATCH_FIELDNAMES = [
    "id", "description", "length", "type",
    "gc_percent", "at_percent",
    "molecular_weight_kda", "isoelectric_point",
    "gravy", "instability_index", "instability_label",
    "error",
]


@celery.task(bind=True, name="process_batch_task", max_retries=0,
             time_limit=600, soft_time_limit=540)
def process_batch_task(self, fasta_content):
    """
    Parse a multi-FASTA string and compute per-sequence metrics.

    Progress is reported via Celery task state so the Flask status
    endpoint can relay it to the frontend polling loop.

    Returns a dict with keys:
        status  – 'complete'
        total   – number of records processed
        failed  – number of records that raised an exception
        csv     – full CSV text (stored in Redis result backend)
    """
    if not _BIOPYTHON_AVAILABLE:
        return {"status": "error", "message": "BioPython is not available on this worker."}

    from io import StringIO

    # ── Parse FASTA ────────────────────────────────────────────
    records = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
    if not records:
        return {
            "status": "error",
            "message": (
                "No sequences found. Ensure the file is valid FASTA format "
                "(each record must start with a '>' header line)."
            ),
        }

    total   = len(records)
    rows    = []
    failed  = 0

    for idx, record in enumerate(records):
        raw_seq  = str(record.seq).upper().strip()
        seq_type = _detect_seq_type(raw_seq)

        row = {
            "id":                   record.id,
            "description":          record.description,
            "length":               len(raw_seq),
            "type":                 seq_type.upper(),
            "gc_percent":           "",
            "at_percent":           "",
            "molecular_weight_kda": "",
            "isoelectric_point":    "",
            "gravy":                "",
            "instability_index":    "",
            "instability_label":    "",
            "error":                "",
        }

        try:
            if seq_type in ("dna", "rna"):
                if len(raw_seq) > 0:
                    gc_count = raw_seq.count("G") + raw_seq.count("C")
                    at_count = raw_seq.count("A") + raw_seq.count("T") + raw_seq.count("U")
                    n        = len(raw_seq)
                    row["gc_percent"] = "{:.1f}".format(gc_count / n * 100)
                    row["at_percent"] = "{:.1f}".format(at_count / n * 100)
                    try:
                        seq_obj = Seq(raw_seq)
                        mol_type = "DNA" if seq_type == "dna" else "RNA"
                        mw = _bio_mw(seq_obj, seq_type=mol_type)
                        row["molecular_weight_kda"] = "{:.3f}".format(mw / 1000.0)
                    except Exception:
                        pass

            elif seq_type == "protein":
                # Truncate at first internal stop codon; map selenocysteine to Cys
                clean = raw_seq.split("*")[0].replace("U", "C")
                clean = "".join(c for c in clean if c in _PROTEIN_RESIDUES)
                if len(clean) >= 4:
                    pa = _ProteinAnalysis(clean)
                    row["molecular_weight_kda"] = "{:.3f}".format(pa.molecular_weight() / 1000.0)
                    row["isoelectric_point"]    = "{:.2f}".format(pa.isoelectric_point())
                    row["gravy"]                = "{:.3f}".format(pa.gravy())
                    inst                        = pa.instability_index()
                    row["instability_index"]    = "{:.2f}".format(inst)
                    row["instability_label"]    = "Stable" if inst < 40 else "Unstable"

        except Exception as exc:
            row["error"] = str(exc)
            failed += 1

        rows.append(row)

        # Report progress every 10 sequences and on final sequence
        processed = idx + 1
        if processed % 10 == 0 or processed == total:
            pct = int(processed / total * 100)
            self.update_state(
                state="PROGRESS",
                meta={"processed": processed, "total": total, "progress": pct},
            )

    # ── Serialize to CSV ───────────────────────────────────────
    csv_buf = _io.StringIO()
    writer  = csv.DictWriter(csv_buf, fieldnames=_BATCH_FIELDNAMES, lineterminator="\r\n")
    writer.writeheader()
    writer.writerows(rows)

    return {
        "status": "complete",
        "total":  total,
        "failed": failed,
        "csv":    csv_buf.getvalue(),
    }


if __name__ == "__main__":
    celery.start()
