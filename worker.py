import os
from celery import Celery

# Celery Configuration
CELERY_BROKER_URL = os.environ.get("REDIS_URL", "redis://localhost:6379/0")
CELERY_RESULT_BACKEND = os.environ.get("REDIS_URL", "redis://localhost:6379/0")

celery = Celery(
    "biotoolkit-worker",
    broker=CELERY_BROKER_URL,
    backend=CELERY_RESULT_BACKEND
)

# Standard Celery settings for reliability
celery.conf.update(
    task_serializer='json',
    accept_content=['json'],
    result_serializer='json',
    timezone='UTC',
    enable_utc=True,
    task_track_started=True,
    task_time_limit=300, # 5 minute execution limit
)

try:
    from Bio import Align, motifs, SeqIO, Restriction
    from Bio.Seq import Seq
    from Bio.SeqUtils import GC_skew, MeltingTemp as mt
    _BIOPYTHON_AVAILABLE = True
except ImportError:
    _BIOPYTHON_AVAILABLE = False

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

if __name__ == "__main__":
    celery.start()
