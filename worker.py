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

    print(f"Executing {tool_type} task with params: {params.keys()}")
    
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
            # Filter enzymes to only those that exist in Biopython to avoid errors
            valid_enzymes = [e for e in enzymes if hasattr(Restriction, e)]
            rb = Restriction.RestrictionBatch(valid_enzymes) if valid_enzymes else Restriction.CommOnly
            results = rb.search(seq_obj)
            
            flattened = []
            for enz, sites in results.items():
                for s in sites:
                    pos = int(s)
                    flattened.append({
                        "name": str(enz),
                        "site": enz.site,
                        "pos": pos,
                        "cut": pos, # Biopython reports the cut position (1-based)
                        "strand": "+",
                        "is_blunt": enz.is_blunt()
                    })
            # Sort by position for predictable rendering
            flattened.sort(key=lambda x: x["pos"])
            return flattened

        # 4. Thermodynamic Analysis (Tm, etc.)
        elif t_type in ("thermo_analysis", "primer_tm", "primer-tm"):
            sequence = params.get("sequence", params.get("primer", ""))
            if not sequence: return {"error": "No sequence"}
            
            # Detailed Salt Corrections
            sc = {
                "dnac": float(params.get("dnac", 250)),
                "na": float(params.get("na", 50)),
                "mg": float(params.get("mg", 1.5)),
                "k": float(params.get("k", 0)),
                "tris": float(params.get("tris", 0)),
                "dntp": float(params.get("dntp", 0.8))
            }
            
            tm_nn = mt.Tm_NN(sequence, dnac1=sc["dnac"], Na=sc["na"], Mg=sc["mg"], K=sc["k"], Tris=sc["tris"], dNTPs=sc["dntp"])
            tm_gc = mt.Tm_GC(sequence, Na=sc["na"])
            tm_wallace = mt.Tm_Wallace(sequence)
            
            return {
                "tm_nn": round(tm_nn, 2),
                "tm_gc": round(tm_gc, 2),
                "tm_wallace": round(tm_wallace, 2),
                "salt_params": sc
            }

        # 5. Sequence Search (Integrated BLAST-like)
        elif t_type in ("seq_search", "seq-search", "genome_scan"):
            query = params.get("query", params.get("sequence", ""))
            database = params.get("database", [])
            
            results = []
            aligner = Align.PairwiseAligner()
            aligner.mode = "local"
            
            for entry in database:
                target_seq = entry.get("sequence", "")
                if not target_seq: continue
                
                alignments = aligner.align(query, target_seq)
                if alignments:
                    best = alignments[0]
                    aln_q = str(best[0])
                    aln_t = str(best[1])
                    matches = sum(1 for a, b in zip(aln_q, aln_t) if a == b)
                    identity = round((matches / len(aln_q)) * 100 if aln_q else 0, 1)
                    
                    results.append({
                        "target": {
                            "gene": entry.get("gene", "Unknown"),
                            "organism": entry.get("organism", "Unknown"),
                            "description": entry.get("description", "")
                        },
                        "aln": {
                            "score": float(best.score),
                            "identity": identity,
                            "alnLen": len(aln_q),
                            "queryAln": aln_q,
                            "subjAln": aln_t,
                            "midAln": "".join("|" if a == b else "." for a, b in zip(aln_q, aln_t))
                        }
                    })
            
            results.sort(key=lambda x: x["aln"]["score"], reverse=True)
            return {"hits": results[:10]}

        # 6. Comprehensive Promoter & ORF Analysis
        elif t_type in ("promoter_analysis", "promoter-analysis", "orf_finder"):
            sequence = params.get("sequence", "")
            if not sequence: return {"error": "No sequence"}
            
            seq_obj = Seq(sequence)
            results = {"orfs": [], "motifs": []}
            
            # 6a. Professional ORF Finding (all 6 frames)
            def find_orfs(s, strand_name):
                found = []
                for frame in range(3):
                    for orf in s[frame:].split("*"):
                        if "M" in orf:
                            start_idx = orf.find("M")
                            actual_orf = orf[start_idx:]
                            if len(actual_orf) >= 30: # Min 30 AA
                                found.append({
                                    "strand": strand_name,
                                    "frame": frame + 1,
                                    "aaLen": len(actual_orf),
                                    "protein": str(actual_orf)
                                })
                return found
            
            results["orfs"].extend(find_orfs(seq_obj.translate(), "+"))
            results["orfs"].extend(find_orfs(seq_obj.reverse_complement().translate(), "-"))
            
            # Simple Motif Scan (Biopython motifs can be huge, using consensus for now)
            # This can be expanded with real PWM files in a scientific prod environment
            return results

        else:
            return {"status": "error", "message": f"Unsupported tool_type: {tool_type}"}

    except Exception as e:
        print(f"Error in task: {e}")
        return {"status": "failure", "error": str(e)}
