import os
import io as _io
import csv
import math
from celery import Celery

# Celery Configuration
CELERY_BROKER_URL     = os.environ.get("REDIS_URL", "redis://localhost:6379/0")
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

# Standard 20 amino acids (selenocysteine U is mapped to C before analysis)
_PROTEIN_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Full IUPAC extended nucleotide alphabets.
# Without this, sequences with R/Y/S/W/K/M/B/D/H/V (common in PCR products,
# primer consensus, and population variant data) misclassify as protein.
_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_IUPAC_RNA = frozenset("ACGURYSWKMBDHVN")


def _detect_seq_type(seq_str):
    """Return 'dna', 'rna', or 'protein' using IUPAC-aware composition."""
    s   = seq_str.upper()
    tot = len(s)
    if tot == 0:
        return "unknown"

    dna_hits = sum(1 for c in s if c in _IUPAC_DNA)
    rna_hits = sum(1 for c in s if c in _IUPAC_RNA)

    if "U" in s and rna_hits / tot > 0.85:
        return "rna"
    if dna_hits / tot > 0.85:
        return "dna"
    return "protein"


def _cpg_obs_exp(seq_str):
    """CpG observed/expected ratio (Gardiner-Garden & Frommer, J Mol Biol 1987).

    Formula: (count_CpG * length) / (count_C * count_G)
    Values >= 0.60 are indicative of CpG islands (unmethylated promoter regions).
    Returns None when C or G count is zero (undefined).
    """
    s = seq_str.upper()
    n = len(s)
    if n < 2:
        return None
    c_count = s.count("C")
    g_count = s.count("G")
    if c_count == 0 or g_count == 0:
        return None
    cg_count = sum(1 for i in range(n - 1) if s[i] == "C" and s[i + 1] == "G")
    return round((cg_count * n) / (c_count * g_count), 4)


def _shannon_entropy(seq_str):
    """Normalized Shannon entropy of sequence composition (0-1 scale).

    Normalized by log2(alphabet_size) so 1.0 = maximally complex.
    Values < 0.70 flag low-complexity or repetitive regions (DUST/SEG territory).
    """
    s = seq_str.upper()
    n = len(s)
    if n == 0:
        return 0.0
    counts = {}
    for c in s:
        counts[c] = counts.get(c, 0) + 1
    alphabet_size = len(counts)
    if alphabet_size <= 1:
        return 0.0
    ent = -sum((v / n) * math.log2(v / n) for v in counts.values() if v > 0)
    return round(ent / math.log2(alphabet_size), 4)


# ─────────────────────────────────────────────────────────────
# General tool task (alignment, GC-skew, restriction map)
# ─────────────────────────────────────────────────────────────

@celery.task(bind=True, max_retries=3)
def process_sequence_task(self, tool_type, params):
    """Scientific-grade bioinformatics worker task."""
    if not _BIOPYTHON_AVAILABLE:
        return {"status": "error", "message": "Biopython is not available on this worker."}

    print("Executing {} task with params: {}".format(tool_type, list(params.keys())))

    try:
        t_type = tool_type.lower().replace('-', '_')

        # 1. Alignment (Global/Local — Needleman-Wunsch / Smith-Waterman)
        if t_type in ("alignment", "align_global", "align_local"):
            seq1 = params.get("seq1", params.get("a", ""))
            seq2 = params.get("seq2", params.get("b", ""))
            mode = "local" if t_type == "align_local" else "global"

            aligner = Align.PairwiseAligner()
            aligner.mode = mode
            if "match"    in params: aligner.match_score    = float(params["match"])
            if "mismatch" in params: aligner.mismatch_score = float(params["mismatch"])
            if "gap"      in params: aligner.open_gap_score  = float(params["gap"])

            alignments = aligner.align(seq1, seq2)
            if not alignments:
                return {"score": 0, "alnA": "", "alnB": "", "mid": ""}

            best = alignments[0]
            alnA = str(best[0])
            alnB = str(best[1])
            mid  = "".join(
                "|" if a == b else (" " if a == "-" or b == "-" else "·")
                for a, b in zip(alnA, alnB)
            )
            return {"score": float(best.score), "alnA": alnA, "alnB": alnB, "mid": mid}

        # 2. GC-Skew Analysis
        elif t_type in ("gc_skew", "gc-skew"):
            sequence = params.get("sequence", "")
            window   = int(params.get("window", 100))
            if not sequence:
                return {"error": "No sequence"}
            skew = GC_skew(sequence, window)
            return {"skew": skew}

        # 3. Restriction Enzyme Mapping
        elif t_type in ("restriction_search", "restriction-search", "restriction_map"):
            sequence = params.get("sequence", "")
            enzymes  = params.get("enzymes", [])
            if not sequence:
                return {"error": "No sequence"}
            seq_obj       = Seq(sequence)
            valid_enzymes = [e for e in enzymes if hasattr(Restriction, e)]
            rb            = Restriction.RestrictionBatch(valid_enzymes) if valid_enzymes else Restriction.CommOnly
            results       = rb.search(seq_obj)

            flattened = []
            for enz, sites in results.items():
                for s in sites:
                    flattened.append({"name": str(enz), "pos": int(s), "site": str(enz.site)})
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

# Column order is fixed for reproducibility; downstream scripts should refer
# to column names, not positions.
_BATCH_FIELDNAMES = [
    # Identity
    "id", "description", "length", "type",
    # Nucleotide base composition (DNA/RNA only; canonical bases)
    "gc_percent", "at_percent",
    "a_percent", "t_percent", "g_percent", "c_percent",
    # DNA-specific methylation context
    "cpg_obs_exp",
    # Thermodynamics
    "melting_temp_c",
    # Sequence complexity
    "shannon_entropy",
    # Molecular weight (shared)
    "molecular_weight_kda",
    # Protein physicochemical properties (protein only)
    "isoelectric_point",
    "gravy",
    "instability_index", "instability_label",
    "aromaticity",
    "aliphatic_index",
    "net_charge_ph74",
    "extinction_coeff_reduced",
    "extinction_coeff_oxidized",
    # Per-record quality annotations
    "warning", "error",
]

_VALID_MODES = frozenset(("all", "nucleotide", "protein", "quick"))


@celery.task(bind=True, name="process_batch_task", max_retries=0,
             time_limit=600, soft_time_limit=540)
def process_batch_task(self, fasta_content, mode="all"):
    """
    Parse a multi-FASTA string and compute per-sequence metrics.

    Parameters
    ----------
    fasta_content : str
        Raw FASTA text (already validated by the Flask route).
    mode : str
        'all'        -- full metrics for nucleotides and proteins
        'nucleotide' -- nucleotide composition + Tm + entropy only
        'protein'    -- ProtParam metrics only
        'quick'      -- id, description, length, type only (fast scan)

    Progress is reported via Celery task state so the Flask status
    endpoint can relay it to the frontend polling loop.

    Returns a dict with keys:
        status  - 'complete'
        total   - number of records processed
        failed  - number of records that raised an exception
        mode    - analysis mode used
        csv     - full CSV text (stored in Redis result backend)
    """
    if not _BIOPYTHON_AVAILABLE:
        return {"status": "error", "message": "BioPython is not available on this worker."}

    if mode not in _VALID_MODES:
        mode = "all"

    from io import StringIO

    # ── Parse FASTA ──────────────────────────────────────────────
    records = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
    if not records:
        return {
            "status": "error",
            "message": (
                "No sequences could be parsed. Verify the file uses standard FASTA format: "
                "each record must begin with a '>' header line immediately followed by "
                "sequence data on the next line(s)."
            ),
        }

    # Duplicate ID detection — flag per record but do not abort the job.
    # Duplicate IDs cause ambiguity in downstream analyses (BLAST, alignment tools).
    _id_positions = {}
    for i, rec in enumerate(records):
        _id_positions.setdefault(rec.id, []).append(i)
    duplicate_ids = {rid for rid, positions in _id_positions.items() if len(positions) > 1}

    do_nuc  = mode in ("all", "nucleotide")
    do_prot = mode in ("all", "protein")

    total  = len(records)
    rows   = []
    failed = 0

    for idx, record in enumerate(records):
        raw_seq  = str(record.seq).upper().strip()
        seq_type = _detect_seq_type(raw_seq)
        warnings = []

        if record.id in duplicate_ids:
            warnings.append("Duplicate sequence ID — results may be ambiguous in downstream tools")

        # Initialise every column to empty string so the CSV is always rectangular
        row = {f: "" for f in _BATCH_FIELDNAMES}
        row["id"]          = record.id
        row["description"] = record.description
        row["length"]      = len(raw_seq)
        row["type"]        = seq_type.upper()

        try:
            if seq_type in ("dna", "rna"):
                n = len(raw_seq)

                if do_nuc and n > 0:
                    # Canonical base counts only; degenerate codes are excluded
                    # from per-base percentages to prevent misleading fractions.
                    a = raw_seq.count("A")
                    g = raw_seq.count("G")
                    c = raw_seq.count("C")
                    # T and U both represent the same pairing partner in their alphabets
                    t = raw_seq.count("T") + raw_seq.count("U")
                    gc_count = g + c
                    at_count = a + t

                    row["gc_percent"] = "{:.1f}".format(gc_count / n * 100)
                    row["at_percent"] = "{:.1f}".format(at_count / n * 100)
                    row["a_percent"]  = "{:.1f}".format(a / n * 100)
                    row["t_percent"]  = "{:.1f}".format(t / n * 100)
                    row["g_percent"]  = "{:.1f}".format(g / n * 100)
                    row["c_percent"]  = "{:.1f}".format(c / n * 100)

                    # CpG O/E ratio — biologically meaningful for DNA only
                    # (methylation-sensitive promoter/island context)
                    if seq_type == "dna":
                        cpg = _cpg_obs_exp(raw_seq)
                        if cpg is not None:
                            row["cpg_obs_exp"] = "{:.4f}".format(cpg)

                    # Tm_Wallace: 4*(G+C) + 2*(A+T) using canonical base counts.
                    # Degenerate positions are excluded from the tally (conservative).
                    # Calibrated for short oligonucleotides; flagged for n > 50.
                    if n >= 4:
                        tm_val = float(4 * gc_count + 2 * at_count)
                        row["melting_temp_c"] = "{:.1f}".format(tm_val)
                        if n > 50:
                            warnings.append(
                                "Tm_Wallace is calibrated for oligonucleotides "
                                "(<= 50 nt) -- value is a rough estimate for longer sequences"
                            )

                    row["shannon_entropy"] = str(_shannon_entropy(raw_seq))

                    # Molecular weight via BioPython
                    try:
                        mol_type = "DNA" if seq_type == "dna" else "RNA"
                        mw = _bio_mw(Seq(raw_seq), seq_type=mol_type)
                        row["molecular_weight_kda"] = "{:.3f}".format(mw / 1000.0)
                    except Exception:
                        pass

            elif seq_type == "protein":
                # Truncate at first internal stop codon; map selenocysteine to Cys
                clean = raw_seq.split("*")[0].replace("U", "C")

                # Non-standard residues are removed with a warning rather than
                # aborting, preserving traceability for ambiguous sequences.
                non_standard = set(clean) - _PROTEIN_RESIDUES
                if non_standard:
                    warnings.append(
                        "Non-standard residues removed before ProtParam: {}".format(
                            ", ".join(sorted(non_standard))
                        )
                    )
                clean = "".join(ch for ch in clean if ch in _PROTEIN_RESIDUES)

                if len(clean) < 4:
                    warnings.append(
                        "Sequence too short ({} aa after cleaning) for reliable "
                        "physicochemical calculations -- BioPython ProtParam requires >= 4 aa".format(
                            len(clean)
                        )
                    )
                elif do_prot:
                    pa = _ProteinAnalysis(clean)

                    row["molecular_weight_kda"]     = "{:.3f}".format(pa.molecular_weight() / 1000.0)
                    row["isoelectric_point"]         = "{:.2f}".format(pa.isoelectric_point())
                    row["gravy"]                     = "{:.3f}".format(pa.gravy())
                    inst                             = pa.instability_index()
                    row["instability_index"]         = "{:.2f}".format(inst)
                    row["instability_label"]         = "Stable" if inst < 40 else "Unstable"
                    row["aromaticity"]               = "{:.4f}".format(pa.aromaticity())
                    row["aliphatic_index"]           = "{:.2f}".format(pa.aliphatic_index())
                    row["net_charge_ph74"]           = "{:.2f}".format(pa.charge_at_pH(7.4))
                    ec_red, ec_ox                    = pa.molar_extinction_coefficient()
                    row["extinction_coeff_reduced"]  = str(ec_red)
                    row["extinction_coeff_oxidized"] = str(ec_ox)
                    row["shannon_entropy"]           = str(_shannon_entropy(clean))

        except Exception as exc:
            row["error"] = str(exc)
            failed += 1

        row["warning"] = "; ".join(warnings)
        rows.append(row)

        # Report progress every 10 sequences and on the final sequence
        processed = idx + 1
        if processed % 10 == 0 or processed == total:
            pct = int(processed / total * 100)
            self.update_state(
                state="PROGRESS",
                meta={"processed": processed, "total": total, "progress": pct},
            )

    # ── Serialize to CSV ─────────────────────────────────────────
    csv_buf = _io.StringIO()
    writer  = csv.DictWriter(csv_buf, fieldnames=_BATCH_FIELDNAMES, lineterminator="\r\n")
    writer.writeheader()
    writer.writerows(rows)

    return {
        "status": "complete",
        "total":  total,
        "failed": failed,
        "mode":   mode,
        "csv":    csv_buf.getvalue(),
    }


if __name__ == "__main__":
    celery.start()
