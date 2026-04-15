/* ============================================================
   fetch-api.js — Database Proxy Client
   All external requests (NCBI, Ensembl, UniProt, KEGG) now go
   through our Flask backend at /api/fetch/* to avoid CORS and
   securely proxy without exposing API keys in the browser.
   ============================================================ */

(function () {
  'use strict';

  // ─────────────────────────────────────────────────────────────
  // Custom error class for sequence data integrity failures
  // ─────────────────────────────────────────────────────────────
  class SequenceIntegrityError extends Error {
    constructor(message, metadata = {}) {
      super(message);
      this.name = 'SequenceIntegrityError';
      this.metadata = metadata; // carries invalidChars, source, etc.
    }
  }

  // ─────────────────────────────────────────────────────────────
  // API Route Map (all routed through our server proxy)
  // ─────────────────────────────────────────────────────────────
  const API_BASE = "/api/fetch";

  function buildUrl(endpoint, params) {
    const url = new URL(API_BASE + endpoint, window.location.origin);
    Object.entries(params).forEach(([k, v]) => {
      if (v !== undefined && v !== null && v !== "") url.searchParams.set(k, v);
    });
    return url.toString();
  }

  const ROUTES = {
    ncbi:        (id) => buildUrl("", { id, db: "ncbi" }),
    ncbiprotein: (id) => buildUrl("", { id, db: "ncbiprotein" }),
    ncbisymbol:  (symbol, db) => buildUrl("", { id: symbol, db: "ncbisymbol", species: "Homo sapiens" }),
    ensembl:     (id) => buildUrl("", { id, db: "ensembl" }),
    uniprot:     (id) => buildUrl("", { id, db: "uniprot" }),
    kegg:        (id) => buildUrl("", { id, db: "kegg", type: "ntseq" }),
    keggprotein: (id) => buildUrl("", { id, db: "kegg", type: "aaseq" }),
  };

  // ─────────────────────────────────────────────────────────────
  // Placeholder hints per database
  // ─────────────────────────────────────────────────────────────
  const PLACEHOLDERS = {
    ncbi:        "e.g. NM_007294  (BRCA1)",
    ncbiprotein: "e.g. NP_000537.3  (TP53)",
    uniprot:     "e.g. P04637  (TP53)",
    ensembl:     "e.g. ENSG00000139618 or BRCA2",
    kegg:        "e.g. hsa:7157  (TP53)",
  };

  // ─────────────────────────────────────────────────────────────
  // FASTA Parser
  // ─────────────────────────────────────────────────────────────
  function parseFasta(text) {
    const lines = text.trim().split("\n");
    let header = "";
    const seqLines = [];
    lines.forEach(line => {
      if (line.startsWith(">")) { header = line.slice(1).trim(); }
      else if (line.trim())     { seqLines.push(line.trim()); }
    });
    return { header, sequence: seqLines.join("") };
  }

  // ─────────────────────────────────────────────────────────────
  // Loading State
  // ─────────────────────────────────────────────────────────────
  function setLoading(btn, on) {
    if (!btn) return;
    if (on) {
      btn.dataset.origText = btn.innerHTML;
      btn.innerHTML = '<span class="fetch-spinner"></span> Fetching…';
      btn.disabled = true;
    } else {
      btn.innerHTML = btn.dataset.origText || "Fetch";
      btn.disabled = false;
    }
  }

  // ─────────────────────────────────────────────────────────────
  // Client-side rate limit (backup; server also enforces)
  // ─────────────────────────────────────────────────────────────
  let _lastFetch = 0;
  const MIN_INTERVAL_MS = 600;

  // ─────────────────────────────────────────────────────────────
  // Alphabet Validation (Issue #2 — character-level integrity)
  // ─────────────────────────────────────────────────────────────
  const IUPAC_DNA_CANONICAL = new Set(['A','T','G','C']);
  const IUPAC_RNA_CANONICAL = new Set(['A','U','G','C']);
  const IUPAC_AMBIGUOUS     = new Set(['R','Y','S','W','K','M','B','D','H','V','N','-']);

  /**
   * Validates raw sequence alphabet after NCBI fetch.
   * Pure function — no side effects, no DOM access, no window references.
   * @param {string} sequence - Raw sequence string (post-FASTA header strip)
   * @param {string} sequenceType - 'DNA' | 'RNA' | 'AUTO'
   * @returns {object} ValidationResult
   */
  function validateSequenceAlphabet(sequence, sequenceType = 'AUTO') {
    // Gap 3C: Replace standard and non-standard whitespace before scanning
    // Gap 3B: Multi-line FASTA header stripping
    const clean = sequence
      .replace(/^>.*$/mg, '')                 // strip all FASTA header lines
      .replace(/[\s\u00A0\u200B\uFEFF]/g, '') // strip standard and non-standard whitespace
      .toUpperCase();

    // Gap 3A: Empty string guard (prevents NaN in canonicalPercent)
    if (clean.length === 0) {
      return {
        isValid: false,
        sequenceType: sequenceType,
        hasAmbiguous: false,
        ambiguousCount: 0,
        ambiguousTypes: [],
        invalidChars: [],
        invalidCount: 0,
        canonicalPercent: 0,
        message: 'Sequence rejected: empty sequence after header and whitespace removal.'
      };
    }

    // AUTO-detect sequence type
    const hasT = /T/.test(clean);
    const hasU = /U/.test(clean);
    let detectedType = sequenceType;
    if (sequenceType === 'AUTO') {
      if (hasT && hasU) {
        return {
          isValid: false,
          sequenceType: 'UNKNOWN',
          hasAmbiguous: false,
          ambiguousCount: 0,
          ambiguousTypes: [],
          invalidChars: ['T', 'U'],
          invalidCount: 2,
          canonicalPercent: 0,
          message: 'Sequence rejected: contains both T and U — cannot determine nucleotide type.'
        };
      }
      detectedType = hasU ? 'RNA' : 'DNA';
    }

    const canonical = detectedType === 'RNA' ? IUPAC_RNA_CANONICAL : IUPAC_DNA_CANONICAL;

    // Single-pass scanner
    let canonicalCount = 0;
    let ambiguousCount = 0;
    let invalidCount = 0;
    const ambiguousTypesSet = new Set();
    const invalidCharsSet = new Set();

    for (let i = 0; i < clean.length; i++) {
      const ch = clean[i];
      if (canonical.has(ch)) {
        canonicalCount++;
      } else if (IUPAC_AMBIGUOUS.has(ch)) {
        ambiguousTypesSet.add(ch);
        ambiguousCount++;
      } else {
        invalidCharsSet.add(ch);
        invalidCount++;
      }
    }

    const isValid = invalidCount === 0;
    const hasAmbiguous = ambiguousCount > 0;
    const ambiguousTypes = Array.from(ambiguousTypesSet).sort();
    const invalidChars  = Array.from(invalidCharsSet).sort();
    const canonicalPercent = parseFloat(((canonicalCount / clean.length) * 100).toFixed(1));

    let message;
    if (!isValid) {
      message = `Sequence rejected: ${invalidCount} invalid character(s) found [${invalidChars.join(', ')}]. Fetch aborted.`;
    } else if (hasAmbiguous) {
      message = `Sequence contains ${ambiguousCount} ambiguous IUPAC base(s) [${ambiguousTypes.join(', ')}]. Tm and MW results may be approximate.`;
    } else {
      message = `Sequence validated: ${canonicalPercent}% canonical ${detectedType} bases.`;
    }

    return { isValid, sequenceType: detectedType, hasAmbiguous, ambiguousCount,
             ambiguousTypes, invalidChars, invalidCount, canonicalPercent, message };
  }

  // ─────────────────────────────────────────────────────────────
  // User-Friendly Error Handling (Issue #5)
  // ─────────────────────────────────────────────────────────────
  const API_ERROR_MESSAGES = {
      'NETWORK_TIMEOUT':     'Database server not responding. Check your connection and try again.',
      'INVALID_ACCESSION':   'Invalid accession format. Use a standard ID like NM_007294 or P04637.',
      'ACCESSION_NOT_FOUND': 'Record not found (404). Verify the accession ID exists in the selected database.',
      'DB_SERVER_ERROR':     'Database server error (500). Retrying once…',
      'MALFORMED_DATA':      'Sequence data is corrupted or in an unsupported format.',
      'RATE_LIMITED':        'Rate limit exceeded. Wait a moment and try again.',
      'GENERIC_API_ERROR':   'NCBI/Database service error. Please try again shortly.'
  };

  /**
   * Internal error classifier to prevent raw JSON/stack traces in UI.
   */
  function handleAPIError(err, context = '', db = '') {
      let errorType = 'GENERIC_API_ERROR';
      const msg = (typeof err === 'string' ? err : (err.message || '')).toLowerCase();

      if (msg.includes('timeout') || msg.includes('network error') || msg.includes('aborted')) {
          errorType = 'NETWORK_TIMEOUT';
      } else if (msg.includes('404') || msg.includes('not found') || msg.includes('incorrect') || msg.includes('no sequence')) {
          errorType = 'ACCESSION_NOT_FOUND';
      } else if (msg.includes('500') || msg.includes('server error') || msg.includes('internal')) {
          errorType = 'DB_SERVER_ERROR';
      } else if (msg.includes('429') || msg.includes('rate limit') || msg.includes('too many')) {
          errorType = 'RATE_LIMITED';
      } else if (msg.includes('400') || msg.includes('invalid') || msg.includes('accession')) {
          errorType = 'INVALID_ACCESSION';
      } else if (msg.includes('fasta') || msg.includes('parse')) {
          errorType = 'MALFORMED_DATA';
      }

      const userMessage = API_ERROR_MESSAGES[errorType];

      console.error(`[BioKit API Error] DB: ${db} | Type: ${errorType} | Context: ${context}`, {
          originalError: err,
          timestamp: new Date().toISOString()
      });

      window.showToast(`❌ ${userMessage}`, 5000);
      return { success: false, type: errorType, message: userMessage };
  }

  // ─────────────────────────────────────────────────────────────
  // Core Fetch Function
  // ─────────────────────────────────────────────────────────────
  window.fetchSequence = async function (db, id, targetInputId, metaId, btn) {
    // Gap 1: Reset warnings at start and initialize recursive fetch guard
    const fetchId = Symbol(); 
    window.BioKit._activeFetchId = fetchId;
    window.BioKit.pendingWarnings = [];

    const cleanId = (id || "").trim();
    if (!cleanId) {
      window.showToast("Please enter a gene / accession ID.");
      return;
    }

    const now = Date.now();
    if (now - _lastFetch < MIN_INTERVAL_MS) {
      window.showToast("⏳ Please wait a moment before fetching again.", 2500);
      return;
    }
    _lastFetch = now;

    // Context-aware protein mode: if the caller DB is a protein source,
    // skip the nucleotide-based AUTO alphabet detection entirely.
    const isProteinDB = ['uniprot', 'ncbiprotein', 'keggprotein'].includes(db);

    setLoading(btn, true);

    // Internal single-retry helper for 500-class server errors.
    async function attemptFetch(attempt) {
      const urlFn = ROUTES[db];
      if (!urlFn) throw new Error('Unknown database: ' + db);
      const url = urlFn(cleanId);
      const headers = {};
      const userKey = localStorage.getItem('btk_ncbi_api_key');
      if (userKey) headers['X-NCBI-API-Key'] = userKey;

      const res = await fetch(url, { headers });
      const text = await res.text();

      if (!res.ok) {
        if (res.status >= 500 && attempt === 1) {
          window.showToast('⚠️ Database server error (500). Retrying in 3\u202fs…', 3500);
          await new Promise(r => setTimeout(r, 3000));
          return attemptFetch(2);       // single retry
        }
        if (res.status === 404) {
          throw new Error(`404 Not found: the ID '${cleanId}' does not exist in ${db.toUpperCase()}.`);
        }
        if (res.status >= 500) {
          throw new Error(`500 Server error from ${db.toUpperCase()} — database may be temporarily unavailable.`);
        }
        throw new Error(`${db.toUpperCase()} response: ${res.status} ${text}`);
      }
      return { res, text };
    }

    try {
      const { text } = await attemptFetch(1);

      // Parse JSON or fallback to FASTA
      let header = "";
      let sequence = "";
      let metadata = null;
      let transcripts = [];
      let uniprotName = null;   // confirmed UniProt canonical name

      try {
        const jsonData = JSON.parse(text);
        if (jsonData.sequence) {
          sequence = jsonData.sequence;
          metadata = jsonData.metadata;
          transcripts = jsonData.transcripts || [metadata];
          header = metadata ? (metadata.description || "") : "";

          // UniProt canonical label extraction (addresses isoform confusion)
          if (db === 'uniprot' && jsonData.proteinDescription?.recommendedName?.fullName?.value) {
            uniprotName = jsonData.proteinDescription.recommendedName.fullName.value;
          }
        }
      } catch (e) {
        // Fallback: Validate legacy FASTA
        if (!text.trim() || text.trim().startsWith('<')) {
          window.showToast('❌ Unexpected response from server. ID may be incorrect.', 4000);
          return;
        }
        const parsed = parseFasta(text);
        header = parsed.header;
        sequence = parsed.sequence;
      }

      if (!sequence) {
        window.showToast("❌ Could not parse a sequence from the response.", 4000);
        return;
      }

      // Extract expected length from header (if available) e.g., "(7088 bp)"
      const lengthMatch = header.match(/(\d+)\s*bp/);
      const expectedLength = lengthMatch ? parseInt(lengthMatch[1]) : null;

      if (expectedLength && sequence.length < expectedLength * 0.95) {
          const errMsg = `Sequence incomplete: received ${sequence.length.toLocaleString()} bp, expected ~${expectedLength.toLocaleString()} bp.`;
          window.showToast(`❌ ${errMsg}`, 6000);
          console.error(`[CRITICAL] Sequence integrity: ${errMsg}`);
          const activePanel = document.querySelector('.tool-panel.active');
          if (activePanel && window.showValidationWarning) {
              window.showValidationWarning(activePanel.id, {
                  valid: false, msg: "Data Integrity Error: Partial Sequence Loaded", suggestion: errMsg
              });
          }
          return;
      }

      if (expectedLength && sequence.length < expectedLength) {
          const cap = ((sequence.length / expectedLength) * 100).toFixed(1);
          window.showToast(`⚠️ Sequence loaded at ${cap}% capacity`, 5000);
          console.warn(`[Data Integrity] Minor truncation: ${sequence.length}/${expectedLength} bp`);
      }

      // ── Context-aware alphabet integrity check ────────────────────────────
      // Protein databases bypass nucleotide AUTO detection, which previously
      // hard-rejected amino-acid sequences containing L, E, F, etc. as "invalid DNA".
      let alphaResult;
      if (isProteinDB) {
        const proteinValid = /^[ACDEFGHIKLMNPQRSTVWYXUZ*\s\r\n]+$/i.test(sequence);
        alphaResult = {
          isValid: proteinValid,
          sequenceType: 'PROTEIN',
          hasAmbiguous: false,
          ambiguousCount: 0,
          ambiguousTypes: [],
          invalidChars: [],
          canonicalPercent: 100,
          message: proteinValid
            ? 'Protein sequence validated.'
            : 'Sequence contains characters outside the IUPAC protein alphabet.'
        };
      } else {
        alphaResult = validateSequenceAlphabet(sequence, 'AUTO');
      }

      if (!alphaResult.isValid) {
        throw new SequenceIntegrityError(
          `Alphabet validation failed: ${alphaResult.message}`,
          { invalidChars: alphaResult.invalidChars, source: 'NCBI_FETCH' }
        );
      }

      if (alphaResult.hasAmbiguous) {
        if (window.BioKit._activeFetchId === fetchId) {
          window.BioKit.pendingWarnings = window.BioKit.pendingWarnings || [];
          window.BioKit.pendingWarnings.push({
            type: 'AMBIGUOUS_BASES',
            message: alphaResult.message,
            ambiguousCount: alphaResult.ambiguousCount,
            ambiguousTypes: alphaResult.ambiguousTypes,
            canonicalPercent: alphaResult.canonicalPercent
          });
        }
      }

      // Race-condition storage guard
      if (window.BioKit._activeFetchId !== fetchId) { return; }

      const normalizedSequence = sequence.toUpperCase().replace(/[\s\r\n]/g, '');
      window.BioKit.activeAnalysis = {
        sequence: normalizedSequence,
        sequenceType: alphaResult.sequenceType,
        hasAmbiguousBases: alphaResult.hasAmbiguous,
        metadata,
        transcripts
      };

      // Populate the target textarea
      const input = document.getElementById(targetInputId);
      if (input) {
        input.value = sequence;
        input.dispatchEvent(new Event("input"));
      }

      // Show confirmed UniProt canonical name label (isoform verification)
      if (uniprotName) {
        const labelEl = document.getElementById('uniprotNameLabel');
        if (labelEl) {
          labelEl.textContent = `✔ Canonical: ${uniprotName}`;
          labelEl.classList.remove('hidden');
        }
      }

      // Update meta label
      const meta = metaId ? document.getElementById(metaId) : null;
      if (meta) {
        const unit = ["uniprot", "ncbiprotein", "keggprotein"].includes(db) ? "aa" : "bp";
        const truncHeader = header.length > 70 ? header.slice(0, 70) + "…" : header;
        let info = `${sequence.length.toLocaleString()} ${unit} · ${truncHeader}`;
        if (metadata && metadata.strand) info += ` [Strand: ${metadata.strand}]`;
        meta.textContent = info;
      }

      window.showToast(`✅ ${sequence.length.toLocaleString()} residues loaded from ${db.toUpperCase()}`);

    } catch (err) {
      if (err instanceof SequenceIntegrityError) {
        window.showToast(`❌ ${err.message}`, 6000);
        console.error('[SequenceIntegrityError]', err.metadata);
        return;
      }
      handleAPIError(err, 'Fetch Exception', db);
    } finally {
      setLoading(btn, false);
    }
  };


  // ─────────────────────────────────────────────────────────────
  // Wire up all fetch widgets on DOM ready
  // ─────────────────────────────────────────────────────────────
  document.addEventListener("DOMContentLoaded", () => {
    // Fetch buttons
    document.querySelectorAll(".fetch-btn").forEach(btn => {
      btn.addEventListener("click", () => {
        const widget = btn.closest(".fetch-widget");
        if (!widget) return;

        let db = widget.querySelector(".fetch-db-select")?.value;
        const id = widget.querySelector(".fetch-id-input")?.value || "";
        const targetId = widget.dataset.target;
        const metaId = widget.dataset.meta || null;

        // Remap KEGG to protein variant if inside the protein panel
        if (db === "kegg") {
          const panel = widget.closest(".tool-panel");
          if (panel && (panel.id === "panel-protein" || targetId === "proteinInput")) {
            db = "keggprotein";
          }
        }

        window.fetchSequence(db, id, targetId, metaId, btn);
      });
    });

    // Update placeholder hints when database changes
    document.querySelectorAll(".fetch-db-select").forEach(sel => {
      sel.addEventListener("change", () => {
        const widget = sel.closest(".fetch-widget");
        const idInput = widget?.querySelector(".fetch-id-input");
        if (idInput) {
          idInput.placeholder = PLACEHOLDERS[sel.value] || "Enter ID";
        }
      });
    });
  });

})();
