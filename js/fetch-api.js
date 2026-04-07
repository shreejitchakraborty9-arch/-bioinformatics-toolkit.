/* ============================================================
   fetch-api.js — Database Proxy Client
   All external requests (NCBI, Ensembl, UniProt, KEGG) now go
   through our Flask backend at /api/fetch/* to avoid CORS and
   securely proxy without exposing API keys in the browser.
   ============================================================ */

(function () {
  'use strict';

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
    ncbi:        (id) => buildUrl("/ncbi",    { id, db: "nucleotide" }),
    ncbiprotein: (id) => buildUrl("/ncbi",    { id, db: "protein" }),
    ncbisymbol:  (symbol, db) => buildUrl("/ncbi/search", { symbol, db: db || "nucleotide" }),
    ensembl:     (id) => buildUrl("/ensembl", { id }),
    uniprot:     (id) => buildUrl("/uniprot", { id }),
    kegg:        (id) => buildUrl("/kegg",    { id, type: "ntseq" }),
    keggprotein: (id) => buildUrl("/kegg",    { id, type: "aaseq" }),
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
  // Core Fetch Function
  // ─────────────────────────────────────────────────────────────
  window.fetchSequence = async function (db, id, targetInputId, metaId, btn) {
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

    setLoading(btn, true);

    try {
      const urlFn = ROUTES[db];
      if (!urlFn) {
        window.showToast("Unknown database selected.");
        return;
      }

      const url = urlFn(cleanId);
      const headers = {};
      const userKey = localStorage.getItem('btk_ncbi_api_key');
      if (userKey) headers['X-NCBI-API-Key'] = userKey;

      const res = await fetch(url, { headers });
      const text = await res.text();

      if (!res.ok) {
        // Server returns JSON error objects
        let errMsg = `${db.toUpperCase()} error (${res.status}).`;
        try {
          const errData = JSON.parse(text);
          if (errData.error) errMsg = errData.error;
        } catch (_) { /* plain text fallback */ }
        window.showToast(`❌ ${errMsg}`, 4500);
        return;
      }

      // Validate: must look like FASTA
      if (!text.trim() || text.trim().startsWith("<")) {
        window.showToast("❌ Unexpected response from server. ID may be incorrect.", 4000);
        return;
      }

      const { header, sequence } = parseFasta(text);

      if (!sequence) {
        window.showToast("❌ Could not parse a sequence from the response.", 4000);
        return;
      }

      // Populate the target textarea
      const input = document.getElementById(targetInputId);
      if (input) {
        input.value = sequence;
        input.dispatchEvent(new Event("input"));
      }

      // Update meta label
      const meta = metaId ? document.getElementById(metaId) : null;
      if (meta) {
        const unit = ["uniprot", "ncbiprotein", "keggprotein"].includes(db) ? "aa" : "bp";
        const truncHeader = header.length > 70 ? header.slice(0, 70) + "…" : header;
        meta.textContent = `${sequence.length.toLocaleString()} ${unit} · ${truncHeader}`;
      }

      window.showToast(`✅ ${sequence.length.toLocaleString()} residues loaded from ${db.toUpperCase()}`);

    } catch (err) {
      console.error("[BioToolkit fetch-api]", err);
      window.showToast("❌ Network error. Check your connection and try again.", 4000);
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
