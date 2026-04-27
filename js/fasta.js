/* ============================================================
   fasta.js — FASTA Parser (Fixed DOM IDs + Drag-Drop)

   ARCHITECTURE CHANGE (Refactor):
   - parseFASTA() has been REMOVED from the main thread.
   - Parsing is now delegated to worker.js (PARSE_FASTA task)
     via window.WorkerManager, keeping the UI thread free.
   - An unmount() hook is registered with window.ViewManager
     so that navigating away from this panel physically removes
     heavy DOM nodes and nulls the result dataset.
   - All biological output logic (displayResults, stat cards,
     sequence cards) is byte-for-byte identical to the original.
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // Track the last parsed dataset so unmount() can null it
  let _lastRecords = null;

  // ── Unmount Hook (Garbage Collection) ──────────────────
  // Registered with ViewManager. Executes when user navigates
  // away from the FASTA panel. Physically removes DOM nodes and
  // releases the reference to parsed record objects.
  function unmount() {
    // 1. Destroy heavy DOM — sequences container can hold thousands of nodes
    const seqContainer = document.getElementById('fastaSequences');
    if (seqContainer) seqContainer.innerHTML = '';

    const statRow = document.getElementById('fastaStatRow');
    if (statRow) statRow.innerHTML = '';

    // 2. Hide results panel
    document.getElementById('fastaResults')?.classList.add('hidden');

    // 3. Null the JS reference so GC can reclaim the record objects
    _lastRecords = null;
  }

  // Register with ViewManager as soon as module loads
  if (window.ViewManager) {
    window.ViewManager.registerUnmount('fasta', unmount);
  } else {
    // ViewManager not yet available — defer to next tick (app.js loads first, so this is a safety net)
    window.addEventListener('load', () => {
      window.ViewManager?.registerUnmount('fasta', unmount);
    });
  }

  // ── Run (Web Worker Dispatch) ──────────────────────────
  // parseFASTA logic has moved to worker.js:runFastaParse().
  // This function now only orchestrates the worker call and
  // hands the result to displayResults() for rendering.
  function run() {
    const input = document.getElementById('fastaText');
    if (!input) return;

    const raw = input.value.trim();
    if (!raw) {
      window.showToast('Please enter FASTA data.');
      return;
    }

    window.withLoading('panel-fasta', (records) => {
      _lastRecords = records;
      displayResults(records);
    }, raw.length, {
      type: 'PARSE_FASTA',
      payload: { rawText: raw }
    });
  }

  // ── Display ───────────────────────────────────────────
  // Not changed — identical to original fasta.js
  function displayResults(records) {
    const container = document.getElementById('fastaResults');
    if (!container) return;
    container.classList.remove('hidden');

    // Stat row
    window.buildStatCards('fastaStatRow', [
      { val: records.length, label: 'Sequences' },
      { val: records.reduce((a, r) => a + r.sequence.length, 0).toLocaleString(), label: 'Total Residues' },
    ]);

    const seqContainer = document.getElementById('fastaSequences');
    if (!seqContainer) return;

    const html = records.map((r, i) => {
      const format = window.detectFormat(r.sequence);
      const gc = format === 'dna' ? ((r.sequence.match(/[GC]/g)||[]).length / r.sequence.length * 100).toFixed(1) : null;
      return `
        <div class="fasta-seq-card">
          <div class="fasta-seq-header">
            <div>
              <span class="fasta-seq-id">${window.escapeHTML(r.header.split(' ')[0] || 'Seq' + (i+1))}</span>
              <span class="fasta-seq-desc">${window.escapeHTML(r.header.split(' ').slice(1).join(' '))}</span>
            </div>
            <div class="fasta-seq-meta">
              <span class="fasta-seq-badge">${r.sequence.length} ${format === 'dna' ? 'bp' : 'aa'}</span>
              <span class="fasta-seq-badge">${format.toUpperCase()}</span>
              ${gc ? `<span class="fasta-seq-badge">GC: ${gc}%</span>` : ''}
            </div>
          </div>
          <div class="fasta-seq-body">${window.escapeHTML(r.sequence.substring(0, 200))}${r.sequence.length > 200 ? '…' : ''}</div>
        </div>`;
    }).join('');
    seqContainer.innerHTML = window.sanitizeHTML(html);
  }

  // ── Events ────────────────────────────────────────────
  document.getElementById('fastaParse')?.addEventListener('click', run);

  document.getElementById('fastaClear')?.addEventListener('click', () => {
    const input = document.getElementById('fastaText');
    if (input) input.value = '';
    unmount(); // Reuse unmount to clear state + DOM
  });

  document.getElementById('fastaSample')?.addEventListener('click', () => {
    const input = document.getElementById('fastaText');
    if (input) {
      input.value = '>Seq1 GFP protein\nMSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK\n>Seq2 BRCA1 fragment\nATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAGTACAATAAAGAATTCATTCATTGATTCTTAGCTTGATTGAACATTCAATTCACAGCTTTCAATCA';
      input.dispatchEvent(new Event('input'));
    }
  });

  // ── Drag and Drop ─────────────────────────────────────
  const dropZone  = document.getElementById('fastaDropZone');
  const fileInput = document.getElementById('fastaFileInput');

  if (dropZone) {
    dropZone.addEventListener('dragover', (e) => {
      e.preventDefault();
      dropZone.classList.add('drag-over');
    });
    dropZone.addEventListener('dragleave', () => {
      dropZone.classList.remove('drag-over');
    });
    dropZone.addEventListener('drop', (e) => {
      e.preventDefault();
      dropZone.classList.remove('drag-over');
      const file = e.dataTransfer.files[0];
      if (file) readFile(file);
    });
  }

  if (fileInput) {
    fileInput.addEventListener('change', (e) => {
      const file = e.target.files[0];
      if (file) readFile(file);
    });
  }

  function readFile(file) {
    if (file.size > 5 * 1024 * 1024) { window.showToast?.('Error: File exceeds 5MB limit. Please split massive FASTA files for Free-Tier processing.', 'error'); return; }
    const reader = new FileReader();
    reader.onload = (e) => {
      const input = document.getElementById('fastaText');
      if (input) {
        input.value = e.target.result;
        run();
      }
    };
    reader.readAsText(file);
  }

})();
