/* ============================================================
   fasta.js — FASTA Parser (Fixed DOM IDs + Drag-Drop)
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── FASTA Parser ──────────────────────────────────────
  function parseFASTA(text) {
    const records = [];
    const lines = text.split('\n');
    let current = null;

    lines.forEach(line => {
      const trimmed = line.trim();
      if (trimmed.startsWith('>')) {
        if (current) records.push(current);
        current = { header: trimmed.substring(1), sequence: '' };
      } else if (current && trimmed) {
        current.sequence += trimmed.toUpperCase().replace(/\s/g, '');
      }
    });
    if (current) records.push(current);
    return records;
  }

  // ── Run ───────────────────────────────────────────────
  function run() {
    const input = document.getElementById('fastaText');
    if (!input) return;

    const raw = input.value.trim();
    if (!raw) {
      window.showToast('Please enter FASTA data.');
      return;
    }

    window.withLoading('panel-fasta', () => {
      const records = parseFASTA(raw);
      displayResults(records);
    });
  }

  // ── Display ───────────────────────────────────────────
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
    document.getElementById('fastaResults')?.classList.add('hidden');
  });

  document.getElementById('fastaSample')?.addEventListener('click', () => {
    const input = document.getElementById('fastaText');
    if (input) {
      input.value = '>Seq1 GFP protein\nMSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK\n>Seq2 BRCA1 fragment\nATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAGTACAATAAAGAATTCATTCATTGATTCTTAGCTTGATTGAACATTCAATTCACAGCTTTCAATCA';
      input.dispatchEvent(new Event('input'));
    }
  });

  // ── Drag and Drop ─────────────────────────────────────
  const dropZone = document.getElementById('fastaDropZone');
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
