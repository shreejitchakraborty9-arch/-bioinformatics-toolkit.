/* ============================================================
   seq-analyzer.js — DNA · RNA Sequence Analyzer
   Phase 1 (crash fix) + Phase 2 (feature wiring)
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── Helpers ─────────────────────────────────────────────
  const complement = { A: 'T', T: 'A', G: 'C', C: 'G', U: 'A' };

  function reverseComplement(seq) {
    return seq.split('').reverse().map(b => complement[b] || 'N').join('');
  }

  function transcribe(seq) {
    return seq.replace(/T/g, 'U');
  }

  function codonUsage(seq) {
    const freq = {};
    for (let i = 0; i <= seq.length - 3; i += 3) {
      const codon = seq.substring(i, i + 3);
      if (/^[ATGC]{3}$/.test(codon)) {
        freq[codon] = (freq[codon] || 0) + 1;
      }
    }
    return freq;
  }

  // ── Main Run ───────────────────────────────────────────
  function run() {
    const input = document.getElementById('dnaInput');
    if (!input) return;

    const raw = input.value;
    const validation = window.validateSequence(raw, 'dna');

    if (!validation.valid) {
      window.showValidationWarning('panel-dna', validation);
      return;
    }

    const seq = validation.clean;
    const BioMath = window.BioKit.core.BioMath;
    const readThrough = document.getElementById('dnaReadThrough')?.checked !== false;

    window.withLoading('panel-dna', () => {
      const gc = BioMath.calculateGC(seq);
      const tmRes = BioMath.calculateTmNN(seq);
      const ss_mw = BioMath.calculateDNA_MW(seq, false);
      const ds_mw = BioMath.calculateDNA_MW(seq, true);
      const ext = BioMath.calculateExtinctionCoefficient(seq);

      // Stat cards
      window.buildStatCards('dnaStatRow', [
        { val: seq.length.toLocaleString(), label: 'Length (bp)' },
        { val: gc + '%', label: 'GC Content' },
        { val: tmRes.tm + '°C', label: 'Melting Temp' },
        { val: ext.toLocaleString(), label: 'Ext. Coeff (L/mol·cm)' },
        { val: ss_mw + ' kDa', label: 'ssDNA MW' },
        { val: ds_mw + ' kDa', label: 'dsDNA MW' }
      ]);

      // Reverse complement
      const rc = reverseComplement(seq);
      const rcEl = document.getElementById('rcDisplay');
      if (rcEl) rcEl.textContent = rc;

      // RNA transcript
      const rna = transcribe(seq);
      const rnaEl = document.getElementById('rnaDisplay');
      if (rnaEl) rnaEl.textContent = rna;

      // Protein translation (3 frames)
      const framesEl = document.getElementById('proteinFrames');
      if (framesEl) {
        let framesHtml = '';
        for (let f = 0; f < 3; f++) {
          const frameSeq = seq.substring(f);
          const protein = BioMath.translateDNA(frameSeq, { stopAtFirst: !readThrough });
          framesHtml += `
            <div class="frame-block">
              <div class="frame-label">Frame +${f + 1}</div>
              <div class="seq-display" style="color:var(--violet);max-height:120px;">${window.escapeHTML(protein)}</div>
            </div>`;
        }
        framesEl.innerHTML = window.sanitizeHTML(framesHtml);
      }

      // Codon usage chart
      const freq = codonUsage(seq);
      renderCodonChart(freq);

      // Show results, hide empty state
      const empty = document.getElementById('dnaEmptyState');
      const content = document.querySelector('#dnaResults .results-content');
      if (empty) empty.classList.add('hidden');
      if (content) content.classList.remove('hidden');

      // Visual feedback
      window.showToast('Analysis Complete');
    }, seq.length);
  }

  // ── Codon Usage Bar Chart ──────────────────────────────
  function renderCodonChart(freq) {
    const canvas = document.getElementById('codonCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');

    const entries = Object.entries(freq).sort((a, b) => b[1] - a[1]).slice(0, 20);
    if (entries.length === 0) return;

    canvas.width  = canvas.offsetWidth || 700;
    canvas.height = 300;
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    const W      = canvas.width;
    const H      = canvas.height;
    const pad    = { top: 20, right: 20, bottom: 50, left: 40 };
    const plotW  = W - pad.left - pad.right;
    const plotH  = H - pad.top - pad.bottom;
    const maxVal = Math.max(...entries.map(e => e[1]));
    const barW   = Math.floor(plotW / entries.length) - 4;

    // Background
    ctx.fillStyle = '#0a1020';
    if (ctx.roundRect) { ctx.roundRect(0, 0, W, H, 8); ctx.fill(); }
    else { ctx.fillRect(0, 0, W, H); }

    // Grid
    ctx.strokeStyle = '#1e2d4a';
    ctx.lineWidth   = 1;
    for (let i = 0; i <= 4; i++) {
      const y = pad.top + (plotH / 4) * i;
      ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(W - pad.right, y); ctx.stroke();
    }

    // Bars
    entries.forEach(([codon, cnt], i) => {
      const x  = pad.left + i * (barW + 4) + 2;
      const h  = (cnt / maxVal) * plotH;
      const y  = pad.top + plotH - h;

      const grad = ctx.createLinearGradient(x, y, x, y + h);
      grad.addColorStop(0, '#00d4aa');
      grad.addColorStop(1, '#7c3aed');
      ctx.fillStyle = grad;
      ctx.beginPath();
      if (ctx.roundRect) { ctx.roundRect(x, y, barW, h, [4, 4, 0, 0]); ctx.fill(); }
      else { ctx.fillRect(x, y, barW, h); }

      // Count label
      ctx.fillStyle = '#e8f0fe';
      ctx.font      = '10px "JetBrains Mono"';
      ctx.textAlign = 'center';
      if (h > 14) ctx.fillText(cnt, x + barW / 2, y - 3);

      // Codon label
      ctx.fillStyle = '#8899bb';
      ctx.font      = '9px "JetBrains Mono"';
      ctx.save();
      ctx.translate(x + barW / 2, H - pad.bottom + 6);
      ctx.rotate(-Math.PI / 3);
      ctx.fillText(codon, 0, 0);
      ctx.restore();
    });
  }

  // ── Events ─────────────────────────────────────────────
  document.getElementById('dnaAnalyzeBtn')?.addEventListener('click', window.debounce(run, 300));

  document.getElementById('dnaSampleBtn')?.addEventListener('click', () => {
    const input = document.getElementById('dnaInput');
    if (input) {
      input.value = 'ATGGCTATCAAGCAGAAGTTTGATGCCATCAAGAAGCTGGAGAAGATCGGCATGAAGCTGGACATCGCCAAGGAGCTGAACTGA';
      input.dispatchEvent(new Event('input'));
    }
  });

  document.getElementById('dnaClearBtn')?.addEventListener('click', () => {
    const input = document.getElementById('dnaInput');
    if (input) { input.value = ''; input.dispatchEvent(new Event('input')); }
    document.getElementById('dnaResults')?.classList.add('hidden');
  });

  const input = document.getElementById('dnaInput');
  const meta = document.getElementById('dnaInputMeta');
  if (input && meta) {
    input.addEventListener('input', () => {
      const seq = window.cleanSeq(input.value);
      meta.textContent = `Length: ${seq.length} bp`;
    });
  }

})();
