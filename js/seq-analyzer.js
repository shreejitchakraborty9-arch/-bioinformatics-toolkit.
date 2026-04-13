/* ============================================================
   seq-analyzer.js — DNA · RNA Sequence Analyzer
   Phase 1 (crash fix) + Phase 2 (feature wiring)
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── Helpers ─────────────────────────────────────────────
  const complement = { 
    A: 'T', T: 'A', G: 'C', C: 'G', U: 'A',
    R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
    B: 'V', V: 'B', D: 'H', H: 'D', N: 'N'
  };

  function reverseComplement(seq) {
    return seq.toUpperCase().split('').reverse().map(b => complement[b] || 'N').join('');
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

    // Non-blocking soft warning for short sequences (2–9 bp)
    window.clearValidationAlert('panel-dna');
    if (validation.warning) {
      window.showValidationWarning('panel-dna', {
        valid: true,
        msg: validation.warning,
        suggestion: 'Analysis will proceed, but treat Tm values with caution.'
      });
    }

    // ── Surface any pending fetch warnings (e.g., ambiguous IUPAC bases) ──
    if (window.BioKit.pendingWarnings?.length > 0) {
      window.BioKit.pendingWarnings.forEach(w => {
        if (w.type === 'AMBIGUOUS_BASES') {
          window.showToast(
            `⚠️ Ambiguous bases detected: ${w.ambiguousTypes.join(', ')} ` +
            `(${w.ambiguousCount} total, ${w.canonicalPercent}% canonical). ` +
            `Tm and MW calculations may be approximate.`,
            5000
          );
        }
      });
      window.BioKit.pendingWarnings = []; // clear after display — do not repeat on next run
    }

    const seq = validation.clean;
    const BioMath = window.BioKit.core.BioMath;
    const readThrough = document.getElementById('dnaReadThrough')?.checked !== false;
    const saltFormVal = document.getElementById('dnaSaltForm')?.value || 'Na_50';

    // ── Salt Condition Map (affects both Tm and MW) ───────
    const SALT_MAP = {
      'Na_50':   { naConc_mM: 50,   saltFormKey: 'Na', label: '50 mM Na⁺' },
      'Na_100':  { naConc_mM: 100,  saltFormKey: 'Na', label: '100 mM Na⁺' },
      'Na_200':  { naConc_mM: 200,  saltFormKey: 'Na', label: '200 mM Na⁺' },
      'Na_1000': { naConc_mM: 1000, saltFormKey: 'Na', label: '1000 mM Na⁺ (Saturated)' },
      'K_50':    { naConc_mM: 50,   saltFormKey: 'K',  label: '50 mM K⁺' }, // K+ uses same eqNa correction in NN model
      'Free':    { naConc_mM: 1,    saltFormKey: 'Free', label: 'Free Acid (~1 mM)' },
      'custom':  { naConc_mM: parseFloat(document.getElementById('dnaCustomSaltConc')?.value) || 50, saltFormKey: 'Na', label: 'Custom' }
    };
    const saltConfig = SALT_MAP[saltFormVal] || SALT_MAP['Na_50'];
    const tmOptions = { naConc_mM: saltConfig.naConc_mM };

    window.withLoading('panel-dna', () => {
      const gc = BioMath.calculateGC(seq);
      const tmRes = BioMath.calculateTm(seq, tmOptions);
      const ss_mw = BioMath.calculateDNA_MW(seq, false, { saltForm: saltConfig.saltFormKey });
      const ds_mw = BioMath.calculateDNA_MW(seq, true,  { saltForm: saltConfig.saltFormKey });
      const ext = BioMath.calculateExtinctionCoefficient(seq);
      const ratioRes = BioMath.calculateA260_A280_Theoretical(seq);

      // Stat cards
      window.buildStatCards('dnaStatRow', [
        { val: seq.length.toLocaleString(), label: 'Length (bp)' },
        { val: gc + '%', label: 'GC Content' },
        { 
          val: tmRes.isValid ? tmRes.tm + '°C' : '⚠ ' + tmRes.tm + '°C', 
          label: `Tm at ${saltConfig.label}` 
        },
        { val: ss_mw + ' kDa', label: 'ssDNA MW' },
        { val: ds_mw + ' kDa', label: 'dsDNA MW' }
      ]);

      // A260/A280 Theoretical Disclaimer
      const ratioDisplay = document.getElementById('dnaA260A280Display');
      if (ratioDisplay) {
        let htmlContent = `
          <div style="background: rgba(255,170,0,0.05); border: 1px solid rgba(255,170,0,0.2); padding: 12px; border-radius: 6px; font-size: 0.8rem; color: var(--text-secondary); margin-bottom: 12px;">
            <strong style="color: var(--text-primary);">A260/A280 Ratio (Theoretical):</strong> ${ratioRes.ratio}
            <div style="color: #f59e0b; margin-top: 4px;">⚠️ ${ratioRes.note}</div>
            <div style="margin-top: 4px; font-size: 0.75rem;">Expected range for pure DNA: ${ratioRes.expectedRange.pureDNA || '1.8-1.9'}</div>
          </div>
        `;
        
        // Add Tm Model Warning if present (e.g. low salt extrapolation)
        if (tmRes.warning) {
            htmlContent += `
              <div style="background: rgba(245, 158, 11, 0.08); border: 1px solid rgba(245, 158, 11, 0.35); padding: 12px; border-radius: 6px; font-size: 0.8rem; color: var(--text-secondary);">
                <strong style="color: #f59e0b;">⚠️ Thermodynamic Model Warning</strong>
                <div style="margin-top: 4px;">${tmRes.warning}</div>
              </div>
            `;
            window.showToast("⚠️ Tm Model Warning: Salt condition extrapolates beyond mathematical validation.", 5000);
        }
        
        ratioDisplay.innerHTML = window.sanitizeHTML(htmlContent);
      }

      // Reverse complement
      const rc = reverseComplement(seq);
      const rcEl = document.getElementById('rcDisplay');
      if (rcEl) rcEl.textContent = rc;

      // RNA transcript
      const rna = transcribe(seq);
      const rnaEl = document.getElementById('rnaDisplay');
      if (rnaEl) rnaEl.textContent = rna;

      // Protein translation (6 frames)
      const framesEl = document.getElementById('proteinFrames');
      if (framesEl) {
        let framesHtml = '';
        
        // Forward Frames (+1, +2, +3)
        for (let f = 0; f < 3; f++) {
          const frameSeq = seq.substring(f);
          const protein = BioMath.translateDNA(frameSeq, { stopAtFirst: !readThrough });
          framesHtml += `
            <div class="frame-block">
              <div class="frame-label">Frame +${f + 1}</div>
              <div class="seq-display" style="color:var(--violet);max-height:120px;">${window.escapeHTML(protein)}</div>
            </div>`;
        }

        // Reverse Frames (-1, -2, -3)
        const rcSeq = rc; // Already computed above at line 139
        for (let f = 0; f < 3; f++) {
          const frameSeq = rcSeq.substring(f);
          const protein = BioMath.translateDNA(frameSeq, { stopAtFirst: !readThrough });
          framesHtml += `
            <div class="frame-block">
              <div class="frame-label">Frame -${f + 1} (Reverse Strand)</div>
              <div class="seq-display" style="color:var(--blue);max-height:120px;">${window.escapeHTML(protein)}</div>
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

  // ── Custom Salt Input Toggle ────────────────────────────
  const saltFormEl = document.getElementById('dnaSaltForm');
  const customSaltBlock = document.getElementById('dnaCustomSaltBlock');
  if (saltFormEl && customSaltBlock) {
    saltFormEl.addEventListener('change', () => {
      customSaltBlock.style.display = saltFormEl.value === 'custom' ? 'block' : 'none';
    });
  }

})();
