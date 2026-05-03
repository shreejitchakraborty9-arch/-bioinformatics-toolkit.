/* ============================================================
   seq-analyzer.js — DNA · RNA Sequence Analyzer (Scientific Grade)
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── Scientific Initialization ──────────────────────────
  function initScientificUI() {
    const codeSelect = document.getElementById('dnaGeneticCode');
    if (codeSelect && window.BioKit.data.GENETIC_CODES) {
      const codes = window.BioKit.data.GENETIC_CODES;
      Object.entries(codes).forEach(([id, info]) => {
        const opt = document.createElement('option');
        opt.value = id;
        opt.textContent = `[${id}] ${info.name}`;
        codeSelect.appendChild(opt);
      });
    }

    // Toggle Collapsible Settings
    const toggle = document.getElementById('dnaScientificToggle');
    const content = document.getElementById('dnaScientificContent');
    if (toggle && content) {
      toggle.addEventListener('click', () => {
        content.classList.toggle('hidden');
        toggle.querySelector('.toggle-icon')?.classList.toggle('rotated');
      });
    }
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

    window.clearValidationAlert('panel-dna');
    const seq = validation.clean;
    const BioMath = window.BioKit.core.BioMath;

    // ── Auto-detect sequence type (DNA vs RNA) ────────────
    const seqUpper = seq.toUpperCase();
    const seqType = seqUpper.includes('U') ? 'rna' : 'dna';
    const tmType  = seqType === 'rna' ? 'rna-rna' : 'dna-dna';
    const tmModel = seqType === 'rna' ? 'Turner 2004' : 'SantaLucia 1998';

    // ── Retrieve Scientific Parameters ────────────────────
    const saltFormVal = document.getElementById('dnaSaltForm')?.value || 'Na_50';
    const SALT_MAP = {
      'Na_50':   0.05, 'Na_100':  0.1, 'Na_200':  0.2, 'Na_1000': 1.0,
      'K_50':    0.05, 'Free':    0.001,
      'custom':  (parseFloat(document.getElementById('dnaCustomSaltConc')?.value) || 50) / 1000
    };
    const NaConc = SALT_MAP[saltFormVal];

    const mgVal  = parseFloat(document.getElementById('dnaMgConc')?.value) || 0;
    const Mg     = mgVal / 1000; // mM → M

    const ctVal = parseFloat(document.getElementById('dnaCt')?.value) || 0.25;
    const ctUnit = document.getElementById('dnaCtUnit')?.value || 'uM';
    const Ct = ctUnit === 'uM' ? ctVal * 1e-6 : ctVal * 1e-9;
    
    const codeId = parseInt(document.getElementById('dnaGeneticCode')?.value) || 1;
    const topology = document.getElementById('dnaTopology')?.value || 'linear';
    const terminal = document.getElementById('dnaTerminal')?.value || 'OH';
    const minORF = parseInt(document.getElementById('dnaMinORF')?.value) || 30;

    window.withLoading('panel-dna', () => {
      // 1. DNA Analysis (Backward Compatibility / Baseline)
      const gc = BioMath.calculateGC(seq);
      
      // 2. Scientific RNA Engine (UPGRADE)
      const rnaMass = BioMath.calculateScientificRNA_MW(seq, { topology, terminal });
      const rnaTm = BioMath.calculateScientificTm(seq, { type: tmType, Ct, Na: NaConc, Mg });
      const indMetrics = BioMath.calculateIndustrialMetrics(seq, { topology, terminal });
      
      // 3. Advanced Translation (ORF Analysis UPGRADE)
      const orfs = BioMath.analyzeORFs(seq, { codeId, minLen: minORF });

      // 4. Synthesis & Viability QC (GC, homopolymer, thermodynamic self-dimer, hairpin)
      // Pass actual [Na⁺] and free [Mg²⁺] so the ΔG engine uses the same salt
      // conditions as the Tm calculation. rnaTm.freeMg is already Von Ahsen-corrected.
      const qcFlags = BioMath.runQualityControl(
        seq,
        NaConc,
        rnaTm && rnaTm.freeMg !== null ? rnaTm.freeMg : 0
      );
      const qcHtml = qcFlags.length > 0
        ? qcFlags.map(f => `<div class="metric-note">\u26A0\uFE0F ${f}</div>`).join('')
        : '<div class="metric-note">\u2713 All synthesis &amp; viability QC checks passed.</div>';

      // Pre-compute free Mg note for the Tm card
      const freeMgNote = (rnaTm && rnaTm.freeMg !== null)
        ? ` | Free $[Mg^{2+}] = ${(rnaTm.freeMg * 1000).toFixed(2)}\\ mM$ (Von Ahsen 2001)`
        : '';

      // ── Update Stat Cards ───────────────────────────────
      window.buildStatCards('dnaStatRow', [
        { val: seq.length.toLocaleString(), label: 'Length (bp)' },
        { val: gc + '%', label: 'GC Content' },
        { val: rnaTm.tm + '°C', label: `${seqType.toUpperCase()} Tm (${tmModel})` },
        { val: rnaMass.avg.toLocaleString(), label: 'RNA Mass (Avg Da)' },
        { val: orfs.length, label: 'Potential ORFs' }
      ]);

      // ── Molar & Industrial Summary Area ────────────────
      const ratioDisplay = document.getElementById('dnaA260A280Display');
      if (ratioDisplay) {
        ratioDisplay.innerHTML = window.sanitizeHTML(`
          <div class="scientific-metrics-grid">
            <div class="metric-card scientific-card">
              <div class="metric-header">RNA Molecular Mass Engine (Exact)</div>
              <div class="metric-body">
                <div>Average Mass: $${rnaMass.avg}$ Da</div>
                <div class="mass-range">IUPAC Range: [$${rnaMass.min}$, $${rnaMass.max}$] Da</div>
                <div class="metric-note">Topology: ${topology.toUpperCase()} | 5' Group: ${terminal}</div>
              </div>
            </div>
            <div class="metric-card scientific-card">
              <div class="metric-header">Thermodynamic Stability Module (NN)</div>
              <div class="metric-body">
                <div class="formula-block">$${rnaTm.formula}$</div>
                <div class="result-block">Predicted $${rnaTm.latex}$</div>
                <div class="metric-note">Model: ${tmModel} | $C_t = ${ctVal}\ ${ctUnit}$ | $[Na^+] = ${NaConc}\ M$${Mg > 0 ? ` | Total $[Mg^{2+}] = ${mgVal}\\ mM$${freeMgNote}` : ''}${rnaTm && rnaTm.degenerateNote ? `<br>${rnaTm.degenerateNote}` : ''}</div>
              </div>
            </div>
            <div class="metric-card scientific-card">
              <div class="metric-header">Industrial Metrics (ϵ260)</div>
              <div class="metric-body">
                <div>Extinction Coefficient: $${indMetrics.latex.e260}$</div>
                <div>Copy Estimate: $${indMetrics.latex.copies}$</div>
                <div class="metric-note">Standard Mass Conc (RNA): $1 A_{260} = 40 \mu g/ml$</div>
              </div>
            </div>
            <div class="metric-card scientific-card">
              <div class="metric-header">Synthesis &amp; Viability QC</div>
              <div class="metric-body">${qcHtml}</div>
            </div>
          </div>
        `);
        BioMath.formatLaTeX(ratioDisplay);
      }

      // ── Tab Content: Translation (ORF Highlight) ───────
      const framesEl = document.getElementById('proteinFrames');
      if (framesEl) {
        let orfHtml = '<div class="orf-results-header">Systematic ORF Map (Top Hits)</div>';
        if (orfs.length === 0) {
            orfHtml += '<div class="no-results">No ORFs found above threshold.</div>';
        } else {
            orfs.slice(0, 10).forEach((orf, i) => {
                orfHtml += `
                  <div class="orf-block">
                    <div class="orf-badge-row">
                      <span class="orf-id">ORF ${i+1}</span>
                      <span class="orf-meta">${orf.strand}${orf.frame} | ${orf.length} nt | ${orf.translation.length} aa</span>
                      <span class="kozak-badge ${orf.kozak.toLowerCase()}">Kozak: ${orf.kozak}</span>
                    </div>
                    <div class="seq-display orf-seq" style="font-family:var(--mono);word-break:break-all;">${window.escapeHTML(orf.translation)}</div>
                    <div class="orf-indices">Coordinates: ${orf.start + 1} ... ${orf.end}</div>
                  </div>`;
            });
        }
        framesEl.innerHTML = window.sanitizeHTML(orfHtml);
      }

      // Tab Content: RNA
      const rnaEl = document.getElementById('rnaDisplay');
      if (rnaEl) rnaEl.textContent = seq.toUpperCase().replace(/T/g, 'U');

      // Tab Content: RevComp (Full IUPAC)
      const rcEl = document.getElementById('rcDisplay');
      if (rcEl) rcEl.textContent = BioMath.reverseComplement(seq);

      const empty = document.getElementById('dnaEmptyState');
      const content = document.querySelector('#dnaResults .results-content');
      if (empty) empty.classList.add('hidden');
      if (content) content.classList.remove('hidden');

      window.showToast('Scientific Analysis Complete');
    }, seq.length);
  }

  // ── Events ─────────────────────────────────────────────
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initScientificUI);
  } else {
    initScientificUI();
  }
  
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
      const detectedType = seq.toUpperCase().includes('U') ? 'RNA' : 'DNA';
      const unit = detectedType === 'RNA' ? 'nt' : 'bp';
      meta.textContent = `Length: ${seq.length} ${unit} (${detectedType} detected)`;
    });
  }

  const saltFormEl = document.getElementById('dnaSaltForm');
  const customSaltBlock = document.getElementById('dnaCustomSaltBlock');
  if (saltFormEl && customSaltBlock) {
    saltFormEl.addEventListener('change', () => {
      customSaltBlock.style.display = saltFormEl.value === 'custom' ? 'block' : 'none';
    });
  }

  // ── Quick-Action Buttons ───────────────────────────────
  function _showQuickResult() {
    const empty = document.getElementById('dnaEmptyState');
    const content = document.querySelector('#dnaResults .results-content');
    if (empty) empty.classList.add('hidden');
    if (content) content.classList.remove('hidden');
  }

  document.getElementById('dnaRevCompBtn')?.addEventListener('click', () => {
    const raw = document.getElementById('dnaInput')?.value || '';
    const v = window.validateSequence(raw, 'dna');
    if (!v.valid) { window.showValidationWarning('panel-dna', v); return; }
    window.clearValidationAlert('panel-dna');
    const rcEl = document.getElementById('rcDisplay');
    if (rcEl) rcEl.textContent = window.BioKit.core.BioMath.reverseComplement(v.clean);
    _showQuickResult();
  });

  document.getElementById('dnaToRnaBtn')?.addEventListener('click', () => {
    const raw = document.getElementById('dnaInput')?.value || '';
    const v = window.validateSequence(raw, 'dna');
    if (!v.valid) { window.showValidationWarning('panel-dna', v); return; }
    window.clearValidationAlert('panel-dna');
    const rnaEl = document.getElementById('rnaDisplay');
    if (rnaEl) rnaEl.textContent = v.clean.toUpperCase().replace(/T/g, 'U');
    _showQuickResult();
  });

  document.getElementById('dnaTranslateBtn')?.addEventListener('click', window.debounce(run, 300));

})();
