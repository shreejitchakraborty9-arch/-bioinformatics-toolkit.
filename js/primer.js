/* ============================================================
   primer.js — Primer Design & Analysis
   ============================================================ */

(function () {
  'use strict';

  // ── Helpers ────────────────────────────────────────────────
  function calcTm(seq, options) {
    const BioMath = window.BioKit && window.BioKit.core && window.BioKit.core.BioMath;
    if (BioMath) {
      return BioMath.calculateTmNN(seq, options);
    }
    return { tm: "0.0", confidence: "", warning: null };
  }

  function tmDiff(tm1, tm2) {
    return Math.abs(parseFloat(tm1) - parseFloat(tm2)).toFixed(1);
  }

  function updateTemplateUI(data) {
    const meta = document.getElementById('templateMeta');
    if (meta) meta.textContent = `Length: ${data.sequence.length} bp | Source: ${data.accession || 'Database'}`;
    // Optional: map features to genome viewer or summary
  }

  // ── State (PROMPT 2) ──────────────────────────────────────
  const state = {
    template: {
      sequence: '',
      organism: 'human',
      type: 'genomic',
      features: []
    },
    assayType: 'pcr',
    constraints: {}, // Populated by preset
    biologicalOptions: {
      avoidExonJunction: false,
      intronSpanning: false,
      excludeSnp: false
    },
    lastAnalysis: null
  };

  // ── Helpers (Prompt 4) ──────────────────────────────────────
  function loadAssayPreset(type) {
    const BioMath = window.BioKit && window.BioKit.core && window.BioKit.core.BioMath;
    if (!BioMath || !BioMath.ASSAY_PRESETS) return;
    const presets = BioMath.ASSAY_PRESETS;
    if (presets[type]) {
      state.constraints = JSON.parse(JSON.stringify(presets[type]));
      renderConstraints();
    }
  }

  function renderConstraints() {
    const display = document.getElementById('constraintsDisplay');
    const c = state.constraints;
    if (!c || !c.length) return;

    display.innerHTML = `
      <div class="summary-item"><span>Length</span><span>${c.length.min}-${c.length.max} bp</span></div>
      <div class="summary-item"><span>Tm</span><span>${c.tm.min}-${c.tm.max} °C</span></div>
      <div class="summary-item"><span>Tm Diff</span><span>≤ ${c.tmDifference} °C</span></div>
      <div class="summary-item"><span>GC Content</span><span>${c.gc.min}-${c.gc.max}%</span></div>
      <div class="summary-item"><span>3' GC Clamp</span><span>${c.gc3prime ? 'Required' : 'Optional'}</span></div>
      <div class="summary-item"><span>Polymer Run</span><span>≤ ${c.homopolymerLimit} bp</span></div>
      <div class="summary-item"><span>Amplicon</span><span>${c.productSize.min}-${c.productSize.max} bp</span></div>
    `;
    
    // Sync to edit form
    document.getElementById('conMinLen').value = c.length.min;
    document.getElementById('conMaxLen').value = c.length.max;
    document.getElementById('conMinTm').value = c.tm.min;
    document.getElementById('conMaxTm').value = c.tm.max;
    document.getElementById('conTmDiff').value = c.tmDifference;
    document.getElementById('conMinGC').value = c.gc.min;
    document.getElementById('conMaxGC').value = c.gc.max;
    document.getElementById('conHomoLimit').value = c.homopolymerLimit;
    document.getElementById('conMinProd').value = c.productSize.min;
    document.getElementById('conMaxProd').value = c.productSize.max;
    document.getElementById('conGC3').checked = c.gc3prime;
  }

  // ── Input live meta ────────────────────────────────────────
  const templateInput = document.getElementById('primerTemplate');
  const templateMeta = document.getElementById('templateMeta');
  const fwdInput  = document.getElementById('fwdPrimer');
  const revInput  = document.getElementById('revPrimer');
  const fwdMeta   = document.getElementById('fwdMeta');
  const revMeta   = document.getElementById('revMeta');

  // ── Collapsible Logic ──────────────────────────────────────
  const togglePanel = (headerId, contentId) => {
    const hdr = document.getElementById(headerId);
    const content = document.getElementById(contentId);
    if (!hdr || !content) return;
    hdr.addEventListener('click', (e) => {
      if (e.target.tagName === 'BUTTON') return;
      content.classList.toggle('hidden');
      hdr.classList.toggle('open');
      const icon = hdr.querySelector('.toggle-icon');
      if (icon) {
        if (content.classList.contains('hidden')) {
          icon.style.transform = 'rotate(0deg)';
        } else {
          icon.style.transform = 'rotate(90deg)';
        }
      }
    });
  };
  togglePanel('templateSettingsToggle', 'templateSettingsContent');
  togglePanel('constraintsToggle', 'constraintsContent');

  // ── Constraint Persistence & Editing ───────────────────────
  const editBtn = document.getElementById('editConstraintsBtn');
  const resetBtn = document.getElementById('resetConstraintsBtn');
  const saveBtn = document.getElementById('saveConstraintsBtn');
  const cancelBtn = document.getElementById('cancelConstraintsBtn');
  const editForm = document.getElementById('constraintsEditForm');
  const displayArea = document.getElementById('constraintsDisplay');

  if (editBtn) editBtn.addEventListener('click', () => {
    editForm.classList.remove('hidden');
    displayArea.classList.add('hidden');
  });
  if (cancelBtn) cancelBtn.addEventListener('click', () => {
    editForm.classList.add('hidden');
    displayArea.classList.remove('hidden');
  });
  if (resetBtn) resetBtn.addEventListener('click', () => {
    loadAssayPreset(state.assayType);
    window.showToast(`✓ Reset to ${state.assayType.toUpperCase()} defaults`);
  });
  if (saveBtn) saveBtn.addEventListener('click', () => {
    state.constraints.length.min = parseInt(document.getElementById('conMinLen').value);
    state.constraints.length.max = parseInt(document.getElementById('conMaxLen').value);
    state.constraints.tm.min = parseFloat(document.getElementById('conMinTm').value);
    state.constraints.tm.max = parseFloat(document.getElementById('conMaxTm').value);
    state.constraints.tmDifference = parseFloat(document.getElementById('conTmDiff').value);
    state.constraints.gc.min = parseFloat(document.getElementById('conMinGC').value);
    state.constraints.gc.max = parseFloat(document.getElementById('conMaxGC').value);
    state.constraints.homopolymerLimit = parseInt(document.getElementById('conHomoLimit').value);
    state.constraints.productSize.min = parseInt(document.getElementById('conMinProd').value);
    state.constraints.productSize.max = parseInt(document.getElementById('conMaxProd').value);
    state.constraints.gc3prime = document.getElementById('conGC3').checked;
    
    renderConstraints();
    editForm.classList.add('hidden');
    displayArea.classList.remove('hidden');
    window.showToast('✓ Custom constraints saved');
  });

  // ── Source Switching ──────────────────────────────────────
  const sourceRadios = document.querySelectorAll('input[name="templateSource"]');
  const pasteGroup = document.getElementById('pasteTemplateGroup');
  const fetchGroup = document.getElementById('fetchTemplateGroup');
  
  sourceRadios.forEach(radio => {
    radio.addEventListener('change', () => {
      if (radio.value === 'paste') {
        pasteGroup.classList.remove('hidden');
        fetchGroup.classList.add('hidden');
      } else {
        pasteGroup.classList.add('hidden');
        fetchGroup.classList.remove('hidden');
      }
    });
  });

  // ── State Listeners ────────────────────────────────────────
  if (templateInput) templateInput.addEventListener('input', () => {
    const s = (templateInput.value || '').toUpperCase().replace(/[^A-Z]/gi, '');
    state.template.sequence = s;
    if (templateMeta) templateMeta.textContent = `Length: ${s.length} bp`;
    
    // Auto-toggle to paste view if sequence suddenly appears (from fetch)
    if (s.length > 0) {
      const pasteRadio = document.querySelector('input[name="templateSource"][value="paste"]');
      if (pasteRadio && !pasteRadio.checked) {
        pasteRadio.click();
      }
    }
  });

  const organismEl = document.getElementById('organism');
  if (organismEl) organismEl.addEventListener('change', (e) => {
    state.template.organism = e.target.value;
  });

  document.querySelectorAll('input[name="templateType"]').forEach(r => {
    r.addEventListener('change', () => state.template.type = r.value);
  });

  document.querySelectorAll('input[name="assayType"]').forEach(r => {
    r.addEventListener('change', () => {
      state.assayType = r.value;
      loadAssayPreset(r.value);
    });
  });
  
  const avoidExonJunctionEl = document.getElementById('avoidExonJunction');
  if (avoidExonJunctionEl) avoidExonJunctionEl.addEventListener('change', (e) => {
    state.biologicalOptions.avoidExonJunction = e.target.checked;
  });
  const intronSpanningEl = document.getElementById('intronSpanning');
  if (intronSpanningEl) intronSpanningEl.addEventListener('change', (e) => {
    state.biologicalOptions.intronSpanning = e.target.checked;
  });
  const excludeSnpEl = document.getElementById('excludeSnp');
  if (excludeSnpEl) excludeSnpEl.addEventListener('change', (e) => {
    state.biologicalOptions.excludeSnp = e.target.checked;
  });
  
  // Initialize Defaults
  loadAssayPreset('pcr');


  if (fwdInput) fwdInput.addEventListener('input', () => {
    const s = window.cleanSeq(fwdInput.value);
    const na = parseFloat(document.getElementById('saltConc').value) || 50;
    const k = parseFloat(document.getElementById('kConc').value) || 0;
    const options = { naConc_mM: na + k, oligoConc_nM: 250 };
    const tmRes = calcTm(s, options);
    const validation = window.validateSequence(s, 'dna');
    const gc = window.BioKit.core.BioMath.calculateGC(s);
    let metaText = `Length: ${s.length} nt  |  GC: ${gc}%  |  Tm ≈ ${tmRes.tm} °C`;
    if (!validation.valid || tmRes.warning) {
      metaText = `⚠ ${metaText} ${tmRes.warning ? '(Warning)' : '(Limit Exceeded)'}`;
      if (fwdMeta) fwdMeta.style.color = 'var(--rose)';
    } else {
      if (fwdMeta) fwdMeta.style.color = 'var(--text-muted)';
    }
    if (fwdMeta) fwdMeta.textContent = metaText;
  });
  if (revInput) revInput.addEventListener('input', () => {
    const s = window.cleanSeq(revInput.value);
    const na = parseFloat(document.getElementById('saltConc').value) || 50;
    const k = parseFloat(document.getElementById('kConc').value) || 0;
    const options = { naConc_mM: na + k, oligoConc_nM: 250 };
    const tmRes = calcTm(s, options);
    const validation = window.validateSequence(s, 'dna');
    const gc = window.BioKit.core.BioMath.calculateGC(s);
    let metaText = `Length: ${s.length} nt  |  GC: ${gc}%  |  Tm ≈ ${tmRes.tm} °C`;
    if (!validation.valid || tmRes.warning) {
      metaText = `⚠ ${metaText} ${tmRes.warning ? '(Warning)' : '(Limit Exceeded)'}`;
      if (revMeta) revMeta.style.color = 'var(--rose)';
    } else {
      if (revMeta) revMeta.style.color = 'var(--text-muted)';
    }
    if (revMeta) revMeta.textContent = metaText;
  });

  // ── Controls ───────────────────────────────────────────────
  const analyzeBtn = document.getElementById('primerAnalyzeBtn');
  if (analyzeBtn) analyzeBtn.addEventListener('click', window.debounce(analyze, 300));
  
  const sampleBtn = document.getElementById('primerSampleBtn');
  if (sampleBtn) sampleBtn.addEventListener('click', () => {
    templateInput.value = "ATGGCTATCAAGCAGAAGTTTGATGCCATCAAGAAGCTGGAGGAGCAGCTGACCAAGGACATCCAGTACAACATGGGCCTGGCCGACATGGCCGGCATCGTGGTGCACGGCCACCACATCAAGAAGCTGTGA";
    fwdInput.value = 'ATGGCTATCAAGCAGAAGTTTG';
    revInput.value = 'TCACAGCTTCTTGATGTGGTG';
    templateInput.dispatchEvent(new Event('input'));
    fwdInput.dispatchEvent(new Event('input'));
    revInput.dispatchEvent(new Event('input'));
  });

  // ── Source & Assay Toggles (PROMPT 2 & 4) ──────────────────
  document.querySelectorAll('input[name="templateSource"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
      const isFetch = e.target.value === 'fetch';
      document.getElementById('pasteTemplateGroup').classList.toggle('hidden', isFetch);
      document.getElementById('fetchTemplateGroup').classList.toggle('hidden', !isFetch);
    });
  });

  document.querySelectorAll('input[name="assayType"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
      state.assayType = e.target.value;
      loadAssayPreset(e.target.value);
    });
  });

  // Bind biological checkboxes
  ['avoidExonJunction', 'intronSpanning', 'excludeSnp'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.addEventListener('change', (e) => {
      state.biologicalOptions[id] = e.target.checked;
    });
  });

  // ── Database Fetching (Unified) ──────────────────────────
  // Handled globally by fetch-api.js via .fetch-btn class in index.html

  // ── Auto-Design Integration ────────────────────────────────
  const suggestBtn = document.getElementById('primerSuggestBtn');
  if (suggestBtn) suggestBtn.addEventListener('click', window.debounce(suggest, 300));

  async function suggest() {
    const template = window.cleanSeq(templateInput.value);
    if (!template || template.length < 50) {
      window.showToast('⚠ Provide a template sequence (min 50bp) for design');
      return;
    }

    const BioMath = window.BioKit.core.BioMath;
    const suggestionsGrid = document.getElementById('suggestionGrid');
    const suggestionsArea = document.getElementById('primerSuggestionsArea');

    // Gather buffer options
    const options = {
      naConc_mM: (parseFloat(document.getElementById('saltConc').value) || 50) + (parseFloat(document.getElementById('kConc').value) || 0),
      mgConc_mM: parseFloat(document.getElementById('mgConc').value) || 1.5,
      dntpConc_mM: parseFloat(document.getElementById('dntpConc').value) || 0.8,
      oligoConc_nM: parseFloat(document.getElementById('primerConc').value) || 250
    };

    window.withLoading('panel-primer', async () => {
      suggestionsArea.classList.remove('hidden');
      suggestionsGrid.innerHTML = '<div style="grid-column: 1/-1; padding: 40px; text-align:center; color:var(--text-muted);"><div class="spinner-sm" style="margin-bottom:12px;"></div>Scanning template for optimal pairs...</div>';

      // Small delay to ensure spinner renders
      await new Promise(r => setTimeout(r, 100));

      const pairs = await BioMath.suggestPrimerPairs(template, state.constraints, options);

      if (!pairs || pairs.length === 0) {
        suggestionsGrid.innerHTML = '<div style="grid-column: 1/-1; padding: 40px; text-align:center; color:var(--rose);">No pairs found matching constraints. Try relaxing Tm or Product Size limits.</div>';
        return;
      }

      suggestionsGrid.innerHTML = pairs.map((p, idx) => `
        <div class="suggestion-card" style="background:var(--card-bg); border:1px solid var(--border-color); padding:16px; border-radius:12px; transition:all 0.2s hover:border-teal;">
          <div style="display:flex; justify-content:space-between; margin-bottom:12px;">
            <strong style="color:var(--teal);">Pair #${idx + 1}</strong>
            <span style="font-size:0.75rem; color:var(--text-muted);">Score: ${p.score.toFixed(1)}</span>
          </div>
          
          <div style="font-size:0.85rem; margin-bottom:14px; color:var(--text-main);">
            <div style="margin-bottom:4px;"><strong>Fwd:</strong> <span style="font-family:var(--mono);">...${p.fwd.seq.slice(-15)}</span></div>
            <div><strong>Rev:</strong> <span style="font-family:var(--mono);">...${p.rev.seq.slice(-15)}</span></div>
          </div>

          <div style="display:grid; grid-template-columns: 1fr 1fr; gap:8px; font-size:0.75rem; color:var(--text-muted); margin-bottom:16px;">
            <div>Tm: ${p.fwd.tm}°/${p.rev.tm}°</div>
            <div>Avg GC: ${((parseFloat(p.fwd.gc) + parseFloat(p.rev.gc))/2).toFixed(1)}%</div>
            <div>Size: ${p.size} bp</div>
            <div>Hetero ΔG: ${p.heteroDg}</div>
          </div>

          <button class="btn-primary btn-sm suggest-apply-btn" style="width:100%;" 
                  data-fwd="${p.fwd.seq}" data-rev="${p.rev.seq}">
            Use This Pair
          </button>
        </div>
      `).join('');

      // Bind Apply buttons
      suggestionsGrid.querySelectorAll('.suggest-apply-btn').forEach(btn => {
        btn.addEventListener('click', () => {
          fwdInput.value = btn.dataset.fwd;
          revInput.value = btn.dataset.rev;
          fwdInput.dispatchEvent(new Event('input'));
          revInput.dispatchEvent(new Event('input'));
          analyze(); // Trigger full audit
          window.showToast('✅ Primers applied to audit engine');
          window.location.hash = 'primerResults';
        });
      });
    });
  }

  function analyze() {
    try {
      const BioMath = window.BioKit.core.BioMath;
      const na = parseFloat(document.getElementById('saltConc').value) || 50;
      const k = parseFloat(document.getElementById('kConc').value) || 0;
      const mg = parseFloat(document.getElementById('mgConc').value) || 1.5;
      const dntp = parseFloat(document.getElementById('dntpConc').value) || 0.8;
      const conc = parseFloat(document.getElementById('primerConc').value) || 250;
      const dmso = parseFloat(document.getElementById('dmsoConc').value) || 0;
      const formamide = parseFloat(document.getElementById('formamideConc').value) || 0;
      
      const assay = state.assayType;
      const organism = state.template.organism;
      const template = state.template.sequence;

      const options = {
        naConc_mM: na + k,
        mgConc_mM: mg,
        dntpConc_mM: dntp,
        oligoConc_nM: conc,
        dmso_pct: dmso,
        formamide_m: formamide
      };

      const fwdStr = window.cleanSeq(fwdInput.value);
      const revStr = window.cleanSeq(revInput.value);

      if (!fwdStr && !revStr) { window.showToast('⚠ Enter at least one primer'); return; }

      const valFwd = window.validateSequence(fwdInput.value, 'dna');
      const valRev = window.validateSequence(revInput.value, 'dna');

      if (!valFwd.valid || !valRev.valid) {
        window.showValidationWarning('panel-primer', !valFwd.valid ? valFwd : valRev);
        return;
      }

      const thresholds = state.constraints;
      const grid = document.getElementById('primerResultGrid');
      const primers = [];
      if (fwdStr.length >= 8) primers.push({ label: 'Forward Primer', seq: fwdStr });
      if (revStr.length >= 8) primers.push({ label: 'Reverse Primer', seq: revStr });

      window.withLoading('panel-primer', async () => {
        grid.textContent = ''; 
        const pairEl = document.getElementById('primerPairResults');
        pairEl.style.display = 'none';

        const results = primers.map(({ label, seq }) => {
          const tmRes = calcTm(seq, options);
          const tm = tmRes.tm;
          const gc = BioMath.calculateGC(seq);
          const len = seq.length;
          const hairpin = BioMath.detectHairpin(seq);
          const dimer = BioMath.calculateDimerThermodynamics(seq, seq);
          const terminal = BioMath.assess3PrimeStability(seq);
          const homoRes = BioMath.detectHomopolymerRuns(seq, thresholds.homopolymerLimit);

          // Validation using core engine
          const lenOk = BioMath.validateLength(seq, thresholds);
          const tmOk  = BioMath.validateTm(tm, thresholds);
          const gcOk  = BioMath.validateGC(gc, thresholds);
          const clampOk = BioMath.validateGCClamp(seq, thresholds);
          const homoOk  = !homoRes.hasViolation;
          const dgOk = parseFloat(dimer.deltaG) > (thresholds.dimerDeltaGLimit || -6.0);

          let failReasons = [];
          if (!lenOk) failReasons.push(`Length (${len}nt) outside ${thresholds.length.min}-${thresholds.length.max} range`);
          if (!gcOk) failReasons.push(`GC (${gc}%) outside ${thresholds.gc.min}-${thresholds.gc.max}% range`);
          if (!tmOk) failReasons.push(`Tm (${tm}°C) outside ${thresholds.tm.min}-${thresholds.tm.max}°C range`);
          if (hairpin.found) failReasons.push(`Hairpin risk detected`);
          if (!homoOk) failReasons.push(`Homopolymer detected: ${homoRes.summary}`);
          if (!clampOk) failReasons.push(`Missing 3' GC clamp (G or C within terminal 3nt)`);
          if (terminal.verdict === 'CRITICAL') failReasons.push(terminal.message);

          const cardHtml = `
            <div class="primer-card ${failReasons.length > 2 ? 'fail-border' : failReasons.length > 0 ? 'warn-border' : ''}">
              <div class="primer-card-header">
                <h4>${window.escapeHTML(label)}</h4>
                <span class="qc-status ${failReasons.length === 0 ? 'pass' : 'fail'}">${failReasons.length === 0 ? 'QC PASS' : 'QC FLAG'}</span>
              </div>
              <div class="primer-seq">${window.escapeHTML(seq)}</div>
              
              <div class="primer-stat-grid">
                <div class="primer-stat-item"><span class="primer-stat-key">Length</span>
                  <span class="primer-stat-val ${lenOk?'pass':'warn'}">${len} nt</span></div>
                <div class="primer-stat-item"><span class="primer-stat-key">Tm (SantaLucia NN)</span>
                  <span class="primer-stat-val ${tmOk?'pass':'warn'}">${tm} °C</span></div>
                <div class="primer-stat-item"><span class="primer-stat-key">GC Content</span>
                  <span class="primer-stat-val ${gcOk?'pass':'warn'}">${gc}%</span></div>
                <div class="primer-stat-item"><span class="primer-stat-key">3′ Integrity</span>
                  <span class="primer-stat-val ${terminal.verdict === 'STABLE' ? 'pass' : 'warn'}">${terminal.verdict}</span></div>
              </div>

              <div class="structural-report">
                <div class="struct-stat ${hairpin.found?'fail':'pass'}">Hairpin: ${hairpin.found?`Fail` : 'Pass'}</div>
                <div class="struct-stat ${dgOk?'pass':'warn'}">Self-Dimer: ${dimer.deltaG} kcal/mol</div>
                <div class="struct-stat ${homoOk?'pass':'fail'}">Homopolymer: ${homoRes.hasViolation ? homoRes.summary : 'None'}</div>
              </div>

              ${failReasons.length > 0 ? `
              <div class="fail-list">
                <strong>Reasons for Flag:</strong>
                <ul>${failReasons.map(r => `<li>${window.escapeHTML(r)}</li>`).join('')}</ul>
              </div>` : ''}

              <div class="thermo-transparency-section">
                <div class="thermo-header">Thermodynamic Transparency</div>
                <div class="thermo-grid">
                  <div class="thermo-item"><span>Model</span><span>SantaLucia 1998 NN</span></div>
                  <div class="thermo-item"><span>Salt [Na⁺/K⁺]</span><span>${na + k} mM</span></div>
                  <div class="thermo-item"><span>Mg²⁺</span><span>${mg} mM</span></div>
                  <div class="thermo-item"><span>dNTPs</span><span>${dntp} mM</span></div>
                  <div class="thermo-item"><span>Primer</span><span>${(conc/1000).toFixed(2)} µM</span></div>
                  <div class="thermo-item"><span>Additives</span><span>${(dmso > 0 || formamide > 0) ? `${dmso}% DMSO / ${formamide}M Form` : 'None'}</span></div>
                </div>
              </div>
            </div>`;
          grid.innerHTML += window.sanitizeHTML(cardHtml);
          return { label, seq, tm, selfDg: dimer.deltaG };
        });

        if (results.length === 2) {
          const diff = tmDiff(results[0].tm, results[1].tm);
          const diffOk = BioMath.validateTmDifference(results[0].tm, results[1].tm, thresholds);
          const hetero = BioMath.calculateDimerThermodynamics(results[0].seq, results[1].seq);
          const heteroDg = hetero.deltaG;

          const pairFailures = [];
          const pairWarnings = [];
          
          if (!diffOk) pairWarnings.push(`Tm difference (${diff}°C) exceeds assay tolerance (max ${thresholds.tmDifference}°C)`);
          const dLimit = thresholds.dimerDeltaGLimit || -9.0;
          if (parseFloat(heteroDg) <= dLimit) pairFailures.push(`High heterodimer risk (ΔG: ${heteroDg} kcal/mol)`);
          else if (parseFloat(heteroDg) <= (dLimit + 3.0)) pairWarnings.push(`Moderate heterodimer risk (ΔG: ${heteroDg} kcal/mol)`);

          // ── Specificity Section ──────────────────────
          let specificityHtml = '';
          if (template) {
            const fwdSpec = BioMath.checkTemplateSpecificity(fwdStr, template, false);
            const revSpec = BioMath.checkTemplateSpecificity(revStr, template, true);
            const specStatus = (fwdSpec.specificity === 'UNIQUE' && revSpec.specificity === 'UNIQUE') ? 'UNIQUE' : 'MULTIPLE';
            const specClass = specStatus === 'UNIQUE' ? 'pass' : 'warn';
            
            if (specStatus === 'MULTIPLE') pairFailures.push("Non-specific binding detected in template");

            specificityHtml = `
              <div class="verdict-section-title" style="margin-top:20px;">Specificity (Template-Local)</div>
              <div class="verdict-panel spec-panel ${specClass}">
                <div class="verdict-header">
                  <span>Target Specificity:</span>
                  <strong class="status-${specClass}">${specStatus}</strong>
                </div>
                <div class="verdict-row"><span>Forward Primer</span><span>${fwdSpec.primarySite.count} binding site(s)</span></div>
                <div class="verdict-row"><span>Reverse Primer</span><span>${revSpec.primarySite.count} binding site(s)</span></div>
                <div class="blast-note">
                  ⚠️ Genome-wide specificity screening not available.<br/>
                  (Future: Integrate NCBI BLAST)
                </div>
              </div>
            `;
          }

          // ── Amplicon Section ─────────────────────────────────────
          let ampliconHtml = '';
          let ampliconVerdict = 'PASS';
          let ampInfo = null;
          if (template) {
            ampInfo = BioMath.calculateAmpliconInfo(template, fwdStr, revStr, thresholds);
            if (ampInfo.amplicon.size > 0) {
              const size = ampInfo.amplicon.size;
              const sizeWarning = ampInfo.warnings.find(w => w.includes('size'));
              if (sizeWarning) {
                ampliconVerdict = 'CAUTION';
                pairWarnings.push(sizeWarning);
              }
              
              const seqDisp = size > 100 ? `${ampInfo.amplicon.sequence.slice(0, 50)}...${ampInfo.amplicon.sequence.slice(-50)}` : ampInfo.amplicon.sequence;
              
              ampliconHtml = `
                <div class="verdict-section-title">AMPLICON ANALYSIS</div>
                <div class="verdict-panel amplicon-panel">
                  <div class="amplicon-seq-box">${seqDisp} <span class="badge-size">${size} bp</span></div>
                  <div class="verdict-row"><span>Forward Binding</span><span>1–${fwdStr.length} bp</span></div>
                  <div class="verdict-row"><span>Reverse Binding</span><span>${size - revStr.length + 1}–${size} bp</span></div>
                  <div class="verdict-row"><span>Product Size</span><span>${size} bp</span></div>
                  <div class="verdict-row" style="margin-top:8px; border-top:1px solid rgba(255,255,255,0.1); padding-top:8px;">
                    <span>Size Constraint</span>
                    <span class="status-${ampliconVerdict.toLowerCase()}">${ampliconVerdict}</span>
                  </div>
                  <div style="font-size:0.75rem; color:var(--text-muted); margin-top:4px;">(Assay Range: ${thresholds.productSize.min}–${thresholds.productSize.max} bp)</div>
                </div>
              `;
            } else {
              pairFailures.push("Amplicon mapping failed (primers not found or oriented incorrectly)");
            }
          }

          // ── Exon-Junction Spanning ─────────────────
          let exonHtml = '';
          if (template && state.template.features && state.template.features.length > 0 && ampInfo && ampInfo.amplicon.size > 0) {
              const exonRes = BioMath.checkExonSpanning(
                ampInfo.amplicon.start, ampInfo.amplicon.start + fwdStr.length - 1,
                ampInfo.amplicon.end - revStr.length + 1, ampInfo.amplicon.end,
                state.template.features,
                { 
                  intronSpanningRequired: state.biologicalOptions.intronSpanning,
                  avoidExonJunction: state.biologicalOptions.avoidExonJunction
                }
              );

              if (state.biologicalOptions.intronSpanning || state.biologicalOptions.avoidExonJunction) {
                if (exonRes.verdict === 'FAIL') {
                  pairFailures.push(exonRes.reason);
                }
                const exonClass = exonRes.verdict === 'FAIL' ? 'warn' : 'pass';
                
                exonHtml = `
                  <div class="verdict-section-title" style="margin-top:20px;">Exon-Junction Spanning</div>
                  <div class="verdict-panel junction-panel ${exonClass}">
                    <div class="verdict-row"><span>Preference</span><span>${state.biologicalOptions.intronSpanning ? 'Intron-Spanning Required' : 'Avoid Junctions'}</span></div>
                    <div class="verdict-row"><span>Forward Primer</span><span class="${exonRes.fwdSpans?'pass':'text-muted'}">${exonRes.fwdSpans?'Spans Junction':'Within Exon'}</span></div>
                    <div class="verdict-row"><span>Reverse Primer</span><span class="${exonRes.revSpans?'pass':'text-muted'}">${exonRes.revSpans?'Spans Junction':'Within Exon'}</span></div>
                    <div class="verdict-row"><span>Amplicon</span><span class="${exonRes.ampliconSpans?'pass':'text-muted'}">${exonRes.ampliconSpans?'Spans Intron':'Single Exon'}</span></div>
                  </div>
                `;
              }
          }

          // ── Final Verdict Engine ───────────────────
          const finalVerdict = BioMath.generatePrimerPairVerdict(
            [results[0], results[1]],
            { tmDiff: diff, heteroDg: heteroDg },
            thresholds,
            template,
            state.template.features,
            state.biologicalOptions
          );

          const pairHtml = `
            ${specificityHtml}
            ${ampliconHtml}
            ${exonHtml}

            <div class="verdict-section-title" style="margin-top:24px;">OVERALL PAIR VERDICT</div>
            <div class="verdict-panel verdict-summary ${finalVerdict.class}">
              <div class="verdict-header">
                <span>Overall Status:</span>
                <strong>${finalVerdict.status}</strong>
              </div>
              
              <ul class="verdict-checklist">
                ${finalVerdict.checklist.map(c => `<li class="${c.pass?'check-pass':'check-fail'}"><span>${c.pass?'✓':'✗'}</span> ${c.label}</li>`).join('')}
              </ul>

              ${(finalVerdict.failures.length > 0 || finalVerdict.warnings.length > 0) ? `
                <div class="verdict-action-box">
                  <strong>Findings &amp; Actions:</strong>
                  <ul>
                    ${[...finalVerdict.failures, ...finalVerdict.warnings].map(err => `<li>${window.escapeHTML(err)}</li>`).join('')}
                  </ul>
                  <p style="margin-top:10px; font-weight:600; font-size:0.8rem;">
                    Action: ${finalVerdict.actionItems.length > 0 ? finalVerdict.actionItems.join(' ') : 'Review primer design and target site'}
                  </p>
                </div>
              ` : '<p style="margin-top:12px; font-size:0.85rem; color:var(--teal);">✓ All laboratory validation criteria met.</p>'}
            </div>

            <div class="export-container" style="margin-top:20px;">
              <button id="exportReportBtn" class="premium-btn-secondary" style="width:100%; justify-content:center;">
                📥 Download QC Report (JSON)
              </button>
            </div>`;

          state.lastAnalysis = {
            forward: results[0],
            reverse: results[1],
            verdict: finalVerdict,
            assay: assay,
            organism: organism
          };

          pairEl.innerHTML = window.sanitizeHTML(pairHtml);
          pairEl.style.display = 'block';

          // Bind export button
          const exportBtn = document.getElementById('exportReportBtn');
          if (exportBtn) exportBtn.addEventListener('click', exportPrimerQCReport);
        }
      });
    } catch (err) {
      window.handleError(err, 'Primer Analysis');
    }
  }

  function exportPrimerQCReport() {
    if (!state.lastAnalysis) return;

    const report = {
      timestamp: new Date().toISOString(),
      template: {
        accession: state.template.accession || "User Pasted Sequence",
        organism: state.lastAnalysis.organism || "Unknown",
        length: state.template.sequence?.length || 0,
        type: state.template.type
      },
      assayType: state.lastAnalysis.assay,
      forwardPrimer: state.lastAnalysis.forward,
      reversePrimer: state.lastAnalysis.reverse,
      verdict: state.lastAnalysis.verdict,
      constraints: state.constraints,
      meta: {
        model: "SantaLucia 1998 NN",
        buffer: {
          na: document.getElementById('saltConc').value,
          k: document.getElementById('kConc').value,
          mg: document.getElementById('mgConc').value,
          dntp: document.getElementById('dntpConc').value
        }
      }
    };

    const dataStr = JSON.stringify(report, null, 2);
    const dataBlob = new Blob([dataStr], { type: 'application/json' });
    const url = URL.createObjectURL(dataBlob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `primer_qc_report_${new Date().toISOString().slice(0, 10)}.json`;
    link.click();
    URL.revokeObjectURL(url);
  }

})();
