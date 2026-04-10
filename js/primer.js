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

  function gcFrac(seq) {
    return (seq.match(/[GC]/g) || []).length / seq.length;
  }
  function gcPct(seq) { return (gcFrac(seq) * 100).toFixed(1); }

  function hasGCClamp(seq) {
    const last3 = seq.slice(-3);
    const gcCount = (last3.match(/[GC]/g) || []).length;
    return gcCount >= 1 && gcCount <= 3;
  }

  function detectHairpin(seq) {
    // Look for self-complementary stem ≥ 4 bp with loop ≥ 3
    const comp = { A:'T', T:'A', G:'C', C:'G' };
    const minStem = 4, minLoop = 3;
    for (let stemLen = minStem; stemLen <= Math.floor((seq.length - minLoop) / 2); stemLen++) {
      const stem5 = seq.slice(0, stemLen);
      for (let loop = minLoop; loop <= seq.length - 2 * stemLen; loop++) {
        const stem3 = seq.slice(stemLen + loop, stemLen + loop + stemLen);
        const rc3   = stem3.split('').reverse().map(b => comp[b] || 'N').join('');
        if (stem5 === rc3) return { found: true, stemLen, loop };
      }
    }
    return { found: false };
  }

  function detectSelfDimer(seq) {
    const comp = { A:'T', T:'A', G:'C', C:'G' };
    const rc   = seq.split('').reverse().map(b => comp[b] || 'N').join('');
    let maxRun = 0, run = 0;
    for (let i = 0; i < seq.length; i++) {
      if (seq[i] === rc[i]) { run++; maxRun = Math.max(maxRun, run); } else run = 0;
    }
    return { found: maxRun >= 4, maxRun };
  }



  function complement(seq) {
    const comp = { A:'T', T:'A', G:'C', C:'G' };
    return seq.split('').map(b => comp[b] || 'N').join('');
  }

  function reverseComplement(seq) {
    return complement(seq).split('').reverse().join('');
  }

  function tmDiff(tm1, tm2) {
    return Math.abs(parseFloat(tm1) - parseFloat(tm2)).toFixed(1);
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
    const presets = window.BioKit.core.BioMath.ASSAY_PRESETS;
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
    hdr.addEventListener('click', (e) => {
      if (e.target.tagName === 'BUTTON') return;
      content.classList.toggle('hidden');
      hdr.classList.toggle('open');
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

  editBtn.addEventListener('click', () => {
    editForm.classList.remove('hidden');
    displayArea.classList.add('hidden');
  });
  cancelBtn.addEventListener('click', () => {
    editForm.classList.add('hidden');
    displayArea.classList.remove('hidden');
  });
  resetBtn.addEventListener('click', () => {
    loadAssayPreset(state.assayType);
    window.showToast(`✓ Reset to ${state.assayType.toUpperCase()} defaults`);
  });
  saveBtn.addEventListener('click', () => {
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
  templateInput.addEventListener('input', () => {
    const s = (templateInput.value || '').toUpperCase().replace(/[^A-Z]/gi, '');
    state.template.sequence = s;
    templateMeta.textContent = `Length: ${s.length} bp`;
  });

  document.getElementById('organism').addEventListener('change', (e) => {
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
  
  document.getElementById('avoidExonJunction').addEventListener('change', (e) => {
    state.biologicalOptions.avoidExonJunction = e.target.checked;
  });
  document.getElementById('intronSpanning').addEventListener('change', (e) => {
    state.biologicalOptions.intronSpanning = e.target.checked;
  });
  document.getElementById('excludeSnp').addEventListener('change', (e) => {
    state.biologicalOptions.excludeSnp = e.target.checked;
  });
  
  // Initialize Defaults
  loadAssayPreset('pcr');

  ['avoidExonJunction', 'intronSpanning', 'excludeSnp'].forEach(id => {
    const cb = document.getElementById(id);
    cb.addEventListener('change', () => state.biologicalOptions[id] = cb.checked);
  });

  // ── Fetch Handler (PROMPT 3) ─────────────────────────────
  const fetchBtn = document.getElementById('primerFetchBtn');
  fetchBtn.addEventListener('click', async () => {
    const db = document.getElementById('primerTemplateDb').value;
    const id = document.getElementById('primerTemplateId').value;
    if (!id) { window.showToast('⚠ Enter an Accession ID'); return; }

    fetchBtn.disabled = true;
    fetchBtn.textContent = 'Fetching...';
    try {
      const resp = await fetch('/api/fetch', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ 
          db: db,
          accession: id, 
          organism: state.template.organism 
        })
      });
      
      if (!resp.ok) {
        throw new Error('Accession not found or database error');
      }

      const data = await resp.json();
      if (data.sequence) {
        state.template.sequence = data.sequence.toUpperCase();
        state.template.features = data.features || [];
        
        templateInput.value = state.template.sequence;
        updateTemplateUI(data);

        // Automatically switch back to paste view to show result
        document.querySelector('input[name="templateSource"][value="paste"]').click();
        window.showToast('✓ Template sequence and features loaded');
      } else {
        window.showToast('⚠ Accession found but no sequence data returned');
      }
    } catch (err) {
      window.handleError(err, 'Fetch Template');
      window.showToast(`⚠ Error: ${err.message}`);
    } finally {
      fetchBtn.disabled = false;
      fetchBtn.textContent = 'Fetch';
    }
  });

  fwdInput.addEventListener('input', () => {
    const s = cleanSeq(fwdInput.value);
    const tmRes = calcTm(s, 50, 250);
    const validation = window.validateSequence(s, 'dna');
    let metaText = `Length: ${s.length} nt  |  GC: ${gcPct(s)}%  |  Tm ≈ ${tmRes.tm} °C`;
    if (!validation.valid || tmRes.warning) {
      metaText = `⚠ ${metaText} ${tmRes.warning ? '(Warning)' : '(Limit Exceeded)'}`;
      fwdMeta.style.color = 'var(--rose)';
    } else {
      fwdMeta.style.color = 'var(--text-muted)';
    }
    fwdMeta.textContent = metaText;
  });
  revInput.addEventListener('input', () => {
    const s = cleanSeq(revInput.value);
    const tmRes = calcTm(s, 50, 250);
    const validation = window.validateSequence(s, 'dna');
    let metaText = `Length: ${s.length} nt  |  GC: ${gcPct(s)}%  |  Tm ≈ ${tmRes.tm} °C`;
    if (!validation.valid || tmRes.warning) {
      metaText = `⚠ ${metaText} ${tmRes.warning ? '(Warning)' : '(Limit Exceeded)'}`;
      revMeta.style.color = 'var(--rose)';
    } else {
      revMeta.style.color = 'var(--text-muted)';
    }
    revMeta.textContent = metaText;
  });

  // ── Controls ───────────────────────────────────────────────
  document.getElementById('primerAnalyzeBtn').addEventListener('click', window.debounce(analyze, 300));
  document.getElementById('primerSampleBtn').addEventListener('click', () => {
    templateInput.value = "ATGGCTATCAAGCAGAAGTTTGATGCCATCAAGAAGCTGGAGGAGCAGCTGACCAAGGACATCCAGTACAACATGGGCCTGGCCGACATGGCCGGCATCGTGGTGCACGGCCACCACATCAAGAAGCTGTGA";
    fwdInput.value = 'ATGGCTATCAAGCAGAAGTTTG';
    revInput.value = 'TCACAGCTTCTTGATGTGGTG';
    templateInput.dispatchEvent(new Event('input'));
    fwdInput.dispatchEvent(new Event('input'));
    revInput.dispatchEvent(new Event('input'));
  });

  function analyze() {
    try {
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

      const fwd = window.cleanSeq(fwdInput.value);
      const rev = window.cleanSeq(revInput.value);

      if (!fwd && !rev) { window.showToast('⚠ Enter at least one primer'); return; }

      const valFwd = window.validateSequence(fwdInput.value, 'dna');
      const valRev = window.validateSequence(revInput.value, 'dna');

      if (!valFwd.valid || !valRev.valid) {
        window.showValidationWarning('panel-primer', !valFwd.valid ? valFwd : valRev);
        return;
      }

      const thresholds = state.constraints;

      const grid = document.getElementById('primerResultGrid');
      const primers = [];
      if (fwd.length >= 8) primers.push({ label: 'Forward Primer', seq: fwd });
      if (rev.length >= 8) primers.push({ label: 'Reverse Primer', seq: rev });

      window.withLoading('panel-primer', async () => {
        grid.textContent = ''; 
        const pairEl = document.getElementById('primerPairResults');
        pairEl.style.display = 'none';

        const results = await Promise.all(primers.map(async ({ label, seq }) => {
          const workerRes = await window.BioKit.core.WorkerManager.runTask({
            type: 'THERMO_ANALYSIS',
            strategy: 'auto',
            payload: { sequence: seq, ...options }
          });

          const tm = workerRes.tm_nn || calcTm(seq, options).tm;
          const gc = gcPct(seq);
          const len = seq.length;
          const clamp = hasGCClamp(seq);
          const hairpin = detectHairpin(seq);
          const dimer = detectSelfDimer(seq);
          const homoRes = BioMath.detectHomopolymerRuns(seq, thresholds.homopolymerLimit);

          const BioMath = window.BioKit.core.BioMath;

          // Validation using core engine
          const lenOk = BioMath.validateLength(seq, thresholds);
          const tmOk  = BioMath.validateTm(tm, thresholds);
          const gcOk  = BioMath.validateGC(gc, thresholds);
          const clampOk = BioMath.validateGCClamp(seq, thresholds);
          const homoOk  = !homoRes.hasViolation;

          let failReasons = [];
          if (!lenOk) failReasons.push(`Length (${len}nt) outside ${thresholds.length.min}-${thresholds.length.max} range`);
          if (!gcOk) failReasons.push(`GC (${gc}%) outside ${thresholds.gc.min}-${thresholds.gc.max}% range`);
          if (!tmOk) failReasons.push(`Tm (${tm}°C) outside ${thresholds.tm.min}-${thresholds.tm.max}°C range`);
          if (hairpin.found) failReasons.push(`Hairpin risk detected`);
          if (!homoOk) failReasons.push(`Homopolymer detected: ${homoRes.summary}`);
          if (!clampOk) failReasons.push(`Missing 3' GC clamp (G or C at end)`);

          const cardHtml = `
            <div class="primer-card ${failReasons.length > 2 ? 'fail-border' : failReasons.length > 0 ? 'warn-border' : ''}">
              <div class="primer-card-header">
                <h4>${escapeHTML(label)}</h4>
                <span class="qc-status ${failReasons.length === 0 ? 'pass' : 'fail'}">${failReasons.length === 0 ? 'QC PASS' : 'QC FLAG'}</span>
              </div>
              <div class="primer-seq">${escapeHTML(seq)}</div>
              
              <div class="primer-stat-grid">
                <div class="primer-stat-item"><span class="primer-stat-key">Length</span>
                  <span class="primer-stat-val ${lenOk?'pass':'warn'}">${len} nt</span></div>
                <div class="primer-stat-item"><span class="primer-stat-key">Tm (SantaLucia NN)</span>
                  <span class="primer-stat-val ${tmOk?'pass':'warn'}">${tm} °C</span></div>
                <div class="primer-stat-item"><span class="primer-stat-key">GC Content</span>
                  <span class="primer-stat-val ${gcOk?'pass':'warn'}">${gc}%</span></div>
                <div class="primer-stat-item"><span class="primer-stat-key">3′ GC Clamp</span>
                  <span class="primer-stat-val ${clamp?'pass':'warn'}">${clamp?'Yes':'No'}</span></div>
              </div>

              <div class="structural-report">
                <div class="struct-stat ${hairpin.found?'fail':'pass'}">Hairpin: ${hairpin.found?`Fail` : 'Pass'}</div>
                <div class="struct-stat ${dimer.maxRun >= 4?'warn':'pass'}">Self-Dimer: ${dimer.maxRun}nt</div>
                <div class="struct-stat ${homoOk?'pass':'fail'}">Homopolymer: ${homoRes.hasViolation ? homoRes.summary : 'None'}</div>
              </div>

              ${failReasons.length > 0 ? `
              <div class="fail-list">
                <strong>Reasons for Flag:</strong>
                <ul>${failReasons.map(r => `<li>${escapeHTML(r)}</li>`).join('')}</ul>
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
          return { label, seq, tm };
        }));

          // ── Pair Assessment & Verdict (PROMPT 8) ──────────────────
          const pairFailures = [];
          const pairWarnings = [];
          
          if (!diffOk) pairWarnings.push(`Tm difference (${diff}°C) exceeds assay tolerance (max ${thresholds.tmDifference}°C)`);
          if (maxHetero >= 5) pairFailures.push(`High heterodimer risk (${maxHetero}bp complementarity)`);
          else if (maxHetero >= 4) pairWarnings.push(`Moderate heterodimer risk (${maxHetero}bp)`);

          const fwdResults = results[0];
          const revResults = results[1];
          const fwdStats = { 
            tm: parseFloat(fwdResults.tm), 
            len: fwdResults.seq.length,
            homo: BioMath.validateHomopolymer(fwdResults.seq, thresholds)
          };
          const revStats = { 
            tm: parseFloat(revResults.tm), 
            len: revResults.seq.length,
            homo: BioMath.validateHomopolymer(revResults.seq, thresholds)
          };

          // ── Specificity Section (PROMPT 10) ──────────────────────
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
          if (template) {
            const ampInfo = BioMath.calculateAmpliconInfo(template, fwdStr, revStr, thresholds);
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

          // ── Exon-Junction Spanning (PROMPT 11) ─────────────────
          let exonHtml = '';
          if (template && state.template.features && state.template.features.length > 0) {
            const ampInfo = BioMath.calculateAmpliconInfo(template, fwdStr, revStr, thresholds);
            if (ampInfo.amplicon.size > 0) {
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
          }

          // ── Final Verdict Engine (PROMPT 12) ───────────────────
          const finalVerdict = BioMath.generatePrimerPairVerdict(
            [fwdResults, revResults],
            { tmDiff: diff, maxHetero: maxHetero },
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
                    ${[...finalVerdict.failures, ...finalVerdict.warnings].map(err => `<li>${escapeHTML(err)}</li>`).join('')}
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
            forward: fwdResults,
            reverse: revResults,
            verdict: finalVerdict,
            assay: assay,
            organism: organism
          };

          pairEl.innerHTML = window.sanitizeHTML(pairHtml);
          pairEl.style.display = 'block';

          // Bind export button
          document.getElementById('exportReportBtn').addEventListener('click', exportPrimerQCReport);
        }
      }, Math.max(fwd.length, rev.length));
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
          na: document.getElementById('na_conc').value,
          mg: document.getElementById('mg_conc').value,
          dntp: document.getElementById('dntp_conc').value
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
