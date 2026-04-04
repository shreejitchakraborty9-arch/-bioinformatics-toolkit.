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

  function tmDiff(tm1, tm2) {
    return Math.abs(parseFloat(tm1) - parseFloat(tm2)).toFixed(1);
  }

  function rating(value, pass, warn) {
    if (pass) return 'pass';
    if (warn) return 'warn';
    return 'fail';
  }

  // ── Input live meta ────────────────────────────────────────
  const fwdInput  = document.getElementById('fwdPrimer');
  const revInput  = document.getElementById('revPrimer');
  const fwdMeta   = document.getElementById('fwdMeta');
  const revMeta   = document.getElementById('revMeta');

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
    fwdInput.value = 'GCTTAAGCGCTATCAAGCCATG';
    revInput.value = 'TTCAGCTTGATGGCATCAAAGCC';
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

      const grid = document.getElementById('primerResultGrid');
      const primers = [];
      if (fwd.length >= 10) primers.push({ label: 'Forward Primer', seq: fwd });
      if (rev.length >= 10) primers.push({ label: 'Reverse Primer', seq: rev });

      window.withLoading('panel-primer', () => {
        grid.textContent = ''; 
        primers.forEach(({ label, seq }) => {
          const tmRes = calcTm(seq, options);
          const tm = tmRes.tm;
          const gc = gcPct(seq);
          const len = seq.length;
          const clamp = hasGCClamp(seq);
          const hairpin = detectHairpin(seq);
          const dimer = detectSelfDimer(seq);

          const lenOk = len >= 18 && len <= 25;
          const gcOk = parseFloat(gc) >= 40 && parseFloat(gc) <= 60;
          const tmOk = parseFloat(tm) >= 52 && parseFloat(tm) <= 68;

          const cardHtml = `
            <div class="primer-card">
              <h4>${escapeHTML(label)}</h4>
              <div class="primer-seq">${escapeHTML(seq)}</div>
              <div class="primer-stat-item"><span class="primer-stat-key">Length</span>
                <span class="primer-stat-val ${lenOk?'pass':len<15||len>30?'fail':'warn'}">${escapeHTML(len)} nt ${lenOk?'✓':len<18?'(too short)':'(too long)'}</span></div>
              <div class="primer-stat-item"><span class="primer-stat-key">Tm (NN model)</span>
                <span class="primer-stat-val ${tmOk?'pass':Math.abs(parseFloat(tm)-60)>10?'fail':'warn'}">${escapeHTML(tm)} °C ${tmRes.confidence ? `<small style="opacity:0.6"> ${escapeHTML(tmRes.confidence)}</small>` : ''}</span></div>
              ${tmRes.warning ? `<div style="font-size:0.75rem; color:var(--rose); margin-top:-4px; margin-bottom:8px;">⚠ ${escapeHTML(tmRes.warning)}</div>` : ''}
              <div class="primer-stat-item"><span class="primer-stat-key">GC Content</span>
                <span class="primer-stat-val ${gcOk?'pass':'warn'}">${escapeHTML(gc)}% ${gcOk?'✓':'(ideal: 40–60%)'}</span></div>
              <div class="primer-stat-item"><span class="primer-stat-key">GC Clamp (3′ end)</span>
                <span class="primer-stat-val ${clamp?'pass':'warn'}">${clamp?'✓ Good':'⚠ Weak 3′ end'}</span></div>
              <div class="primer-stat-item"><span class="primer-stat-key">Hairpin Risk</span>
                <span class="primer-stat-val ${hairpin.found?'fail':'pass'}">${hairpin.found?`⚠ Detected (stem ${escapeHTML(hairpin.stemLen)}bp)`:'✓ None'}</span></div>
              <div class="primer-stat-item"><span class="primer-stat-key">Self-Dimer Risk</span>
                <span class="primer-stat-val ${dimer.found?'warn':'pass'}">${dimer.found?`⚠ ${escapeHTML(dimer.maxRun)}nt run`:'✓ Low'}</span></div>
            </div>`;
          grid.innerHTML += window.sanitizeHTML(cardHtml);
        });

        const pairEl = document.getElementById('primerPairResults');
        if (primers.length === 2) {
          const tm1Res = calcTm(primers[0].seq, options);
          const tm2Res = calcTm(primers[1].seq, options);
          const diff = tmDiff(tm1Res.tm, tm2Res.tm);
          const diffOk = parseFloat(diff) <= 5;

          const comp = { A:'T', T:'A', G:'C', C:'G' };
          const rc2  = primers[1].seq.split('').reverse().map(b=>comp[b]||'N').join('');
          let maxHetero = 0, run = 0;
          for (let i = 0; i < Math.min(primers[0].seq.length, rc2.length); i++) {
            if (primers[0].seq[i] === rc2[i]) { run++; maxHetero = Math.max(maxHetero, run); } else run = 0;
          }
          const heteroRisk = maxHetero >= 4;

          const pairHtml = `
            <h4>🔗 Primer Pair Analysis</h4>
            <div class="primer-stat-item"><span class="primer-stat-key">Tm Difference</span>
              <span class="primer-stat-val ${diffOk?'pass':'warn'}">${escapeHTML(diff)} °C ${diffOk?'✓':'(ideal ≤ 5°C)'}</span></div>
            <div class="primer-stat-item"><span class="primer-stat-key">Heterodimer Risk</span>
              <span class="primer-stat-val ${heteroRisk?'warn':'pass'}">${heteroRisk?`⚠ ${escapeHTML(maxHetero)}nt complementarity`:'✓ Low'}</span></div>
            <div class="primer-stat-item"><span class="primer-stat-key">Forward Tm</span>
              <span class="primer-stat-val">${escapeHTML(tm1Res.tm)} °C</span></div>
            <div class="primer-stat-item"><span class="primer-stat-key">Reverse Tm</span>
              <span class="primer-stat-val">${escapeHTML(tm2Res.tm)} °C</span></div>`;
          pairEl.innerHTML = window.sanitizeHTML(pairHtml);
          pairEl.style.display = 'block';
        } else {
          pairEl.style.display = 'none';
        }
      }, Math.max(fwd.length, rev.length));
    } catch (err) {
      window.handleError(err, 'Primer Analysis');
    }
  }
})();
