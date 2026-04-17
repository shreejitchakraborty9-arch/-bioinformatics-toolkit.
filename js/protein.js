/* ============================================================
   protein.js — Protein Analyzer
   Phase 1 (crash fix) + Phase 2 (feature wiring) + Phase 3 (pI)
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };
  const BioMath = window.BioKit.core.BioMath;

  let currentSequence = "";
  let lastMetrics = null;

  function getThemeColors() {
    const style = getComputedStyle(document.documentElement);
    return {
      bg: style.getPropertyValue('--bg-card').trim() || '#17191c',
      text: style.getPropertyValue('--text-primary').trim() || '#f8f9fa',
      muted: style.getPropertyValue('--text-muted').trim() || '#868e96',
      border: style.getPropertyValue('--border').trim() || '#2c3035',
      teal: style.getPropertyValue('--teal').trim() || '#14b8a6'
    };
  }

  // ── Amino Acid Data ────────────────────────────────────
  const AA_NAMES = {
    A:'Ala',R:'Arg',N:'Asn',D:'Asp',C:'Cys',E:'Glu',Q:'Gln',G:'Gly',
    H:'His',I:'Ile',L:'Leu',K:'Lys',M:'Met',F:'Phe',P:'Pro',S:'Ser',
    T:'Thr',W:'Trp',Y:'Tyr',V:'Val'
  };

  const AA_COLORS = {
    A:'#94a3b8', G:'#94a3b8', V:'#94a3b8', L:'#94a3b8', I:'#94a3b8', M:'#94a3b8', P:'#94a3b8', // Non-polar (Slate)
    S:'#4ade80', T:'#4ade80', C:'#4ade80', N:'#4ade80', Q:'#4ade80', // Polar (Green)
    D:'#f87171', E:'#f87171', // Acidic (Red)
    K:'#60a5fa', R:'#60a5fa', H:'#60a5fa', // Basic (Blue)
    F:'#c084fc', W:'#c084fc', Y:'#c084fc'  // Aromatic (Purple)
  };

  const AA_GROUPS = {
    A: 'Non-polar', G: 'Non-polar', V: 'Non-polar', L: 'Non-polar', I: 'Non-polar', M: 'Non-polar', P: 'Non-polar',
    S: 'Polar', T: 'Polar', C: 'Polar', N: 'Polar', Q: 'Polar',
    D: 'Acidic', E: 'Acidic',
    K: 'Basic', R: 'Basic', H: 'Basic',
    F: 'Aromatic', W: 'Aromatic', Y: 'Aromatic'
  };

  // ── AA Composition ────────────────────────────────────
  function getAAComposition(seq) {
    const comp = {};
    for (const aa of seq) {
      if (AA_NAMES[aa]) comp[aa] = (comp[aa] || 0) + 1;
    }
    return Object.entries(comp).sort((a, b) => b[1] - a[1]);
  }

  // ── Main Run ──────────────────────────────────────────
  function run(silent) {
    const input = document.getElementById('proteinInput');
    if (!input) return;

    const raw = input.value;

    if (silent && raw.replace(/[^A-Za-z*]/g, '').length < 4) return;

    const validation = window.validateSequence(raw, 'protein');

    if (!validation.valid) {
      window.showValidationWarning('panel-protein', validation);
      return;
    }

    // Read scientific parameters from UI toggles
    const massType  = document.getElementById('proteinMassType')?.value  || 'average';
    const redoxState = document.getElementById('proteinRedoxState')?.value || 'reduced';

    const seq = validation.clean.replace(/\*/g, '');
    currentSequence = seq;

    window.withLoading('panel-protein', () => {
      const mwResult  = BioMath.calculateProtein_MW(seq, massType);
      const hydro     = BioMath.calculateHydrophobicity(seq);
      const pI        = BioMath.calculatePI(seq);
      const extCoeff  = BioMath.calculateExtinctionCoefficient(seq, 'protein', redoxState);
      const aliphatic   = BioMath.calculateAliphaticIndex(seq);
      const instability = BioMath.calculateInstabilityIndex(seq);
      const instabLabel = instability < 40 ? ' (Stable)' : ' (Unstable)';
      const netCharge74 = BioMath.calculateNetChargeAtPhysiological(seq);

      lastMetrics = {
        sequence: seq,
        mwResult,
        massType,
        pI,
        extCoeff,
        redoxState,
        netCharge74,
        aaComposition: getAAComposition(seq)
      };

      // Internal stop-codon warning: report the exact truncated length
      if (mwResult.internalStop) {
        const warn = { valid: true, msg: `Internal stop codon (*) detected. Sequence truncated to ${mwResult.truncatedLength} aa for mass calculation.` };
        window.showValidationWarning('panel-protein', warn);
      }

      // Ambiguity warning banner
      const ambiguityBanner = document.getElementById('proteinAmbiguityWarning');
      if (ambiguityBanner) {
        if (mwResult.hasAmbiguousResidues) {
          ambiguityBanner.textContent = `⚠ Contains ambiguous residues (X/B/Z); MW precision reduced (±${mwResult.errorMarginKDa} kDa).`;
          ambiguityBanner.classList.remove('hidden');
        } else {
          ambiguityBanner.classList.add('hidden');
        }
      }

      const massLabel = massType === 'monoisotopic' ? 'MW (Mono, kDa)' : 'Mol. Weight (kDa)';
      const ecLabel   = redoxState === 'oxidized' ? 'Ext. Coeff [Ox] (M⁻¹cm⁻¹)' : 'Ext. Coeff [Red] (M⁻¹cm⁻¹)';
      const mwDisplay = mwResult.hasAmbiguousResidues
        ? `${mwResult.kDa} ±${mwResult.errorMarginKDa} kDa`
        : `${mwResult.kDa} kDa`;

      // Stat cards
      window.buildStatCards('proteinStatRow', [
        { val: seq.length.toLocaleString(),      label: 'Length (aa)' },
        { val: mwDisplay,                        label: massLabel },
        { val: pI,                               label: 'Isoelectric Pt (pI)' },
        { val: extCoeff.toLocaleString(),        label: ecLabel },
        { val: hydro,                            label: 'GRAVY Index' },
        { val: aliphatic,                        label: 'Aliphatic Index' },
        { val: instability + instabLabel,        label: 'Instability Index' },
        { val: netCharge74,                      label: 'Net Charge (pH 7.4)' }
      ]);

      // AA Composition chart
      renderAAChart(seq);

      // AA Table
      renderAATable(seq);

      // Hydrophobicity plot — window size from UI selector (default 7)
      const winEl = document.getElementById('hydroWindowSize');
      const winSize = winEl ? parseInt(winEl.value, 10) : 7;
      renderHydroPlot(seq, winSize);

      // Titration curve
      renderTitrationCurve(seq);

      // Show results, hide empty state
      const empty = document.getElementById('proteinEmptyState');
      const content = document.querySelector('#proteinResults .results-content');
      if (empty) empty.classList.add('hidden');
      if (content) content.classList.remove('hidden');

      if (!silent) window.showToast('Analysis Complete');
    }, seq.length);
  }

  // ── AA Composition Bar Chart ──────────────────────────
  function renderAAChart(seq) {
    const canvas = document.getElementById('aaCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const comp = getAAComposition(seq);
    if (comp.length === 0) return;

    canvas.width = canvas.offsetWidth || 700;
    canvas.height = 320;
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    const W = canvas.width, H = canvas.height;
    const pad = { top: 20, right: 20, bottom: 50, left: 50 };
    const plotW = W - pad.left - pad.right;
    const plotH = H - pad.top - pad.bottom;
    const colors = getThemeColors();
    const barW = Math.floor(plotW / comp.length) - 4;

    ctx.fillStyle = colors.bg;
    ctx.fillRect(0, 0, W, H);

    const maxVal = Math.max(...comp.map(c => c[1]));

    comp.forEach(([aa, cnt], i) => {
      const x = pad.left + i * (barW + 4) + 2;
      const h = (cnt / maxVal) * plotH;
      const y = pad.top + plotH - h;

      ctx.fillStyle = AA_COLORS[aa] || '#94a3b8';
      ctx.fillRect(x, y, barW, h);

      ctx.fillStyle = colors.text;
      ctx.font = '10px "JetBrains Mono"';
      ctx.textAlign = 'center';
      if (h > 14) ctx.fillText(cnt, x + barW / 2, y - 3);

      ctx.fillStyle = colors.muted;
      ctx.font = '11px "JetBrains Mono"';
      ctx.fillText(aa, x + barW / 2, H - pad.bottom + 14);
    });
  }

  // ── AA Composition Table ──────────────────────────────
  function renderAATable(seq) {
    const wrap = document.getElementById('aaTableWrap');
    if (!wrap) return;
    const comp = getAAComposition(seq);
    const total = seq.length;

    let html = `<table class="aa-table"><thead><tr>
      <th>AA</th><th>Name</th><th>Group</th><th>Count</th><th>%</th>
    </tr></thead><tbody>`;
    comp.forEach(([aa, cnt]) => {
      const pct = ((cnt / total) * 100).toFixed(1);
      const group = AA_GROUPS[aa] || 'Other';
      const color = AA_COLORS[aa] || '#fff';
      html += `<tr>
        <td style="color:${color};font-weight:700;">${aa}</td>
        <td>${AA_NAMES[aa]||aa}</td>
        <td><span style="font-size:0.7rem;opacity:0.8;">${group}</span></td>
        <td>${cnt}</td>
        <td>${pct}%</td>
      </tr>`;
    });
    html += '</tbody></table>';
    wrap.innerHTML = window.sanitizeHTML(html);
  }

  // ── Hydrophobicity Plot ───────────────────────────────
  function renderHydroPlot(seq, windowSize) {
    windowSize = (windowSize >= 5 && windowSize <= 21) ? windowSize : 7;
    const canvas = document.getElementById('hydroCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');

    canvas.width = canvas.offsetWidth || 700;
    canvas.height = 220;
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    const W = canvas.width, H = canvas.height;
    const pad = { top: 28, right: 20, bottom: 30, left: 40 };
    const plotW = W - pad.left - pad.right;
    const plotH = H - pad.top - pad.bottom;

    // Compute sliding window
    const values = [];
    for (let i = 0; i <= seq.length - windowSize; i++) {
      let sum = 0;
      for (let j = i; j < i + windowSize; j++) {
        const val = BioMath.KD_SCALE[seq[j]];
        sum += val || 0;
      }
      values.push(sum / windowSize);
    }

    if (values.length === 0) return;

    const minVal = Math.min(...values, -2);
    const maxVal = Math.max(...values, 2);
    const colors = getThemeColors();
    const range = maxVal - minVal || 1;

    // Background
    ctx.fillStyle = colors.bg;
    ctx.fillRect(0, 0, W, H);

    // Window label
    ctx.fillStyle = colors.muted;
    ctx.font = '10px Inter';
    ctx.textAlign = 'left';
    ctx.fillText(`Kyte-Doolittle  window = ${windowSize}`, pad.left, 14);

    // Zero line
    const zeroY = pad.top + plotH - ((-minVal) / range) * plotH;
    ctx.strokeStyle = colors.border;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(pad.left, zeroY);
    ctx.lineTo(W - pad.right, zeroY);
    ctx.stroke();

    ctx.fillStyle = colors.muted;
    ctx.font = '10px Inter';
    ctx.textAlign = 'right';
    ctx.fillText('0', pad.left - 5, zeroY + 3);

    // Draw line
    ctx.beginPath();
    ctx.strokeStyle = colors.teal;
    ctx.lineWidth = 1.5;
    values.forEach((v, i) => {
      const x = pad.left + (i / values.length) * plotW;
      const y = pad.top + plotH - ((v - minVal) / range) * plotH;
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();

    // Fill hydrophobic region
    ctx.globalAlpha = 0.15;
    ctx.beginPath();
    ctx.moveTo(pad.left, zeroY);
    values.forEach((v, i) => {
      if (v > 0) {
        const x = pad.left + (i / values.length) * plotW;
        const y = pad.top + plotH - ((v - minVal) / range) * plotH;
        ctx.lineTo(x, y);
      }
    });
    ctx.lineTo(W - pad.right, zeroY);
    ctx.fillStyle = '#f59e0b';
    ctx.fill();
    ctx.globalAlpha = 1;
  }

  // ── Titration Curve (Charge vs. pH) ──────────────────
  function renderTitrationCurve(seq) {
    const canvas = document.getElementById('titrationCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');

    canvas.width = canvas.offsetWidth || 700;
    canvas.height = 220;
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    const W = canvas.width, H = canvas.height;
    const pad = { top: 20, right: 20, bottom: 36, left: 52 };
    const plotW = W - pad.left - pad.right;
    const plotH = H - pad.top - pad.bottom;
    const colors = getThemeColors();

    // Generate data: pH 0 to 14 at 0.5 intervals
    const points = [];
    for (let ph = 0; ph <= 14; ph += 0.5) {
      points.push({ ph, charge: BioMath.calculateNetChargeAtPH(seq, ph) });
    }

    const charges = points.map(p => p.charge);
    const minC = Math.min(...charges);
    const maxC = Math.max(...charges);
    const rangeC = (maxC - minC) || 1;

    const toX = (ph) => pad.left + (ph / 14) * plotW;
    const toY = (c)  => pad.top + plotH - ((c - minC) / rangeC) * plotH;

    // Background
    ctx.fillStyle = colors.bg;
    ctx.fillRect(0, 0, W, H);

    // Zero charge line
    const zeroY = toY(0);
    ctx.strokeStyle = colors.border;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(pad.left, zeroY);
    ctx.lineTo(W - pad.right, zeroY);
    ctx.stroke();
    ctx.setLineDash([]);

    // pH 7.4 marker
    const x74 = toX(7.4);
    ctx.strokeStyle = '#f59e0b';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(x74, pad.top);
    ctx.lineTo(x74, pad.top + plotH);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#f59e0b';
    ctx.font = '9px Inter';
    ctx.textAlign = 'center';
    ctx.fillText('pH 7.4', x74, pad.top + plotH + 20);

    // Curve
    ctx.beginPath();
    ctx.strokeStyle = colors.teal;
    ctx.lineWidth = 2;
    points.forEach((p, i) => {
      const x = toX(p.ph);
      const y = toY(p.charge);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();

    // Axes
    ctx.fillStyle = colors.muted;
    ctx.font = '10px Inter';
    ctx.textAlign = 'center';

    // X axis labels (pH)
    [0, 2, 4, 6, 8, 10, 12, 14].forEach(ph => {
      ctx.fillText(ph, toX(ph), pad.top + plotH + 14);
    });
    ctx.fillText('pH', pad.left + plotW / 2, H - 2);

    // Y axis labels (charge)
    ctx.textAlign = 'right';
    [minC, 0, maxC].forEach(c => {
      ctx.fillText(c.toFixed(1), pad.left - 5, toY(c) + 3);
    });
  }

  // ── CSV Export ─────────────────────────────────────────
  function exportMetricsToCSV() {
    if (!lastMetrics) return;

    const { sequence, mwResult, massType, pI, extCoeff, redoxState, netCharge74, aaComposition } = lastMetrics;
    const total = sequence.length;
    const ts = new Date().toISOString().replace('T', ' ').slice(0, 19);

    const rows = [];
    rows.push(['# BioToolkit Protein Analysis Export']);
    rows.push([`# Date,${ts}`]);
    rows.push([`# Sequence Length,${total} aa`]);
    rows.push(['']);
    rows.push(['Metric', 'Value', 'Unit']);
    rows.push([`Molecular Weight (${massType === 'monoisotopic' ? 'Monoisotopic' : 'Average'})`, mwResult.kDa, 'kDa']);
    rows.push(['Isoelectric Point (pI)', pI, '']);
    rows.push([`Extinction Coefficient (${redoxState === 'oxidized' ? 'Oxidized' : 'Reduced'})`, extCoeff, 'M-1cm-1']);
    rows.push(['Net Charge (pH 7.4)', netCharge74, '']);
    rows.push(['Sequence Length', total, 'aa']);
    rows.push(['']);
    rows.push(['AA', 'Name', 'Group', 'Count', 'Percentage']);
    aaComposition.forEach(([aa, cnt]) => {
      const pct = ((cnt / total) * 100).toFixed(2);
      rows.push([aa, AA_NAMES[aa] || aa, AA_GROUPS[aa] || 'Other', cnt, `${pct}%`]);
    });

    const csvContent = rows.map(r => r.join(',')).join('\r\n');
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `protein_analysis_${ts.replace(/[: ]/g, '-')}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  // ── UniProt Biological Context (Step 6) ───────────────

  function clearUniProtCard() {
    const card = document.getElementById('uniprotMetaCard');
    if (card) card.classList.add('hidden');
  }

  function renderUniProtCard(data) {
    const card = document.getElementById('uniprotMetaCard');
    if (!card) return;

    const setText = (id, val) => {
      const el = document.getElementById(id);
      if (el) el.textContent = val || '';
    };
    const setDisplay = (id, show) => {
      const el = document.getElementById(id);
      if (el) el.style.display = show ? '' : 'none';
    };

    setText('uniprotMetaAccessionBadge', data.accession || '');
    setText('uniprotMetaProteinName', data.protein_name || '');
    setText('uniprotMetaOrganism', data.organism || '');

    // Subcellular locations
    const locTags = document.getElementById('uniprotMetaLocTags');
    if (locTags && data.subcellular_locations && data.subcellular_locations.length) {
      locTags.innerHTML = data.subcellular_locations.map(loc =>
        `<span style="background:rgba(20,184,166,0.12);border:1px solid rgba(20,184,166,0.3);color:var(--teal);border-radius:3px;padding:2px 8px;font-size:0.7rem;">${window.sanitizeHTML ? window.sanitizeHTML(loc) : loc}</span>`
      ).join('');
      setDisplay('uniprotMetaLocBlock', true);
    } else {
      setDisplay('uniprotMetaLocBlock', false);
    }

    // PTM descriptions
    const ptmList = document.getElementById('uniprotMetaPTMList');
    if (ptmList && data.ptms && data.ptms.length) {
      ptmList.innerHTML = data.ptms.map(p =>
        `<div style="padding:2px 0 2px 8px; border-left:2px solid var(--border);">• ${window.sanitizeHTML ? window.sanitizeHTML(p) : p}</div>`
      ).join('');
      setDisplay('uniprotMetaPTMBlock', true);
    } else {
      setDisplay('uniprotMetaPTMBlock', false);
    }

    // Structural features (TM helices, signal peptides, etc.)
    const featTags = document.getElementById('uniprotMetaFeatTags');
    if (featTags && data.features && data.features.length) {
      const FEAT_COLORS = {
        'Transmembrane':    { bg: 'rgba(139,92,246,0.12)', border: 'rgba(139,92,246,0.35)', text: '#a78bfa' },
        'Signal':           { bg: 'rgba(245,158,11,0.12)', border: 'rgba(245,158,11,0.35)', text: '#f59e0b' },
        'Modified residue': { bg: 'rgba(248,113,113,0.12)',border: 'rgba(248,113,113,0.35)', text: '#f87171' },
        'Glycosylation':    { bg: 'rgba(52,211,153,0.12)', border: 'rgba(52,211,153,0.35)', text: '#34d399' },
        'Disulfide bond':   { bg: 'rgba(251,191,36,0.12)', border: 'rgba(251,191,36,0.35)', text: '#fbbf24' },
        'Lipidation':       { bg: 'rgba(96,165,250,0.12)', border: 'rgba(96,165,250,0.35)', text: '#60a5fa' },
      };
      featTags.innerHTML = data.features.map(f => {
        const c = FEAT_COLORS[f.type] || { bg: 'rgba(148,163,184,0.1)', border: 'rgba(148,163,184,0.3)', text: '#94a3b8' };
        const range = (f.start != null && f.end != null) ? ` ${f.start}–${f.end}` : '';
        const desc = f.description ? `: ${f.description}` : '';
        const label = `${f.type}${range}${desc}`;
        return `<span style="background:${c.bg};border:1px solid ${c.border};color:${c.text};border-radius:3px;padding:2px 7px;font-size:0.68rem;" title="${f.type}">${window.sanitizeHTML ? window.sanitizeHTML(label) : label}</span>`;
      }).join('');
      setDisplay('uniprotMetaFeatBlock', true);
    } else {
      setDisplay('uniprotMetaFeatBlock', false);
    }

    card.classList.remove('hidden');
  }

  function fetchAndDisplayUniProtMeta(accessionId) {
    fetch('/api/uniprot/' + encodeURIComponent(accessionId))
      .then(r => r.json())
      .then(data => {
        if (data.error) { clearUniProtCard(); return; }
        renderUniProtCard(data);
      })
      .catch(() => clearUniProtCard());
  }

  // Hook into the protein panel fetch button — intercept UniProt fetches only
  (function () {
    const panel = document.getElementById('panel-protein');
    if (!panel) return;
    const fetchBtn  = panel.querySelector('.fetch-btn');
    const dbSelect  = panel.querySelector('.fetch-db-select');
    const idInput   = panel.querySelector('.fetch-id-input');
    if (!fetchBtn) return;

    fetchBtn.addEventListener('click', () => {
      const db = dbSelect?.value;
      const id = idInput?.value.trim();
      if (db === 'uniprot' && id) {
        fetchAndDisplayUniProtMeta(id);
      } else {
        clearUniProtCard();
      }
    });
  })();

  // ── Events ─────────────────────────────────────────────
  document.getElementById('proteinAnalyzeBtn')?.addEventListener('click', () => run());

  document.getElementById('proteinSampleBtn')?.addEventListener('click', () => {
    const input = document.getElementById('proteinInput');
    if (input) {
      input.value = 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK';
      input.dispatchEvent(new Event('input'));
    }
  });

  document.getElementById('proteinClearBtn')?.addEventListener('click', () => {
    const input = document.getElementById('proteinInput');
    if (input) { input.value = ''; input.dispatchEvent(new Event('input')); }
    const label = document.getElementById('uniprotNameLabel');
    if (label) { label.textContent = ''; label.classList.add('hidden'); }
    clearUniProtCard();
    const empty = document.getElementById('proteinEmptyState');
    const content = document.querySelector('#proteinResults .results-content');
    if (empty) empty.classList.remove('hidden');
    if (content) content.classList.add('hidden');
  });

  // ── Clean Sequence Button ──────────────────────────────
  // Triggers the FASTA header stripper and collapses all whitespace/numbers on demand.
  document.getElementById('proteinCleanBtn')?.addEventListener('click', () => {
    const input = document.getElementById('proteinInput');
    if (!input) return;
    const cleaned = (window.stripFastaHeader ? window.stripFastaHeader(input.value) : input.value)
      .replace(/[^ACDEFGHIKLMNPQRSTVWYXUZ*]/gi, '');
    input.value = cleaned;
    input.dispatchEvent(new Event('input'));
    window.showToast('✨ Sequence cleaned');
  });

  // ── Paste: auto-strip FASTA headers ───────────────────
  // Intercepts paste events on the protein input and strips FASTA header lines
  // before they land in the textarea, preventing false "invalid characters" rejections.
  document.getElementById('proteinInput')?.addEventListener('paste', (e) => {
    if (!window.stripFastaHeader) return; // guard — app.js not yet loaded
    const raw = e.clipboardData?.getData('text') || '';
    if (!raw.trimStart().startsWith('>')) return; // plain sequence — don't interfere
    e.preventDefault();
    const stripped = window.stripFastaHeader(raw);
    const input = e.target;
    // Insert at cursor position for correct undo/redo behaviour
    const start = input.selectionStart ?? 0;
    const end   = input.selectionEnd   ?? 0;
    const current = input.value;
    input.value = current.slice(0, start) + stripped + current.slice(end);
    input.selectionStart = input.selectionEnd = start + stripped.length;
    input.dispatchEvent(new Event('input'));
  });

  // ── Live meta counter ──────────────────────────────────
  // Bound strictly to the sanitized string length so "Length: N aa" always
  // reflects what the solver will actually receive — eliminates ghost-lengths.
  const pInput = document.getElementById('proteinInput');
  const pMeta  = document.getElementById('proteinInputMeta');
  if (pInput && pMeta) {
    pInput.addEventListener('input', () => {
      const stripped = window.stripFastaHeader ? window.stripFastaHeader(pInput.value) : pInput.value;
      const seq = stripped.replace(/[^ACDEFGHIKLMNPQRSTVWYXUZ*]/gi, '');
      pMeta.textContent = `Length: ${seq.length} aa`;
    });
  }

  // ── Reactive debounced input (500 ms) ─────────────────
  (function () {
    const el = document.getElementById('proteinInput');
    if (!el) return;
    let timer;
    el.addEventListener('input', () => {
      clearTimeout(timer);
      timer = setTimeout(() => run(true), 500);
    });
  })();

  document.getElementById('proteinExportCSVBtn')?.addEventListener('click', exportMetricsToCSV);

  // ── Hydropathy window selector ─────────────────────────
  document.getElementById('hydroWindowSize')?.addEventListener('change', () => {
    if (!currentSequence) return;
    const winEl = document.getElementById('hydroWindowSize');
    const winSize = winEl ? parseInt(winEl.value, 10) : 7;
    renderHydroPlot(currentSequence, winSize);
  });

  // ── Theme Sync ──────────────────────────────────────────
  window.addEventListener('biokit-theme-change', () => {
    const isVisible = !document.getElementById('panel-protein').classList.contains('hidden') &&
                      !document.querySelector('#proteinResults .results-content').classList.contains('hidden');
    if (isVisible && currentSequence) {
      const winEl = document.getElementById('hydroWindowSize');
      const winSize = winEl ? parseInt(winEl.value, 10) : 7;
      renderAAChart(currentSequence);
      renderHydroPlot(currentSequence, winSize);
      renderTitrationCurve(currentSequence);
    }
  });

})();
