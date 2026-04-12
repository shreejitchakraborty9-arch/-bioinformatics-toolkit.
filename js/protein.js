/* ============================================================
   protein.js — Protein Analyzer
   Phase 1 (crash fix) + Phase 2 (feature wiring) + Phase 3 (pI)
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };
  const BioMath = window.BioKit.core.BioMath;

  let currentSequence = "";

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
  function run() {
    const input = document.getElementById('proteinInput');
    if (!input) return;

    const raw = input.value;
    const validation = window.validateSequence(raw, 'protein');

    if (!validation.valid) {
      window.showValidationWarning('panel-protein', validation);
      return;
    }

    const seq = validation.clean.replace(/\*/g, '');
    currentSequence = seq;

    window.withLoading('panel-protein', () => {
      const hydro = BioMath.calculateHydrophobicity(seq);
      const mw = BioMath.calculateProtein_MW(seq);
      const pI = BioMath.calculatePI(seq);
      const extMatch = BioMath.calculateExtinctionCoefficient(seq, 'protein'); 
      
      const aliphatic = BioMath.calculateAliphaticIndex(seq);
      const instability = BioMath.calculateInstabilityIndex(seq);
      
      const instabLabel = instability < 40 ? ' (Stable)' : ' (Unstable)';

      // Stat cards
      window.buildStatCards('proteinStatRow', [
        { val: seq.length.toLocaleString(), label: 'Length (aa)' },
        { val: mw + ' kDa', label: 'Mol. Weight' },
        { val: pI, label: 'Isoelectric Pt' },
        { val: extMatch.toLocaleString(), label: 'Ext. Coeff (M⁻¹cm⁻¹)' },
        { val: hydro, label: 'GRAVY Index' },
        { val: instability + instabLabel, label: 'Instability Index' }
      ]);

      // AA Composition chart
      renderAAChart(seq);

      // AA Table
      renderAATable(seq);

      // Hydrophobicity plot
      renderHydroPlot(seq);

      // Show results, hide empty state
      const empty = document.getElementById('proteinEmptyState');
      const content = document.querySelector('#proteinResults .results-content');
      if (empty) empty.classList.add('hidden');
      if (content) content.classList.remove('hidden');

      window.showToast('Analysis Complete');
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
  function renderHydroPlot(seq) {
    const canvas = document.getElementById('hydroCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');

    canvas.width = canvas.offsetWidth || 700;
    canvas.height = 220;
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    const W = canvas.width, H = canvas.height;
    const pad = { top: 20, right: 20, bottom: 30, left: 40 };
    const plotW = W - pad.left - pad.right;
    const plotH = H - pad.top - pad.bottom;
    const windowSize = 7;

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

  // ── Events ─────────────────────────────────────────────
  document.getElementById('proteinAnalyzeBtn')?.addEventListener('click', window.debounce(run, 300));

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
    document.getElementById('proteinResults')?.classList.add('hidden');
  });

  // Live meta
  const pInput = document.getElementById('proteinInput');
  const pMeta = document.getElementById('proteinInputMeta');
  if (pInput && pMeta) {
    pInput.addEventListener('input', () => {
      const seq = window.cleanSeq(pInput.value);
      pMeta.textContent = `Length: ${seq.length} aa`;
    });
  }

  // ── Theme Sync ──────────────────────────────────────────
  window.addEventListener('biokit-theme-change', () => {
    const isVisible = !document.getElementById('panel-protein').classList.contains('hidden') && 
                      !document.querySelector('#proteinResults .results-content').classList.contains('hidden');
    if (isVisible && currentSequence) {
      renderAAChart(currentSequence);
      renderHydroPlot(currentSequence);
    }
  });

})();
