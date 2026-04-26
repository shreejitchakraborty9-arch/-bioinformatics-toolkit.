/* ============================================================
   genome-viewer.js — IGV-style Canvas Genome Viewer
   ============================================================ */

(function () {
  'use strict';

  // State
  let canvas, ctx;
  let sequence = "";
  let cdsData = [];
  let motifData = [];
  let storedPromoterLen = 1000;
  let themeColorCache = null;
  let rafPending = false;

  function getThemeColors() {
    if (themeColorCache !== null) return themeColorCache;
    const style = getComputedStyle(document.documentElement);
    themeColorCache = {
      bg:      style.getPropertyValue('--bg-base').trim()    || '#111827',
      surface: style.getPropertyValue('--bg-sidebar').trim() || '#1e293b',
      text:    style.getPropertyValue('--text-primary').trim()|| '#f8f9fa',
      muted:   style.getPropertyValue('--text-muted').trim() || '#94a3b8',
      border:  style.getPropertyValue('--border').trim()     || '#334155',
      grid:    style.getPropertyValue('--border-muted').trim()|| '#1e293b'
    };
    return themeColorCache;
  }

  // Viewport (0-based internal coordinates, matching the 0-based CDS/motif data fed by promoter.js)
  let bpStart = 0;
  let bpEnd   = 1000;
  let isDragging  = false;
  let dragStartX  = 0;
  let dragStartBp = 0;
  let dragMoved   = false; // true when the mousedown→mouseup movement exceeded the click threshold

  // Track Layout
  const TRACK_HEIGHT  = 40;
  const HEADER_HEIGHT = 30;
  const MARGIN_TOP    = 20;
  const TRACK_SPACING = 30;

  // Stable frame→color mapping so a given reading frame always gets the same color
  const FRAME_COLORS = {
    '+1': '#14b8a6', '+2': '#8b5cf6', '+3': '#ec4899',
    '-1': '#0f766e', '-2': '#6d28d9', '-3': '#be185d'
  };

  function frameColor(frame) {
    const key = (frame > 0 ? '+' : '') + frame;
    return FRAME_COLORS[key] || '#64748b';
  }

  const TRACKS = [
    { id: 'seq',      label: 'DNA Sequence',    y: 0 },
    { id: 'cds',      label: 'CDS / ORFs',      y: 0 },
    { id: 'promoter', label: 'Promoter Region', y: 0 },
    { id: 'motif',    label: 'Cis-Elements',    y: 0 }
  ];

  function initViewer() {
    canvas = document.getElementById('genomeCanvas');
    if (!canvas) return;

    ctx = canvas.getContext('2d', { alpha: false });

    window.addEventListener('resize', resizeCanvas);
    resizeCanvas();

    canvas.addEventListener('mousedown', onMouseDown);
    window.addEventListener('mousemove', onMouseMove);
    window.addEventListener('mouseup',   onMouseUp);
    canvas.addEventListener('click',     onClick);

    requestAnimationFrame(render);
  }

  function resizeCanvas() {
    const parent = canvas.parentElement;
    if (!parent) return;

    const rect = parent.getBoundingClientRect();
    const dpr  = window.devicePixelRatio || 1;

    canvas.width  = rect.width * dpr;
    canvas.height = 400 * dpr;

    canvas.style.width  = `${rect.width}px`;
    canvas.style.height = '400px';

    // Reset before scaling to prevent cumulative multiplication
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.scale(dpr, dpr);

    let currentY = MARGIN_TOP + HEADER_HEIGHT;
    for (const t of TRACKS) {
      t.y = currentY;
      currentY += TRACK_HEIGHT + TRACK_SPACING;
    }

    requestAnimationFrame(render);
  }

  // ── Data Loading ─────────────────────────────────────────────
  window.loadGenomeData = function (seq, cds, promoterLen, motifs) {
    sequence  = seq || '';
    cdsData   = cds || [];
    cdsData.sort((a, b) => a.start - b.start);
    motifData = motifs || [];
    motifData.sort((a, b) => a.position - b.position);
    storedPromoterLen = (typeof promoterLen === 'number' && promoterLen > 0) ? promoterLen : 1000;

    bpStart = 0;
    bpEnd   = Math.min(100, sequence.length || 100);

    resizeCanvas();
  };

  window.setZoom = function (level) {
    if (!sequence) return;
    const centerBp = bpStart + (bpEnd - bpStart) / 2;
    let windowSize = 100;

    if (level === '1bp')    windowSize = 50;
    if (level === '10bp')   windowSize = 250;
    if (level === '100bp')  windowSize = 2500;
    if (level === '1000bp') windowSize = 25000;

    windowSize = Math.min(windowSize, sequence.length);

    bpStart = Math.max(0, centerBp - windowSize / 2);
    bpEnd   = Math.min(sequence.length, centerBp + windowSize / 2);

    if (bpStart === 0)              bpEnd   = Math.min(sequence.length, windowSize);
    if (bpEnd   === sequence.length) bpStart = Math.max(0, sequence.length - windowSize);

    requestAnimationFrame(render);
  };

  // ── Interaction ──────────────────────────────────────────────
  function onMouseDown(e) {
    isDragging  = true;
    dragMoved   = false;
    dragStartX  = e.clientX;
    dragStartBp = bpStart;
    canvas.style.cursor = 'grabbing';
  }

  function onMouseMove(e) {
    if (!isDragging || !sequence) return;

    const dx = e.clientX - dragStartX;
    if (Math.abs(dx) > 3) dragMoved = true; // 3 px threshold: drag vs. click

    const cw      = canvas.clientWidth;
    const bpWidth = bpEnd - bpStart;
    const bpShift = (dx / cw) * bpWidth;

    bpStart = dragStartBp - bpShift;
    bpEnd   = bpStart + bpWidth;

    if (bpStart < 0) {
      bpStart = 0;
      bpEnd   = bpWidth;
    }
    if (bpEnd > sequence.length) {
      bpEnd   = sequence.length;
      bpStart = Math.max(0, sequence.length - bpWidth);
    }

    if (!rafPending) {
      rafPending = true;
      requestAnimationFrame(() => { rafPending = false; render(); });
    }
  }

  function onMouseUp() {
    isDragging = false;
    if (canvas) canvas.style.cursor = 'grab';
  }

  function onClick(e) {
    if (!sequence || dragMoved) return; // suppress tooltip after a drag

    const rect    = canvas.getBoundingClientRect();
    const x       = e.clientX - rect.left;
    const y       = e.clientY - rect.top;
    const cw      = rect.width;
    const bpWidth = bpEnd - bpStart;
    const clickBp = bpStart + (x / cw) * bpWidth;

    // CDS track
    const cdsTrack = TRACKS.find(t => t.id === 'cds');
    if (y >= cdsTrack.y && y <= cdsTrack.y + TRACK_HEIGHT) {
      for (const cds of cdsData) {
        if (clickBp >= cds.start && clickBp <= cds.end) {
          const strand       = cds.frame > 0 ? '+' : '−';
          const displayStart = cds.start + 1; // convert to 1-based for display
          showCanvasTooltip(
            `CDS Frame ${strand}${Math.abs(cds.frame)}\nPos: ${displayStart}–${cds.end}\nLen: ${cds.length} bp`,
            e.clientX, e.clientY
          );
          return;
        }
      }
    }

    // Motif track
    const motifTrack = TRACKS.find(t => t.id === 'motif');
    if (y >= motifTrack.y && y <= motifTrack.y + TRACK_HEIGHT) {
      for (const m of motifData) {
        const mMatchLen = m.matched?.length || (m.end != null ? m.end - m.position : 6);
        const mEnd      = m.position + mMatchLen;
        if (clickBp >= m.position && clickBp <= mEnd) {
          const displayPos = m.position + 1; // 1-based for display
          showCanvasTooltip(
            `Motif: ${m.name}\nPos: ${displayPos}\nSeq: ${m.matched || '—'}`,
            e.clientX, e.clientY
          );
          // Snap viewport to centre on the motif
          const center = m.position + mMatchLen / 2;
          const w      = bpEnd - bpStart;
          bpStart = Math.max(0, center - w / 2);
          bpEnd   = Math.min(sequence.length, center + w / 2);
          requestAnimationFrame(render);
          return;
        }
      }
    }

    hideCanvasTooltip();
  }

  // ── Tooltip ──────────────────────────────────────────────────
  function showCanvasTooltip(text, x, y) {
    let tt = document.getElementById('gvTooltip');
    if (!tt) {
      tt = document.createElement('div');
      tt.id = 'gvTooltip';
      // Use CSS variables that actually exist in style.css
      tt.style.cssText = [
        'position:fixed',
        'background:var(--bg-surface,#1e293b)',
        'color:var(--text-primary,#f8f9fa)',
        'border:1px solid var(--border,#334155)',
        'padding:8px 12px',
        'border-radius:6px',
        'box-shadow:0 10px 15px -3px rgba(0,0,0,0.5)',
        'pointer-events:none',
        'z-index:1000',
        'font-family:var(--font,"Inter",system-ui,sans-serif)',
        'font-size:0.85rem',
        'white-space:pre-line',
        'max-width:260px'
      ].join(';');
      document.body.appendChild(tt);
    }
    tt.textContent  = text;
    tt.style.left   = (x + 15) + 'px';
    tt.style.top    = (y + 15) + 'px';
    tt.style.display = 'block';
  }

  function hideCanvasTooltip() {
    const tt = document.getElementById('gvTooltip');
    if (tt) tt.style.display = 'none';
  }

  // ── Rendering Engine ─────────────────────────────────────────
  function render() {
    if (!ctx || !canvas) return;

    const cw     = canvas.clientWidth;
    const ch     = canvas.clientHeight;
    const colors = getThemeColors();

    ctx.clearRect(0, 0, cw, ch);
    ctx.fillStyle = colors.bg;
    ctx.fillRect(0, 0, cw, ch);

    if (!sequence) {
      ctx.fillStyle    = colors.muted;
      ctx.font         = '14px "Inter", system-ui, sans-serif';
      ctx.textAlign    = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('No sequence loaded. Analyze a sequence to view.', cw / 2, ch / 2);
      return;
    }

    drawRuler(cw);
    drawSequenceTrack(cw);
    drawCDSTrack(cw);
    drawPromoterTrack(cw);
    drawMotifTrack(cw);
  }

  function xForBp(bp, cw) {
    return ((bp - bpStart) / (bpEnd - bpStart)) * cw;
  }

  // Binary search: first index where arr[i][startKey] >= minBp.
  // Steps back one to catch features whose start is before minBp but end extends into view.
  function findFirstVisibleIndex(arr, startKey, minBp) {
    let low = 0, high = arr.length - 1, res = arr.length;
    while (low <= high) {
      const mid = (low + high) >> 1;
      if (arr[mid][startKey] >= minBp) { res = mid; high = mid - 1; }
      else                               low = mid + 1;
    }
    return Math.max(0, res - 1);
  }

  function drawRuler(cw) {
    const colors = getThemeColors();
    ctx.fillStyle = colors.surface;
    ctx.fillRect(0, 0, cw, HEADER_HEIGHT);

    ctx.fillStyle    = colors.muted;
    ctx.font         = '11px "JetBrains Mono","Courier New",monospace';
    ctx.textAlign    = 'center';
    ctx.textBaseline = 'middle';

    const bpWidth = bpEnd - bpStart;

    let tickInterval = 10;
    if (bpWidth >  200)  tickInterval = 50;
    if (bpWidth > 1000)  tickInterval = 200;
    if (bpWidth > 5000)  tickInterval = 1000;
    if (bpWidth > 20000) tickInterval = 5000;

    const firstTick = Math.ceil(bpStart / tickInterval) * tickInterval;

    ctx.strokeStyle = colors.border;
    ctx.lineWidth   = 1;
    ctx.beginPath();
    for (let bp = firstTick; bp <= bpEnd; bp += tickInterval) {
      const x = xForBp(bp, cw);
      ctx.moveTo(x, HEADER_HEIGHT - 6);
      ctx.lineTo(x, HEADER_HEIGHT);
    }
    ctx.stroke();

    // Display 1-based positions (internal 0-based + 1) for biological convention
    for (let bp = firstTick; bp <= bpEnd; bp += tickInterval) {
      ctx.fillText((bp + 1).toLocaleString(), xForBp(bp, cw), HEADER_HEIGHT / 2);
    }

    ctx.beginPath();
    ctx.moveTo(0,  HEADER_HEIGHT);
    ctx.lineTo(cw, HEADER_HEIGHT);
    ctx.strokeStyle = colors.border;
    ctx.stroke();
  }

  function drawTrackBackground(t, cw) {
    const colors = getThemeColors();
    ctx.fillStyle    = colors.muted;
    ctx.font         = '12px "Inter", system-ui, sans-serif';
    ctx.textAlign    = 'left';
    ctx.textBaseline = 'bottom';
    ctx.fillText(t.label, 10, t.y - 4);

    ctx.lineWidth   = 1;
    ctx.strokeStyle = colors.grid;
    ctx.beginPath();
    ctx.moveTo(0,  t.y + TRACK_HEIGHT);
    ctx.lineTo(cw, t.y + TRACK_HEIGHT);
    ctx.stroke();
  }

  function drawSequenceTrack(cw) {
    const t       = TRACKS.find(t => t.id === 'seq');
    drawTrackBackground(t, cw);

    const bpWidth = bpEnd - bpStart;

    if (bpWidth < 150) {
      ctx.font         = '14px "JetBrains Mono","Courier New",monospace';
      ctx.textAlign    = 'center';
      ctx.textBaseline = 'middle';

      const startIdx = Math.max(0, Math.floor(bpStart));
      const endIdx   = Math.min(sequence.length, Math.ceil(bpEnd));

      for (let i = startIdx; i < endIdx; i++) {
        const char = sequence[i].toUpperCase();
        const x    = xForBp(i + 0.5, cw);

        if      (char === 'A')              ctx.fillStyle = '#22c55e';
        else if (char === 'T' || char === 'U') ctx.fillStyle = '#ef4444';
        else if (char === 'G')              ctx.fillStyle = '#f97316';
        else if (char === 'C')              ctx.fillStyle = '#3b82f6';
        else                               ctx.fillStyle = '#94a3b8';

        ctx.fillText(char, x, t.y + TRACK_HEIGHT / 2);
      }
    } else {
      // Zoomed-out: show a labelled placeholder until GC-density histogram is approved
      const colors = getThemeColors();
      ctx.fillStyle = colors.border;
      ctx.fillRect(0, t.y + 10, cw, TRACK_HEIGHT - 20);

      ctx.fillStyle    = colors.muted;
      ctx.font         = '11px "Inter", system-ui, sans-serif';
      ctx.textAlign    = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('Zoom in below 150 bp to see individual bases', cw / 2, t.y + TRACK_HEIGHT / 2);
    }
  }

  function drawCDSTrack(cw) {
    const t       = TRACKS.find(t => t.id === 'cds');
    drawTrackBackground(t, cw);

    const bpWidth  = bpEnd - bpStart;
    const firstIdx = findFirstVisibleIndex(cdsData, 'start', bpStart);

    for (let i = firstIdx; i < cdsData.length; i++) {
      const cds = cdsData[i];
      if (cds.end   < bpStart) continue; // entirely before viewport
      if (cds.start > bpEnd)   break;    // entirely after viewport

      const x1        = Math.max(0, xForBp(cds.start, cw));
      const x2raw     = xForBp(cds.end, cw);
      const x2        = Math.min(cw, x2raw);
      const w         = Math.max(2, x2 - x1);
      const isReverse = cds.frame < 0;
      const arrowSize = Math.min(10, w * 0.25);
      const y         = t.y + 10;
      const h         = 20;

      ctx.fillStyle = frameColor(cds.frame);
      ctx.beginPath();

      if (!isReverse) {
        // Forward strand: arrowhead on right (→)
        ctx.moveTo(x1, y);
        if (cds.end <= bpEnd && w > arrowSize) {
          ctx.lineTo(x1 + w - arrowSize, y);
          ctx.lineTo(x1 + w,             y + h / 2);
          ctx.lineTo(x1 + w - arrowSize, y + h);
        } else {
          ctx.lineTo(x1 + w, y);
          ctx.lineTo(x1 + w, y + h);
        }
        ctx.lineTo(x1, y + h);
      } else {
        // Reverse strand: arrowhead on left (←)
        ctx.moveTo(x1 + w, y);
        if (cds.start >= bpStart && w > arrowSize) {
          ctx.lineTo(x1 + arrowSize, y);
          ctx.lineTo(x1,             y + h / 2);
          ctx.lineTo(x1 + arrowSize, y + h);
        } else {
          ctx.lineTo(x1, y);
          ctx.lineTo(x1, y + h);
        }
        ctx.lineTo(x1 + w, y + h);
      }

      ctx.closePath();
      ctx.fill();

      if (w > 40 && bpWidth < 5000) {
        ctx.fillStyle    = '#ffffff';
        ctx.font         = '10px "Inter", system-ui, sans-serif';
        ctx.textAlign    = isReverse ? 'right' : 'left';
        ctx.textBaseline = 'middle';
        const strand     = cds.frame > 0 ? '+' : '−';
        const labelX     = isReverse ? x1 + w - 5 : x1 + 5;
        ctx.fillText(`ORF ${strand}${Math.abs(cds.frame)}`, labelX, y + h / 2);
      }
    }
  }

  function drawPromoterTrack(cw) {
    const t = TRACKS.find(t => t.id === 'promoter');
    drawTrackBackground(t, cw);

    if (cdsData.length === 0) return;

    // Anchor to the first forward-strand ORF; fall back to the first ORF overall
    const fwdOrfs = cdsData.filter(o => o.frame > 0);
    const anchor  = fwdOrfs.length > 0 ? fwdOrfs[0] : cdsData[0];

    const pStart = Math.max(0, anchor.start - storedPromoterLen);
    const pEnd   = anchor.start;

    if (pEnd > bpStart && pStart < bpEnd) {
      const x1 = Math.max(0, xForBp(pStart, cw));
      const x2 = Math.min(cw, xForBp(pEnd,   cw));
      const w  = x2 - x1;

      ctx.fillStyle = 'rgba(245, 158, 11, 0.2)';
      ctx.fillRect(x1, t.y + 10, w, 20);

      ctx.strokeStyle = '#f59e0b';
      ctx.lineWidth   = 1;
      ctx.strokeRect(x1, t.y + 10, w, 20);
    }
  }

  function drawMotifTrack(cw) {
    const t        = TRACKS.find(t => t.id === 'motif');
    drawTrackBackground(t, cw);

    const firstIdx = findFirstVisibleIndex(motifData, 'position', bpStart);

    for (let i = firstIdx; i < motifData.length; i++) {
      const m        = motifData[i];
      const motifLen = m.matched?.length || (m.end != null ? m.end - m.position : 6);
      const mEnd     = m.position + motifLen;

      if (mEnd   < bpStart) continue;
      if (m.position > bpEnd) break;

      const x1 = Math.max(0, xForBp(m.position, cw));
      const x2 = Math.min(cw, xForBp(mEnd, cw));
      const w  = Math.max(2, x2 - x1);

      // Color matches the MOTIF_META legend in promoter.js
      if      (m.name.includes('TATA')) ctx.fillStyle = '#ef4444';
      else if (m.name.includes('CAAT')) ctx.fillStyle = '#3b82f6';
      else if (m.name.includes('GC'))   ctx.fillStyle = '#8b5cf6';
      else if (m.name.includes('-35'))  ctx.fillStyle = '#22c55e';
      else if (m.name.includes('-10'))  ctx.fillStyle = '#f97316';
      else                              ctx.fillStyle = '#22d3ee';

      ctx.fillRect(x1, t.y + 5, w, TRACK_HEIGHT - 10);

      if (w > 30) {
        ctx.fillStyle    = '#ffffff';
        ctx.font         = '9px "JetBrains Mono","Courier New",monospace';
        ctx.textAlign    = 'left';
        ctx.textBaseline = 'middle';
        ctx.fillText(m.name, x1 + 3, t.y + TRACK_HEIGHT / 2);
      }
    }
  }

  // pUC19 lacZ-alpha region: lac promoter, Shine-Dalgarno, ATG start codon (GenBank L09137)
  const DEMO_SEQUENCE =
    'GCGGATAACAATTTCACACAGGAAACAGCTATGACCATGATTACGCCAAGCTTGCATGCCTGCAG' +
    'GTCGACGGATCCCCGGGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCA' +
    'TGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATT' +
    'CCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACT' +
    'CACATTAATTGCGTTACGCTGCGGTTTTCATGAGAATGTTTTTCTTTTCATGAGAAAAGCCCGGC' +
    'TTATAGTTTGCTTTTTATTTGATTTGAGTAATTTTGTTTTTTATACTATTTTTTTTTGAGCTTTT' +
    'GTTCGTTCAGAGTTTATTCGCTTCATTTAAATGGTATGAAATTTACTGATAATGATATTTTTAT';

  document.addEventListener('DOMContentLoaded', () => {
    initViewer();

    const btn = document.getElementById('gvLoadDemoBtn');
    if (btn) {
      btn.addEventListener('click', () => {
        const input      = document.getElementById('proInput');
        const analyzeBtn = document.getElementById('proAnalyzeBtn');
        if (!input || !analyzeBtn) return;

        input.value = DEMO_SEQUENCE;
        input.dispatchEvent(new Event('input'));
        analyzeBtn.click();

        const panel = document.getElementById('panel-genome-viewer');
        if (panel && !panel.querySelector('.demo-banner')) {
          const banner = document.createElement('div');
          banner.className = 'demo-banner';
          banner.style.cssText =
            'display:flex;align-items:center;justify-content:space-between;gap:8px;' +
            'padding:8px 12px;margin-bottom:10px;' +
            'background:var(--warning-bg,#fef3c7);border:1px solid var(--warning-border,#f59e0b);' +
            'border-radius:6px;font-size:0.85rem;color:var(--warning-text,#92400e);';
          banner.innerHTML =
            '<span>⚠ Viewing Educational Demo. Paste your own sequence to calculate custom coordinates.</span>' +
            '<button style="background:none;border:none;cursor:pointer;font-size:1rem;line-height:1;' +
            'padding:0 2px;color:inherit;" aria-label="Dismiss">&times;</button>';
          banner.querySelector('button').addEventListener('click', () => banner.remove());
          panel.insertBefore(banner, panel.firstChild);
        }
      });
    }
  });

  // ── Theme Sync ───────────────────────────────────────────────
  window.addEventListener('biokit-theme-change', () => {
    themeColorCache = null;
    const panel     = document.getElementById('panel-genome-viewer');
    const isVisible = panel && !panel.classList.contains('hidden');
    if (isVisible) requestAnimationFrame(render);
  });

})();
