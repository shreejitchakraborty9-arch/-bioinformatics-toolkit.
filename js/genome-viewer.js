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

    ctx = canvas.getContext('2d', { alpha: true });

    const canvasContainer = canvas.closest('.gv-canvas-container') || canvas.parentElement;
    if (canvasContainer) {
      new ResizeObserver(() => {
        if (canvasContainer.clientWidth > 0) {
          resizeCanvas();
          render();
        }
      }).observe(canvasContainer);
    }

    window.addEventListener('resize', resizeCanvas);
    resizeCanvas();

    canvas.addEventListener('mousedown', onMouseDown);
    window.addEventListener('mousemove', onMouseMove);
    window.addEventListener('mouseup',   onMouseUp);
    canvas.addEventListener('click',     onClick);

    const zoomInBtn  = document.getElementById('gvZoomIn');
    const zoomOutBtn = document.getElementById('gvZoomOut');

    if (zoomInBtn) {
      zoomInBtn.addEventListener('click', () => {
        windowSize = Math.max(50, windowSize / 1.5);
        updateZoom();
      });
    }
    if (zoomOutBtn) {
      zoomOutBtn.addEventListener('click', () => {
        windowSize = Math.min(sequence.length || 20000, windowSize * 1.5);
        updateZoom();
      });
    }

    canvas.addEventListener('wheel', (e) => {
      e.preventDefault();
      if (e.deltaY < 0) {
        windowSize = Math.max(50, windowSize / 1.5);
      } else {
        windowSize = Math.min(sequence.length || 20000, windowSize * 1.5);
      }
      updateZoom();
    }, { passive: false });

    // Keyboard navigation: arrow keys to pan, +/- to zoom, Home/End to jump (UX-01)
    document.addEventListener('keydown', (e) => {
      const panel = document.getElementById('panel-genome-viewer');
      if (!panel || panel.classList.contains('hidden')) return;
      if (!sequence) return;

      // Ignore when focus is inside an input/textarea/select
      const tag = document.activeElement?.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;

      const bpWidth   = bpEnd - bpStart;
      const panAmount = Math.max(1, bpWidth * 0.1);

      if (e.key === 'ArrowLeft') {
        e.preventDefault();
        bpStart = Math.max(0, bpStart - panAmount);
        bpEnd   = bpStart + bpWidth;
        if (bpEnd > sequence.length) { bpEnd = sequence.length; bpStart = Math.max(0, bpEnd - bpWidth); }
        requestAnimationFrame(render);
      } else if (e.key === 'ArrowRight') {
        e.preventDefault();
        bpEnd   = Math.min(sequence.length, bpEnd + panAmount);
        bpStart = bpEnd - bpWidth;
        if (bpStart < 0) { bpStart = 0; bpEnd = Math.min(sequence.length, bpWidth); }
        requestAnimationFrame(render);
      } else if (e.key === '+' || e.key === '=') {
        windowSize = Math.max(50, windowSize / 1.5);
        updateZoom();
      } else if (e.key === '-') {
        windowSize = Math.min(sequence.length || 20000, windowSize * 1.5);
        updateZoom();
      } else if (e.key === 'Home') {
        e.preventDefault();
        bpStart = 0;
        bpEnd   = bpWidth;
        requestAnimationFrame(render);
      } else if (e.key === 'End') {
        e.preventDefault();
        bpEnd   = sequence.length;
        bpStart = Math.max(0, sequence.length - bpWidth);
        requestAnimationFrame(render);
      }
    });

    // Fix 5: re-render when the panel is re-shown after navigation away
    const gvPanel = document.getElementById('panel-genome-viewer');
    if (gvPanel) {
      gvPanel.addEventListener('panel-shown', () => {
        resizeCanvas();
        requestAnimationFrame(render);
      });
    }

    // Fix 7: wire custom sequence accordion
    const seqTextarea = document.getElementById('gvSeqTextarea');
    const seqCounter  = document.getElementById('gvSeqCounter');
    const analyzeBtn  = document.getElementById('gvAnalyzeBtn');
    const accordion   = document.getElementById('gvCustomSeqAccordion');

    if (seqTextarea && seqCounter) {
      seqTextarea.addEventListener('input', () => {
        const clean = seqTextarea.value.replace(/^>.*$/mg, '').replace(/\s+/g, '');
        const count = clean.length;
        const warn  = count > 45000;
        seqCounter.textContent = `${count.toLocaleString()} bp${warn ? ' — approaching 50,000 bp limit' : ''}`;
        seqCounter.style.color = warn ? '#f59e0b' : 'var(--text-muted)';
      });
    }

    if (analyzeBtn) {
      analyzeBtn.addEventListener('click', async () => {
        if (!seqTextarea) return;
        const raw = seqTextarea.value.trim();
        if (!raw) { alert('Please paste a sequence first.'); return; }

        const clean = raw.replace(/^>.*$/mg, '').replace(/\s+/g, '').toUpperCase();
        if (clean.length === 0) { alert('No sequence bases found.'); return; }
        if (clean.length > 50000) {
          alert(`Sequence is ${clean.length.toLocaleString()} bp. Maximum allowed is 50,000 bp.`);
          return;
        }
        const invalid = [...new Set(clean.split('').filter(c => !/^[ACGTN]$/.test(c)))];
        if (invalid.length > 0) {
          alert(`Invalid characters found: ${invalid.join(', ')}\nOnly A, C, G, T, N are accepted.`);
          return;
        }

        const organism = document.getElementById('gvOrganismSelect')?.value || 'prokaryote';
        analyzeBtn.disabled    = true;
        analyzeBtn.textContent = 'Analyzing…';

        try {
          const resp = await fetch('/api/analyze/promoter', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sequence: clean, organism_type: organism })
          });
          if (!resp.ok) {
            const err = await resp.json().catch(() => ({}));
            throw new Error(err.error || `Server error ${resp.status}`);
          }
          const data = await resp.json();

          const cds = (data.orfs || []).map(o => ({
            start:   o.start - 1,
            end:     o.end,
            frame:   o.frame,
            length:  o.length_nt || o.length,
            name:    o.name || `ORF frame ${o.frame}`,
            partial: o.partial || false
          }));
          const motifs = (data.motifs || []).map(m => ({
            position: m.position - 1,
            end:      m.end != null ? m.end - 1 : undefined,
            matched:  m.matched,
            name:     m.name
          }));
          const promoterLen = data.promoter_length || data.promoterLen || 1000;

          window.loadGenomeData(clean, cds, promoterLen, motifs);
          window.setZoom('10bp');
          if (accordion) accordion.open = false;

        } catch (err) {
          alert(`Analysis failed: ${err.message}`);
        } finally {
          analyzeBtn.disabled    = false;
          analyzeBtn.textContent = 'Analyze & Visualize';
        }
      });
    }

    requestAnimationFrame(render);
  }

  function resizeCanvas() {
    if (!canvas) return;
    const parent = canvas.parentElement;
    if (!parent) return;

    // Obliterate any rogue CSS constraints before measuring — inline styles
    // override non-!important stylesheet rules regardless of specificity.
    canvas.style.maxWidth = 'none';
    canvas.style.minWidth = '100%';
    canvas.style.width    = '100%';
    canvas.style.display  = 'block';

    const rect = parent.getBoundingClientRect();
    const dpr  = window.devicePixelRatio || 1;

    // Compute canvas height from track layout so adding a 5th track
    // requires no magic-number changes here (CQ-03).
    const canvasHeight = MARGIN_TOP + HEADER_HEIGHT + TRACKS.length * (TRACK_HEIGHT + TRACK_SPACING);

    canvas.width  = rect.width * dpr;
    canvas.height = canvasHeight * dpr;

    canvas.style.height = canvasHeight + 'px';

    // Reset before scaling to prevent cumulative DPR multiplication.
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.scale(dpr, dpr);

    // Fill immediately so no black frame reaches the browser before RAF.
    ctx.fillStyle = getThemeColors().bg;
    ctx.fillRect(0, 0, rect.width, canvasHeight);

    let currentY = MARGIN_TOP + HEADER_HEIGHT;
    for (const t of TRACKS) {
      t.y = currentY;
      currentY += TRACK_HEIGHT + TRACK_SPACING;
    }

    requestAnimationFrame(render);
  }

  // ── Data Loading ─────────────────────────────────────────────
  //
  // Coordinate contract (canonical, 0-based half-open):
  //   CDS:   { start: 0-based inclusive, end: 0-based exclusive, frame, length, name }
  //   Motif: { position: 0-based inclusive, end: 0-based exclusive, matched, name }
  //
  // Both paths that call this function MUST convert API 1-based coordinates
  // before calling (Path A: promoter.js; Path B: analyzeBtn handler above).
  window.loadGenomeData = function (seq, cds, promoterLen, motifs) {
    sequence  = seq || '';
    cdsData   = cds || [];
    cdsData.sort((a, b) => a.start - b.start);
    motifData = motifs || [];
    motifData.sort((a, b) => a.position - b.position);
    storedPromoterLen = (typeof promoterLen === 'number' && promoterLen > 0) ? promoterLen : 1000;

    bpStart    = 0;
    windowSize = Math.min(100, sequence.length || 100);
    bpEnd      = windowSize;

    const label = document.getElementById('gvZoomLabel');
    if (label) label.textContent = Math.round(windowSize) + ' bp';

    resizeCanvas();
  };

  // Current window size in bp; kept in module scope so buttons and wheel share state
  let windowSize = 1000;

  function updateZoom() {
    const label = document.getElementById('gvZoomLabel');
    if (label) label.textContent = Math.round(windowSize) + ' bp';

    const centerBp = bpStart + (bpEnd - bpStart) / 2;
    bpStart = Math.max(0, centerBp - windowSize / 2);
    bpEnd   = Math.min(sequence.length || windowSize, centerBp + windowSize / 2);

    if (bpStart === 0)                       bpEnd   = Math.min(sequence.length || windowSize, windowSize);
    if (bpEnd   === (sequence.length || 0))  bpStart = Math.max(0, (sequence.length || 0) - windowSize);

    requestAnimationFrame(render);
  }

  window.setZoom = function (level) {
    if (!sequence) return;
    if (level === '1bp')    windowSize = 50;
    else if (level === '10bp')   windowSize = 250;
    else if (level === '100bp')  windowSize = 2500;
    else if (level === '1000bp') windowSize = 25000;
    windowSize = Math.min(windowSize, sequence.length);
    updateZoom();
  };

  // Namespaced API — preferred over raw window.* globals (CQ-02)
  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };
  window.BioKit.genomeViewer = {
    loadData: window.loadGenomeData,
    setZoom:  window.setZoom
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
    if (!sequence) return;

    // ── Drag panning ────────────────────────────────────────────
    if (isDragging) {
      const dx = e.clientX - dragStartX;
      if (Math.abs(dx) > 3) dragMoved = true;

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
      hideCanvasTooltip();
      return;
    }

    // ── Hover hit-detection ──────────────────────────────────────
    const rect    = canvas.getBoundingClientRect();
    const x       = e.clientX - rect.left;
    const y       = e.clientY - rect.top;
    const cw      = rect.width;
    const bpWidth = bpEnd - bpStart;
    const hoverBp = bpStart + (x / cw) * bpWidth;

    // CDS / ORF track
    const cdsTrack = TRACKS.find(t => t.id === 'cds');
    if (y >= cdsTrack.y && y <= cdsTrack.y + TRACK_HEIGHT) {
      for (const cds of cdsData) {
        if (hoverBp >= cds.start && hoverBp <= cds.end) {
          const strand = cds.frame > 0 ? '+' : '−';
          showTooltip(
            `<strong>${cds.name || 'ORF'}</strong><br>` +
            `Frame: ${strand}${Math.abs(cds.frame)}<br>` +
            `Start: ${cds.start + 1} | End: ${cds.end}<br>` +
            `Length: ${cds.length ?? '—'} bp` +
            (cds.partial ? '<br><em style="color:#f59e0b">partial (no stop codon)</em>' : ''),
            e.clientX, e.clientY,
            { isHTML: true }
          );
          return;
        }
      }
    }

    // Motif / Cis-element track
    const motifTrack = TRACKS.find(t => t.id === 'motif');
    if (y >= motifTrack.y && y <= motifTrack.y + TRACK_HEIGHT) {
      for (const m of motifData) {
        const mMatchLen = m.matched?.length || (m.end != null ? m.end - m.position : 6);
        const mEnd      = m.position + mMatchLen;
        if (hoverBp >= m.position && hoverBp <= mEnd) {
          showTooltip(
            `<strong>${m.name || 'Motif'}</strong><br>` +
            `Start: ${m.position + 1}<br>` +
            `Length: ${mMatchLen} bp<br>` +
            `Seq: <code>${m.matched || '—'}</code>`,
            e.clientX, e.clientY,
            { isHTML: true }
          );
          return;
        }
      }
    }

    hideCanvasTooltip();
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
          showTooltip(
            `CDS Frame ${strand}${Math.abs(cds.frame)}\nPos: ${displayStart}–${cds.end}\nLen: ${cds.length ?? '—'} bp` +
            (cds.partial ? '\n(partial)' : ''),
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
          showTooltip(
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
  // Single entry-point for all tooltip content.
  // isHTML=true: content is sanitized via DOMPurify before insertion (SEC-01).
  // Viewport clamping prevents clips at right/bottom screen edges (UX-02).
  function showTooltip(content, x, y, { isHTML = false } = {}) {
    let tt = document.getElementById('gvTooltip');
    if (!tt) {
      tt = document.createElement('div');
      tt.id = 'gvTooltip';
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
        'line-height:1.55',
        'max-width:260px'
      ].join(';');
      document.body.appendChild(tt);
    }

    if (isHTML) {
      tt.innerHTML = (window.DOMPurify ? DOMPurify.sanitize(content) : content);
    } else {
      tt.textContent = content;
    }

    // Clamp to viewport so the tooltip never clips off-screen (UX-02)
    const TT_W = 280;
    const TT_H = 110;
    const left = Math.min(x + 15, window.innerWidth  - TT_W);
    const top  = Math.min(y + 15, window.innerHeight - TT_H);
    tt.style.left    = Math.max(0, left) + 'px';
    tt.style.top     = Math.max(0, top)  + 'px';
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

    ctx.clearRect(0, 0, canvas.width, canvas.height);

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

    // Ticks at 0-based positions whose 1-based equivalent (bp+1) is a multiple
    // of tickInterval — so the ruler reads "10, 20, 30..." not "11, 21, 31..."
    const firstTick1based = Math.ceil((bpStart + 1) / tickInterval) * tickInterval;
    const firstTick       = firstTick1based - 1; // back to 0-based

    ctx.strokeStyle = colors.border;
    ctx.lineWidth   = 1;
    ctx.beginPath();
    for (let bp = firstTick; bp <= bpEnd; bp += tickInterval) {
      const x = xForBp(bp, cw);
      ctx.moveTo(x, HEADER_HEIGHT - 6);
      ctx.lineTo(x, HEADER_HEIGHT);
    }
    ctx.stroke();

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
    ctx.save();
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
    ctx.restore();
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
      // GC-density histogram (IGV coverage track style)
      const colors   = getThemeColors();
      const trackY   = t.y;
      const barAreaH = TRACK_HEIGHT - 6;

      // 50% threshold guide line
      const threshY = trackY + 3 + barAreaH * 0.5;
      ctx.save();
      ctx.strokeStyle = colors.border;
      ctx.lineWidth   = 1;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(0,  threshY);
      ctx.lineTo(cw, threshY);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.restore();

      ctx.save();
      ctx.fillStyle    = colors.muted;
      ctx.font         = '9px "Inter", system-ui, sans-serif';
      ctx.textAlign    = 'right';
      ctx.textBaseline = 'middle';
      ctx.fillText('50%', cw - 2, threshY);
      ctx.restore();

      // Bin size: target ~200 bins across the viewport, minimum 1 bp
      const bpWidth  = bpEnd - bpStart;
      const BIN_SIZE = Math.max(1, Math.ceil(bpWidth / 200));
      const seqS     = Math.max(0, Math.floor(bpStart));
      const seqE     = Math.min(sequence.length, Math.ceil(bpEnd));

      for (let binS = seqS; binS < seqE; binS += BIN_SIZE) {
        const binE = Math.min(seqE, binS + BIN_SIZE);
        let gc = 0;
        for (let i = binS; i < binE; i++) {
          const b = sequence[i].toUpperCase();
          if (b === 'G' || b === 'C') gc++;
        }
        const gcPct = (binE > binS) ? gc / (binE - binS) : 0;

        const x1 = xForBp(binS, cw);
        const x2 = xForBp(binE, cw);
        const w  = Math.max(1, x2 - x1);

        const barH = gcPct * barAreaH;
        const barY = trackY + 3 + (barAreaH - barH);

        // Linear interpolation blue (#3b82f6) → teal/green (#14b8a6) by GC fraction
        const r = Math.round(59  + (20  - 59)  * gcPct);
        const g = Math.round(130 + (184 - 130) * gcPct);
        const b = Math.round(246 + (166 - 246) * gcPct);
        ctx.fillStyle = `rgb(${r},${g},${b})`;
        ctx.fillRect(x1, barY, w, barH);
      }
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

      // Partial (run-off) ORFs: dashed stroke overlay to signal no stop codon found
      if (cds.partial) {
        ctx.save();
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth   = 1.5;
        ctx.setLineDash([4, 3]);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.restore();
      }

      if (w > 40 && bpWidth < 5000) {
        ctx.save();
        ctx.fillStyle    = '#ffffff';
        ctx.font         = '10px "Inter", system-ui, sans-serif';
        ctx.textAlign    = isReverse ? 'right' : 'left';
        ctx.textBaseline = 'middle';
        const strand     = cds.frame > 0 ? '+' : '−';
        const labelX     = isReverse ? x1 + w - 5 : x1 + 5;
        ctx.fillText(`ORF ${strand}${Math.abs(cds.frame)}`, labelX, y + h / 2);
        ctx.restore();
      }
    }
  }

  function drawPromoterTrack(cw) {
    const t = TRACKS.find(t => t.id === 'promoter');
    drawTrackBackground(t, cw);

    if (cdsData.length === 0) return;

    // Draw one upstream promoter window per ORF on each strand (SCI-01).
    // Forward ORFs: upstream region is immediately 5' of the ATG (lower coordinates).
    // Reverse ORFs: upstream region is immediately 5' of the ATG in reverse direction
    //               (higher coordinates, upstream on the minus strand).
    const seqLen = sequence.length;
    const drawn  = new Set();

    for (const orf of cdsData) {
      let pStart, pEnd;
      if (orf.frame > 0) {
        pStart = Math.max(0, orf.start - storedPromoterLen);
        pEnd   = orf.start;
      } else {
        pStart = orf.end;
        pEnd   = Math.min(seqLen, orf.end + storedPromoterLen);
      }

      const key = `${pStart}-${pEnd}`;
      if (drawn.has(key)) continue;
      drawn.add(key);

      if (pEnd <= bpStart || pStart >= bpEnd) continue;

      const x1 = Math.max(0, xForBp(pStart, cw));
      const x2 = Math.min(cw, xForBp(pEnd, cw));
      const w  = x2 - x1;
      if (w <= 0) continue;

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
        ctx.save();
        ctx.fillStyle    = '#ffffff';
        ctx.font         = '9px "JetBrains Mono","Courier New",monospace';
        ctx.textAlign    = 'left';
        ctx.textBaseline = 'middle';
        ctx.fillText(m.name, x1 + 3, t.y + TRACK_HEIGHT / 2);
        ctx.restore();
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

  // Fix 1: script is lazy-injected after DOMContentLoaded has fired — call directly.
  initViewer();

  const btn = document.getElementById('gvLoadDemoBtn');
  if (btn) {
    btn.addEventListener('click', () => {
      // Fix 2: use the defined DEMO_SEQUENCE (pUC19 lacZ-alpha, GenBank L09137)
      // with biologically correct features derived from the actual sequence.
      const demoCDS = [
        // lacZ-alpha ORF: ATG at position 30 (0-based), frame +1
        // The ORF extends beyond this fragment; end is set to sequence end.
        {
          start:  30,
          end:    DEMO_SEQUENCE.length,
          frame:  1,
          length: DEMO_SEQUENCE.length - 30,
          name:   'lacZ-alpha (partial)'
        }
      ];

      // 30 bp upstream of the ATG encompasses the Shine-Dalgarno site
      const demoPromoterLen = 30;

      const demoMotifs = [
        // Shine-Dalgarno: AGGAAA at positions 19–24 (0-based)
        { position: 19, end: 25, matched: 'AGGAAA', name: 'Shine-Dalgarno' },
        // HindIII: AAGCTT at positions 47–52 (0-based)
        { position: 47, end: 53, matched: 'AAGCTT', name: 'HindIII'        },
        // SalI: GTCGAC at positions 65–70 (0-based, start of MCS)
        { position: 65, end: 71, matched: 'GTCGAC', name: 'SalI (MCS)'    },
        // BamHI: GGATCC at positions 71–76 (0-based)
        { position: 71, end: 77, matched: 'GGATCC', name: 'BamHI (MCS)'   },
        // EcoRI: GAATTC at positions 81–86 (0-based)
        { position: 81, end: 87, matched: 'GAATTC', name: 'EcoRI (MCS)'   }
      ];

      window.loadGenomeData(DEMO_SEQUENCE, demoCDS, demoPromoterLen, demoMotifs);
      window.setZoom('10bp');
    });
  }

  // ── Theme Sync ───────────────────────────────────────────────
  window.addEventListener('biokit-theme-change', () => {
    themeColorCache = null;
    const panel     = document.getElementById('panel-genome-viewer');
    const isVisible = panel && !panel.classList.contains('hidden');
    if (isVisible) requestAnimationFrame(render);
  });

})();
