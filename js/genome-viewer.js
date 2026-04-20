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

  function getThemeColors() {
    const style = getComputedStyle(document.documentElement);
    return {
      bg: style.getPropertyValue('--bg-base').trim() || '#111827',
      surface: style.getPropertyValue('--bg-sidebar').trim() || '#1e293b',
      text: style.getPropertyValue('--text-primary').trim() || '#f8f9fa',
      muted: style.getPropertyValue('--text-muted').trim() || '#94a3b8',
      border: style.getPropertyValue('--border').trim() || '#334155',
      grid: style.getPropertyValue('--border-muted').trim() || '#1e293b'
    };
  }
  
  // Viewport
  let bpStart = 0;       // Leftmost base pair visible
  let bpEnd = 1000;      // Rightmost base pair visible
  let zoomLevel = 1;     // scale multiplier
  let isDragging = false;
  let dragStartX = 0;
  let dragStartBp = 0;

  // Track Layout Settings
  const TRACK_HEIGHT = 40;
  const HEADER_HEIGHT = 30;
  const MARGIN_TOP = 20;
  const TRACK_SPACING = 30;

  const TRACKS = [
    { id: 'seq', label: 'DNA Sequence', y: 0 },
    { id: 'cds', label: 'CDS / ORFs', y: 0 },
    { id: 'promoter', label: 'Promoter Region', y: 0 },
    { id: 'motif', label: 'Cis-Elements', y: 0 }
  ];

  function initViewer() {
    canvas = document.getElementById('genomeCanvas');
    if (!canvas) return;
    
    ctx = canvas.getContext('2d', { alpha: false });
    
    // Handle resizing
    window.addEventListener('resize', resizeCanvas);
    resizeCanvas();

    // Mouse Events for Panning
    canvas.addEventListener('mousedown', onMouseDown);
    window.addEventListener('mousemove', onMouseMove);
    window.addEventListener('mouseup', onMouseUp);

    // Mouse Events for Clicking (Interactivity)
    canvas.addEventListener('click', onClick);

    // Initial render
    requestAnimationFrame(render);
  }

  function resizeCanvas() {
    const parent = canvas.parentElement;
    if (!parent) return;

    const rect = parent.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;

    canvas.width = rect.width * dpr;
    canvas.height = 400 * dpr;
    
    canvas.style.width = `${rect.width}px`;
    canvas.style.height = `400px`;

    // Reset transform THEN scale to prevent cumulative multiplication
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.scale(dpr, dpr);

    let currentY = MARGIN_TOP + HEADER_HEIGHT;
    for (let t of TRACKS) {
      t.y = currentY;
      currentY += TRACK_HEIGHT + TRACK_SPACING;
    }

    requestAnimationFrame(render);
  }

  // ── Data Loading ────────────────────────────────────────────
  window.loadGenomeData = function(seq, cds, promoterLen, motifs) {
    sequence = seq || "";
    cdsData = cds || [];
    motifData = motifs || [];
    
    // Set view to start of sequence
    bpStart = 0;
    bpEnd = Math.min(100, sequence.length || 100); // start zoomed in to see bases
    
    resizeCanvas(); // Trigger a render
  }

  window.setZoom = function(level) {
    // level: '1bp', '10bp', '100bp', '1000bp'
    if (!sequence) return;
    const centerBp = bpStart + (bpEnd - bpStart) / 2;
    let windowSize = 100;

    if (level === '1bp') windowSize = 50;
    if (level === '10bp') windowSize = 250;
    if (level === '100bp') windowSize = 2500;
    if (level === '1000bp') windowSize = 25000;

    windowSize = Math.min(windowSize, sequence.length);

    bpStart = Math.max(0, centerBp - windowSize / 2);
    bpEnd = Math.min(sequence.length, centerBp + windowSize / 2);
    
    // adjust if hitting boundaries
    if (bpStart === 0) bpEnd = Math.min(sequence.length, windowSize);
    if (bpEnd === sequence.length) bpStart = Math.max(0, sequence.length - windowSize);

    requestAnimationFrame(render);
  }

  // ── Interaction ─────────────────────────────────────────────
  function onMouseDown(e) {
    isDragging = true;
    dragStartX = e.clientX;
    dragStartBp = bpStart;
    canvas.style.cursor = 'grabbing';
  }

  function onMouseMove(e) {
    if (!isDragging || !sequence) return;
    
    const dx = e.clientX - dragStartX;
    const cw = canvas.clientWidth;
    const bpWidth = bpEnd - bpStart;
    
    // Calculate how many base pairs the mouse movement represents
    const bpShift = (dx / cw) * bpWidth;
    
    bpStart = dragStartBp - bpShift;
    bpEnd = bpStart + bpWidth;

    // Clamp to boundaries
    if (bpStart < 0) {
      bpStart = 0;
      bpEnd = bpWidth;
    }
    if (bpEnd > sequence.length) {
      bpEnd = sequence.length;
      bpStart = Math.max(0, sequence.length - bpWidth);
    }

    requestAnimationFrame(render);
  }

  function onMouseUp() {
    isDragging = false;
    if (canvas) canvas.style.cursor = 'grab';
  }

  function onClick(e) {
      if (!sequence) return;
      const rect = canvas.getBoundingClientRect();
      const x = e.clientX - rect.left;
      const y = e.clientY - rect.top;

      const cw = rect.width;
      const bpWidth = bpEnd - bpStart;
      const clickBp = bpStart + (x / cw) * bpWidth;

      // Check CDS clicks
      const cdsTrack = TRACKS.find(t => t.id === 'cds');
      if (y >= cdsTrack.y && y <= cdsTrack.y + TRACK_HEIGHT) {
          for (const cds of cdsData) {
              if (clickBp >= cds.start && clickBp <= cds.end) {
                  showCanvasTooltip(`CDS Frame +${cds.frame}\nPos: ${cds.start}-${cds.end}\nLen: ${cds.length}bp`, e.clientX, e.clientY);
                  return;
              }
          }
      }

      // Check Motif clicks
      const motifTrack = TRACKS.find(t => t.id === 'motif');
      if (y >= motifTrack.y && y <= motifTrack.y + TRACK_HEIGHT) {
          for (const m of motifData) {
              const mEnd = m.position + m.matched.length;
              if (clickBp >= m.position && clickBp <= mEnd) {
                   showCanvasTooltip(`Motif: ${m.name}\nPos: ${m.position}\nSeq: ${m.matched}`, e.clientX, e.clientY);
                   
                   // Highlight by snapping to it briefly
                   const center = m.position + (m.matched.length/2);
                   const w = bpEnd - bpStart;
                   bpStart = Math.max(0, center - w/2);
                   bpEnd = Math.min(sequence.length, center + w/2);
                   requestAnimationFrame(render);
                   return;
              }
          }
      }
      hideCanvasTooltip();
  }

  // Tooltip DOM
  function showCanvasTooltip(text, x, y) {
      let tt = document.getElementById('gvTooltip');
      if (!tt) {
          tt = document.createElement('div');
          tt.id = 'gvTooltip';
          tt.style.position = 'fixed';
          tt.style.background = 'var(--panel-bg)';
          tt.style.color = '#fff';
          tt.style.border = '1px solid var(--border-color)';
          tt.style.padding = '8px 12px';
          tt.style.borderRadius = '6px';
          tt.style.boxShadow = '0 10px 15px -3px rgba(0,0,0,0.5)';
          tt.style.pointerEvents = 'none';
          tt.style.zIndex = '1000';
          tt.style.fontFamily = 'var(--ui-font)';
          tt.style.fontSize = '0.85rem';
          tt.style.whiteSpace = 'pre-line';
          document.body.appendChild(tt);
      }
      tt.textContent = text;
      tt.style.left = (x + 15) + 'px';
      tt.style.top = (y + 15) + 'px';
      tt.style.display = 'block';
  }
  function hideCanvasTooltip() {
      const tt = document.getElementById('gvTooltip');
      if (tt) tt.style.display = 'none';
  }

  // ── Rendering Engine ────────────────────────────────────────
  function render() {
    if (!ctx || !canvas) return;

    const cw = canvas.clientWidth;
    const ch = canvas.clientHeight;
    const colors = getThemeColors();

    // 1. Clear stale frame, then draw background
    ctx.clearRect(0, 0, cw, ch);
    ctx.fillStyle = colors.bg;
    ctx.fillRect(0, 0, cw, ch);

    if (!sequence) {
      ctx.fillStyle = colors.muted;
      ctx.font = '14px Inter';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('No sequence loaded. Analyze a sequence to view.', cw / 2, ch / 2);
      return;
    }

    // 2. Draw Ruler (Header)
    drawRuler(cw);

    // 3. Draw Tracks
    drawSequenceTrack(cw);
    drawCDSTrack(cw);
    drawPromoterTrack(cw);
    drawMotifTrack(cw);
  }

  function xForBp(bp, cw) {
    return ((bp - bpStart) / (bpEnd - bpStart)) * cw;
  }

  // Binary search to find the first index of an array where item.pos/start >= bpStart
  function findFirstVisibleIndex(arr, startKey, minBp) {
    let low = 0, high = arr.length - 1;
    let res = arr.length;
    while (low <= high) {
      let mid = (low + high) >> 1;
      if (arr[mid][startKey] >= minBp) {
          res = mid;
          high = mid - 1;
      } else {
          low = mid + 1;
      }
    }
    // Step back once to catch items that might span into the view from the left
    return Math.max(0, res - 1);
  }

  function drawRuler(cw) {
    const colors = getThemeColors();
    ctx.fillStyle = colors.surface;
    ctx.fillRect(0, 0, cw, HEADER_HEIGHT);
    
    ctx.fillStyle = colors.muted;
    ctx.font = '11px JetBrains Mono';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    const bpWidth = bpEnd - bpStart;
    
    // Determine tick interval based on zoom level
    let tickInterval = 10;
    if (bpWidth > 200) tickInterval = 50;
    if (bpWidth > 1000) tickInterval = 200;
    if (bpWidth > 5000) tickInterval = 1000;
    if (bpWidth > 20000) tickInterval = 5000;

    const firstTick = Math.ceil(bpStart / tickInterval) * tickInterval;

    for (let bp = firstTick; bp <= bpEnd; bp += tickInterval) {
      const x = xForBp(bp, cw);
      
      // Draw tick mark
      ctx.beginPath();
      ctx.moveTo(x, HEADER_HEIGHT - 6);
      ctx.lineTo(x, HEADER_HEIGHT);
      ctx.strokeStyle = colors.border;
      ctx.stroke();

      // Draw label
      ctx.fillText(bp.toLocaleString(), x, HEADER_HEIGHT / 2);
    }
    
    ctx.beginPath();
    ctx.moveTo(0, HEADER_HEIGHT);
    ctx.lineTo(cw, HEADER_HEIGHT);
    ctx.strokeStyle = colors.border;
    ctx.stroke();
  }

  function drawTrackBackground(t, cw) {
    const colors = getThemeColors();
    // Track Label
    ctx.fillStyle = colors.muted;
    ctx.font = '12px Inter';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'bottom';
    ctx.fillText(t.label, 10, t.y - 4);
    
    // Track bounds (light alternating backgrounds could go here, but keep it clean)
    ctx.beginPath();
    ctx.moveTo(0, t.y + TRACK_HEIGHT);
    ctx.lineTo(cw, t.y + TRACK_HEIGHT);
    ctx.strokeStyle = colors.grid;
    ctx.stroke();
  }

  function drawSequenceTrack(cw) {
    const t = TRACKS.find(t => t.id === 'seq');
    drawTrackBackground(t, cw);

    const bpWidth = bpEnd - bpStart;
    
    // Only render actual letters if zoomed in enough (< 150bp)
    if (bpWidth < 150) {
        ctx.font = '14px JetBrains Mono';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        
        const startIdx = Math.max(0, Math.floor(bpStart));
        const endIdx = Math.min(sequence.length, Math.ceil(bpEnd));
        
        for (let i = startIdx; i < endIdx; i++) {
           const char = sequence[i];
           const x = xForBp(i + 0.5, cw); // center of the bp
           
           // Color map
           if (char === 'A') ctx.fillStyle = '#22c55e';
           else if (char === 'T' || char === 'U') ctx.fillStyle = '#ef4444';
           else if (char === 'G') ctx.fillStyle = '#f97316';
           else if (char === 'C') ctx.fillStyle = '#3b82f6';
           else ctx.fillStyle = '#94a3b8';
           
           ctx.fillText(char, x, t.y + TRACK_HEIGHT / 2);
        }
    } else {
        // Render a density map or solid block representing sequence presence
        const colors = getThemeColors();
        ctx.fillStyle = colors.border;
        ctx.fillRect(0, t.y + 10, cw, TRACK_HEIGHT - 20);
    }
  }

  function drawCDSTrack(cw) {
    const t = TRACKS.find(t => t.id === 'cds');
    drawTrackBackground(t, cw);

    // Draw ORFs as distinct colored arrows
    const bpWidth = bpEnd - bpStart;
    
    // Sort if not already sorted (should be sorted by results logic)
    cdsData.sort((a, b) => a.start - b.start);
    
    const firstIdx = findFirstVisibleIndex(cdsData, 'start', bpStart);

    for (let i = firstIdx; i < cdsData.length; i++) {
       const cds = cdsData[i];
       if (cds.start > bpEnd) break; // Finished visible range
       
       const x1 = Math.max(0, xForBp(cds.start, cw));
       let x2 = Math.min(cw, xForBp(cds.end, cw));
       let w = x2 - x1;
       if (w < 2) w = 2; // minimum visibility
       
       // Cycle colors for diff frames
       const colors = ['#14b8a6', '#8b5cf6', '#ec4899'];
       ctx.fillStyle = colors[i % 3];
       
       const y = t.y + 10;
       const h = 20;
       
       // Arrow shape logic
       ctx.beginPath();
       ctx.moveTo(x1, y);
       
       // if end is visible, draw arrowhead
       if (cds.end <= bpEnd && w > 10) {
           ctx.lineTo(x1 + w - 10, y);
           ctx.lineTo(x1 + w, y + h/2);
           ctx.lineTo(x1 + w - 10, y + h);
       } else {
           ctx.lineTo(x1 + w, y);
           ctx.lineTo(x1 + w, y + h);
       }
       ctx.lineTo(x1, y + h);
       ctx.fill();
       
       // Label if wide enough
       if (w > 40 && bpWidth < 5000) {
           ctx.fillStyle = '#ffffff';
           ctx.font = '10px Inter';
           ctx.textAlign = 'left';
           ctx.textBaseline = 'middle';
           ctx.fillText(`ORF ${cds.frame}`, x1 + 5, y + h/2);
       }
    }
  }

  function drawPromoterTrack(cw) {
    const t = TRACKS.find(t => t.id === 'promoter');
    drawTrackBackground(t, cw);
    
    // In this toolkit, promoter regions are often defined relatively early.
    // For now we will highlight [0 to first CDS] if it exists, or just a mock region if passed.
    // Let's rely on global context or a specific variable passed to `loadGenomeData`.
    // If no explicit promoter data, we skip. But we have motifs within the promoter. Let's assume the promoter spans the first 1000bp or up to first CDS.
    if (cdsData.length > 0) {
        const pStart = Math.max(0, cdsData[0].start - 1000);
        const pEnd = cdsData[0].start;
        
        if (pEnd > bpStart && pStart < bpEnd) {
             const x1 = Math.max(0, xForBp(pStart, cw));
             const x2 = Math.min(cw, xForBp(pEnd, cw));
             
             ctx.fillStyle = 'rgba(245, 158, 11, 0.2)'; // amber transp
             ctx.fillRect(x1, t.y + 10, (x2 - x1), 20);
             
             ctx.strokeStyle = '#f59e0b';
             ctx.lineWidth = 1;
             ctx.strokeRect(x1, t.y + 10, (x2 - x1), 20);
        }
    }
  }

  function drawMotifTrack(cw) {
    const t = TRACKS.find(t => t.id === 'motif');
    drawTrackBackground(t, cw);

    motifData.sort((a, b) => a.position - b.position);
    const firstIdx = findFirstVisibleIndex(motifData, 'position', bpStart);

    for (let i = firstIdx; i < motifData.length; i++) {
       const m = motifData[i];
       if (m.position > bpEnd) break;
       
       const mEnd = m.position + m.matched.length;
       const x1 = Math.max(0, xForBp(m.position, cw));
       const x2 = Math.min(cw, xForBp(mEnd, cw));
       let w = x2 - x1;
       if (w < 2) w = 2; // min visible
       
       ctx.fillStyle = '#22d3ee'; // cyan
       if (m.name.includes("TATA")) ctx.fillStyle = '#ef4444';
       if (m.name.includes("CAAT")) ctx.fillStyle = '#3b82f6';
       
       ctx.fillRect(x1, t.y + 5, w, TRACK_HEIGHT - 10);
    }
  }

  // pUC19 lacZ-alpha region: lac promoter (-35/-10 boxes), Shine-Dalgarno, ATG start codon (GenBank L09137)
  const DEMO_SEQUENCE =
    'GCGGATAACAATTTCACACAGGAAACAGCTATGACCATGATTACGCCAAGCTTGCATGCCTGCAG' +
    'GTCGACGGATCCCCGGGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCA' +
    'TGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATT' +
    'CCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACT' +
    'CACATTAATTGCGTTACGCTGCGGTTTTCATGAGAATGTTTTTCTTTTCATGAGAAAAGCCCGGC' +
    'TTATAGTTTGCTTTTTATTTGATTTGAGTAATTTTGTTTTTTATACTATTTTTTTTTGAGCTTTT' +
    'GTTCGTTCAGAGTTTATTCGCTTCATTTAAATGGTATGAAATTTACTGATAATGATATTTTTAT';

  // Initialize on load
  document.addEventListener('DOMContentLoaded', () => {
     initViewer();

     const btn = document.getElementById('gvLoadDemoBtn');
     if (btn) {
         btn.addEventListener('click', () => {
             const input = document.getElementById('proInput');
             const analyzeBtn = document.getElementById('proAnalyzeBtn');
             if (!input || !analyzeBtn) return;

             input.value = DEMO_SEQUENCE;
             input.dispatchEvent(new Event('input'));
             analyzeBtn.click();

             const panel = document.getElementById('panel-genome-viewer');
             if (panel && !panel.querySelector('.demo-banner')) {
               const banner = document.createElement('div');
               banner.className = 'demo-banner';
               banner.style.cssText = 'display:flex;align-items:center;justify-content:space-between;gap:8px;padding:8px 12px;margin-bottom:10px;background:var(--warning-bg,#fef3c7);border:1px solid var(--warning-border,#f59e0b);border-radius:6px;font-size:0.85rem;color:var(--warning-text,#92400e);';
               banner.innerHTML = '<span>⚠ Viewing Educational Demo. Paste your own sequence to calculate custom coordinates.</span><button style="background:none;border:none;cursor:pointer;font-size:1rem;line-height:1;padding:0 2px;color:inherit;" aria-label="Dismiss">&times;</button>';
               banner.querySelector('button').addEventListener('click', () => banner.remove());
               panel.insertBefore(banner, panel.firstChild);
             }
         });
     }
  });

  // ── Theme Sync ──────────────────────────────────────────
  window.addEventListener('biokit-theme-change', () => {
    const isVisible = !document.getElementById('panel-genome-viewer').classList.contains('hidden');
    if (isVisible) {
      requestAnimationFrame(render);
    }
  });

})();
