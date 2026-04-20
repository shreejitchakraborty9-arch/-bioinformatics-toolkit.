/* ============================================================
   promoter.js — Promoter / CDS / ORF Analyzer
   All heavy computation delegated to /api/analyze/promoter (Python/Biopython PSSM).
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── Motif Display Metadata ────────────────────────────────
  // Maps backend motif names → CSS class and description for rendering.
  // Source: Bucher 1990 (NAR 18:6299), Breathnach & Chambon 1981.
  const MOTIF_META = {
    'TATA-box': {
      cssClass: 'motif-red',
      desc:     'Core promoter element (Bucher 1990); peak ~−30 from TSS'
    },
    'CAAT-box': {
      cssClass: 'motif-yellow',
      desc:     'Common cis-acting element in promoter and enhancer regions (Breathnach & Chambon 1981); peak ~−80 from TSS'
    },
    'GC-box': {
      cssClass: 'motif-violet',
      desc:     'SP1 binding site; enhancer element (Bucher 1990); peak ~−90 from TSS'
    },
    '-35 box': {
      cssClass: 'motif-green',
      desc:     'σ⁷⁰ promoter −35 hexamer (Harley & Reynolds 1987, NAR 15:2343); consensus TTGACA'
    },
    '-10 box': {
      cssClass: 'motif-amber',
      desc:     'σ⁷⁰ Pribnow box (Pribnow 1975; Harley & Reynolds 1987, NAR 15:2343); consensus TATAAT'
    }
  };

  // ── Results Cache (for Export) ────────────────────────────
  let _lastResults = { seq: '', motifs: [], orfs: [], upstreamLen: 1000 };

  // ── API Response Normalizers ──────────────────────────────
  // Adapt backend schema to the contract expected by displayResults().

  function normalizeORFs(apiOrfs) {
    return apiOrfs
      .map(o => ({
        frame:         o.frame,
        start:         o.start,
        end:           o.end,
        length:        o.length_nt,
        protein:       o.protein,
        proteinLength: o.length_aa,
        partial:       false   // backend only emits stop-codon-terminated ORFs
      }))
      .sort((a, b) => b.length - a.length);
  }

  function normalizeMotifs(apiMotifs) {
    return apiMotifs
      .map(m => {
        const meta = MOTIF_META[m.motif] || { cssClass: 'motif-slate', desc: '' };
        return {
          name:       m.motif,
          position:   m.start,
          matched:    m.sequence,
          cssClass:   meta.cssClass,
          desc:       meta.desc,
          score:      m.score_pct.toFixed(1) + '%',
          strand:     m.strand || '+',
          spacerToBp: m.spacer_to_partner_bp ?? null
        };
      })
      .sort((a, b) => a.position - b.position);
  }

  // ── Loading State Helpers ─────────────────────────────────

  function showLoadingState(panel, btn) {
    const resultsArea = panel.querySelector('.results-area');
    let progContainer = panel.querySelector('.progress-container');

    if (!progContainer) {
      progContainer = document.createElement('div');
      progContainer.className = 'progress-container';
      progContainer.innerHTML = `
        <div class="progress-bar-wrap">
          <div class="progress-bar"></div>
        </div>
        <div class="progress-meta">
          <span class="progress-text">Running PSSM motif scan &amp; 6-frame ORF detection…</span>
        </div>`;
      const btnRow = panel.querySelector('.btn-row');
      if (btnRow) btnRow.parentNode.insertBefore(progContainer, btnRow.nextSibling);
    }

    progContainer.classList.add('active');
    const bar = progContainer.querySelector('.progress-bar');
    bar.style.transition = 'none';
    bar.style.width = '0%';
    // After next paint, animate to show indeterminate progress
    requestAnimationFrame(() => {
      bar.style.transition = 'width 1.8s cubic-bezier(0.4, 0, 0.2, 1)';
      bar.style.width = '75%';
    });

    if (resultsArea) resultsArea.classList.add('hidden');
    if (btn) btn.disabled = true;

    return { progContainer, bar, resultsArea };
  }

  function hideLoadingState({ progContainer, bar, resultsArea }, btn) {
    bar.style.transition = 'width 0.2s ease-out';
    bar.style.width = '100%';
    setTimeout(() => {
      progContainer.classList.remove('active');
      if (resultsArea) resultsArea.classList.remove('hidden');
      if (btn) btn.disabled = false;
      if (window.lucide) window.lucide.createIcons();
    }, 220);
  }

  function abortLoadingState({ progContainer, bar, resultsArea }, btn) {
    progContainer.classList.remove('active');
    if (resultsArea) resultsArea.classList.remove('hidden');
    if (btn) btn.disabled = false;
  }

  // ── Main Run ──────────────────────────────────────────────

  async function run() {
    const input = document.getElementById('proInput');
    if (!input) return;

    const raw = input.value;
    const validation = window.validateSequence(raw, 'dna');

    if (!validation.valid) {
      window.showValidationWarning('panel-promoter', validation);
      return;
    }

    const seq          = validation.clean;
    const upstreamLen  = parseInt(document.getElementById('proUpstream').value) || 1000;
    const organismType = document.getElementById('proOrganismType')?.value || 'eukaryote';
    const btn          = document.getElementById('proAnalyzeBtn');
    const panel        = document.getElementById('panel-promoter');

    if (!panel) return;

    const loadState = showLoadingState(panel, btn);

    try {
      const res = await fetch('/api/analyze/promoter', {
        method:  'POST',
        headers: { 'Content-Type': 'application/json' },
        body:    JSON.stringify({ sequence: seq, upstream_length: upstreamLen, organism_type: organismType })
      });

      if (!res.ok) {
        const errData = await res.json().catch(() => ({}));
        throw new Error(errData.error || `Server responded with status ${res.status}`);
      }

      const data = await res.json();
      if (data.error) throw new Error(data.error);

      const motifs = normalizeMotifs(data.motifs || []);
      const orfs   = normalizeORFs(data.orfs || []);

      // Cache for export
      _lastResults = { seq, motifs, orfs, upstreamLen };

      hideLoadingState(loadState, btn);

      // Small delay matches hideLoadingState's DOM reset timing
      setTimeout(() => displayResults(seq, motifs, orfs, upstreamLen), 240);

    } catch (err) {
      abortLoadingState(loadState, btn);
      window.showToast?.(`Analysis failed: ${err.message}`, 'error');
      console.error('[Promoter] /api/analyze/promoter error:', err);
    }
  }

  // ── Helper: Extract Spliced CDS (Annotation Mode) ────────
  function extractCDS(seq, regions, strand) {
    if (!regions || regions.length === 0) return '';
    const BioMath = window.BioKit.core.BioMath;
    const sorted = [...regions].sort((a, b) => strand === '-' ? b[0] - a[0] : a[0] - b[0]);
    let full = '';
    sorted.forEach(([start, end]) => {
      let chunk = seq.substring(Math.max(0, start - 1), Math.min(seq.length, end));
      if (strand === '-') chunk = BioMath.reverseComplement(chunk);
      full += chunk;
    });
    return full;
  }

  // Sliding-window GC-density scan for promoter-region annotation
  function detectPromoterMotifs(sequence) {
    const hits = [];
    const windowSize = 50;
    const threshold  = 0.6;

    if (!sequence || sequence.length < windowSize) return hits;

    for (let i = 0; i <= sequence.length - windowSize; i++) {
      const chunk = sequence.substring(i, i + windowSize);
      let gcCount = 0;
      for (let j = 0; j < chunk.length; j++) {
        const b = chunk[j].toUpperCase();
        if (b === 'G' || b === 'C') gcCount++;
      }
      const gcContent = gcCount / windowSize;
      if (gcContent >= threshold) {
        hits.push({
          name:     'GC-Rich Region',
          position: i + 1,
          matched:  chunk.substring(0, 8) + '...',
          cssClass: 'motif-sky',
          desc:     `High GC density region (${(gcContent * 100).toFixed(0)}%)`,
          score:    (gcContent * 100).toFixed(1) + '%'
        });
        i += 15;
      }
    }
    return hits;
  }

  // ── Display ───────────────────────────────────────────────

  function displayResults(seq, motifs, orfs, upstreamLen) {
    const container = document.getElementById('proResults');
    if (!container) return;
    container.classList.remove('hidden');

    // Identify Anchor ORF
    let anchor = null;
    if (orfs.length > 0) {
      const metaText = document.getElementById('proInputMeta')?.textContent || '';
      const hintMinus = metaText.includes('[strand: -]');
      const hintPlus  = metaText.includes('[strand: +]');

      if (hintMinus)      anchor = orfs.find(o => o.frame < 0) || orfs[0];
      else if (hintPlus)  anchor = orfs.find(o => o.frame > 0) || orfs[0];
      else                anchor = orfs[0];
    }

    window.buildStatCards('proStatRow', [
      { val: seq.length.toLocaleString(), label: 'Length (bp)' },
      { val: motifs.length,               label: 'Cis-Elements' },
      { val: orfs.length,                 label: 'ORFs Found'   },
      { val: anchor ? anchor.proteinLength + ' aa' : '—', label: 'Largest ORF' },
    ]);

    renderORFTable(orfs);
    renderMotifTable(motifs);
    renderLegend(motifs);

    const activeAnalysis = window.BioKit.activeAnalysis;
    const BioMath = window.BioKit.core.BioMath;

    let meta = activeAnalysis?.metadata;
    if (activeAnalysis?.transcripts?.length > 0) {
      meta = activeAnalysis.transcripts[0];
    }

    if (meta && (meta.cds_regions?.length > 0 || meta.tss)) {
      // ── A) DATABASE ANNOTATION MODE ─────────────────────────
      const strand  = meta.strand || '+';
      const isMinus = strand === '-';
      const tss     = meta.tss || 1;

      const cdsSeq = extractCDS(seq, meta.cds_regions, strand);
      const txInfo = meta.id && meta.id !== 'undefined' ? `Isoform ${meta.id}` : 'Annotated Transcript';
      document.getElementById('proCDSSeqRange').textContent = `${txInfo} · ${cdsSeq.length} bp · Strand [${strand}]`;
      renderPlainSequence(cdsSeq, 'proCDSSeqViewer');

      const protein = BioMath.translateDNA(cdsSeq);
      document.getElementById('proProteinSeqRange').textContent = `Annotated Translation · ${protein.length} aa`;
      renderPlainSequence(protein, 'proProteinSeqViewer', true);

      let proStart, proEnd, promoterSeq;
      if (!isMinus) {
        proStart     = Math.max(1, tss - upstreamLen);
        proEnd       = tss - 1;
        promoterSeq  = seq.substring(proStart - 1, proEnd);
      } else {
        proStart     = tss + 1;
        proEnd       = Math.min(seq.length, tss + upstreamLen);
        promoterSeq  = BioMath.reverseComplement(seq.substring(proStart - 1, proEnd));
      }

      const deepMotifs   = detectPromoterMotifs(promoterSeq);
      const mappedMotifs = motifs
        .filter(m => m.position >= proStart && m.position <= proEnd)
        .map(m => ({ ...m, position: m.position - (proStart - 1) }));
      const finalPromoterMotifs = [...mappedMotifs, ...deepMotifs].sort((a, b) => a.position - b.position);

      document.getElementById('proSeqRange').textContent =
        `Promoter [TSS ${tss}] · ${promoterSeq.length} bp ${isMinus ? '(RC)' : '(FWD)'}`;
      if (promoterSeq) {
        renderHighlightedSequence(promoterSeq, finalPromoterMotifs, 'proSeqViewer');
      } else {
        document.getElementById('proSeqViewer').innerHTML =
          '<p class="text-muted">Insufficient upstream sequence for promoter extraction.</p>';
      }

    } else if (anchor) {
      // ── B) PREDICTION MODE (ORF Anchor) ─────────────────────
      const isMinus = anchor.frame < 0;

      let cdsSeq = seq.substring(anchor.start - 1, anchor.end);
      if (isMinus) cdsSeq = BioMath.reverseComplement(cdsSeq);

      document.getElementById('proCDSSeqRange').textContent =
        `Predicted CDS [${anchor.start} … ${anchor.end}] · ${anchor.length} bp · Frame ${anchor.frame > 0 ? '+' : ''}${anchor.frame}`;
      renderPlainSequence(cdsSeq, 'proCDSSeqViewer');

      document.getElementById('proProteinSeqRange').textContent =
        `Predicted Translation · ${anchor.proteinLength} aa`;
      renderPlainSequence(anchor.protein, 'proProteinSeqViewer', true);

      let proStart, proEnd, promoterSeq;
      if (!isMinus) {
        proStart    = Math.max(1, anchor.start - upstreamLen);
        proEnd      = anchor.start - 1;
        promoterSeq = seq.substring(proStart - 1, proEnd);
      } else {
        proStart    = anchor.end + 1;
        proEnd      = Math.min(seq.length, anchor.end + upstreamLen);
        promoterSeq = BioMath.reverseComplement(seq.substring(proStart - 1, proEnd));
      }

      const proMotifs = motifs
        .filter(m => m.position >= proStart && m.position <= proEnd)
        .map(m => ({ ...m, position: m.position - (proStart - 1) }));

      document.getElementById('proSeqRange').textContent =
        `Predicted Promoter [${proStart} … ${proEnd}] · ${promoterSeq.length} bp ${isMinus ? '(RC)' : '(FWD)'}`;
      if (promoterSeq) {
        renderHighlightedSequence(promoterSeq, proMotifs, 'proSeqViewer');
      } else {
        document.getElementById('proSeqViewer').innerHTML =
          '<p class="text-muted">Insufficient sequence for promoter extraction.</p>';
      }
    } else {
      ['proCDSSeqViewer', 'proProteinSeqViewer', 'proSeqViewer'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.innerHTML =
          '<p style="color:var(--text-muted)">No annotation or ORF found to anchor sequence extraction.</p>';
      });
      document.querySelectorAll('.pro-seq-meta').forEach(el => { el.textContent = '—'; });
    }

    // Feed genome viewer
    if (window.loadGenomeData) {
      const cdsData   = orfs.map(o => ({ start: o.start - 1, end: o.end, frame: o.frame, length: o.length }));
      const motifData = motifs.map(m => ({ position: m.position - 1, matched: m.matched, name: m.name }));
      window.loadGenomeData(seq, cdsData, 200, motifData);
    }
  }

  function renderPlainSequence(seq, viewerId, isProtein = false) {
    const viewer = document.getElementById(viewerId);
    if (!viewer || !seq) return;
    viewer.innerHTML = `<div class="seq-display" style="font-family:var(--mono); line-height:1.6; word-break:break-all; color:${isProtein ? 'var(--violet)' : 'var(--teal)'}">${window.escapeHTML(seq)}</div>`;
  }

  function renderHighlightedSequence(seq, motifs, viewerId) {
    const viewer = document.getElementById(viewerId);
    if (!viewer) return;

    const hlMap = new Array(seq.length).fill(null);
    motifs.forEach(m => {
      const end = Math.min(seq.length, m.position - 1 + m.matched.replace('...', '').length);
      for (let i = m.position - 1; i < end; i++) {
        if (!hlMap[i]) hlMap[i] = m.cssClass;
      }
    });

    let html = '';
    for (let i = 0; i < seq.length; i++) {
      const c = window.escapeHTML(seq[i]);
      html += hlMap[i] ? `<span class="motif-highlight ${hlMap[i]}">${c}</span>` : c;
    }
    viewer.innerHTML = window.sanitizeHTML(html);
  }

  function renderLegend(motifs) {
    const legend = document.getElementById('proCisLegend');
    if (!legend) return;
    const seen = new Set();
    let html = '';
    motifs.forEach(m => {
      if (seen.has(m.name)) return;
      seen.add(m.name);
      html += `<span class="motif-legend-item"><span class="motif-dot ${m.cssClass}"></span>${window.escapeHTML(m.name)}</span>`;
    });
    legend.innerHTML = window.sanitizeHTML(html);
  }

  // Pair prokaryotic -35/-10 box hits by exact spacer-derived position.
  // +strand: -10 starts at (-35 start) + 6 + spacer
  // -strand: -10 starts at (-35 start) - spacer - 6
  // Source: Harley & Reynolds 1987 spacer convention.
  function _pairProkaryoticMotifs(motifs) {
    const boxes35    = motifs.filter(m => m.name === '-35 box');
    const boxes10    = motifs.filter(m => m.name === '-10 box');
    const pairs      = [];
    const usedIdx10  = new Set();

    for (const b35 of boxes35) {
      if (b35.spacerToBp === null) continue;
      const expFwd    = b35.position + 6 + b35.spacerToBp;
      const expRev    = b35.position - b35.spacerToBp - 6;
      const partnerIdx = boxes10.findIndex((b10, i) =>
        !usedIdx10.has(i) && (b10.position === expFwd || b10.position === expRev)
      );
      if (partnerIdx !== -1) {
        usedIdx10.add(partnerIdx);
        pairs.push({ b35, b10: boxes10[partnerIdx] });
      }
    }
    return pairs;
  }

  function renderMotifTable(motifs) {
    const table = document.getElementById('proCisTable');
    if (!table) return;

    if (motifs.length === 0) {
      table.innerHTML = '<p style="color:var(--text-muted)">No cis-elements found.</p>';
      return;
    }

    const hasProkMotifs = motifs.some(m => m.name === '-35 box' || m.name === '-10 box');

    let html = `<table class="data-table"><thead><tr>
      <th>Element</th><th>Position</th><th>Sequence</th><th>PSSM Score</th><th>Description</th>
    </tr></thead><tbody>`;

    if (hasProkMotifs) {
      const pairs       = _pairProkaryoticMotifs(motifs);
      const pairedKeys  = new Set();
      pairs.forEach(({ b35, b10 }) => {
        pairedKeys.add(`${b35.name}|${b35.position}`);
        pairedKeys.add(`${b10.name}|${b10.position}`);
      });

      pairs.forEach(({ b35, b10 }) => {
        const spacerLabel = b35.spacerToBp !== null ? `Spacer: ${b35.spacerToBp} bp` : '';
        html += `<tr style="background:rgba(0,212,170,0.04)">
          <td>
            <div style="display:flex;flex-direction:column;gap:4px;">
              <span class="motif-badge motif-sky" style="white-space:nowrap">−35/−10 Complex</span>
              ${spacerLabel ? `<span style="font-size:0.7rem;color:var(--text-muted);padding-left:2px">${window.escapeHTML(spacerLabel)}</span>` : ''}
            </div>
          </td>
          <td>
            <div style="display:flex;flex-direction:column;gap:3px;font-size:0.85rem;">
              <span><span class="motif-badge motif-green" style="font-size:0.65rem;padding:1px 6px;margin-right:4px">−35</span>${b35.position}</span>
              <span><span class="motif-badge motif-amber" style="font-size:0.65rem;padding:1px 6px;margin-right:4px">−10</span>${b10.position}</span>
            </div>
          </td>
          <td style="font-family:var(--mono);color:var(--teal)">
            <div style="display:flex;flex-direction:column;gap:3px;">
              <span>${window.escapeHTML(b35.matched)}</span>
              <span>${window.escapeHTML(b10.matched)}</span>
            </div>
          </td>
          <td>
            <div style="display:flex;flex-direction:column;gap:3px;">
              <span style="font-weight:500">${window.escapeHTML(b35.score)}</span>
              <span style="font-weight:500">${window.escapeHTML(b10.score)}</span>
            </div>
          </td>
          <td style="color:var(--text-muted)">σ⁷⁰ promoter complex; validates 15–19 bp spacer (Harley &amp; Reynolds 1987, NAR 15:2343)</td>
        </tr>`;
      });

      // Render any unpaired hits individually (defensive — backend guarantees pairing)
      motifs.forEach(m => {
        if (!pairedKeys.has(`${m.name}|${m.position}`)) {
          html += `<tr>
            <td><span class="motif-badge ${m.cssClass}">${window.escapeHTML(m.name)}</span></td>
            <td>${m.position}</td>
            <td style="font-family:var(--mono);color:var(--teal)">${window.escapeHTML(m.matched)}</td>
            <td style="font-weight:500">${window.escapeHTML(m.score || '—')}</td>
            <td style="color:var(--text-muted)">${window.escapeHTML(m.desc)}</td>
          </tr>`;
        }
      });

    } else {
      motifs.forEach(m => {
        html += `<tr>
          <td><span class="motif-badge ${m.cssClass}">${window.escapeHTML(m.name)}</span></td>
          <td>${m.position}</td>
          <td style="font-family:var(--mono);color:var(--teal)">${window.escapeHTML(m.matched)}</td>
          <td style="font-weight:500">${window.escapeHTML(m.score || '—')}</td>
          <td style="color:var(--text-muted)">${window.escapeHTML(m.desc)}</td>
        </tr>`;
      });
    }

    html += '</tbody></table>';
    table.innerHTML = window.sanitizeHTML(html);
  }

  function renderORFTable(orfs) {
    const table = document.getElementById('proOrfTable');
    if (!table) return;

    if (orfs.length === 0) {
      table.innerHTML = '<p style="color:var(--text-muted)">No ORFs found.</p>';
      return;
    }

    let html = `<table class="data-table"><thead><tr>
      <th>Frame</th><th>Start</th><th>End</th><th>Length (nt)</th><th>Protein (aa)</th><th>Note</th>
    </tr></thead><tbody>`;
    orfs.forEach(o => {
      html += `<tr>
        <td>${o.frame > 0 ? '+' : ''}${o.frame}</td>
        <td>${o.start}</td>
        <td>${o.end}</td>
        <td>${o.length}</td>
        <td>${o.proteinLength}</td>
        <td style="color:var(--text-muted)">${o.partial ? '⚠ Partial (no stop)' : '✓ Complete'}</td>
      </tr>`;
    });
    html += '</tbody></table>';
    table.innerHTML = window.sanitizeHTML(html);
  }

  // ── Events ────────────────────────────────────────────────

  document.getElementById('proAnalyzeBtn')?.addEventListener('click', window.debounce(run, 300));

  function showDemoBanner(panelId) {
    const panel = document.getElementById(panelId);
    if (!panel || panel.querySelector('.demo-banner')) return;
    const banner = document.createElement('div');
    banner.className = 'demo-banner';
    banner.style.cssText = 'display:flex;align-items:center;justify-content:space-between;gap:8px;padding:8px 12px;margin-bottom:10px;background:var(--warning-bg,#fef3c7);border:1px solid var(--warning-border,#f59e0b);border-radius:6px;font-size:0.85rem;color:var(--warning-text,#92400e);';
    banner.innerHTML = '<span>⚠ Viewing Educational Demo. Paste your own sequence to calculate custom coordinates.</span><button style="background:none;border:none;cursor:pointer;font-size:1rem;line-height:1;padding:0 2px;color:inherit;" aria-label="Dismiss">&times;</button>';
    banner.querySelector('button').addEventListener('click', () => banner.remove());
    panel.insertBefore(banner, panel.firstChild);
  }

  document.getElementById('proSampleBtn')?.addEventListener('click', () => {
    const input = document.getElementById('proInput');
    if (input) {
      // E. coli lacZ promoter region + start of gene
      input.value = 'GCTTATTATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATTATAAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTAACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGTTGGAGTGACGGCAGTTATCTGGAAGATCAGG';
      input.dispatchEvent(new Event('input'));
      showDemoBanner('panel-promoter');
    }
  });

  document.getElementById('proClearBtn')?.addEventListener('click', () => {
    const input = document.getElementById('proInput');
    if (input) input.value = '';
    document.getElementById('proResults')?.classList.add('hidden');
    input?.dispatchEvent(new Event('input'));
  });

  document.getElementById('proExportBtn')?.addEventListener('click', () => {
    const { seq, orfs, upstreamLen } = _lastResults;
    if (!seq) return;

    const anchor = orfs.length > 0 ? orfs[0] : null;

    let text = `BioToolkit Promoter/CDS Analysis Report\n`;
    text += `Generated: ${new Date().toLocaleString()}\n`;
    text += `Source Sequence Length: ${seq.length} bp\n`;
    text += `Motif Engine: PSSM (Biopython/Bucher 1990)\n`;
    text += `ORF Engine: 6-frame scan (min 30 aa)\n`;
    text += `----------------------------------------\n\n`;

    if (anchor) {
      const cds      = seq.substring(anchor.start - 1, anchor.end);
      const proStart = Math.max(1, anchor.start - upstreamLen);
      const promoter = seq.substring(proStart - 1, anchor.start - 1);

      text += `>Primary_Anchor_CDS [Range: ${anchor.start}-${anchor.end}] [Frame: ${anchor.frame}]\n${cds}\n\n`;
      text += `>Primary_Anchor_Translation [${anchor.proteinLength} aa]\n${anchor.protein}\n\n`;
      text += `>Promoter_Region [Range: ${proStart}-${anchor.start - 1}] [Upstream: ${upstreamLen}bp]\n${promoter}\n\n`;
    } else {
      text += `No Anchor ORF found in query sequence.\n`;
    }

    window.BioKit.utils.downloadText(`biokit_analysis_${Date.now()}.txt`, text);
  });

})();
