/* ============================================================
   promoter.js — Promoter / CDS / ORF Analyzer (Phase 3 Upgrade)
   Full cis-element scanning + ORF finder
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── Cis-Element Patterns (Expanded PlantCARE / Common DB) ────────
  // Note: IUPAC codes are supported in regex forms (e.g. [GT] for K)
  const CIS_ELEMENTS = [
    // Core Promoter Elements
    { name: 'TATA-box',    regex: /TATA[AT]A[AT]/g,     category: 'Core', desc: 'Core promoter element around -30 of transcription start' },
    { name: 'CAAT-box',    regex: /GG[CT]CAATCT/g,      category: 'Core', desc: 'Common cis-acting element in promoter and enhancer regions' },
    { name: 'GC-motif',    regex: /GGGCGG/g,            category: 'Core', desc: 'Enhancer binding site for Sp1' },
    { name: 'Inr',         regex: /[CT][CT]AN[TA][CT][CT]/g, category: 'Core', desc: 'Initiator element, surrounds the transcription start site' },

    // Light Responsive Elements
    { name: 'G-box',       regex: /CACGTG/g,            category: 'Light Responsive', desc: 'Cis-acting regulatory element involved in light responsiveness' },
    { name: 'GT1-motif',   regex: /GGTTAA[AT][AT][AT]/g, category: 'Light Responsive', desc: 'Light responsive element' },
    { name: 'Sp1',         regex: /GGGCGG/g,            category: 'Light Responsive', desc: 'Light responsive element' },
    { name: 'I-box',       regex: /GATAA[GA]/g,         category: 'Light Responsive', desc: 'Part of a light responsive element' },
    { name: 'AE-box',      regex: /AGAAACAA/g,          category: 'Light Responsive', desc: 'Part of a module for light response' },
    { name: 'AAAC-motif',  regex: /AAAC[AT]CC/g,        category: 'Light Responsive', desc: 'Light responsive element' },
    { name: 'ACE',         regex: /AAAACGTT/g,          category: 'Light Responsive', desc: 'Cis-acting element involved in light responsiveness' },
    { name: 'ATCT-motif',  regex: /AATCGATG/g,          category: 'Light Responsive', desc: 'Part of a conserved DNA module involved in light responsiveness' },
    { name: 'Box 4',       regex: /ATTAAT/g,            category: 'Light Responsive', desc: 'Part of a conserved DNA module involved in light responsiveness' },
    { name: 'Box I',       regex: /TTTCAAA/g,           category: 'Light Responsive', desc: 'Light responsive element' },
    { name: 'Box II',      regex: /TGGTAAAT/g,          category: 'Light Responsive', desc: 'Part of a light responsive element' },
    { name: 'CAG-motif',   regex: /CAGAC[AT]C/g,        category: 'Light Responsive', desc: 'Light responsive element' },
    { name: 'Chs-CMA1a',   regex: /TTACTTAA/g,          category: 'Light Responsive', desc: 'Light responsive element' },
    { name: 'Chs-CMA2a',   regex: /TAGTAC/g,            category: 'Light Responsive', desc: 'Light responsive element' },
    { name: 'GA-motif',    regex: /ATAGATC[AT]C/g,      category: 'Light Responsive', desc: 'Part of a light responsive element' },
    { name: 'Gap-box',     regex: /CAAATGAA[AG]A/g,     category: 'Light Responsive', desc: 'Part of a light responsive element' },
    { name: 'TCCC-motif',  regex: /TCTCCCT/g,           category: 'Light Responsive', desc: 'Light responsive element' },
    { name: 'ATC-motif',   regex: /AATC[AT]CC/g,        category: 'Light Responsive', desc: 'Part of a conserved DNA module involved in light responsiveness' },
    { name: 'MRE',         regex: /AACCTAA/g,           category: 'Light Responsive', desc: 'MYB binding site involved in light responsiveness' },

    // Phytohormone Responsive Elements
    { name: 'ABRE',        regex: /ACGTG[GT]C/g,        category: 'Hormone Responsive', desc: 'Cis-acting element involved in the abscisic acid (ABA) responsiveness' },
    { name: 'CGTCA-motif', regex: /CGTCA/g,             category: 'Hormone Responsive', desc: 'Cis-acting regulatory element involved in the MeJA-responsiveness' },
    { name: 'TGACG-motif', regex: /TGACG/g,             category: 'Hormone Responsive', desc: 'Cis-acting regulatory element involved in the MeJA-responsiveness' },
    { name: 'TCA-element', regex: /CCATCTTTTT/g,        category: 'Hormone Responsive', desc: 'Cis-acting element involved in salicylic acid responsiveness' },
    { name: 'AuxRR-core',  regex: /GGTCCAT/g,           category: 'Hormone Responsive', desc: 'Cis-acting regulatory element involved in auxin responsiveness' },
    { name: 'TGA-element', regex: /AACGAC/g,            category: 'Hormone Responsive', desc: 'Auxin-responsive element' },
    { name: 'GARE-motif',  regex: /TCTGTTG/g,           category: 'Hormone Responsive', desc: 'Gibberellin-responsive element' },
    { name: 'P-box',       regex: /CCTTTTG/g,           category: 'Hormone Responsive', desc: 'Gibberellin-responsive element' },
    { name: 'TATC-box',    regex: /TATCCCA/g,           category: 'Hormone Responsive', desc: 'Cis-acting element involved in gibberellin-responsiveness' },
    { name: 'ERE',         regex: /ATTTCAAA/g,          category: 'Hormone Responsive', desc: 'Ethylene-responsive element' },

    // Environmental Stress Responsive Elements
    { name: 'ARE',         regex: /AAACCA/g,            category: 'Stress Responsive', desc: 'Cis-acting regulatory element essential for the anaerobic induction' },
    { name: 'MBS',         regex: /CAACTG/g,            category: 'Stress Responsive', desc: 'MYB binding site involved in drought-inducibility' },
    { name: 'LTR',         regex: /CCGAAA/g,            category: 'Stress Responsive', desc: 'Cis-acting element involved in low-temperature responsiveness' },
    { name: 'TC-rich repeats', regex: /GTTTTCTTAC/g,    category: 'Stress Responsive', desc: 'Cis-acting element involved in defense and stress responsiveness' },
    { name: 'HSE',         regex: /AGAANNATTCT/g,       category: 'Stress Responsive', desc: 'Cis-acting element involved in heat stress responsiveness' },
    { name: 'GC-motif',    regex: /CCCCCG/g,            category: 'Stress Responsive', desc: 'Enhancer-like element involved in anoxic specific inducibility' },
    { name: 'WUN-motif',   regex: /AAATTTCCT/g,         category: 'Stress Responsive', desc: 'Wound-responsive element' },
    { name: 'DRE',         regex: /TACCGACAT/g,         category: 'Stress Responsive', desc: 'Dehydration, low-temp, salt stress responsive' },
    { name: 'STRE',        regex: /AGGGG/g,             category: 'Stress Responsive', desc: 'Stress response element (heat, osmotic, oxidative)' },

    // Development & Tissue Specific Elements
    { name: 'CAT-box',     regex: /GCCACT/g,            category: 'Development', desc: 'Cis-acting element related to meristem expression' },
    { name: 'CCGTCC-box',  regex: /CCGTCC/g,            category: 'Development', desc: 'Meristem specific activation' },
    { name: 'O2-site',     regex: /GATGACATG/g,         category: 'Development', desc: 'Cis-acting regulatory element involved in zein metabolism regulation' },
    { name: 'Skn-1_motif', regex: /GTCAT/g,             category: 'Development', desc: 'Cis-acting regulatory element required for endosperm expression' },
    { name: 'GCN4_motif',  regex: /TGTGTCA/g,           category: 'Development', desc: 'Cis-regulatory element involved in endosperm expression' },
    { name: 'RY-element',  regex: /CATGCATG/g,          category: 'Development', desc: 'Cis-acting regulatory element involved in seed-specific regulation' },
    { name: 'HD-Zip 1',    regex: /CAAT[AT]ATTG/g,      category: 'Development', desc: 'Involved in differentiation of the palisade mesophyll cells' },
    { name: 'MSA-like',    regex: /T[CT]CAA[AC]AA/g,    category: 'Development', desc: 'Cis-acting element involved in cell cycle regulation' },
    { name: 'TGA-box',     regex: /TGACG/g,             category: 'Development', desc: 'Part of a cis-acting element involved in seed-specific regulation' },

    // Other / Binding Sites
    { name: 'W-box',       regex: /TTGAC[CT]/g,         category: 'Binding Site', desc: 'WRKY binding site involved in defense response' },
    { name: 'NON-box',     regex: /AGAT[AT]TAGA/g,      category: 'Binding Site', desc: 'Cis-acting regulatory element related to meristem specific activation' },
    { name: 'CE3',         regex: /CACGCG/g,            category: 'Binding Site', desc: 'Cis-acting element involved in ABA and VP1 responsiveness' },
    { name: 'MYB-core',    regex: /CAGGTTA/g,           category: 'Binding Site', desc: 'MYB transcription factor binding site' },
    { name: 'MYC-core',    regex: /CACGTG/g,            category: 'Binding Site', desc: 'bHLH / MYC transcription factor binding site' },
    
    // Circadian Control
    { name: 'circadian',   regex: /CAANNNNATC/g,        category: 'Circadian', desc: 'Cis-acting regulatory element involved in circadian control' },

    // General Animal/Bacterial (for broad toolkit use)
    { name: 'NF-κB',       regex: /GGGRNNYYCC/g,        category: 'Animal', desc: 'NF-κB transcription factor binding site' },
    { name: 'AP-1',        regex: /TGANTCA/g,           category: 'Animal', desc: 'Activator protein 1 binding site' },
    { name: 'E-box',       regex: /CANNTG/g,            category: 'Animal', desc: 'Enhancer box (bHLH binding)' },
    { name: 'CREB',        regex: /TGACGTCA/g,          category: 'Animal', desc: 'cAMP response element' },
    { name: 'Pribnow Box', regex: /TATAAT/g,            category: 'Bacterial', desc: 'Bacterial -10 promoter element' },
    { name: '-35 Element', regex: /TTGACA/g,            category: 'Bacterial', desc: 'Bacterial -35 sigma factor binding site' }
  ];

  // Helper mapping categories to CSS colors for visual consistency in UI
  const CATEGORY_COLORS = {
    'Core': 'motif-red',
    'Light Responsive': 'motif-yellow',
    'Hormone Responsive': 'motif-amber',
    'Stress Responsive': 'motif-sky',
    'Development': 'motif-green',
    'Binding Site': 'motif-violet',
    'Circadian': 'motif-blue',
    'Animal': 'motif-slate',
    'Bacterial': 'motif-slate'
  };

  // Assign CSS classes automatically based on category
  CIS_ELEMENTS.forEach(el => {
    el.cssClass = CATEGORY_COLORS[el.category] || 'motif-slate';
  });

  // ── ORF Finder ────────────────────────────────────────
  function getReverseComplement(seq) {
    const comp = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N' };
    return seq.split('').reverse().map(b => comp[b] || 'N').join('');
  }

  function findORFs(seq, minLength = 30) {
    const orfs = [];
    const BioMath = window.BioKit.core.BioMath;
    const rcSeq = getReverseComplement(seq);

    function scanStrand(strandSeq, isRev) {
      for (let frame = 0; frame < 3; frame++) {
        let inORF = false;
        let orfStart = -1;

        for (let i = frame; i <= strandSeq.length - 3; i += 3) {
          const codon = strandSeq.substring(i, i + 3);

          if (!inORF && codon === 'ATG') {
            inORF = true;
            orfStart = i;
          }

          if (inORF && (codon === 'TAA' || codon === 'TAG' || codon === 'TGA')) {
            const length = i + 3 - orfStart;
            if (length >= minLength) {
              const orfSeq = strandSeq.substring(orfStart, i + 3);
              const protein = BioMath ? BioMath.translateDNA(orfSeq) : '';
              
              let actualStart, actualEnd;
              if (isRev) {
                // Reverse strand logic: 
                // orfStart is 0-indexed on the rcSeq.
                // It ends at i + 2 on rcSeq.
                // Map to 1-indexed forward strand:
                actualEnd = seq.length - orfStart;
                actualStart = seq.length - (i + 2);
              } else {
                actualStart = orfStart + 1;
                actualEnd = i + 3;
              }

              orfs.push({
                frame: isRev ? -(frame + 1) : (frame + 1),
                start: actualStart,
                end: actualEnd,
                length: length,
                protein: protein,
                proteinLength: Math.floor(length / 3)
              });
            }
            inORF = false;
          }
        }

        // Handle unterminated ORF at end
        if (inORF) {
          const length = strandSeq.length - orfStart;
          if (length >= minLength) {
            const orfSeq = strandSeq.substring(orfStart);
            const protein = BioMath ? BioMath.translateDNA(orfSeq) : '';
            
            let actualStart, actualEnd;
            if (isRev) {
              actualEnd = seq.length - orfStart;
              actualStart = 1;
            } else {
              actualStart = orfStart + 1;
              actualEnd = seq.length;
            }

            orfs.push({
              frame: isRev ? -(frame + 1) : (frame + 1),
              start: actualStart,
              end: actualEnd,
              length: length,
              protein: protein,
              proteinLength: Math.floor(length / 3),
              partial: true
            });
          }
        }
      }
    }

    scanStrand(seq, false); // +1, +2, +3
    scanStrand(rcSeq, true); // -1, -2, -3

    return orfs.sort((a, b) => b.length - a.length);
  }

  // ── Cis-Element Scanner ───────────────────────────────
  function scanCisElements(seq) {
    const found = [];
    const BioMath = window.BioKit.core.BioMath;

    // Use PWM engine for high-precision motifs
    if (BioMath && BioMath.scoreSequenceWithPWM && BioMath.PWM_DATABASE) {
      Object.entries(BioMath.PWM_DATABASE).forEach(([name, pwm]) => {
        const hits = BioMath.scoreSequenceWithPWM(seq, pwm, 0.8);
        const elInfo = CIS_ELEMENTS.find(e => e.name === name);
        hits.forEach(hit => {
          found.push({
            name: name,
            position: hit.pos,
            matched: hit.matched,
            cssClass: elInfo ? elInfo.cssClass : 'motif-slate',
            desc: elInfo ? elInfo.desc : 'PWM detected consensus',
            score: (hit.score * 100).toFixed(1) + '%'
          });
        });
      });
    }

    // Fallback/Legacy Regex scanning for elements without PWMs
    CIS_ELEMENTS.forEach(el => {
      // Don't duplicate if already handled by PWM
      if (BioMath && BioMath.PWM_DATABASE && BioMath.PWM_DATABASE[el.name]) return;

      const re = new RegExp(el.regex.source, 'g');
      let match;
      while ((match = re.exec(seq)) !== null) {
        found.push({
          name: el.name,
          position: match.index + 1,
          matched: match[0],
          cssClass: el.cssClass,
          desc: el.desc,
          score: '100% (Exact)'
        });
      }
    });

    return found.sort((a, b) => a.position - b.position);
  }

  // ── Main Run ──────────────────────────────────────────
  function run() {
    const input = document.getElementById('proInput');
    if (!input) return;

    const raw = input.value;
    const validation = window.validateSequence(raw, 'dna');

    if (!validation.valid) {
      window.showValidationWarning('panel-promoter', validation);
      return;
    }

    const seq = validation.clean;
    const upstreamLen = parseInt(document.getElementById('proUpstream').value) || 1000;

    window.withLoading('panel-promoter', () => {
      const motifs = scanCisElements(seq);
      const orfs = findORFs(seq);
      displayResults(seq, motifs, orfs, upstreamLen);
    }, seq.length);
  }

  // ── Display ───────────────────────────────────────────
  function displayResults(seq, motifs, orfs, upstreamLen) {
    const container = document.getElementById('proResults');
    if (!container) return;
    container.classList.remove('hidden');

    // Identify Anchor ORF (Longest)
    const anchor = orfs.length > 0 ? orfs[0] : null;

    // Stats
    window.buildStatCards('proStatRow', [
      { val: seq.length.toLocaleString(), label: 'Length (bp)' },
      { val: motifs.length, label: 'Cis-Elements' },
      { val: orfs.length, label: 'ORFs Found' },
      { val: anchor ? anchor.proteinLength + ' aa' : '—', label: 'Largest ORF' },
    ]);

    // Populate CDS/ORF Regions Table
    renderORFTable(orfs);

    // Populate Cis-Elements Table & Legend
    renderMotifTable(motifs);
    renderLegend(motifs);

    // Handle Primary Anchor Extraction (CDS, Promoter, Protein)
    if (anchor) {
      // 1. Extract CDS
      const cdsSeq = seq.substring(anchor.start - 1, anchor.end);
      document.getElementById('proCDSSeqRange').textContent = `CDS [${anchor.start} ... ${anchor.end}] · ${anchor.length} bp · Frame ${anchor.frame > 0 ? '+' : ''}${anchor.frame}`;
      renderPlainSequence(cdsSeq, 'proCDSSeqViewer');

      // 2. Extract Protein
      document.getElementById('proProteinSeqRange').textContent = `Translation Product · ${anchor.proteinLength} aa`;
      renderPlainSequence(anchor.protein, 'proProteinSeqViewer', true);

      // 3. Extract Promoter (Upstream)
      const proStart = Math.max(1, anchor.start - upstreamLen);
      const proEnd = anchor.start - 1;
      const promoterSeq = seq.substring(proStart - 1, proEnd);
      
      const proMotifs = motifs.filter(m => m.position >= proStart && m.position <= proEnd)
                             .map(m => ({ ...m, position: m.position - (proStart - 1) }));

      document.getElementById('proSeqRange').textContent = `Promoter Region [${proStart} ... ${proEnd}] · ${promoterSeq.length} bp`;
      renderHighlightedSequence(promoterSeq, proMotifs, 'proSeqViewer');
    } else {
      ['proCDSSeqViewer', 'proProteinSeqViewer', 'proSeqViewer'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.innerHTML = '<p style="color:var(--text-muted)">No suitable ORF found to anchor sequence extraction.</p>';
      });
      document.querySelectorAll('.pro-seq-meta').forEach(el => el.textContent = '—');
    }

    // Feed genome viewer
    if (window.loadGenomeData) {
      const cdsData = orfs.map(o => ({
        start: o.start - 1, end: o.end, frame: o.frame, length: o.length
      }));
      const motifData = motifs.map(m => ({
        position: m.position - 1, matched: m.matched, name: m.name
      }));
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

    // Build highlight map
    const hlMap = new Array(seq.length).fill(null);
    motifs.forEach(m => {
      for (let i = m.position - 1; i < m.position - 1 + m.matched.length && i < seq.length; i++) {
        if (!hlMap[i]) hlMap[i] = m.cssClass;
      }
    });

    let html = '';
    for (let i = 0; i < seq.length; i++) {
      const c = window.escapeHTML(seq[i]);
      if (hlMap[i]) {
        html += `<span class="motif-highlight ${hlMap[i]}">${c}</span>`;
      } else {
        html += c;
      }
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

  function renderMotifTable(motifs) {
    const table = document.getElementById('proCisTable');
    if (!table) return;

    if (motifs.length === 0) {
      table.innerHTML = '<p style="color:var(--text-muted)">No cis-elements found.</p>';
      return;
    }

    let html = `<table class="data-table"><thead><tr>
      <th>Element</th><th>Position</th><th>Sequence</th><th>Confidence/Score</th><th>Description</th>
    </tr></thead><tbody>`;
    motifs.forEach(m => {
      html += `<tr>
        <td><span class="motif-badge ${m.cssClass}">${window.escapeHTML(m.name)}</span></td>
        <td>${m.position}</td>
        <td style="font-family:var(--mono);color:var(--teal)">${window.escapeHTML(m.matched)}</td>
        <td style="font-weight: 500;">${window.escapeHTML(m.score || '—')}</td>
        <td style="color:var(--text-muted)">${window.escapeHTML(m.desc)}</td>
      </tr>`;
    });
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
      <th>Frame</th><th>Start</th><th>End</th><th>Length (bp)</th><th>Protein (aa)</th><th>Note</th>
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

  // ── Events ────────────────────────────────────────────
  document.getElementById('proAnalyzeBtn')?.addEventListener('click', window.debounce(run, 300));

  document.getElementById('proSampleBtn')?.addEventListener('click', () => {
    const input = document.getElementById('proInput');
    if (input) {
      // E. coli lacZ promoter region + start of gene
      input.value = 'GCTTATTATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATTATAAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTAACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGTTGGAGTGACGGCAGTTATCTGGAAGATCAGG';
      input.dispatchEvent(new Event('input'));
    }
  });

  document.getElementById('proClearBtn')?.addEventListener('click', () => {
    const input = document.getElementById('proInput');
    if (input) input.value = '';
    document.getElementById('proResults')?.classList.add('hidden');
    input.dispatchEvent(new Event('input'));
  });

  document.getElementById('proExportBtn')?.addEventListener('click', () => {
    const input = document.getElementById('proInput');
    if (!input || !input.value) return;
    
    const seq = window.cleanSeq(input.value);
    const orfs = findORFs(seq);
    const anchor = orfs.length > 0 ? orfs[0] : null;
    const upstream = parseInt(document.getElementById('proUpstream').value) || 1000;

    let text = `BioToolkit Promoter/CDS Analysis Report\n`;
    text += `Generated: ${new Date().toLocaleString()}\n`;
    text += `Source Sequence Length: ${seq.length} bp\n`;
    text += `----------------------------------------\n\n`;

    if (anchor) {
      const cds = seq.substring(anchor.start - 1, anchor.end);
      const proStart = Math.max(1, anchor.start - upstream);
      const promoter = seq.substring(proStart - 1, anchor.start - 1);

      text += `>Primary_Anchor_CDS [Range: ${anchor.start}-${anchor.end}] [Frame: ${anchor.frame}]\n${cds}\n\n`;
      text += `>Primary_Anchor_Translation [${anchor.proteinLength} aa]\n${anchor.protein}\n\n`;
      text += `>Promoter_Region [Range: ${proStart}-${anchor.start - 1}] [Upstream: ${upstream}bp]\n${promoter}\n\n`;
    } else {
      text += `No Anchor ORF found in query sequence.\n`;
    }

    window.BioKit.utils.downloadText(`biokit_analysis_${Date.now()}.txt`, text);
  });

})();
