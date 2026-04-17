/* ============================================================
   science-algorithms.js — Scientific Core Engine for BioKit
   ============================================================
   This file implements scientifically validated algorithms.
   References:
   1. SantaLucia Jr, J. (1998). PNAS, 95(4), 1460-1465. (Nearest-Neighbor Tm)
   2. Kyte, J., & Doolittle, R. F. (1982). JMB, 157(1), 105-132. (Hydrophobicity)
   3. NCBI Standard Genetic Code (Transl_table 1).
*/

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };
  window.BioKit.core.BioMath = window.BioKit.core.BioMath || {};
  const BioMath = window.BioKit.core.BioMath;

  // ==========================================================
  // ASSAY PRESETS (PROMPT 4)
  // ==========================================================
  BioMath.ASSAY_PRESETS = {
    "pcr": {
      name: "Standard PCR",
      length: { min: 18, max: 25 },
      tm: { min: 58, max: 65 },
      tmDifference: 5,
      gc: { min: 40, max: 60 },
      gc3prime: true,
      homopolymerLimit: 5,
      productSize: { min: 200, max: 5000 },
      dimerDeltaGLimit: -6.0
    },
    "qpcr": {
      name: "qPCR/RT-PCR",
      length: { min: 18, max: 25 },
      tm: { min: 59, max: 61 },
      tmDifference: 2,
      gc: { min: 45, max: 55 },
      gc3prime: true,
      homopolymerLimit: 4,
      productSize: { min: 50, max: 500 },
      dimerDeltaGLimit: -5.0
    },
    "cloning": {
      name: "Cloning/Sequencing",
      length: { min: 20, max: 30 },
      tm: { min: 60, max: 70 },
      tmDifference: 5,
      gc: { min: 35, max: 65 },
      gc3prime: true,
      homopolymerLimit: 5,
      productSize: { min: 500, max: 10000 },
      dimerDeltaGLimit: -7.0
    },
    "multiplex": {
      name: "Multiplex PCR",
      length: { min: 19, max: 24 },
      tm: { min: 60, max: 65 },
      tmDifference: 1,
      gc: { min: 45, max: 55 },
      gc3prime: true,
      homopolymerLimit: 4,
      productSize: { min: 100, max: 1000 },
      dimerDeltaGLimit: -4.0
    }
  };

  // ==========================================================
  // 0. PWM (Position Weight Matrix) Engine
  // ==========================================================
  BioMath.scoreSequenceWithPWM = function(seq, pwm, threshold = 0.8) {
      const L = pwm.length;
      const res = [];
      const bases = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 };

      for (let i = 0; i <= seq.length - L; i++) {
          let score = 0;
          let minScore = 0;
          let maxScore = 0;

          for (let j = 0; j < L; j++) {
              const base = seq[i + j];
              const bIdx = bases[base];
              
              // Calculate specific score relative to min/max possible for normalization
              const posWeights = pwm[j];
              const currentWeight = bIdx !== undefined ? posWeights[bIdx] : -10; // Penalty for N
              score += currentWeight;
              
              minScore += Math.min(...posWeights);
              maxScore += Math.max(...posWeights);
          }

          // Normalize score to 0.0 - 1.0 (RelScore)
          const normalized = (score - minScore) / (maxScore - minScore);
          if (normalized >= threshold) {
              res.push({ pos: i + 1, score: normalized, matched: seq.substring(i, i + L) });
          }
      }
      return res;
  };

  // Common TFBS PWMs (Log-Odds format)
  // Rows: A, C, G, T
  BioMath.PWM_DATABASE = {
      'TATA-box': [
          [ -10, -10, -10,  10 ], // T
          [  10, -10, -10, -10 ], // A
          [ -10, -10, -10,  10 ], // T
          [  10, -10, -10, -10 ], // A
          [  10, -10, -10,  10 ], // A/T
          [  10, -10, -10, -10 ], // A
          [  10, -10, -10,  10 ], // A/T
          [  10, -10, -10, -10 ]  // A
      ],
      'G-box': [
          [ -10,  10, -10, -10 ], // C
          [  10, -10, -10, -10 ], // A
          [ -10,  10, -10, -10 ], // C
          [ -10, -10,  10, -10 ], // G
          [ -10, -10, -10,  10 ], // T
          [ -10, -10,  10, -10 ]  // G
      ],
      'ABRE': [
          [  10, -10, -10, -10 ], // A
          [ -10,  10, -10, -10 ], // C
          [ -10, -10,  10, -10 ], // G
          [ -10, -10, -10,  10 ], // T
          [ -10, -10,  10, -10 ], // G
          [ -10, -10,  10,  10 ]  // G/T
      ]
  };
  
  BioMath.reverseComplement = function(seq) {
      if (typeof seq !== 'string' || !seq) return '';
      const comp = { 
          'A': 'T', 'T': 'A', 'U': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
          'a': 't', 't': 'a', 'u': 'a', 'c': 'g', 'g': 'c', 'n': 'n',
          'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
          'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
          'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm', 'm': 'k',
          'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd'
      };
      return seq.split('').reverse().map(b => comp[b] !== undefined ? comp[b] : b).join('');
  };

  // ==========================================================
  // 1. Accurate GC Calculation
  // ==========================================================
  BioMath.calculateGC = function(sequence) {
      if (!sequence) return 0;
      const s = sequence.toUpperCase().replace(/[^ATGCNRYMKSWHBVD]/g, '');
      if (s.length === 0) return 0;
      
      let gcSum = 0;
      const gcContrib = {
          'G': 1, 'C': 1, 'S': 1,
          'A': 0, 'T': 0, 'W': 0,
          'R': 0.5, 'Y': 0.5, 'K': 0.5, 'M': 0.5, 'N': 0.5,
          'B': 0.6667, 'V': 0.6667,
          'D': 0.3333, 'H': 0.3333
      };
      
      for (let i = 0; i < s.length; i++) {
          gcSum += gcContrib[s[i]];
      }
      return ((gcSum / s.length) * 100).toFixed(1);
  };

  /**
   * Validates primer length against assay constraints.
   * @param {string} primer - The nucleotide sequence.
   * @param {Object} constraints - The active assay constraints.
   * @returns {boolean}
   */
  BioMath.validateLength = function(primer, constraints) {
    if (!primer || !constraints.length) return false;
    const len = primer.length;
    return len >= constraints.length.min && len <= constraints.length.max;
  };

  /**
   * Validates primer Tm against assay constraints.
   * @param {number|string} tm - The computed melting temperature.
   * @param {Object} constraints - The active assay constraints.
   * @returns {boolean}
   */
  BioMath.validateTm = function(tm, constraints) {
    if (tm === undefined || !constraints.tm) return false;
    const val = parseFloat(tm);
    return val >= constraints.tm.min && val <= constraints.tm.max;
  };

  /**
   * Validates primer GC content against assay constraints.
   * @param {number|string} gc - The GC percentage.
   * @param {Object} constraints - The active assay constraints.
   * @returns {boolean}
   */
  BioMath.validateGC = function(gc, constraints) {
    if (gc === undefined || !constraints.gc) return false;
    const val = parseFloat(gc);
    return val >= constraints.gc.min && val <= constraints.gc.max;
  };

  /**
   * Validates 3' GC clamp existence if required by constraints.
   * @param {string} primer - The nucleotide sequence.
   * @param {Object} constraints - The active assay constraints.
   * @returns {boolean}
   */
  BioMath.validateGCClamp = function(primer, constraints) {
    if (!constraints.gc3prime || !primer) return true;
    const last3 = primer.slice(-3).toUpperCase();
    return last3.includes('G') || last3.includes('C');
  };

  /**
   * Detect runs of identical consecutive bases.
   * @param {string} primer - DNA sequence
   * @param {number} maxAllowed - Max consecutive identical bases allowed
   * @returns {object} { hasViolation, violations, summary }
   */
  BioMath.detectHomopolymerRuns = function(primer, maxAllowed) {
    const result = { hasViolation: false, violations: [], summary: 'None' };
    if (!primer) return result;

    const s = primer.toUpperCase();
    let currentRun = 1;

    for (let i = 1; i <= s.length; i++) {
      if (i < s.length && s[i] === s[i-1]) {
        currentRun++;
      } else {
        if (currentRun > maxAllowed) {
          result.hasViolation = true;
          result.violations.push({
            base: s[i-1],
            position: i - currentRun + 1,
            length: currentRun
          });
        }
        currentRun = 1;
      }
    }

    if (result.hasViolation) {
      result.summary = result.violations.map(v => `${v.base.repeat(v.length)} (${v.length} bp at pos ${v.position})`).join(', ');
    }
    return result;
  };

  /**
   * Validates that homopolymer runs do not exceed constraints.
   * @param {string} primer - The nucleotide sequence.
   * @param {Object} constraints - The active assay constraints.
   * @returns {boolean}
   */
  BioMath.validateHomopolymer = function(primer, constraints) {
    const limit = constraints.homopolymerLimit || 5;
    const res = BioMath.detectHomopolymerRuns(primer, limit);
    return !res.hasViolation;
  };

  /**
   * Validates the Tm difference between a primer pair.
   * @param {number|string} tm1 - First primer Tm.
   * @param {number|string} tm2 - Second primer Tm.
   * @param {Object} constraints - The active assay constraints.
   * @returns {boolean}
   */
  BioMath.validateTmDifference = function(tm1, tm2, constraints) {
    if (!constraints.tmDifference) return true;
    const diff = Math.abs(parseFloat(tm1) - parseFloat(tm2));
    return diff <= constraints.tmDifference;
  };

  /**
   * Scans a template for primer binding sites to determine local specificity.
   * Uses mismatch-aware scanning: counts sites where the primer aligns with
   * at most `maxMismatch` mismatches over its full length. This catches
   * realistic mispriming events that pure indexOf would miss.
   * @param {string} primer - The primer sequence.
   * @param {string} template - The target DNA sequence.
   * @param {boolean} isReverse - Whether to search for the reverse complement.
   * @param {number} maxMismatch - Maximum allowed mismatches for a "binding" hit (default 2).
   * @returns {Object} { primarySite, specificity }
   */
  BioMath.checkTemplateSpecificity = function(primer, template, isReverse = false, maxMismatch = 2) {
    const res = {
      primarySite: { position: -1, count: 0 },
      specificity: "UNIQUE"
    };
    if (!primer || !template) return res;

    const p = isReverse ? BioMath.reverseComplement(primer).toUpperCase() : primer.toUpperCase();
    const t = template.toUpperCase();
    const pLen = p.length;

    if (pLen === 0 || t.length < pLen) return res;

    // Slide the primer across the template and count mismatches at each position
    for (let i = 0; i <= t.length - pLen; i++) {
      let mismatches = 0;
      for (let j = 0; j < pLen; j++) {
        if (t[i + j] !== p[j]) {
          mismatches++;
          if (mismatches > maxMismatch) break; // Early exit for this window
        }
      }
      if (mismatches <= maxMismatch) {
        res.primarySite.count++;
        if (res.primarySite.position === -1) res.primarySite.position = i + 1;
      }
    }

    if (res.primarySite.count > 1) res.specificity = "MULTIPLE";
    else if (res.primarySite.count === 0) res.specificity = "NOT_FOUND";

    return res;
  };

  /**
   * Validates if primers or amplicons span exon-exon junctions or introns.
   * @param {number} fwdStart - Forward 1-based start
   * @param {number} fwdEnd - Forward 1-based end
   * @param {number} revStart - Reverse 1-based start
   * @param {number} revEnd - Reverse 1-based end
   * @param {Array} features - Feature list from template metadata.
   * @param {Object} prefs - { intronSpanningRequired, avoidExonJunction }
   * @returns {Object} { fwdSpans, revSpans, ampliconSpans, verdict, reason }
   */
  BioMath.checkExonSpanning = function(fwdStart, fwdEnd, revStart, revEnd, features, prefs) {
    const res = {
      fwdSpans: false, revSpans: false, ampliconSpans: false, verdict: "PASS", reason: ""
    };
    if (!features || features.length === 0) return res;

    const exons = features.filter(f => f.type.toLowerCase() === 'exon').sort((a,b) => a.start - b.start);
    if (exons.length < 2) return res;

    const checkSpans = (s, e) => {
      for (let i = 0; i < exons.length - 1; i++) {
        const boundary = exons[i].end;
        if (s <= boundary && e > boundary) return true;
      }
      return false;
    };

    res.fwdSpans = checkSpans(fwdStart, fwdEnd);
    res.revSpans = checkSpans(revStart, revEnd);

    const getExonIdx = (p) => exons.findIndex(e => p >= e.start && p <= e.end);
    const fEx = getExonIdx(fwdStart);
    const rEx = getExonIdx(revStart);
    if (fEx !== -1 && rEx !== -1 && fEx !== rEx) res.ampliconSpans = true;

    if (prefs.intronSpanningRequired) {
      if (!res.fwdSpans && !res.revSpans && !res.ampliconSpans) {
        res.verdict = "FAIL";
        res.reason = "Primers do not span an exon-junction or intron (qPCR optimization required).";
      }
    } else if (prefs.avoidExonJunction) {
      if (res.fwdSpans || res.revSpans) {
        res.verdict = "FAIL";
        res.reason = "Primer overlaps an exon-junction (Avoid Junction selected).";
      }
    }

    return res;
  };
   
  /**
   * Calculates amplicon information (coordinates, size, sequence) from template and primers.
   *

  /**
   * Calculates amplicon information (coordinates, size, sequence) from template and primers.
   * @param {string} template - The target DNA sequence.
   * @param {string} fwdSeq - Forward primer sequence (5'->3').
   * @param {string} revSeq - Reverse primer sequence (5'->3').
   * @param {Object} constraints - Active assay constraints for size validation.
   * @returns {Object} { amplicon: { sequence, size, start, end }, warnings: [] }
   */
  BioMath.calculateAmpliconInfo = function(template, fwdSeq, revSeq, constraints = {}) {
    const res = {
      amplicon: { sequence: '', size: 0, start: 0, end: 0 },
      warnings: []
    };

    if (!template || !fwdSeq || !revSeq) return res;

    const t = template.toUpperCase();
    const f = fwdSeq.toUpperCase();
    const r = revSeq.toUpperCase();
    const r_rc = BioMath.reverseComplement(r).toUpperCase();

    // 1. Find Primer Sites (Specificity Check)
    const fIdx = [];
    let posF = t.indexOf(f);
    while (posF !== -1) { fIdx.push(posF); posF = t.indexOf(f, posF + 1); }

    const rIdx = [];
    let posR = t.indexOf(r_rc);
    while (posR !== -1) { rIdx.push(posR); posR = t.indexOf(r_rc, posR + 1); }

    if (fIdx.length === 0) res.warnings.push("Forward primer binding site not found in template.");
    if (rIdx.length === 0) res.warnings.push("Reverse primer binding site not found in template.");

    if (fIdx.length > 1) res.warnings.push("FLAG (Multiple binding sites): Forward primer binds at multiple positions.");
    if (rIdx.length > 1) res.warnings.push("FLAG (Multiple binding sites): Reverse primer binds at multiple positions.");

    if (fIdx.length === 0 || rIdx.length === 0) return res;

    // Use the first Fwd and last Rev to determine the maximum span for diagnostic reporting
    const fwdIdx = fIdx[0];
    const revIdx = rIdx[rIdx.length - 1];

    // 2. Map Boundaries
    const ampliconStart = fwdIdx;
    const ampliconEnd = revIdx + r.length;
    const size = ampliconEnd - ampliconStart;

    res.amplicon.start = ampliconStart + 1; // 1-based
    res.amplicon.end = ampliconEnd;         // 1-based
    res.amplicon.size = size;

    if (size <= 0) {
      res.warnings.push(`Primers oriented incorrectly or overlap significantly (computed size: ${size}bp).`);
      return res;
    }

    res.amplicon.sequence = t.substring(ampliconStart, ampliconEnd);

    // 3. Constraint Validation
    return res;
  };

  /**
   * Calculates the thermodynamic stability (Delta G) of a dimerization event.
   * Based on SantaLucia NN parameters at a given temperature.
   * @param {string} seq1 - First primer sequence.
   * @param {string} seq2 - Second primer sequence (or same for self-dimer).
   * @param {number} temp_C - Evaluation temperature in Celsius (default 37.0).
   * @returns {Object} { deltaG, maxRun, is3PrimeTerminal }
   */
  BioMath.calculateDimerThermodynamics = function (seq1, seq2, temp_C = 37.0) {
    const s1 = seq1.toUpperCase();
    const s2 = seq2.toUpperCase();
    const s2_rc = BioMath.reverseComplement(s2).toUpperCase();
    const temp_K = temp_C + 273.15;
    const R = 1.9872;

    let minDeltaG = 0;
    let bestRun = 0;
    let is3Prime = false;

    // Use a sliding alignment to find the most stable binding site
    for (let shift = -s1.length + 1; shift < s2_rc.length; shift++) {
      let dH = 0, dS = 0, currentRun = 0, maxRunInShift = 0;
      let shiftIs3Prime = false;

      // Check complementarity at this shift
      for (let i = 0; i < s1.length; i++) {
        const j = i + shift;
        if (j >= 0 && j < s2_rc.length) {
          if (s1[i] === s2_rc[j]) {
            currentRun++;
            maxRunInShift = Math.max(maxRunInShift, currentRun);
            
            // If we have at least a dimer pair (2 bases), add NN params
            if (currentRun >= 2) {
              const pair = s1[i-1] + s1[i];
              if (NN_PARAMS[pair]) {
                dH += NN_PARAMS[pair][0];
                dS += NN_PARAMS[pair][1];
              }
            }
            
            // Detect if this match involves the 3' end of either primer
            // Primer 1 3' is at indexed s1.length-1
            // Primer 2 3' is at s2_rc index 0 (if seq2=...5' then s2_rc=3'...5'?) 
            // Wait, reverseComplement of 5'-GAT-3' is 3'-CTA-5'. So base index 0 of s2_rc is its 3' end.
            if (i === s1.length - 1 || j === 0) {
              shiftIs3Prime = true;
            }
          } else {
            currentRun = 0;
          }
        }
      }

      if (maxRunInShift >= 4) {
        // Add initiation parameters using the actual terminal bases of the
        // complementary run (SantaLucia 1998 convention: initiation depends
        // on the identity of the terminal base pair, not a fixed G).
        // Find the first and last matched bases in this shift to pick the
        // correct initiation parameters.
        let firstMatchBase = 'G', lastMatchBase = 'G';
        for (let k = 0; k < s1.length; k++) {
          const kj = k + shift;
          if (kj >= 0 && kj < s2_rc.length && s1[k] === s2_rc[kj]) {
            firstMatchBase = s1[k]; break;
          }
        }
        for (let k = s1.length - 1; k >= 0; k--) {
          const kj = k + shift;
          if (kj >= 0 && kj < s2_rc.length && s1[k] === s2_rc[kj]) {
            lastMatchBase = s1[k]; break;
          }
        }
        const initStart = NN_INIT[firstMatchBase] || NN_INIT['G'];
        const initEnd   = NN_INIT[lastMatchBase]  || NN_INIT['G'];
        dH += initStart.dH + initEnd.dH;
        dS += initStart.dS + initEnd.dS;

        const dG = (dH * 1000 - temp_K * dS) / 1000;
        if (dG < minDeltaG) {
          minDeltaG = dG;
          bestRun = maxRunInShift;
          is3Prime = shiftIs3Prime;
        }
      }
    }

    return { 
      deltaG: minDeltaG.toFixed(2), 
      maxRun: bestRun, 
      is3PrimeTerminal: is3Prime 
    };
  };

  /**
   * Assesses the stability and quality of the 3' end.
   * Identifies excessive GC sequences (over-clamping) or loose A/T ends.
   */
  BioMath.assess3PrimeStability = function(seq) {
    const s = seq.toUpperCase();
    const last5 = s.slice(-5);
    const gcCount = (last5.match(/[GC]/g) || []).length;
    
    let verdict = 'STABLE';
    let message = 'Optimal 3\' terminus complexity.';

    if (gcCount >= 4) {
      verdict = 'CRITICAL';
      message = 'Excessive 3\' GC content: High risk of mispriming.';
    } else if (gcCount <= 1) {
      verdict = 'WEAK';
      message = 'Terminal 3\' AT-richness may reduce priming efficiency.';
    }

    return { verdict, message, gcCount };
  };

  // ==========================================================
  // 2. SantaLucia Nearest-Neighbor Primer Tm Formula
  // ==========================================================
  /*
     Tm (melting temperature) calculation based on Nearest-Neighbor thermodynamics.
     SantaLucia, J. Jr. (1998) Proc Natl Acad Sci USA 95, 1460-1465.
     Tm = (deltaH * 1000) / (deltaS + R * ln(C/4)) - 273.15 + 16.6 * log10([Na+])
     Where:
        deltaH = Enthalpy (kcal/mol -> cal/mol)
        deltaS = Entropy (cal/mol*K)
        R = Universal Gas Constant (1.987 cal/K*mol)
        C = Primer Concentration
  */
   // Thermodynamic values for NN pairs: [deltaH (kcal/mol), deltaS (cal/K*mol)]
  const NN_PARAMS = {
    'AA': [-7.9, -22.2], 'TT': [-7.9, -22.2],
    'AT': [-7.2, -20.4], 'TA': [-7.2, -21.3],
    'CA': [-8.5, -22.7], 'TG': [-8.5, -22.7],
    'GT': [-8.4, -22.4], 'AC': [-8.4, -22.4],
    'CT': [-7.8, -21.0], 'AG': [-7.8, -21.0],
    'GA': [-8.2, -22.2], 'TC': [-8.2, -22.2],
    'CG': [-10.6, -27.2], 'GC': [-9.8, -24.4],
    'GG': [-8.0, -19.9], 'CC': [-8.0, -19.9]
  };

  const IUPAC_DNA = {
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
  };

  // RNA-specific IUPAC ambiguity map — uses U, never T (Fix 3B)
  const IUPAC_RNA = {
    'R': ['A', 'G'],
    'Y': ['C', 'U'],   // U not T
    'S': ['G', 'C'],
    'W': ['A', 'U'],   // U not T
    'K': ['G', 'U'],   // U not T
    'M': ['A', 'C'],
    'B': ['C', 'G', 'U'],  // U not T
    'D': ['A', 'G', 'U'],  // U not T
    'H': ['A', 'C', 'U'],  // U not T
    'V': ['A', 'C', 'G'],
    'N': ['A', 'U', 'G', 'C'],  // U not T
  };

  const NN_INIT = {
    'G': { dH: 0.1, dS: -2.8 }, 'C': { dH: 0.1, dS: -2.8 },
    'A': { dH: 2.3, dS: 4.1 }, 'T': { dH: 2.3, dS: 4.1 }
  };

  const RNA_NN_PARAMS = {
    // Source: Freier et al. (1986) PNAS 83:9373-9377, Table 2
    // Canonical Watson-Crick RNA nearest-neighbor parameters (1M NaCl)
    // Format: { dH: kcal/mol, dS: cal/mol·K }
    'AA': { dH: -6.6,  dS: -18.4 },
    'UU': { dH: -6.6,  dS: -18.4 },
    'AU': { dH: -5.7,  dS: -15.5 },
    'UA': { dH: -8.1,  dS: -22.6 },
    'CU': { dH: -8.6,  dS: -22.2 },
    'AG': { dH: -8.6,  dS: -22.2 },
    'CA': { dH: -10.5, dS: -27.8 },
    'UG': { dH: -10.5, dS: -27.8 },
    'CG': { dH: -8.0,  dS: -19.4 },
    'GC': { dH: -14.2, dS: -34.9 },
    'GG': { dH: -12.2, dS: -29.7 },
    'CC': { dH: -12.2, dS: -29.7 },
    'AC': { dH: -10.2, dS: -26.2 },
    'GU': { dH: -10.2, dS: -26.2 }
  };

  const RNA_INIT_PARAMS = {
    // Initiation parameters for RNA duplexes (Freier 1986)
    GC: { dH: 0.0,  dS: -10.1 },
    AU: { dH: 0.0,  dS: -10.1 } // In Freier 1986, initiation was a fixed -10.1 value
  };

  // Helix penalty if self-complementary (approximated)
  const SYM_PENALTY = { dH: 0, dS: -1.4 }; 

  BioMath.calculateTmNN = function (primerSeq, options = {}, sequenceType = 'DNA') {
    const {
      oligoConc_nM = 250,
      naConc_mM = 50,
      mgConc_mM = 1.5,
      dntpConc_mM = 0.8,
      dmso_pct = 0,
      formamide_m = 0
    } = options;

    // Normalize sequence correctly based on strand type
    const seq = (sequenceType === 'RNA')
      ? primerSeq.toUpperCase().replace(/T/g, 'U').replace(/[^AUGCYRSWKMBVHDN]/g, '')
      : primerSeq.toUpperCase().replace(/U/g, 'T').replace(/[^ATGCRYWSKMBDHVN]/g, '');

    const nnParams = (sequenceType === 'RNA') ? RNA_NN_PARAMS : NN_PARAMS;
    const initParams = (sequenceType === 'RNA') ? RNA_INIT_PARAMS : NN_INIT;
    if (seq.length < 2) return { tm: 0, confidence: 'N/A', warning: 'Sequence too short', isValid: false };

    // Buffer validation and Equivalent Na+ calculation (von Ahsen 2001 / Owczarzy 2008)
    // The Owczarzy (2008) and SantaLucia (1998) salt correction models are empirically 
    // validated for [Na+] between 50 mM and 1000 mM. Below this, phosphate-phosphate 
    // electrostatic repulsion destabilizes the duplex in ways simplistic log equations cannot resolve.
    const SALT_MODEL_MIN_MM = 50;
    
    // Prevent Math.log(0) crashes, but do not silently clamp the calculation bounds
    const safeNa = Math.max(1e-5, naConc_mM);
    let eqNa_mM = safeNa;
    
    // Apply Mg2+ equivalent Na+ calculation (von Ahsen 2001)
    if (mgConc_mM > dntpConc_mM) {
        eqNa_mM += 120 * Math.sqrt(mgConc_mM - dntpConc_mM);
    }
    const na_M = eqNa_mM / 1000;
    
    let modelWarning = null;
    if (eqNa_mM < SALT_MODEL_MIN_MM) {
        modelWarning = `Salt concentration (${eqNa_mM.toFixed(1)} mM eq. Na⁺) is below the validated range of the Owczarzy 2008 model (≥ ${SALT_MODEL_MIN_MM} mM). Result is an extrapolation — interpret with caution.`;
    }

    let dH = 0, dS = 0;

    // Resolve IUPAC ambiguity codes — selects the correct alphabet by strand type
    const ambiguityMap = (sequenceType === 'RNA') ? IUPAC_RNA : IUPAC_DNA;
    const resolveBase = (b) => ambiguityMap[b] || [b];
    const avgInit = (base) => {
      const variants = resolveBase(base);
      let th = 0, ts = 0;
      variants.forEach(v => {
        let p;
        if (sequenceType === 'RNA') {
            const rnaKey = (v === 'G' || v === 'C') ? 'GC' : 'AU';
            p = initParams[rnaKey] || initParams['GC'];
        } else {
            p = initParams[v] || initParams['G'];
        }
        th += p.dH; ts += p.dS;
      });
      return { dH: th / variants.length, dS: ts / variants.length };
    };

    const start = avgInit(seq[0]);
    const end = avgInit(seq[seq.length - 1]);
    dH += start.dH + end.dH;
    dS += start.dS + end.dS;

    // NN sums with IUPAC weighting
    for (let i = 0; i < seq.length - 1; i++) {
      const b1s = resolveBase(seq[i]);
      const b2s = resolveBase(seq[i + 1]);
      let pairH = 0, pairS = 0, count = 0;

      b1s.forEach(v1 => {
        b2s.forEach(v2 => {
          const pair = v1 + v2;
          if (nnParams[pair]) {
            pairH += nnParams[pair].dH !== undefined ? nnParams[pair].dH : nnParams[pair][0];
            pairS += nnParams[pair].dS !== undefined ? nnParams[pair].dS : nnParams[pair][1];
            count++;
          }
        });
      });
      if (count > 0) { dH += pairH / count; dS += pairS / count; }
    }

    dH *= 1000; // kcal to cal
    const R = 1.9872;
    const C = (oligoConc_nM * 1e-9) / 4;

    let Tm = (dH / (dS + R * Math.log(C))) - 273.15;
    
    // Physical Corrections
    Tm += (16.6 * Math.log10(na_M)); // Salt
    Tm -= (0.6 * dmso_pct);          // DMSO
    Tm -= (0.65 * formamide_m);     // Formamide

    // Range validity check (SantaLucia 1998)
    let warnings = [];
    let isValid = true;
    if (seq.length < 14 || seq.length > 50) {
      warnings.push(`CRITICAL: Sequence length (${seq.length} bp) outside validated NN model range (14-50 bp).`);
      isValid = false;
    }
    if (modelWarning) {
      warnings.push(modelWarning);
    } // High salt warning
    else if (eqNa_mM > 1000) {
      warnings.push('Salt concentration > 1000mM is outside validated range.');
    }

    return {
      tm: Math.max(0, Tm).toFixed(1),
      confidence: isValid && !modelWarning ? '±1.5°C' : 'UNCERTAIN',
      method: 'SantaLucia (1998) NN',
      warning: warnings.length > 0 ? warnings.join(' ') : null,
      isValid: isValid
    };
  };

  /**
   * Calculates Tm using the GC% formula (Schildkraut & Lifson, 1965).
   * Valid for long DNA sequences (>50 bp).
   * Includes salt corrections (von Ahsen, 2001).
   */
  BioMath.calculateTmGC = function (seq, options = {}) {
    const { naConc_mM = 50, mgConc_mM = 1.5, dntpConc_mM = 0.8 } = options;
    const cleanSeq = seq.toUpperCase().replace(/[^ATGC]/g, '');
    const len = cleanSeq.length;
    if (len === 0) return { tm: 0, confidence: 'N/A', isValid: false };

    const gc = parseFloat(BioMath.calculateGC(cleanSeq));
    
    // Equivalent Sodium calculation for laboratory-grade accuracy
    let eqNa_mM = naConc_mM;
    if (mgConc_mM > dntpConc_mM) {
      eqNa_mM += 120 * Math.sqrt(mgConc_mM - dntpConc_mM);
    }
    const na_M = Math.max(1e-5, eqNa_mM / 1000);

    // Schildkraut & Lifson (1965) formula
    // Tm = 81.5 + 16.6 * log10[Na+] + 0.41 * %GC - 675/length
    const tm = 81.5 + 16.6 * Math.log10(na_M) + 0.41 * gc - (675 / len);

    return {
      tm: Math.max(0, tm).toFixed(1),
      confidence: '±2.0°C',
      method: 'Schildkraut & Lifson (1965) GC%',
      warning: len < 50 ? 'GC% formula is less accurate for short sequences (<50bp).' : null,
      isValid: true
    };
  };

  /**
   * High-level Tm routing engine.
   * Automatically selects the most scientifically accurate model based on sequence length.
   */
  BioMath.calculateTm = function (seq, options = {}, sequenceType = 'DNA') {
    const s = seq.toUpperCase().replace(/[^ATGCNU]/g, '');
    if (s.length < 8) {
      return {
        tm: 'N/A',
        confidence: 'N/A',
        method: 'None',
        warning: 'Sequence too short for reliable Tm calculation.',
        isValid: false
      };
    }

    if (s.length <= 50) {
      return BioMath.calculateTmNN(s, options, sequenceType);
    } else {
      return BioMath.calculateTmGC(s, options); // We will assume calculateTmGC does not need RNA specifics for now, per prompt instructions.
    }
  };

  // ==========================================================
  // 3. Molecular Weight (MW) & Extinction Coefficient
  // ==========================================================
  BioMath.calculateDNA_MW = function(seq, isDoubleStranded = true, options = {}) {
      if (!seq) return 0;
      const s = seq.toUpperCase();
      const { saltForm = 'Na' } = options; // 'Na', 'K', 'Free'
      
      const monoisotopic = {
          'A': 313.21, 'T': 304.20, 'C': 289.18, 'G': 329.21, 'U': 306.20,
          'R': 321.21, 'Y': 296.69, 'S': 309.20, 'W': 308.71,
          'K': 316.71, 'M': 301.20, 'B': 307.53, 'D': 315.54,
          'H': 302.20, 'V': 310.53, 'N': 308.95
      };

      // Named constants — do not inline these
      const WATER_MW = 18.0153; // g/mol — H₂O molecular weight

      // Complement base map for dsDNA strand 2
      const COMPLEMENT = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A',
                           'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
                           'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N' };

      // ── Strand 1 (input sequence) ─────────────────────────────
      let ss_mw = 0;
      let validLen = 0;
      for (let i = 0; i < s.length; i++) {
          if (monoisotopic[s[i]]) {
              ss_mw += monoisotopic[s[i]];
              validLen++;
          }
      }

      if (validLen === 0) return 0;

      const PHOSPHATE_OFFSET = 61.965; // PO3- minus H (offset from OH to phosphate)
      const TERMINAL_WATER = 18.0153; // H on 5', OH on 3'

      ss_mw += TERMINAL_WATER;

      if (options.fivePrimePhosphate) {
          ss_mw += PHOSPHATE_OFFSET;
      }

      let saltAdj = 0;
      if (saltForm === 'Free') saltAdj = -21.99; // Remove Na, add H
      else if (saltForm === 'K') saltAdj = 17.10; // Swap Na (22.99) for K (39.1)

      ss_mw += (validLen * saltAdj); // Adjust for salt form counter-ions

      if (isDoubleStranded) {
          // ── Strand 2 (reverse complement) ───────────────────────
          // MW_dsDNA = MW_strand1 + MW_strand2. No water term at duplex level.
          let comp_mw = 0;
          let compLen = 0;
          for (let i = 0; i < s.length; i++) {
              const cb = COMPLEMENT[s[i]] || 'N';
              if (monoisotopic[cb]) {
                  comp_mw += monoisotopic[cb];
                  compLen++;
              }
          }
          comp_mw += TERMINAL_WATER;
          if (options.fivePrimePhosphate) {
              comp_mw += PHOSPHATE_OFFSET;
          }
          comp_mw += (compLen * saltAdj);

          const ds_mw = ss_mw + comp_mw;
          return (ds_mw / 1000).toFixed(2);
      } else {
          return (ss_mw / 1000).toFixed(2);
      }
  };

  /**
   * Calculates extinction coefficient at 280nm for proteins (Pace 1995)
   * or at 260nm for nucleic acids (Nearest-Neighbor model).
   * @param {string} seq - The sequence.
   * @param {string} type - 'dna' | 'rna' | 'protein'
   * @param {string} redoxState - For proteins: 'reduced' (default) | 'oxidized'
   *   Reduced: all Cys are free thiols, ε(Cys)=0.
   *   Oxidized: Cys form disulfide bonds, ε(half-Cys)=125 M⁻¹cm⁻¹.
   */
  BioMath.calculateExtinctionCoefficient = function(seq, type='dna', redoxState='reduced') {
      // Calculates Extinction Coefficient (e) at 260nm for DNA/RNA
      // Useful for A260 -> Concentration conversions
      if (!seq || seq.length === 0) return 0;
      const s = seq.toUpperCase().replace(/[^ATGCU]/g, '');
      
      let e = 0;
      // Nearest-Neighbor Extinction Coefficient values at 260 nm (L/(mol·cm))
      const NN_EC = {
        'AA': 27400, 'AC': 21200, 'AG': 25000, 'AT': 22800, 'AU': 22800,
        'CA': 21200, 'CC': 14600, 'CG': 18000, 'CT': 15200, 'CU': 15200,
        'GA': 25200, 'GC': 17600, 'GG': 21600, 'GT': 20000, 'GU': 20000,
        'TA': 23400, 'TC': 16200, 'TG': 19000, 'TT': 16800, 'TU': 16800,
        'UA': 24600, 'UC': 17200, 'UG': 19600, 'UT': 16800, 'UU': 19600
      };
      
      // Individual bases (for subtracting overlaps)
      const IND_EC = { 'A': 15400, 'C': 7300, 'G': 11700, 'T': 8700, 'U': 10000 };

      // Protein specific extinction (Pace 1995, Protein Sci. 4:2411-2423)
      // ε(280) = nW * 5500 + nY * 1490 + nC_disulfide * 125
      // Reduced form: free Cys (-SH) contributes ε≈0 at 280nm.
      // Oxidized form: cystine (disulfide, -S-S-) contributes 125 M⁻¹cm⁻¹ per half-cystine.
      if (type === 'protein') {
          const counts = { 'W': 0, 'Y': 0, 'C': 0 };
          for (const aa of seq.toUpperCase()) {
              if (counts[aa] !== undefined) counts[aa]++;
          }
          const cysContribution = (redoxState === 'oxidized') ? (counts['C'] * 125) : 0;
          return (counts['W'] * 5500) + (counts['Y'] * 1490) + cysContribution;
      }

      for (let i = 0; i < s.length - 1; i++) {
          const pair = s.substring(i, i+2);
          if (NN_EC[pair]) e += NN_EC[pair];
      }
      for (let i = 1; i < s.length - 1; i++) {
          if (IND_EC[s[i]]) e -= IND_EC[s[i]];
      }
      return e;
  };

  /**
   * K-mer pre-filter for Smith-Waterman (Larger Scale search optimization)
   * Returns true if sequences share at least one k-mer of length k.
   */
  BioMath.hasKmerMatch = function(seq1, seq2, k = 6) {
      if (seq1.length < k || seq2.length < k) return true; // Fallback to full DP for short seqs
      const set = new Set();
      for (let i = 0; i <= seq1.length - k; i++) {
        set.add(seq1.substring(i, i + k));
      }
      for (let i = 0; i <= seq2.length - k; i++) {
        if (set.has(seq2.substring(i, i + k))) return true;
      }
      return false;
  };

  /**
   * Calculate theoretical A260/A280 ratio based on extinction coefficients.
   * NOTE: This is a THEORETICAL prediction, not a measured value.
   * Actual A260/A280 ratios are measured via UV spectrophotometry (Nanodrop, etc.)
   */
  BioMath.calculateA260_A280_Theoretical = function(seq) {
      if (!seq || seq.length === 0) return { ratio: '1.80', note: '', expectedRange: {} };
      
      // Method: Use Nearest-Neighbor extinction coefficients
      // This gives ε260 and ε280 for the duplex
      const s = seq.toUpperCase();
      
      // NN extinction coefficients (per base pair) at 260 and 280 nm
      const NN_EC_260 = {
          'AA': 27400, 'AC': 21200, 'AG': 25000, 'AT': 22800,
          'CA': 21200, 'CC': 14600, 'CG': 18000, 'CT': 15200,
          'GA': 25200, 'GC': 17600, 'GG': 21600, 'GT': 20000,
          'TA': 23400, 'TC': 16200, 'TG': 19000, 'TT': 16800
      };
      
      const NN_EC_280 = {
          'AA': 3900,  'AC': 2900,  'AG': 3600,  'AT': 2900,
          'CA': 2900,  'CC': 2100,  'CG': 2500,  'CT': 2000,
          'GA': 3600,  'GC': 2500,  'GG': 2900,  'GT': 2500,
          'TA': 2900,  'TC': 2000,  'TG': 2400,  'TT': 2100
      };
      
      let sum260 = 0, sum280 = 0;
      
      // Sum nearest-neighbor pairs
      for (let i = 0; i < s.length - 1; i++) {
          const pair = s.substring(i, i + 2);
          if (NN_EC_260[pair]) {
              sum260 += NN_EC_260[pair];
              sum280 += NN_EC_280[pair];
          }
      }
      
      // Fallback if no valid pairs (e.g., sequence too short)
      if (sum280 === 0) return { ratio: '1.80', note: '', expectedRange: { pureDNA: '1.8–1.9', contaminated: '1.5–1.6' } };
      
      const ratio = (sum260 / sum280).toFixed(2);
      
      return {
          ratio: ratio,
          method: 'Theoretical (Nearest-Neighbor extinction coefficients)',
          note: 'This is NOT a measured A260/A280. ' +
                'For actual measurements, use Nanodrop or UV spectrophotometry.',
          expectedRange: {
              pureDNA: '1.8–1.9',
              contaminated: '1.5–1.6'
          }
      };
  };

  /**
   * Calculates protein molecular weight.
   * @param {string} seq - Canonical protein sequence (no headers, no stop codons).
   * @param {string} massType - 'average' (default) | 'monoisotopic'
   *   Average: standard residue masses for solution-phase estimation.
   *   Monoisotopic: most abundant isotopologue masses for MS applications.
   * @returns {{ kDa: string, internalStop: boolean }}
   */
  BioMath.calculateProtein_MW = function(seq, massType='average') {
      // Average residue masses — Fasman (1989) CRC Handbook of Biochemistry
      const averageMass = {
          'A': 71.0788,  'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
          'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
          'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
          'M': 131.1926, 'F': 147.1766, 'P': 97.1167,  'S': 87.0782,
          'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326,
          // Non-canonical: use best available average masses
          'X': 111.1000, // Unknown — approximate average
          'U': 150.0388, // Selenocysteine (avg)
          'Z': 128.6231  // Glx (Glu/Gln avg)
      };

      // Monoisotopic residue masses — NIST / Roepstorff & Fohlman (1984)
      const monoisotopicMass = {
          'A': 71.03711,  'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
          'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
          'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
          'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
          'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
          // Non-canonical
          'X': 111.00000, // Unknown placeholder
          'U': 150.95363, // Selenocysteine (monoisotopic)
          'Z': 128.05858  // Use Gln mass for Glx
      };

      // Stop codon middleware — strip trailing stops, detect internal stops
      const trailingStripped = seq.replace(/\*+$/, '');
      const internalStop = trailingStripped.includes('*');
      // For mass calculation, truncate at first internal stop
      const workingSeq = internalStop ? trailingStripped.split('*')[0] : trailingStripped;

      const table = (massType === 'monoisotopic') ? monoisotopicMass : averageMass;
      // H₂O: N-terminal H + C-terminal OH = 18.01056 (mono) or 18.01524 (avg)
      let mass = (massType === 'monoisotopic') ? 18.01056 : 18.01524;

      for (let i = 0; i < workingSeq.length; i++) {
          const char = workingSeq[i].toUpperCase();
          if (table[char] !== undefined) mass += table[char];
      }
      return { kDa: (mass / 1000).toFixed(2), internalStop };
  };

  // ==========================================================
  // 4. Kyte-Doolittle Hydrophobicity Scale
  // ==========================================================
  /*
     Kyte, J. and Doolittle, R.F. (1982). J. Mol. Biol. 157:105-132.
  */
  const KD_SCALE = {
      'R': -4.5, 'K': -3.9, 'N': -3.5, 'D': -3.5,
      'Q': -3.5, 'E': -3.5, 'H': -3.2, 'P': -1.6,
      'Y': -1.3, 'W': -0.9, 'S': -0.8, 'T': -0.7,
      'G': -0.4, 'A': 1.8,  'M': 1.9,  'C': 2.5,
  };
  BioMath.KD_SCALE = KD_SCALE;

  BioMath.calculateHydrophobicity = function(seq) {
      if (!seq || seq.length === 0) return 0;
      let sum = 0;
      let count = 0;
      for (let i = 0; i < seq.length; i++) {
          const val = KD_SCALE[seq[i].toUpperCase()];
          if (val !== undefined) {
              sum += val;
              count++;
          }
      }
      return count > 0 ? (sum / count).toFixed(2) : "0.00";
  };

  BioMath.calculateAliphaticIndex = function(seq) {
      // Ikai (1980) J. Biochem. 88, 1895-1898
      // Aliphatic index = X(Ala) + a * X(Val) + b * (X(Ile) + X(Leu))
      if (!seq || seq.length === 0) return 0;
      const s = seq.toUpperCase();
      const len = s.length;
      const A = (s.match(/A/g)||[]).length;
      const V = (s.match(/V/g)||[]).length;
      const I = (s.match(/I/g)||[]).length;
      const L = (s.match(/L/g)||[]).length;
      
      const molPercent = (AA) => (AA / len) * 100;
      // Constants a = 2.9, b = 3.9
      const ai = molPercent(A) + 2.9 * molPercent(V) + 3.9 * (molPercent(I) + molPercent(L));
      return ai.toFixed(2);
  };

  BioMath.calculateInstabilityIndex = function(seq) {
      // Guruprasad et al. (1990) Protein Eng. 4, 155-161
      // Estimates in vivo stability. <40 is stable, >40 is unstable.
      if (!seq || seq.length < 2) return 0;
      const s = seq.toUpperCase();
      
      // DIWV dipeptide instability weight values
      const DIWV = {
          'W':{W:1,V:-7,L:13.34,I:1,M:24.68,F:1,Y:-5.55,C:1,P:1,G:-9,A:-14.03,S:1,T:-14,N:13.34,Q:1,H:1,K:1,R:2,D:1,E:-14},
          'V':{W:-6.54,V:1,L:1,I:1,M:1,F:10,Y:-6.54,C:1,P:20.26,G:-7.49,A:1,S:1,T:7,N:1,Q:1,H:1,K:-1.88,R:1,D:-14,E:1},
          'L':{W:1,V:1,L:1,I:1,M:1,F:1,Y:1,C:1,P:20.26,G:1,A:1,S:1,T:1,N:1,Q:39.1,H:1,K:1,R:30,D:1,E:1},
          'I':{W:1,V:-7,L:20,I:1,M:1,F:1,Y:1,C:1,P:-1.88,G:1,A:1,S:1,T:1,N:1,Q:1,H:13.34,K:1,R:1,D:1,E:44.94},
          'M':{W:1,V:1,L:1,I:1,M:-1.88,F:1,Y:1,C:1,P:44.94,G:1,A:13.34,S:44.94,T:1,N:1,Q:-6.54,H:58.28,K:1,R:-6.54,D:1,E:1},
          'F':{W:1,V:1,L:1,I:1,M:1,F:1,Y:33.6,C:1,P:20.26,G:1,A:1,S:1,T:1,N:1,Q:1,H:1,K:-14,R:1,D:13.34,E:1},
          'Y':{W:-14,V:1,L:1,I:1,M:1,F:1,Y:13.34,C:1,P:1,G:-7.49,A:24.68,S:1,T:-7.49,N:1,Q:1,H:13.34,K:1,R:15,D:24.68,E:-5.55},
          'C':{W:1,V:-6.54,L:20,I:1,M:1,F:1,Y:1,C:1,P:20.26,G:1,A:1,S:1,T:33.6,N:1,Q:-6.54,H:33.6,K:1,R:1,D:20,E:1},
          'P':{W:-1.88,V:20.26,L:1,I:1,M:1,F:20.26,Y:1,C:1,P:20.26,G:1,A:20.26,S:20.26,T:1,N:1,Q:20.26,H:1,K:1,R:1,D:-14,E:18.38},
          'G':{W:-14,V:1,L:1,I:-7.49,M:1,F:1,Y:-7.49,C:1,P:1,G:13.34,A:1,S:1,T:1,N:-7.49,Q:1,H:1,K:1,R:1,D:1,E:-14},
          'A':{W:1,V:1,L:1,I:1,M:1,F:1,Y:1,C:1,P:20.26,G:1,A:1,S:1,T:1,N:1,Q:1,H:1,K:1,R:1,D:1,E:1},
          'S':{W:1,V:1,L:1,I:1,M:1,F:1,Y:1,C:1,P:44.94,G:1,A:1,S:20.26,T:1,N:1,Q:20.26,H:1,K:1,R:20.26,D:1,E:20.26},
          'T':{W:-14,V:1,L:1,I:1,M:1,F:13.34,Y:1,C:1,P:1,G:1,A:1,S:1,T:1,N:-14,Q:-14,H:1,K:1,R:1,D:1,E:20.26},
          'N':{W:1,V:1,L:1,I:44.94,M:1,F:1,Y:1,C:1,P:1,G:1,A:1,S:1,T:1,N:1,Q:-6.54,H:1,K:24.68,R:1,D:1,E:1},
          'Q':{W:1,V:-6.54,L:1,I:1,M:1,F:-6.54,Y:-6.54,C:-47,P:20.26,G:1,A:1,S:44.94,T:1,N:1,Q:20.26,H:1,K:1,R:1,D:20.26,E:20.26},
          'H':{W:-1.88,V:1,L:1,I:44.94,M:1,F:-9,Y:44.94,C:1,P:-1.88,G:1,A:1,S:1,T:-14,N:24.68,Q:1,H:1,K:24.68,R:1,D:1,E:1},
          'K':{W:1,V:1,L:1,I:1,M:44.94,F:1,Y:1,C:1,P:-6.54,G:1,A:1,S:1,T:1,N:1,Q:24.68,H:1,K:1,R:33.6,D:1,E:1},
          'R':{W:58.28,V:1,L:1,I:1,M:1,F:1,Y:-6.54,C:1,P:20.26,G:-7.49,A:1,S:1,T:1,N:13.34,Q:20.26,H:20.26,K:1,R:58.28,D:1,E:1},
          'D':{W:1,V:1,L:1,I:1,M:1,F:-6.54,Y:1,C:1,P:1,G:1,A:1,S:20.26,T:-14,N:1,Q:1,H:1,K:-7.49,R:-14,D:1,E:1},
          'E':{W:-14,V:1,L:1,I:1,M:1,F:1,Y:1,C:44.94,P:20.26,G:1,A:1,S:20.26,T:1,N:1,Q:20.26,H:-6.54,K:1,R:1,D:20.26,E:33.6}
      };
      
      let sum = 0;
      for (let i = 0; i < s.length - 1; i++) {
          const a1 = s[i];
          const a2 = s[i+1];
          // Default to 1 if one of the amino acids is non-standard (e.g. X)
          const weight = (DIWV[a1] && DIWV[a1][a2]) !== undefined ? DIWV[a1][a2] : 1;
          sum += weight;
      }
      
      const ii = (10 / seq.length) * sum;
      return ii.toFixed(2);
  };

  // ==========================================================
  // 5. Isoelectric Point (pI) Solver
  // ==========================================================
  const BJELLQVIST_PKA = {
    'N-term': 7.5, 'C-term': 3.55,
    'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0,
    'H': 5.98, 'K': 10.0, 'R': 12.0
  };

  BioMath.calculatePI = function(seq) {
    if (!seq || seq.length === 0) return 7.0;
    const s = seq.toUpperCase();
    
    // Count titratable groups
    const counts = { D:0, E:0, C:0, Y:0, H:0, K:0, R:0 };
    for (const aa of s) {
      if (counts[aa] !== undefined) counts[aa]++;
    }

    const netCharge = (pH) => {
      let charge = 0;
      // Basic groups (Positive when pH < pKa)
      charge += 1 / (1 + Math.pow(10, pH - BJELLQVIST_PKA['N-term']));
      charge += counts['K'] / (1 + Math.pow(10, pH - BJELLQVIST_PKA['K']));
      charge += counts['R'] / (1 + Math.pow(10, pH - BJELLQVIST_PKA['R']));
      charge += counts['H'] / (1 + Math.pow(10, pH - BJELLQVIST_PKA['H']));
      
      // Acidic groups (Negative when pH > pKa)
      charge -= 1 / (1 + Math.pow(10, BJELLQVIST_PKA['C-term'] - pH));
      charge -= counts['D'] / (1 + Math.pow(10, BJELLQVIST_PKA['D'] - pH));
      charge -= counts['E'] / (1 + Math.pow(10, BJELLQVIST_PKA['E'] - pH));
      charge -= counts['C'] / (1 + Math.pow(10, BJELLQVIST_PKA['C'] - pH));
      charge -= counts['Y'] / (1 + Math.pow(10, BJELLQVIST_PKA['Y'] - pH));
      
      return charge;
    };

    // Iterative binary search for pH where netCharge ≈ 0
    let low = 0, high = 14, pI = 7;
    for (let i = 0; i < 20; i++) {
        pI = (low + high) / 2;
        if (netCharge(pI) > 0) low = pI;
        else high = pI;
    }
    return pI.toFixed(2);
  };

  // ==========================================================
  // 6. Standard Codon Translation
  // ==========================================================

  const GENETIC_CODES = {
    1: { // Standard
      'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
      'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
      'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
      'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
      'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
      'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
      'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
      'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
      'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
      'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
      'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
      'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
      'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
      'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
      'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
      'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
    },
    2: { // Vertebrate Mitochondrial
      'AGA':'*', 'AGG':'*', 'ATA':'M', 'TGA':'W'
    },
    3: { // Yeast Mitochondrial
      'ATA':'M', 'CTT':'T', 'CTC':'T', 'CTA':'T', 'CTG':'T', 'TGA':'W', 'CGA':'X', 'CGC':'X'
    },
    4: { // Mold, Protozoan, Coelenterate Mito
      'TGA':'W'
    },
    5: { // Invertebrate Mitochondrial
      'AGA':'S', 'AGG':'S', 'ATA':'M', 'TGA':'W'
    },
    6: { // Ciliate, Dasycladacean, Hexamita Nuclear
      'TAA':'Q', 'TAG':'Q'
    },
    11: { // Bacterial, Archaeal, Plant Plastid
      // Identical to Standard (1) but with different alternative start codons (not handled fully here)
    }
  };

  BioMath.translateDNA = function (seq, options = {}) {
    const { 
        tableId = 1, 
        stopAtFirst = false, 
        requireStart = false 
    } = options;
    
    const table = Object.assign({}, GENETIC_CODES[1], GENETIC_CODES[tableId] || {});
    const s = seq.toUpperCase().replace(/U/g, 'T').replace(/[^ATCG]/g, 'N');
    
    let protein = '';
    let started = !requireStart;

    for (let i = 0; i < s.length - 2; i += 3) {
      const codon = s.substring(i, i + 3);
      if (!started) {
        if (codon === 'ATG') started = true;
        else continue;
      }

      const aa = table[codon] || 'X';
      if (aa === '*') {
        protein += '*';
        if (stopAtFirst) break;
      } else {
        protein += aa;
      }
    }

    return protein || (requireStart ? 'No ATG start found' : '');
  };
  // ==========================================================
  // 3. SECIENTIFIC UPGRADE: Molecular Weight & Extinction
  // ==========================================================

  /**
   * Scientific-Grade RNA Molecular Mass Engine
   * Uses monoresidue weights: AMP(329.21), UMP(306.17), CMP(305.18), GMP(345.21)
   */
  BioMath.calculateScientificRNA_MW = function(seq, options = {}) {
    if (!seq) return { avg: 0, min: 0, max: 0 };
    const s = seq.toUpperCase().replace(/T/g, 'U');
    const table = window.BioKit.data.RNA_MONORESIDUE_MASS;
    const termTable = window.BioKit.data.RNA_TERMINAL_GROUPS;
    
    let avg = 0, min = 0, max = 0;
    let validLen = 0;

    for (let i = 0; i < s.length; i++) {
        const char = s[i];
        if (table[char]) {
            if (typeof table[char] === 'object') {
                avg += table[char].avg;
                min += table[char].min;
                max += table[char].max;
            } else {
                avg += table[char];
                min += table[char];
                max += table[char];
            }
            validLen++;
        }
    }

    if (validLen === 0) return { avg: 0, min: 0, max: 0 };

    // Apply Topology & Terminal Groups
    const topology = options.topology || 'linear'; // 'linear' or 'circular'
    const terminal = options.terminal || 'OH';     // 'OH', 'MONO', 'TRI'
    const termWeight = termTable[terminal] || 18.015;

    avg += termWeight;
    min += termWeight;
    max += termWeight;

    if (topology === 'circular') {
        avg -= 18.015;
        min -= 18.015;
        max -= 18.015;
    }

    return {
        avg: parseFloat(avg.toFixed(2)),
        min: parseFloat(min.toFixed(2)),
        max: parseFloat(max.toFixed(2)),
        unit: 'Da'
    };
  };

  /**
   * Advanced Thermodynamic Stability (Nearest-Neighbor)
   * Turner 2004 for RNA-RNA, Sugimoto 1995 for RNA-DNA hybrids.
   * Upgrades: IUPAC degenerate base filtering, Von Ahsen (2001) free Mg²⁺
   * correction for dNTP chelation, Owczarzy (2004) threshold-based salt logic.
   */
  BioMath.calculateScientificTm = function(seq, options = {}) {
    if (!seq || seq.length < 2) return null;

    // ── IUPAC Degenerate Base Detection ──────────────────────────
    const IUPAC_DEGENERATE = new Set(['N','R','Y','W','S','M','K','B','V','D','H','I']);
    const rawUpper = seq.toUpperCase();
    const degenerateBases = [...rawUpper].filter(c => IUPAC_DEGENERATE.has(c));
    const hasDegenerateBases = degenerateBases.length > 0;

    // ── Degenerate Permutation Matrix ────────────────────────────
    // Per Biochemical Engine v3: stop stripping degenerate bases. Instead,
    // expand the full permutation pool (up to 64), run the NN Tm engine over
    // every exact variant, and return the [minTm, maxTm] range. The minimum Tm
    // represents the weakest binding variant and governs wet-lab annealing.
    if (hasDegenerateBases) {
      const { Mg = 0, dNTPs = 0.0008 } = options;
      const totalMgDeg = (typeof Mg === 'number' && Mg > 0) ? Mg : 0;
      const freeMgDeg  = Math.max(0, totalMgDeg - dNTPs);

      const expansion = BioMath.expandDegenerateSequence(rawUpper);
      let minTm = Infinity, maxTm = -Infinity;

      for (const exactSeq of expansion.sequences) {
        const r = BioMath.calculateScientificTm(exactSeq, options);
        if (r && r.tm !== 'N/A') {
          const v = parseFloat(r.tm);
          if (!isNaN(v)) {
            if (v < minTm) minTm = v;
            if (v > maxTm) maxTm = v;
          }
        }
      }

      if (!isFinite(minTm)) {
        return {
          tm: 'N/A', hasDegenerateBases: true, degenerateCount: degenerateBases.length,
          freeMg: null, formula: '', latex: `T_m = \\text{N/A}`, degenerateNote: null
        };
      }

      const minStr = minTm.toFixed(2);
      const maxStr = maxTm.toFixed(2);
      const degNote = `Minimum $T_m$ (${minStr}\u00b0C) represents the weakest binding` +
        ` variant in the degenerate pool and should dictate annealing conditions.` +
        (expansion.truncated ? ` ${expansion.warning}` : '');

      return {
        tm:              `${minStr} - ${maxStr}`,
        minTm:           minStr,
        maxTm:           maxStr,
        hasDegenerateBases: true,
        degenerateCount: degenerateBases.length,
        freeMg:          totalMgDeg > 0 ? freeMgDeg : null,
        formula:         `T_m = \\frac{\\Delta H}{\\Delta S + R \\ln(\\frac{C_t}{4})} - 273.15`,
        latex:           `T_m = ${minStr} ^\\circ\\text{C} {-} ${maxStr} ^\\circ\\text{C}`,
        degenerateNote:  degNote
      };
    }

    // ── Exact-sequence NN path (no degenerate bases) ──────────────
    // Convert DNA T → RNA U for the Turner/Sugimoto NN tables.
    const s = rawUpper.replace(/[^ATGCU]/g, '').replace(/T/g, 'U');

    if (s.length < 2) return {
        tm: 'N/A', hasDegenerateBases, degenerateCount: degenerateBases.length,
        freeMg: null, formula: '', latex: `T_m = \\text{N/A}`, degenerateNote: null
    };

    const {
        type  = 'rna-rna', // 'rna-rna' | 'rna-dna'
        Ct    = 0.5e-6,    // strand concentration (M)
        Na    = 0.05,      // [Na⁺] (M)
        Mg    = 0,         // Total [Mg²⁺] (M)
        dNTPs = 0.0008,    // Von Ahsen (2001): standard 4 × 0.2 mM = 0.8 mM total dNTPs
    } = options;

    // ── Von Ahsen (2001) Free Mg²⁺ Correction ────────────────────
    // dNTPs chelate Mg²⁺ 1:1; free Mg²⁺ drives thermodynamics, not total.
    const totalMg = (typeof Mg === 'number' && Mg > 0) ? Mg : 0;
    const freeMg  = Math.max(0, totalMg - dNTPs);

    const data   = window.BioKit.data;
    const params = (type === 'rna-dna') ? data.SUGIMOTO_1995 : data.TURNER_2004;
    const R      = 1.987; // cal / (mol·K)

    let sumH = params.init.dH * 1000; // kcal → cal
    let sumS = params.init.dS;

    for (let i = 0; i < s.length - 1; i++) {
        const pair = s.substring(i, i + 2);
        if (params[pair]) {
            sumH += params[pair].dH * 1000;
            sumS += params[pair].dS;
        } else {
            // Median NN approximation for unrecognised dinucleotides
            sumH += -10.0 * 1000;
            sumS += -25.0;
        }
    }

    // Ct/4 for non-self-complementary strands; Ct for self-complementary
    const cleanForRC   = rawUpper.replace(/[^ATGCU]/g, '').replace(/T/g, 'U');
    const isSelfComp   = (s === BioMath.reverseComplement(cleanForRC));
    const ctAdjustment = isSelfComp ? Ct : (Ct / 4);

    let tm = (sumH / (sumS + R * Math.log(ctAdjustment))) - 273.15;

    // ── Salt / Entropy Correction ─────────────────────────────────
    // Owczarzy (2004): when sqrt([freeMg]) / [Na] > 0.22, Mg²⁺ dominates.
    // Mg²⁺-dominant correction uses Sugimoto (2001) coefficient (12.0).
    let saltAdj;
    if (freeMg > 0 && Na > 0 && (Math.sqrt(freeMg) / Na) > 0.22) {
        saltAdj = 12.0 * Math.log10(freeMg); // Mg²⁺-dominant regime
    } else {
        saltAdj = 16.6 * Math.log10(Na);     // Na⁺-dominant regime (SantaLucia 1998)
    }
    tm += saltAdj;

    return {
        tm: tm.toFixed(2),
        dH: (sumH / 1000).toFixed(2),
        dS: sumS.toFixed(2),
        hasDegenerateBases,
        degenerateCount: degenerateBases.length,
        freeMg: totalMg > 0 ? freeMg : null,
        formula: `T_m = \\frac{\\Delta H}{\\Delta S + R \\ln(\\frac{C_t}{4})} - 273.15`,
        latex:   `T_m = ${tm.toFixed(2)} ^\\circ\\text{C}`,
        degenerateNote: null
    };
  };

  /**
   * Industrial Metrics & Molar Calculations
   */
  BioMath.calculateIndustrialMetrics = function(seq, options = {}) {
    if (!seq) return null;
    const s = seq.toUpperCase();
    const mw = parseFloat(BioMath.calculateScientificRNA_MW(s, options).avg);
    
    // Molar Extinction Coefficient (e260) NN model
    const e260 = BioMath.calculateExtinctionCoefficient(s, 'rna');
    
    // Mass Conc: 1 A260 = 40 ug/ml for RNA
    const massConcPerA260 = 40; // ug/ml
    
    // Molar Conc at 1 A260: (MassConc / MW)
    // 40 ug/ml = 0.04 g/L. Molar = (0.04 / MW) mol/L
    const molarConcPerA260 = (0.04 / mw); 
    
    // Copy Numbers (Avogadro's Number: 6.022e23)
    const Na = 6.02214076e23;
    const copiesPerUg = (1e-6 / mw) * Na;
    
    return {
        mw: mw,
        e260: e260,
        molarConc: molarConcPerA260.toExponential(3),
        copies: copiesPerUg.toExponential(3),
        latex: {
            e260: `\\epsilon_{260} = ${e260.toLocaleString()} \\text{ L mol}^{-1} \\text{ cm}^{-1}`,
            copies: `N = ${copiesPerUg.toExponential(2)} \\text{ copies/\\mu g}`
        }
    };
  };

  /**
   * Expands a degenerate IUPAC DNA sequence into all exact-base permutations.
   * Expansion is capped at 64 sequences. If the true permutation count exceeds
   * this limit, the first 64 are returned with a warning flag set.
   *
   * @param {string} sequence - Input sequence; IUPAC ambiguity codes allowed.
   * @returns {{ sequences: string[], truncated: boolean, warning: string|null }}
   */
  BioMath.expandDegenerateSequence = function(sequence) {
    // Normalise to DNA space (U → T)
    const s = sequence.toUpperCase().replace(/U/g, 'T');

    const IUPAC_EXPAND = {
      'R': ['A', 'G'],       'Y': ['C', 'T'],       'S': ['G', 'C'],
      'W': ['A', 'T'],       'K': ['G', 'T'],       'M': ['A', 'C'],
      'B': ['C', 'G', 'T'],  'D': ['A', 'G', 'T'],  'H': ['A', 'C', 'T'],
      'V': ['A', 'C', 'G'],  'N': ['A', 'C', 'G', 'T']
    };

    const MAX_PERMS = 64;
    let sequences = [''];
    let truncated  = false;

    for (let pos = 0; pos < s.length; pos++) {
      const base       = s[pos];
      const expansions = IUPAC_EXPAND[base] || [base];

      if (expansions.length === 1) {
        // Non-degenerate base: append in-place, no branching
        for (let k = 0; k < sequences.length; k++) sequences[k] += expansions[0];
      } else {
        const newSeqs = [];
        let   limitHit = false;

        for (let k = 0; k < sequences.length && !limitHit; k++) {
          for (let e = 0; e < expansions.length && !limitHit; e++) {
            newSeqs.push(sequences[k] + expansions[e]);
            if (newSeqs.length >= MAX_PERMS) limitHit = true;
          }
        }

        sequences = newSeqs;

        if (limitHit) {
          // Hard limit reached: complete all partial sequences deterministically
          // by appending the first (alphabetically lowest) expansion for every
          // remaining degenerate position. Sequences are fully-formed but the
          // pool is a subset of the true permutation space.
          for (let remPos = pos + 1; remPos < s.length; remPos++) {
            const remBase = s[remPos];
            const remExp  = (IUPAC_EXPAND[remBase] || [remBase])[0];
            for (let k = 0; k < sequences.length; k++) sequences[k] += remExp;
          }
          truncated = true;
          break;
        }
      }
    }

    return {
      sequences,
      truncated,
      warning: truncated ? 'Sequence highly degenerate; calculating partial pool.' : null
    };
  };

  /**
   * Thermodynamic self-dimer ΔG calculation (SantaLucia 1998 Nearest-Neighbor,
   * Owczarzy 2004 / Von Ahsen 2001 salt correction).
   *
   * Algorithm:
   *   1. Compute the reverse complement (RC) of the sequence.
   *   2. Slide the sequence against its RC over all integer offsets to find every
   *      contiguous Watson-Crick matching region (no gaps).
   *   3. For each contiguous run ≥ 4 bp, accumulate ΔH and ΔS via the SantaLucia
   *      1998 Nearest-Neighbor parameter table (NN_PARAMS / NN_INIT).
   *   4. Track the run whose 1 M NaCl ΔG at 310.15 K is most negative.
   *   5. Apply Owczarzy (2004) / Von Ahsen (2001) salt correction:
   *        Tm_1M_K  = ΔH_cal / ΔS_1M
   *        Tm_salt_K = Tm_1M_K + salt_adj   (salt_adj from Owczarzy regime selector)
   *        ΔS_salt  = ΔH_cal / Tm_salt_K    (back-derived from corrected Tm)
   *      This yields:  ΔG(37 °C) = ΔH_kcal · (1 − 310.15 / Tm_salt_K)
   *
   * @param {string} sequence - Primer sequence (U → T normalised internally).
   * @param {number} Na       - [Na⁺] in M.
   * @param {number} freeMg   - Free [Mg²⁺] in M after dNTP chelation (Von Ahsen 2001).
   * @returns {{ deltaG: number, warning: string|null }}
   */
  BioMath.calculateDimerDeltaG = function(sequence, Na, freeMg) {
    const T_K    = 310.15;  // 37 °C in Kelvin (thermodynamic baseline)
    const MIN_BP = 4;       // Minimum contiguous WC base pairs to evaluate

    // Normalise to DNA uppercase; strip non-canonical
    const s  = sequence.toUpperCase().replace(/U/g, 'T').replace(/[^ATGC]/g, '');
    if (s.length < MIN_BP) return { deltaG: 0, warning: null };

    const rc = BioMath.reverseComplement(s);
    const n  = s.length;

    // bestDG_1M tracks the most negative ΔG at 1 M NaCl, used only for ranking
    let bestDG_1M = 0;
    let bestDH    = 0;   // kcal/mol — ΔH of the best run
    let bestDS    = 0;   // cal/mol·K — ΔS of the best run

    // ── Sliding-window alignment: all integer offsets ─────────────────────
    // At each shift, j = i + shift maps a position in s to a position in rc.
    // Checking s[i] === rc[j] is equivalent to asking whether s[i] is the
    // Watson-Crick complement of the base at position (n-1-j) in the original
    // sequence — i.e., the correct antiparallel pairing for a self-dimer.
    for (let shift = -(n - 1); shift <= n - 1; shift++) {
      let runLen = 0, runDH = 0, runDS = 0;
      let firstBase = '', lastBase = '', prevBase = '';

      const evaluateRun = () => {
        if (runLen < MIN_BP) return;
        // Add SantaLucia 1998 initiation parameters for the terminal base pairs
        const iF  = NN_INIT[firstBase] || NN_INIT['G'];
        const iL  = NN_INIT[lastBase]  || NN_INIT['G'];
        const tDH = runDH + iF.dH + iL.dH;   // kcal/mol
        const tDS = runDS + iF.dS + iL.dS;   // cal/mol·K
        const dg  = tDH - T_K * (tDS / 1000);
        if (dg < bestDG_1M) {
          bestDG_1M = dg;
          bestDH    = tDH;
          bestDS    = tDS;
        }
      };

      for (let i = 0; i < n; i++) {
        const j = i + shift;
        if (j >= 0 && j < n && s[i] === rc[j]) {
          if (runLen === 0) firstBase = s[i];
          if (runLen >= 1) {
            const pair = prevBase + s[i];
            if (NN_PARAMS[pair]) { runDH += NN_PARAMS[pair][0]; runDS += NN_PARAMS[pair][1]; }
          }
          prevBase = s[i];
          lastBase = s[i];
          runLen++;
        } else {
          evaluateRun();
          runLen = 0; runDH = 0; runDS = 0;
          firstBase = ''; lastBase = ''; prevBase = '';
        }
      }
      evaluateRun(); // Flush any run reaching the end of the sequence
    }

    // No contiguous WC run of ≥ MIN_BP found at any offset → no significant dimer
    if (bestDH === 0 && bestDS === 0) return { deltaG: 0, warning: null };

    // ── Salt Correction (Owczarzy 2004 / Von Ahsen 2001) ─────────────────
    const na_M = (typeof Na     === 'number' && Na     > 0) ? Na     : 0.05;
    const mg_M = (typeof freeMg === 'number' && freeMg > 0) ? freeMg : 0;

    let saltAdj = 0;
    if (mg_M > 0 && (Math.sqrt(mg_M) / na_M) > 0.22) {
      saltAdj = 12.0 * Math.log10(mg_M);  // Mg²⁺-dominant regime (Owczarzy 2004)
    } else {
      saltAdj = 16.6 * Math.log10(na_M);  // Na⁺-dominant regime (SantaLucia 1998)
    }

    // ── ΔS correction derived from the salt-adjusted Tm ──────────────────
    // At 1 M NaCl (reference): Tm_1M_K = ΔH_cal / ΔS_1M
    // Salt correction shifts Tm by saltAdj (°C = K for differences):
    //   Tm_salt_K = Tm_1M_K + saltAdj
    // At Tm, ΔG = 0, so ΔS_salt = ΔH_cal / Tm_salt_K
    // Therefore: ΔG(310.15 K) = ΔH_kcal · (1 − 310.15 / Tm_salt_K)
    const dH_cal    = bestDH * 1000;
    const Tm_1M_K   = dH_cal / bestDS;
    const Tm_salt_K = Math.max(1, Tm_1M_K + saltAdj);  // clamp: must be > 0
    const finalDG   = bestDH * (1 - T_K / Tm_salt_K);

    // ── QC Threshold ─────────────────────────────────────────────────────
    const CRITICAL_THRESHOLD = -9.0;
    const warning = finalDG <= CRITICAL_THRESHOLD
      ? `CRITICAL: Severe self-dimer potential (\u0394G = ${finalDG.toFixed(2)} kcal/mol` +
        ` \u2264 \u22129.0 kcal/mol). High risk of PCR failure.`
      : null;

    return { deltaG: parseFloat(finalDG.toFixed(2)), warning };
  };

  /**
   * Synthesis & Biological Viability QC
   * Returns an array of human-readable alert strings. Empty array = all clear.
   * Checks: extreme GC content, severe homopolymer runs, secondary structure risk.
   */
  BioMath.runQualityControl = function(sequence, Na, freeMg) {
    const warnings = [];
    // Strip non-canonical characters; treat DNA and RNA uniformly
    const s = sequence.toUpperCase().replace(/[^ATGCU]/g, '');
    if (s.length === 0) return warnings;
    // Salt parameters forwarded to the thermodynamic self-dimer engine.
    // Defaults: 50 mM Na⁺ (standard PCR), 0 free Mg²⁺.
    const na_M  = (typeof Na     === 'number' && Na     > 0) ? Na     : 0.05;
    const fmg_M = (typeof freeMg === 'number' && freeMg >= 0) ? freeMg : 0;

    // ── GC Content ────────────────────────────────────────────────
    const gcCount = (s.match(/[GC]/g) || []).length;
    const gcPct   = (gcCount / s.length) * 100;
    if (gcPct < 30 || gcPct > 70) {
        warnings.push(
            `QC ALERT: Extreme GC Content (${gcPct.toFixed(1)}%) \u2014 Potential amplification failure.`
        );
    }

    // ── Homopolymer Runs ──────────────────────────────────────────
    // > 4 consecutive G's (G-quadruplex / synthesis failure risk)
    if (/G{5,}/i.test(s)) {
        warnings.push(
            'QC ALERT: Severe Homopolymer Run (\u22655 G\'s) \u2014 High synthesis risk / G-quadruplex potential.'
        );
    }
    // > 5 consecutive A, T, or C (synthesis stutter risk)
    if (/A{6,}|T{6,}|C{6,}|U{6,}/i.test(s)) {
        warnings.push(
            'QC ALERT: Severe Homopolymer Run (\u22656 A/T/C) \u2014 High synthesis risk / polymerase slippage.'
        );
    }

    // ── Self-Dimer: Thermodynamic ΔG (SantaLucia 1998 + Owczarzy/Von Ahsen) ──
    // Replaces the heuristic stem-count check with a physics-based ΔG calculation.
    // Warning is injected only when ΔG ≤ −9.0 kcal/mol (see calculateDimerDeltaG).
    const dimerResult = BioMath.calculateDimerDeltaG(s, na_M, fmg_M);
    if (dimerResult.warning) {
        warnings.push(dimerResult.warning);
    }

    // ── Hairpin: sliding inverted-repeat search (stem ≥ 4 bp, loop ≥ 4 nt) ──
    // Retained as a structural geometry screen independent of thermodynamic ΔG.
    // Cap at 200 nt for performance on genomic-length inputs.
    {
        const sDNA      = s.replace(/U/g, 'T');
        const COMP_HP   = { A: 'T', T: 'A', G: 'C', C: 'G' };
        const revCompHP = str => str.split('').reverse().map(c => COMP_HP[c] || 'N').join('');
        const MIN_STEM  = 4;
        const MIN_LOOP  = 4;
        const checkSeq  = sDNA.slice(0, 200);
        const maxStem   = Math.floor((checkSeq.length - MIN_LOOP) / 2);
        let   hairpin   = false;

        for (let stemLen = MIN_STEM; stemLen <= maxStem && !hairpin; stemLen++) {
            for (let i = 0; i <= checkSeq.length - 2 * stemLen - MIN_LOOP && !hairpin; i++) {
                const stem5RC = revCompHP(checkSeq.substr(i, stemLen));
                for (let j = i + stemLen + MIN_LOOP; j <= checkSeq.length - stemLen && !hairpin; j++) {
                    if (checkSeq.substr(j, stemLen) === stem5RC) hairpin = true;
                }
            }
        }

        if (hairpin) {
            warnings.push(
                'QC ALERT: Probable hairpin loop (stem \u22654 bp) \u2014 Secondary structure may inhibit hybridisation.'
            );
        }
    }

    return warnings;
  };

  /**
   * Enhanced 6-Frame Translation & ORF Detection
   * Supports all 33 NCBI genetic codes, alternative starts, and Kozak flagging.
   */
  BioMath.analyzeORFs = function(dnaSeq, options = {}) {
    const { 
        codeId = 1, 
        minLen = 30, // nt
        findAltStarts = true 
    } = options;

    const data = window.BioKit.data;
    const config = data.GENETIC_CODES[codeId] || data.GENETIC_CODES[1];
    const sequence = dnaSeq.toUpperCase();
    const strands = [sequence, BioMath.reverseComplement(sequence)];
    
    const codonMap = {};
    const baseOrder = "TCAG";
    let idx = 0;
    for (const b1 of baseOrder) {
        for (const b2 of baseOrder) {
            for (const b3 of baseOrder) {
                const codon = b1 + b2 + b3;
                codonMap[codon] = { aa: config.table[idx], isStart: config.starts[idx] === 'M' };
                idx++;
            }
        }
    }

    const orfs = [];
    strands.forEach((s, strandIdx) => {
        const strandSign = strandIdx === 0 ? '+' : '-';
        for (let frame = 0; frame < 3; frame++) {
            let currentORF = null;
            for (let i = frame; i <= s.length - 3; i += 3) {
                const codon = s.substring(i, i + 3);
                const info = codonMap[codon] || { aa: 'X', isStart: false };

                if (info.isStart && !currentORF) {
                    currentORF = { start: i, seq: '', frame: frame + 1, strand: strandSign };
                }

                if (currentORF) {
                    currentORF.seq += codon;
                    if (info.aa === '*') {
                        if (currentORF.seq.length >= minLen) {
                            // Kozak consensus check (GCC)RCCAUGG
                            // Pos -3 (index i_start - 3) and +4 (index i_start + 3)
                            const context = s.substring(currentORF.start - 6, currentORF.start + 4);
                            const kozakScore = BioMath.scoreKozak(context);
                            
                            orfs.push({
                                ...currentORF,
                                end: i + 3,
                                length: currentORF.seq.length,
                                translation: BioMath.translateDNA(currentORF.seq, { tableId: codeId }),
                                kozak: kozakScore
                            });
                        }
                        currentORF = null;
                    }
                }
            }
        }
    });

    return orfs.sort((a,b) => b.length - a.length);
  };

  BioMath.scoreKozak = function(context) {
    if (context.length < 10) return 'Weak';
    // Simplified Kozak: R at -3 (index 3) and G at +4 (index 9)
    const rAtMinus3 = context[3] === 'A' || context[3] === 'G';
    const gAtPlus4 = context[9] === 'G';
    if (rAtMinus3 && gAtPlus4) return 'Strong';
    if (rAtMinus3 || gAtPlus4) return 'Moderate';
    return 'Weak';
  };

  /**
   * Full IUPAC Compliance for Reverse Complement
   */
  BioMath.reverseComplement = function(seq) {
    if (!seq) return '';
    const map = {
      'A':'T', 'T':'A', 'G':'C', 'C':'G', 'U':'A',
      'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K',
      'B':'V', 'V':'B', 'D':'H', 'H':'D', 'N':'N'
    };
    return seq.toUpperCase().split('').reverse().map(b => map[b] || 'N').join('');
  };

  // Helper for LaTeX formatting
  BioMath.formatLaTeX = function(element) {
    if (window.renderMathInElement) {
        window.renderMathInElement(element, {
            delimiters: [
                {left: '$$', right: '$$', display: true},
                {left: '$', right: '$', display: false},
                {left: '\\(', right: '\\)', display: false},
                {left: '\\[', right: '\\]', display: true}
            ],
            throwOnError : false
        });
    }
  };

})();
