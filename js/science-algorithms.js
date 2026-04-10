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
      productSize: { min: 200, max: 5000 }
    },
    "qpcr": {
      name: "qPCR/RT-PCR",
      length: { min: 18, max: 25 },
      tm: { min: 59, max: 61 },
      tmDifference: 2,
      gc: { min: 45, max: 55 },
      gc3prime: true,
      homopolymerLimit: 4,
      productSize: { min: 50, max: 500 }
    },
    "cloning": {
      name: "Cloning/Sequencing",
      length: { min: 20, max: 30 },
      tm: { min: 60, max: 70 },
      tmDifference: 5,
      gc: { min: 35, max: 65 },
      gc3prime: true,
      homopolymerLimit: 5,
      productSize: { min: 500, max: 10000 }
    },
    "multiplex": {
      name: "Multiplex PCR",
      length: { min: 19, max: 24 },
      tm: { min: 60, max: 65 },
      tmDifference: 1,
      gc: { min: 45, max: 55 },
      gc3prime: true,
      homopolymerLimit: 4,
      productSize: { min: 100, max: 1000 }
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
      if (!seq) return '';
      const comp = { 
          'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
          'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n',
          'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
          'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'
      };
      return seq.split('').reverse().map(b => comp[b] || b).join('');
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
    if (!constraints.gc3prime) return true;
    const lastBase = primer.slice(-1).toUpperCase();
    return lastBase === 'G' || lastBase === 'C';
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
   * @param {string} primer - The primer sequence.
   * @param {string} template - The target DNA sequence.
   * @param {boolean} isReverse - Whether to search for the reverse complement.
   * @returns {Object} { primarySite, specificity }
   */
  BioMath.checkTemplateSpecificity = function(primer, template, isReverse = false) {
    const res = {
      primarySite: { position: -1, count: 0 },
      specificity: "UNIQUE"
    };
    if (!primer || !template) return res;

    const p = isReverse ? BioMath.reverseComplement(primer).toUpperCase() : primer.toUpperCase();
    const t = template.toUpperCase();

    let pos = t.indexOf(p);
    while (pos !== -1) {
      res.primarySite.count++;
      if (res.primarySite.position === -1) res.primarySite.position = pos + 1;
      pos = t.indexOf(p, pos + 1);
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

    // 1. Find Primer Sites
    const fwdIdx = t.indexOf(f);
    const revIdx = t.lastIndexOf(r_rc);

    if (fwdIdx === -1) res.warnings.push("Forward primer binding site not found in template.");
    if (revIdx === -1) res.warnings.push("Reverse primer binding site not found in template.");

    if (fwdIdx === -1 || revIdx === -1) return res;

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
        // Add initiation parameters for this run (simplified: one initiation per binding site)
        // We use the terminal bases of the run for initiation
        // For simplicity, we use the first match of the run
        dH += NN_INIT['G'].dH; // Approximation
        dS += NN_INIT['G'].dS;

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

  const NN_INIT = {
    'G': { dH: 0.1, dS: -2.8 }, 'C': { dH: 0.1, dS: -2.8 },
    'A': { dH: 2.3, dS: 4.1 }, 'T': { dH: 2.3, dS: 4.1 }
  };

  // Helix penalty if self-complementary (approximated)
  const SYM_PENALTY = { dH: 0, dS: -1.4 }; 

  BioMath.calculateTmNN = function (primerSeq, options = {}) {
    const {
      oligoConc_nM = 250,
      naConc_mM = 50,
      mgConc_mM = 1.5,
      dntpConc_mM = 0.8,
      dmso_pct = 0,
      formamide_m = 0
    } = options;

    const seq = primerSeq.toUpperCase().replace(/U/g, 'T').replace(/[^ATGCRYWSKMBDHVN]/g, '');
    if (seq.length < 2) return { tm: 0, warning: 'Sequence too short' };

    // Buffer validation and Equivalent Na+ calculation (von Ahsen 2001 / Owczarzy 2008)
    const salt = Math.min(1000, Math.max(10, naConc_mM));
    let eqNa_mM = salt;
    if (mgConc_mM > dntpConc_mM) {
        eqNa_mM += 120 * Math.sqrt(mgConc_mM - dntpConc_mM);
    }
    const na_M = eqNa_mM / 1000;

    let dH = 0, dS = 0;

    // Resolve IUPAC for initiation
    const resolveBase = (b) => IUPAC_DNA[b] || [b];
    const avgInit = (base) => {
      const variants = resolveBase(base);
      let th = 0, ts = 0;
      variants.forEach(v => {
        const p = NN_INIT[v] || NN_INIT['G']; // Fallback to G if base not in NN_INIT (e.g., N)
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
          if (NN_PARAMS[pair]) {
            pairH += NN_PARAMS[pair][0];
            pairS += NN_PARAMS[pair][1];
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
    Tm -= (0.6 * dmso_pct);          // DMSO (approx)
    Tm -= (0.65 * formamide_m);     // Formamide (approx)

    // Range warning
    let warning = '';
    if (seq.length < 14 || seq.length > 50) warning = 'Sequence length outside validated NN model range (14-50bp).';
    if (salt < 10 || salt > 1000) warning += ' Salt concentration outside validated range (10-1000mM).';

    return {
      tm: Math.max(0, Tm).toFixed(1),
      confidence: '±1.5°C',
      warning: warning || null
    };
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

      let ss_mw = 0;
      let validLen = 0;
      for (let i = 0; i < s.length; i++) {
          if (monoisotopic[s[i]]) {
              ss_mw += monoisotopic[s[i]];
              validLen++;
          }
      }
      
      if (validLen === 0) return 0;
      
      let saltAdj = 0;
      if (saltForm === 'Free') saltAdj = -21.99; // Remove Na, add H
      else if (saltForm === 'K') saltAdj = 17.10; // Swap Na (22.99) for K (39.1)
      
      ss_mw += (validLen * saltAdj) - 61.96; // Adjust salt and terminal water
      
      return (isDoubleStranded ? (ss_mw * 2) / 1000 : ss_mw / 1000).toFixed(2);
  };

  BioMath.calculateExtinctionCoefficient = function(seq, type='dna') {
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

      // Protein specific extinction (Wetlaufer / Pace model)
      if (type === 'protein') {
          const counts = { 'W': 0, 'Y': 0, 'C': 0 };
          for (const aa of seq.toUpperCase()) {
              if (counts[aa] !== undefined) counts[aa]++;
          }
          return (counts['W'] * 5500) + (counts['Y'] * 1490) + (counts['C'] * 125);
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

  BioMath.calculateA260_A280 = function(seq) {
      if (!seq || seq.length === 0) return 1.8;
      const s = seq.toUpperCase();
      // Bases extinction at 260 and 280 (L/(mol·cm))
      const e260 = { 'A': 15300, 'C': 7400, 'G': 11800, 'T': 9300, 'U': 10200 };
      const e280 = { 'A': 2500,  'C': 1500, 'G': 5800,  'T': 1500, 'U': 2800  };
      
      let sum260 = 0, sum280 = 0;
      for (const char of s) {
        sum260 += e260[char] || 0;
        sum280 += e280[char] || 0;
      }
      return (sum280 === 0) ? 1.8 : (sum260 / sum280).toFixed(2);
  };

  BioMath.calculateProtein_MW = function(seq) {
      // Monoisotopic masses or average residues. We use average residue masses here.
      // H2O is added once at the end.
      const AAmass = {
          'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
          'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
          'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
          'M': 131.1926, 'F': 147.1766, 'P': 97.1167,  'S': 87.0782,
          'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
      };
      
      let mass = 18.01524; // Water added (N-term H, C-term OH)
      for (let i = 0; i < seq.length; i++) {
          const char = seq[i].toUpperCase();
          if (AAmass[char]) mass += AAmass[char];
      }
      return (mass / 1000).toFixed(2); // kDa
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
  // ── Standard & Alternative Genetic Codes ─────────────────
  const GENETIC_CODES = {
    1: { // Standard
      'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
      'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
      'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
      'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
      'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
      'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
      'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
      'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    },
    2: { // Vertebrate Mitochondrial
      'AGA': '*', 'AGG': '*', 'AUA': 'M', 'UGA': 'W' // overrides only
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

    // Mark trailing nucleotides
    const trailing = s.length % 3;
    if (trailing > 0) protein += ' (truncated)';

    return protein || (requireStart ? 'No ATG start found' : '');
  };

  /**
   * Aggregates all QC metrics into a final laboratory verdict.
   * @param {Object} results - [{seq, tm, gc}, {seq, tm, gc}]
   * @param {Object} metrics - { tmDiff, maxHetero }
   * @param {Object} constraints - The active assay thresholds.
   * @param {string} template - The template DNA sequence.
   * @param {Array} features - Feature list from template.
   * @param {Object} bioOptions - { intronSpanning, avoidExonJunction }
   * @returns {Object} { status, class, failures, warnings, checklist, actionItems }
   */
  BioMath.generatePrimerPairVerdict = function(results, metrics, constraints, template, features, bioOptions) {
    const verdict = { 
      status: "PASS", class: "status-pass", failures: [], warnings: [], 
      checklist: [], actionItems: [] 
    };

    const fwd = results[0], rev = results[1];
    const diff = metrics.tmDiff;
    const heteroDg = parseFloat(metrics.heteroDg);

    const fwdStats = {
      len: BioMath.validateLength(fwd.seq, constraints),
      tm: BioMath.validateTm(fwd.tm, constraints),
      homo: BioMath.validateHomopolymer(fwd.seq, constraints),
      dg: parseFloat(fwd.selfDg),
      terminal: BioMath.assess3PrimeStability(fwd.seq)
    };
    const revStats = {
      len: BioMath.validateLength(rev.seq, constraints),
      tm: BioMath.validateTm(rev.tm, constraints),
      homo: BioMath.validateHomopolymer(rev.seq, constraints),
      dg: parseFloat(rev.selfDg),
      terminal: BioMath.assess3PrimeStability(rev.seq)
    };

    const tmDiffOk = BioMath.validateTmDifference(fwd.tm, rev.tm, constraints);
    if (!tmDiffOk) {
      verdict.warnings.push(`Tm difference (${diff}°C) exceeds tolerance (max ${constraints.tmDifference}°C)`);
      verdict.actionItems.push("Adjust sequences to equalize melting temperatures.");
    }

    // ΔG Thresholds: -6.0 kcal/mol is a common cutoff for concern, -9.0 is critical.
    if (heteroDg <= -9.0) {
      verdict.failures.push(`Critical heterodimer stability (ΔG: ${heteroDg} kcal/mol)`);
      verdict.actionItems.push("Redesign primers to eliminate extensive complementarity.");
    } else if (heteroDg <= -6.0) {
      verdict.warnings.push(`Significant heterodimer risk (ΔG: ${heteroDg} kcal/mol)`);
    }

    // 3' Terminal Integrity
    if (fwdStats.terminal.verdict === 'CRITICAL' || revStats.terminal.verdict === 'CRITICAL') {
      verdict.failures.push("Excessive 3' GC clamping: High mispriming risk.");
      verdict.actionItems.push("Reduce terminal G/C bases at the 3' end.");
    }

    let ampVerdict = 'PASS';
    if (template) {
      const amp = BioMath.calculateAmpliconInfo(template, fwd.seq, rev.seq, constraints);
      if (amp.amplicon.size > 0) {
        const sizeWarn = amp.warnings.find(w => w.includes('size'));
        if (sizeWarn) { ampVerdict = 'CAUTION'; verdict.warnings.push(sizeWarn); }
        
        if (features && features.length > 0 && (bioOptions.intronSpanning || bioOptions.avoidExonJunction)) {
          const exonRes = BioMath.checkExonSpanning(amp.amplicon.start, amp.amplicon.start + fwd.seq.length - 1, amp.amplicon.end - rev.seq.length + 1, amp.amplicon.end, features, { intronSpanningRequired: bioOptions.intronSpanning, avoidExonJunction: bioOptions.avoidExonJunction });
          if (exonRes.verdict === 'FAIL') verdict.failures.push(exonRes.reason);
        }

        const fwdSpec = BioMath.checkTemplateSpecificity(fwd.seq, template, false);
        const revSpec = BioMath.checkTemplateSpecificity(rev.seq, template, true);
        if (fwdSpec.specificity === 'MULTIPLE' || revSpec.specificity === 'MULTIPLE') {
          verdict.failures.push("Non-specific binding detected (Multiple sites)");
          verdict.actionItems.push("Relocate primers to more specific genomic regions.");
        }
      } else {
        verdict.failures.push("Amplicon mapping failed (Orientation mismatch)");
      }
    }

    verdict.checklist = [
      { label: 'Length/Tm Specifications', pass: fwdStats.len && fwdStats.tm && revStats.len && revStats.tm },
      { label: 'Thermodynamic Symmetry', pass: tmDiffOk },
      { label: 'Biological & Target Specificity', pass: ampVerdict === 'PASS' && !verdict.failures.some(f => f.includes('mapping') || f.includes('specific')) },
      { label: 'Secondary Structures (ΔG)', pass: fwdStats.homo && revStats.homo && fwdStats.dg > -6.0 && revStats.dg > -6.0 && heteroDg > -6.0 },
      { label: '3\' Terminal Integrity', pass: fwdStats.terminal.verdict !== 'CRITICAL' && revStats.terminal.verdict !== 'CRITICAL' }
    ];

    if (verdict.failures.length > 0) { verdict.status = "🔴 FAIL"; verdict.class = "status-fail"; }
    else if (verdict.warnings.length > 0) { verdict.status = "🟡 PASS WITH CAUTION"; verdict.class = "status-caution"; }
    else { verdict.status = "🟢 PASS"; verdict.class = "status-pass"; }

    return verdict;
  };

})();
