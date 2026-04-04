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
      'F': 2.8,  'L': 3.8,  'V': 4.2,  'I': 4.5
  };

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

})();
