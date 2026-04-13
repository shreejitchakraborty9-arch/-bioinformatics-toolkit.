/**
 * scientific-data.js — High-Fidelity Parameters for Bioinformatics Toolkit
 * 
 * 1. Turner 2004 RNA-RNA Nearest-Neighbor Parameters
 * 2. Sugimoto 1995 RNA-DNA Hybrid Nearest-Neighbor Parameters
 * 3. NCBI Genetic Codes (Transl_table 1-33)
 * 4. Precise Monoresidue Masses (IUPAC)
 */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };
  const data = window.BioKit.data;

  // 1. RNA Molecular Mass Engine Constants (IUPAC High-Fidelity)
  data.RNA_MONORESIDUE_MASS = {
    'A': 329.21,
    'U': 306.17,
    'C': 305.18,
    'G': 345.21,
    // IUPAC Ambiguity Groups (Min, Max, Avg)
    'R': { min: 329.21, max: 345.21, avg: 337.21 }, // A or G
    'Y': { min: 305.18, max: 306.17, avg: 305.675 }, // C or U
    'S': { min: 305.18, max: 345.21, avg: 325.195 }, // G or C
    'W': { min: 306.17, max: 329.21, avg: 317.69 }, // A or U
    'K': { min: 306.17, max: 345.21, avg: 325.69 }, // G or U
    'M': { min: 305.18, max: 329.21, avg: 317.195 }, // A or C
    'B': { min: 305.18, max: 345.21, avg: 318.853 }, // C, G, or U
    'D': { min: 306.17, max: 345.21, avg: 326.863 }, // A, G, or U
    'H': { min: 305.18, max: 329.21, avg: 313.52 }, // A, C, or U
    'V': { min: 305.18, max: 345.21, avg: 326.533 }, // A, C, or G
    'N': { min: 305.18, max: 345.21, avg: 321.4425 } // Any
  };

  data.RNA_TERMINAL_GROUPS = {
    'OH': 18.015,
    'MONO': 79.0,
    'TRI': 159.0
  };

  // 2. Thermodynamic Parameters (Nearest-Neighbor)
  // Turner 2004 RNA-RNA NN Parameters (Mathews et al. 2004)
  // Units: dH (kcal/mol), dS (cal/mol·K)
  data.TURNER_2004 = {
    'AA': { dH: -6.82, dS: -19.0 },
    'UU': { dH: -6.82, dS: -19.0 },
    'AU': { dH: -9.38, dS: -26.7 },
    'UA': { dH: -7.69, dS: -20.5 },
    'CU': { dH: -10.48, dS: -27.1 },
    'AG': { dH: -10.48, dS: -27.1 },
    'CA': { dH: -10.44, dS: -26.9 },
    'UG': { dH: -10.44, dS: -26.9 },
    'GU': { dH: -11.40, dS: -29.5 },
    'AC': { dH: -11.40, dS: -29.5 },
    'GA': { dH: -12.44, dS: -32.5 },
    'UC': { dH: -12.44, dS: -32.5 },
    'CG': { dH: -10.64, dS: -26.7 },
    'GC': { dH: -14.88, dS: -36.9 },
    'GG': { dH: -13.39, dS: -32.7 },
    'CC': { dH: -13.39, dS: -32.7 },
    'init': { dH: 0.0, dS: -10.1 } // Standard RNA initiation penalty
  };

  // Sugimoto 1995 RNA-DNA Hybrid NN Parameters (Sugimoto et al. 1995)
  // Reference: Biochemistry 34:11211-11216
  // Sequence is 5'->3' RNA, matched with 3'<-5' DNA
  data.SUGIMOTO_1995 = {
    'AA': { dH: -7.8, dS: -21.9 }, // rA.dT/rA.dT
    'UU': { dH: -10.4, dS: -29.9 }, // rU.dA/rU.dA
    'AC': { dH: -11.6, dS: -29.9 }, // rA.dT/rC.dG
    'AG': { dH: -10.4, dS: -26.9 }, // rA.dT/rG.dC
    'AU': { dH: -7.2, dS: -18.3 }, // rA.dT/rU.dA
    'CA': { dH: -9.6, dS: -23.5 }, // rC.dG/rA.dT
    'CC': { dH: -14.1, dS: -34.9 }, // rC.dG/rC.dG
    'CG': { dH: -16.3, dS: -47.1 }, // rC.dG/rG.dC
    'CU': { dH: -10.8, dS: -25.2 }, // rC.dG/rU.dA
    'GA': { dH: -14.3, dS: -37.9 }, // rG.dC/rA.dT
    'GC': { dH: -13.3, dS: -31.6 }, // rG.dC/rC.dG
    'GG': { dH: -9.3, dS: -23.2 }, // rG.dC/rG.dC
    'GU': { dH: -13.2, dS: -34.1 }, // rG.dC/rU.dA
    'UA': { dH: -11.5, dS: -36.4 }, // rU.dA/rA.dT
    'UC': { dH: -10.6, dS: -26.7 }, // rU.dA/rC.dG
    'UG': { dH: -13.3, dS: -35.5 }, // rU.dA/rG.dC
    'init': { dH: 1.9, dS: -3.9 }
  };

  // 3. NCBI Genetic Codes (Transl_table 1-33)
  data.GENETIC_CODES = {
    1: { name: "Standard", table: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "---M---------------M---------------M----------------------------" },
    2: { name: "Vertebrate Mitochondrial", table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "--------------------------------MMMM---------------M------------" },
    3: { name: "Yeast Mitochondrial", table: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "----------------------------------M----------------------------" },
    4: { name: "Mold/Protozoan Mitochondrial", table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "--MM---------------M------------MMMM---------------M------------" },
    5: { name: "Invertebrate Mitochondrial", table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "---M----------------------------MMMM---------------M------------" },
    6: { name: "Ciliate Nuclear", table: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    9: { name: "Echinoderm/Flatworm Mitochondrial", table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    10: { name: "Euplotid Nuclear", table: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    11: { name: "Bacterial/Archaeal/Plant Plastid", table: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "---M---------------M------------MMMM---------------M------------" },
    12: { name: "Alternative Yeast Nuclear", table: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-------------------M---------------M----------------------------" },
    13: { name: "Ascidian Mitochondrial", table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    14: { name: "Alternative Flatworm Mitochondrial", table: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    16: { name: "Chlorophycean Mitochondrial", table: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    21: { name: "Trematode Mitochondrial", table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    22: { name: "Scenedesmus obliquus Mitochondrial", table: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    23: { name: "Thraustochytrium Mitochondrial", table: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "--------------------------------M--M---------------M------------" },
    24: { name: "Pterobranchia Mitochondrial", table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", starts: "---M---------------M---------------M----------------------------" },
    25: { name: "Candidate Division SR1/Gracilibacteria", table: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "---M-----------------------------------------------M------------" },
    26: { name: "Pachysolen tannophilus Nuclear", table: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-------------------M---------------M----------------------------" },
    27: { name: "Karyorelict Nuclear", table: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    28: { name: "Condylostoma Nuclear", table: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    29: { name: "Mesodinium Nuclear", table: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    30: { name: "Peritrich Nuclear", table: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    31: { name: "Blastocrithidia Nuclear", table: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    33: { name: "Cephalodiscidae Mitochondrial", table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", starts: "---M---------------M---------------M----------------------------" }
  };

})();
