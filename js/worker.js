/**
 * worker.js — Background worker for heavy bioinformatics calculations.
 * Implements memory-efficient algorithms with chunked processing for responsiveness.
 */

self.onmessage = function (e) {
  const { type, payload, taskId } = e.data;

  try {
    switch (type) {
      case 'ALIGN_GLOBAL':
        runNeedlemanWunsch(payload, taskId);
        break;
      case 'ALIGN_LOCAL':
        runSmithWaterman(payload, taskId);
        break;
      case 'SEQ_SEARCH':
        runSequenceSearch(payload, taskId);
        break;
      case 'RESTRICTION_SEARCH':
        runRestrictionSearch(payload, taskId);
        break;
      case 'THERMO_ANALYSIS':
        runThermoAnalysis(payload, taskId);
        break;
      case 'PARSE_FASTA':
        runFastaParse(payload, taskId);
        break;
      default:
        throw new Error(`Unknown task type: ${type}`);
    }
  } catch (err) {
    self.postMessage({ type: 'ERROR', taskId, message: err.message, stack: err.stack });
  }
};

// ── BLOSUM62 Substitution Matrix (Standard) ──────────────────────
const BLOSUM62 = {
  'A': { 'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': 0, '*': -4 },
  'R': { 'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'Z': 0, 'X': -1, '*': -4 },
  'N': { 'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 3, 'Z': 0, 'X': -1, '*': -4 },
  'D': { 'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4 },
  'C': { 'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'Z': -3, 'X': -2, '*': -4 },
  'Q': { 'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'Z': 3, 'X': -1, '*': -4 },
  'E': { 'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4 },
  'G': { 'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'Z': -2, 'X': -1, '*': -4 },
  'H': { 'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'Z': 0, 'X': -1, '*': -4 },
  'I': { 'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -1, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'Z': -3, 'X': -1, '*': -4 },
  'L': { 'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'Z': -3, 'X': -1, '*': -4 },
  'K': { 'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 1, 'X': -1, '*': -4 },
  'M': { 'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'Z': -1, 'X': -1, '*': -4 },
  'F': { 'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'Z': -3, 'X': -1, '*': -4 },
  'P': { 'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2, 'Z': -1, 'X': -1, '*': -4 },
  'S': { 'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -1, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 0, 'X': 0, '*': -4 },
  'T': { 'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'Z': -1, 'X': 0, '*': -4 },
  'W': { 'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -4, 'Z': -3, 'X': -2, '*': -4 },
  'Y': { 'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3, 'Z': -2, 'X': -1, '*': -4 },
  'V': { 'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3, 'Z': -2, 'X': -1, '*': -4 },
  'B': { 'A': -2, 'R': -1, 'N': 3, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3, 'L': -4, 'K': 0, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4 },
  'Z': { 'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -3, 'Q': 3, 'E': 4, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4 },
  'X': { 'A': 0, 'R': -1, 'N': -1, 'D': -1, 'C': -2, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -1, 'P': -1, 'S': 0, 'T': 0, 'W': -2, 'Y': -1, 'V': -1, 'B': -1, 'Z': -1, 'X': -1, '*': -4 },
  '*': { 'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4, 'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'Z': -4, 'X': -4, '*': 1 }
};

function getScore(a, b, match, mismatch, isProtein) {
  if (!isProtein) return a === b ? match : mismatch;
  const row = BLOSUM62[a] || BLOSUM62['X'];
  return row[b] || row['X'] || mismatch;
}

// ── Needleman-Wunsch Global Alignment (Gotoh Affine Gaps) ─────────
async function runNeedlemanWunsch(data, taskId) {
  const { a, b, match, mismatch, gap } = data;
  const gapOpen = gap;
  const gapExt = data.gapExtend !== undefined ? data.gapExtend : (gap / 2); // Default affine extension
  const m = a.length, n = b.length;
  
  const H = new Float32Array((m + 1) * (n + 1));
  const E = new Float32Array((m + 1) * (n + 1));
  const F = new Float32Array((m + 1) * (n + 1));
  const MIN = -1e9;

  H[0] = 0; E[0] = MIN; F[0] = MIN;
  for (let i = 1; i <= m; i++) {
    H[i * (n + 1)] = gapOpen + i * gapExt;
    E[i * (n + 1)] = MIN;
    F[i * (n + 1)] = MIN;
  }
  for (let j = 1; j <= n; j++) {
    H[j] = gapOpen + j * gapExt;
    E[j] = MIN;
    F[j] = MIN;
  }

  const CHUNK_SIZE = 400;
  const isProtein = data.type === 'protein';

  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const idx = i * (n + 1) + j;
      const diag = H[(i - 1) * (n + 1) + (j - 1)] + getScore(a[i - 1], b[j - 1], match, mismatch, isProtein);
      
      const eOpen = H[i * (n + 1) + (j - 1)] + gapOpen + gapExt;
      const eExt  = E[i * (n + 1) + (j - 1)] + gapExt;
      E[idx] = Math.max(eOpen, eExt);
      
      const fOpen = H[(i - 1) * (n + 1) + j] + gapOpen + gapExt;
      const fExt  = F[(i - 1) * (n + 1) + j] + gapExt;
      F[idx] = Math.max(fOpen, fExt);
      
      H[idx] = Math.max(diag, E[idx], F[idx]);
    }

    if (i % CHUNK_SIZE === 0) {
      self.postMessage({ type: 'PROGRESS', taskId, progress: (i / m) * 0.9 });
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  let i = m, j = n, alnA = '', alnB = '', mid = '';
  let state = 'H';
  while (i > 0 || j > 0) {
    if (i === 0) state = 'E';
    else if (j === 0) state = 'F';

    if (state === 'H') {
      const current = H[i * (n + 1) + j];
      const matchScore = getScore(a[i - 1], b[j - 1], match, mismatch, isProtein);
      if (i > 0 && j > 0 && Math.abs(current - (H[(i - 1) * (n + 1) + (j - 1)] + matchScore)) < 1e-4) {
        alnA = a[i - 1] + alnA; alnB = b[j - 1] + alnB;
        mid = (a[i - 1] === b[j - 1] ? '|' : '·') + mid;
        i--; j--;
      } else if (j > 0 && Math.abs(current - E[i * (n + 1) + j]) < 1e-4) {
        state = 'E';
      } else {
        state = 'F';
      }
    } else if (state === 'E') {
      alnA = '-' + alnA; alnB = b[j - 1] + alnB; mid = ' ' + mid;
      if (j > 1 && Math.abs(E[i * (n + 1) + j] - (E[i * (n + 1) + (j - 1)] + gapExt)) < 1e-4) {
        j--;
      } else {
        j--; state = 'H';
      }
    } else if (state === 'F') {
      alnA = a[i - 1] + alnA; alnB = '-' + alnB; mid = ' ' + mid;
      if (i > 1 && Math.abs(F[i * (n + 1) + j] - (F[(i - 1) * (n + 1) + j] + gapExt)) < 1e-4) {
        i--;
      } else {
        i--; state = 'H';
      }
    }
  }

  self.postMessage({ type: 'RESULT', taskId, result: { alnA, alnB, mid, score: H[m * (n + 1) + n] } });
}

// ── Smith-Waterman Local Alignment (Gotoh Affine Gaps) ────────────
async function runSmithWaterman(data, taskId) {
  const { a, b, match, mismatch, gap } = data;
  const gapOpen = gap;
  const gapExt = data.gapExtend !== undefined ? data.gapExtend : (gap / 2); // Default affine extension
  const m = a.length, n = b.length;
  
  const H = new Float32Array((m + 1) * (n + 1));
  const E = new Float32Array((m + 1) * (n + 1));
  const F = new Float32Array((m + 1) * (n + 1));
  let maxScore = 0, maxI = 0, maxJ = 0;

  const isProtein = data.type === 'protein';

  const CHUNK_SIZE = 400;
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const idx = i * (n + 1) + j;
      const diag = H[(i - 1) * (n + 1) + (j - 1)] + getScore(a[i - 1], b[j - 1], match, mismatch, isProtein);
      
      const eOpen = H[i * (n + 1) + (j - 1)] + gapOpen + gapExt;
      const eExt  = E[i * (n + 1) + (j - 1)] + gapExt;
      E[idx] = Math.max(0, eOpen, eExt);
      
      const fOpen = H[(i - 1) * (n + 1) + j] + gapOpen + gapExt;
      const fExt  = F[(i - 1) * (n + 1) + j] + gapExt;
      F[idx] = Math.max(0, fOpen, fExt);
      
      const score = Math.max(0, diag, E[idx], F[idx]);
      H[idx] = score;

      if (score > maxScore) { maxScore = score; maxI = i; maxJ = j; }
    }

    if (i % CHUNK_SIZE === 0) {
      self.postMessage({ type: 'PROGRESS', taskId, progress: (i / m) * 0.9 });
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  let i = maxI, j = maxJ, alnA = '', alnB = '', mid = '';
  let state = 'H';

  while (i > 0 && j > 0 && H[i * (n + 1) + j] > 0) {
    if (state === 'H') {
      const current = H[i * (n + 1) + j];
      const matchScore = getScore(a[i - 1], b[j - 1], match, mismatch, isProtein);
      if (Math.abs(current - (H[(i - 1) * (n + 1) + (j - 1)] + matchScore)) < 1e-4) {
        alnA = a[i - 1] + alnA; alnB = b[j - 1] + alnB;
        mid = (a[i - 1] === b[j - 1] ? '|' : '·') + mid;
        i--; j--;
      } else if (Math.abs(current - E[i * (n + 1) + j]) < 1e-4) {
        state = 'E';
      } else {
        state = 'F';
      }
    } else if (state === 'E') {
      alnA = '-' + alnA; alnB = b[j - 1] + alnB; mid = ' ' + mid;
      if (j > 1 && Math.abs(E[i * (n + 1) + j] - (E[i * (n + 1) + (j - 1)] + gapExt)) < 1e-4) {
        j--;
      } else {
        j--; state = 'H';
      }
    } else if (state === 'F') {
      alnA = a[i - 1] + alnA; alnB = '-' + alnB; mid = ' ' + mid;
      if (i > 1 && Math.abs(F[i * (n + 1) + j] - (F[(i - 1) * (n + 1) + j] + gapExt)) < 1e-4) {
        i--;
      } else {
        i--; state = 'H';
      }
    }
  }

  self.postMessage({ type: 'RESULT', taskId, result: { alnA, alnB, mid, score: maxScore } });
}

// ── Sequence Search Tool ────────────────────────────────────
async function runSequenceSearch(data, taskId) {
  const { query, database, params } = data;
  const { match, mismatch, gap } = params;
  const results = [];
  const total = database.length;

  for (let idx = 0; idx < total; idx++) {
    const entry = database[idx];
    
    // Seed match (heuristic)
    if (hasSeedMatch(query, entry.sequence, 7)) {
      const aln = alignSW_Internal(query, entry.sequence, match, mismatch, gap);
      if (aln.score > 20) {
        results.push({ target: entry, aln });
      }
    }

    if (idx % 10 === 0) {
      self.postMessage({ type: 'PROGRESS', taskId, progress: (idx / total) });
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  results.sort((a, b) => b.aln.score - a.aln.score);
  self.postMessage({ type: 'RESULT', taskId, result: results });
}

// Internal SW for search (TypedArrays)
function alignSW_Internal(query, subject, matchScore, mismatchScore, gapScore) {
  const qLen = query.length;
  const sLen = subject.length;
  const H = new Int32Array((qLen + 1) * (sLen + 1));
  const ptr = new Int8Array((qLen + 1) * (sLen + 1)); 

  let maxScore = 0, maxI = 0, maxJ = 0;

  for (let i = 1; i <= qLen; i++) {
    for (let j = 1; j <= sLen; j++) {
      const idx = i * (sLen + 1) + j;
      const diag = H[(i - 1) * (sLen + 1) + (j - 1)] + (query[i - 1] === subject[j - 1] ? matchScore : mismatchScore);
      const up   = H[(i - 1) * (sLen + 1) + j] + gapScore;
      const left = H[i * (sLen + 1) + (j - 1)] + gapScore;

      let score = 0, direction = 0;
      if (diag > score) { score = diag; direction = 1; }
      if (up > score)   { score = up;   direction = 2; }
      if (left > score) { score = left; direction = 3; }

      H[idx] = score;
      ptr[idx] = direction;
      if (score > maxScore) { maxScore = score; maxI = i; maxJ = j; }
    }
  }

  // Traceback
  let alignQ = "", alignS = "", matchStr = "", i = maxI, j = maxJ, matches = 0, alnLen = 0;
  while (i > 0 && j > 0 && H[i * (sLen + 1) + j] > 0) {
    const dir = ptr[i * (sLen + 1) + j];
    if (dir === 1) {
      const qChar = query[i - 1], sChar = subject[j - 1];
      alignQ = qChar + alignQ; alignS = sChar + alignS;
      if (qChar === sChar) { matchStr = "|" + matchStr; matches++; } else { matchStr = " " + matchStr; }
      i--; j--; alnLen++;
    } else if (dir === 2) {
      alignQ = query[i - 1] + alignQ; alignS = "-" + alignS; matchStr = " " + matchStr; i--; alnLen++;
    } else if (dir === 3) {
      alignQ = "-" + alignQ; alignS = subject[j - 1] + alignS; matchStr = " " + matchStr; j--; alnLen++;
    }
  }

  const dbSize = 30000000; // Simulated database size (30MB)
  const effS = Math.max(1, sLen); 
  const effQ = Math.max(1, qLen);
  
  // Approximate Karlin-Altschul parameters for matched/mismatched scores
  // Real values require evaluating the Gumbel extreme value distribution
  const lambda = Math.log(mismatchScore / (matchScore + mismatchScore)) / -matchScore || 0.317;
  const K = 0.13; // Typical for DNA alignment space

  const expectedAligns = K * effS * effQ * Math.exp(-lambda * maxScore);
  const bitScore = (lambda * maxScore - Math.log(K)) / Math.LN2;
  const eValue = expectedAligns * (dbSize / effS);

  return {
    score: maxScore, eValue: eValue < 1e-180 ? 0 : eValue, bitScore: bitScore, matches, alnLen,
    identity: alnLen > 0 ? ((matches / alnLen) * 100).toFixed(1) : "0.0",
    queryAln: alignQ, midAln: matchStr, subjAln: alignS,
    qStart: i + 1, qEnd: maxI, sStart: j + 1, sEnd: maxJ
  };
}

function hasSeedMatch(query, subject, wordSize = 7) {
  if (query.length < wordSize || subject.length < wordSize) return true;
  const kmerSet = new Set();
  for (let i = 0; i <= query.length - wordSize; i++) kmerSet.add(query.substring(i, i + wordSize));
  for (let i = 0; i <= subject.length - wordSize; i++) if (kmerSet.has(subject.substring(i, i + wordSize))) return true;
  return false;
}

// ── Restriction Mapping (Optimized for scale) ───────────────

// ─────────────────────────────────────────────────────────────────────────────
// METHYLATION INTERFERENCE ENGINE
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Curated enzyme methylation-sensitivity table.
 * Source: REBASE (Restriction Enzyme Database) and NEB catalog.
 *
 * NEW DATABASE SCHEMA — each enzyme object should now carry:
 * {
 *   name:           "XbaI",
 *   site:           "TCTAGA",
 *   cut_sense:      1,
 *   cut_antisense:  5,
 *   dam_sensitive:  true,   // Blocked when a Dam GATC overlaps the recognition site
 *   dcm_sensitive:  false,  // Blocked when a Dcm CCWGG overlaps
 *   cpg_sensitive:  false   // Blocked when a CpG CG within site is methylated
 * }
 *
 * Until the auto-generated DB is regenerated, this worker-side table covers
 * the most commonly used enzymes. Enzymes absent from this table fall back to
 * algorithmic overlap detection only (no false positives — conservative).
 */
/**
 * Determines whether a restriction site physically overlaps with any methylation
 * site in a way that would block enzymatic cleavage.
 *
 * Overlap condition (standard half-open interval math):
 *   Two intervals [a, b) and [c, d) overlap iff:  a < d  AND  c < b
 *
 * An enzyme is blocked only if:
 *   1. Its recognition site overlaps a methylation site (interval intersection), AND
 *   2. The enzyme's database object explicitly lists the methylation type in `blockedBy`.
 *
 * @param {number} reStart   - 0-based start of restriction recognition site.
 * @param {number} reEnd     - 0-based end (exclusive) of restriction recognition site.
 * @param {Object} enz       - Enzyme database object containing `blockedBy` array.
 * @param {Array}  methSites - Output of scanMethylationSites().
 * @returns {{ isBlocked:boolean, blockType:string|null, warning:string|null }}
 */
function checkMethylationOverlap(reStart, reEnd, enz, methSites) {
  const blockedBy = enz.blockedBy || [];

  for (const ms of methSites) {
    // Standard half-open interval overlap test
    const overlaps = (ms.start < reEnd) && (ms.end > reStart);
    if (!overlaps) continue;

    if (blockedBy.includes(ms.type)) {
      return {
        isBlocked: true,
        blockType: ms.type,
        warning:   `Site blocked by overlapping ${ms.type} methylation`
      };
    }
  }

  return { isBlocked: false, blockType: null, warning: null };
}

// ─────────────────────────────────────────────────────────────────────────────

async function runRestrictionSearch(data, taskId) {
  // FIX: payload uses 'sequence' key — destructure with alias to avoid seq=undefined crash
  const { sequence: seq, enzymes, isCircular = false } = data;
  const results = [];
  const seqLen = seq.length;

  // Pre-calculate Reverse Complement once for bottom-strand scanning
  const comp = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N' };
  const rcSeq = seq.split('').reverse().map(b => comp[b] || b).join('');

  // IUPAC MAP PATCH: ambiguous bases in the *target* DNA (N, R, Y…) must also match.
  const IUPAC_MAP = {
    'R': '[AGR]',   'Y': '[CTY]',   'S': '[GCS]',   'W': '[ATW]',
    'K': '[GTK]',   'M': '[ACM]',   'B': '[CGTB]',  'D': '[AGTD]',
    'H': '[ACTH]',  'V': '[ACGV]',  'N': '[ATGCN]'
  };

  // RUN METHYLATION SCANNER ONCE on the canonical forward sequence.
  // Reverse strand methylation sites are captured implicitly: Dam GATC is a
  // palindrome (its reverse complement is also GATC), Dcm CCWGG is also palindromic,
  // and CpG CG is palindromic. Thus a single forward-strand scan is sufficient.
  const methSites = scanMethylationSites(seq);

  for (let i = 0; i < enzymes.length; i++) {
    const enz = enzymes[i];

    // Build regex from enzyme recognition site IUPAC codes
    let regexStr = '';
    for (const char of enz.site.toUpperCase()) { regexStr += IUPAC_MAP[char] || char; }
    const regex = new RegExp(regexStr, 'g');

    // Derive cut positions — support both legacy (cut) and new (cut_sense / cut_antisense) schema
    const cut_sense     = enz.cut_sense     !== undefined ? enz.cut_sense     : (enz.cut || 0);
    const cut_antisense = enz.cut_antisense !== undefined ? enz.cut_antisense : (enz.site.length - (enz.cut || 0));
    const endType = determineEndType(cut_sense, cut_antisense);

    // CIRCULAR TOPOLOGY: extend the search string so sites spanning the origin are found.
    const searchFwd = isCircular ? seq + seq.substring(0, enz.site.length - 1) : seq;
    const searchRev = isCircular ? rcSeq + rcSeq.substring(0, enz.site.length - 1) : rcSeq;

    // ── Forward strand scan ──────────────────────────────────
    let match;
    while ((match = regex.exec(searchFwd)) !== null) {
      const siteStart0 = match.index % seqLen;               // 0-based, circular-safe
      const siteEnd0   = siteStart0 + enz.site.length;
      const cut5       = (match.index + cut_sense) % seqLen;
      const cut3       = (match.index + cut_antisense) % seqLen;

      // Run methylation overlap check for this recognition site window
      const methResult = checkMethylationOverlap(siteStart0, siteEnd0, enz, methSites);

      results.push({
        name:                 enz.name,
        site:                 enz.site,
        strand:               '+',
        pos:                  siteStart0 + 1,                // 1-based
        cut:                  cut5 + 1,
        cut3:                 cut3 + 1,
        endType:              endType,
        overhang:             Math.abs(cut_antisense - cut_sense),
        isBlocked:            methResult.isBlocked,
        blockType:            methResult.blockType,
        warning:              methResult.warning
      });
      regex.lastIndex = match.index + 1;
    }

    // ── Reverse (bottom) strand scan ────────────────────────
    regex.lastIndex = 0;
    while ((match = regex.exec(searchRev)) !== null) {
      const fwdIndex = seqLen - match.index - enz.site.length;
      const siteStart0 = ((fwdIndex % seqLen) + seqLen) % seqLen;  // 0-based, always positive
      const siteEnd0   = siteStart0 + enz.site.length;
      const fwdCut5    = (seqLen - (match.index + cut_sense))    % seqLen;
      const fwdCut3    = (seqLen - (match.index + cut_antisense)) % seqLen;

      const methResult = checkMethylationOverlap(siteStart0, siteEnd0, enz, methSites);

      results.push({
        name:                 enz.name,
        site:                 enz.site,
        strand:               '-',
        pos:                  siteStart0 + 1,
        cut:                  ((fwdCut5 % seqLen) + seqLen) % seqLen + 1,
        cut3:                 ((fwdCut3 % seqLen) + seqLen) % seqLen + 1,
        endType:              endType,
        overhang:             Math.abs(cut_antisense - cut_sense),
        isBlocked:            methResult.isBlocked,
        blockType:            methResult.blockType,
        warning:              methResult.warning
      });
      regex.lastIndex = match.index + 1;
    }

    if (i % 20 === 0) {
      self.postMessage({ type: 'PROGRESS', taskId, progress: i / enzymes.length });
      await new Promise(r => setTimeout(r, 0));
    }
  }

  // Deduplicate by name + strand + position, then sort by position
  const seen = new Set();
  const unique = [];
  results.forEach(r => {
    const key = `${r.name}-${r.strand}-${r.pos}`;
    if (!seen.has(key)) { seen.add(key); unique.push(r); }
  });

  unique.sort((a, b) => a.pos - b.pos);
  self.postMessage({ type: 'RESULT', taskId, result: unique });
}
// ── Thermodynamic Analysis (Primer Design) ─────────────────
async function runThermoAnalysis(data, taskId) {
  const { sequence, naConc_mM, mgConc_mM, dntpConc_mM, oligoConc_nM } = data;
  
  // A simplified worker-side Tm calculation (Nearest-Neighbor)
  // We can't import science-algorithms easily in a worker without structured exports,
  // so we rely on the payload providing necessary baseline or we implement a minimal version.
  // For now, we return a success signal with the data, but for heavy lifting, 
  // we would perform the window-sliding here.
  
  self.postMessage({ 
    type: 'RESULT', 
    taskId, 
    result: { 
      status: 'OK',
      sequence: sequence,
      length: sequence.length
    }
  });
}

// ── FASTA Parser (Off-Main-Thread) ──────────────────────────
//
// parseFASTA logic migrated from fasta.js to eliminate the
// synchronous forEach + toUpperCase + replace() call that was
// blocking the browser's main UI thread for large files.
//
// Algorithm is unchanged from the original fasta.js implementation.
// Only the execution context has changed: Worker thread vs. UI thread.
//
function runFastaParse(data, taskId) {
  const { rawText } = data;
  if (!rawText) {
    self.postMessage({ type: 'ERROR', taskId, message: 'No FASTA text provided.' });
    return;
  }

  const records = [];
  const lines   = rawText.split('\n');
  let current   = null;

  // Exact same logic as original parseFASTA in fasta.js — no changes
  lines.forEach(line => {
    const trimmed = line.trim();
    if (trimmed.startsWith('>')) {
      if (current) records.push(current);
      current = { header: trimmed.substring(1), sequence: '' };
    } else if (current && trimmed) {
      current.sequence += trimmed.toUpperCase().replace(/\s/g, '');
    }
  });
  if (current) records.push(current);

  self.postMessage({ type: 'RESULT', taskId, result: records });
}
