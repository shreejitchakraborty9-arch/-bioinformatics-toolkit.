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
      default:
        throw new Error(`Unknown task type: ${type}`);
    }
  } catch (err) {
    self.postMessage({ type: 'ERROR', taskId, message: err.message, stack: err.stack });
  }
};

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
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const idx = i * (n + 1) + j;
      const diag = H[(i - 1) * (n + 1) + (j - 1)] + (a[i - 1] === b[j - 1] ? match : mismatch);
      
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
      if (i > 0 && j > 0 && Math.abs(current - (H[(i - 1) * (n + 1) + (j - 1)] + (a[i - 1] === b[j - 1] ? match : mismatch))) < 1e-4) {
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

  const CHUNK_SIZE = 400;
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const idx = i * (n + 1) + j;
      const diag = H[(i - 1) * (n + 1) + (j - 1)] + (a[i - 1] === b[j - 1] ? match : mismatch);
      
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
      if (Math.abs(current - (H[(i - 1) * (n + 1) + (j - 1)] + (a[i - 1] === b[j - 1] ? match : mismatch))) < 1e-4) {
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
