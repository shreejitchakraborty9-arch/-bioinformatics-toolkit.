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

// ── Needleman-Wunsch Global Alignment ───────────────────────
async function runNeedlemanWunsch(data, taskId) {
  const { a, b, match, mismatch, gap } = data;
  const m = a.length, n = b.length;
  
  // Memory optimization: TypedArrays
  // H matrix: (m+1) * (n+1)
  const H = new Int32Array((m + 1) * (n + 1));

  // Initialization
  for (let i = 0; i <= m; i++) H[i * (n + 1)] = i * gap;
  for (let j = 0; j <= n; j++) H[j] = j * gap;

  // Computation in chunks
  const CHUNK_SIZE = 500;
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const diag = H[(i - 1) * (n + 1) + (j - 1)] + (a[i - 1] === b[j - 1] ? match : mismatch);
      const up = H[(i - 1) * (n + 1) + j] + gap;
      const left = H[i * (n + 1) + (j - 1)] + gap;
      H[i * (n + 1) + j] = Math.max(diag, up, left);
    }

    // Yield back to event loop every CHUNK_SIZE rows
    if (i % CHUNK_SIZE === 0) {
      self.postMessage({ type: 'PROGRESS', taskId, progress: (i / m) * 0.9 });
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  // Traceback
  let i = m, j = n, alnA = '', alnB = '', mid = '';
  while (i > 0 || j > 0) {
    const current = H[i * (n + 1) + j];
    if (i > 0 && j > 0 && current === H[(i - 1) * (n + 1) + (j - 1)] + (a[i - 1] === b[j - 1] ? match : mismatch)) {
      alnA = a[i - 1] + alnA; alnB = b[j - 1] + alnB;
      mid = (a[i - 1] === b[j - 1] ? '|' : '·') + mid;
      i--; j--;
    } else if (i > 0 && current === H[(i - 1) * (n + 1) + j] + gap) {
      alnA = a[i - 1] + alnA; alnB = '-' + alnB; mid = ' ' + mid; i--;
    } else {
      alnA = '-' + alnA; alnB = b[j - 1] + alnB; mid = ' ' + mid; j--;
    }
  }

  const result = { alnA, alnB, mid, score: H[m * (n + 1) + n] };
  
  // Cleanup
  self.postMessage({ type: 'RESULT', taskId, result });
}

// ── Smith-Waterman Local Alignment ──────────────────────────
async function runSmithWaterman(data, taskId) {
  const { a, b, match, mismatch, gap } = data;
  const m = a.length, n = b.length;
  const H = new Int32Array((m + 1) * (n + 1));
  let maxScore = 0, maxI = 0, maxJ = 0;

  const CHUNK_SIZE = 500;
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const diag = H[(i - 1) * (n + 1) + (j - 1)] + (a[i - 1] === b[j - 1] ? match : mismatch);
      const up = H[(i - 1) * (n + 1) + j] + gap;
      const left = H[i * (n + 1) + (j - 1)] + gap;
      const score = Math.max(0, diag, up, left);
      H[i * (n + 1) + j] = score;
      if (score > maxScore) { maxScore = score; maxI = i; maxJ = j; }
    }

    if (i % CHUNK_SIZE === 0) {
      self.postMessage({ type: 'PROGRESS', taskId, progress: (i / m) * 0.9 });
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  // Traceback from max
  let i = maxI, j = maxJ, alnA = '', alnB = '', mid = '';
  while (i > 0 && j > 0 && H[i * (n + 1) + j] > 0) {
    const current = H[i * (n + 1) + j];
    if (current === H[(i - 1) * (n + 1) + (j - 1)] + (a[i - 1] === b[j - 1] ? match : mismatch)) {
      alnA = a[i - 1] + alnA; alnB = b[j - 1] + alnB;
      mid = (a[i - 1] === b[j - 1] ? '|' : '·') + mid;
      i--; j--;
    } else if (current === H[(i - 1) * (n + 1) + j] + gap) {
      alnA = a[i - 1] + alnA; alnB = '-' + alnB; mid = ' ' + mid; i--;
    } else {
      alnA = '-' + alnA; alnB = b[j - 1] + alnB; mid = ' ' + mid; j--;
    }
  }

  const result = { alnA, alnB, mid, score: maxScore };
  self.postMessage({ type: 'RESULT', taskId, result });
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

  const K = 0.1, lambda = 0.3;
  const eValueRaw = K * qLen * sLen * Math.exp(-lambda * maxScore);
  const eValue = eValueRaw < 1e-180 ? 0 : eValueRaw;

  return {
    score: maxScore, eValue, matches, alnLen,
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
