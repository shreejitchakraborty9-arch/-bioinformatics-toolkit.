/* ============================================================
   alignment.js — Pairwise Sequence Alignment (Fixed)
   ============================================================ */

(function () {
  'use strict';

  document.getElementById('alignBtn')?.addEventListener('click', runAlignment);
  document.getElementById('alignSampleBtn')?.addEventListener('click', () => {
    document.getElementById('seqA').value = 'ACGTACGTACGT';
    document.getElementById('seqB').value = 'ACGTCCGTACGT';
  });

  function runAlignment() {
    const seqA = window.cleanSeq(document.getElementById('seqA')?.value || '');
    const seqB = window.cleanSeq(document.getElementById('seqB')?.value || '');
    const algo = document.querySelector('input[name="alignAlgo"]:checked')?.value || 'global';
    const match = parseInt(document.getElementById('matchScore')?.value) || 2;
    const mismatch = parseInt(document.getElementById('mismatchScore')?.value) || -1;
    const gap = parseInt(document.getElementById('gapScore')?.value) || -2;

    if (!seqA || !seqB) {
      window.showToast('Please enter both sequences.');
      return;
    }

    const totalLen = seqA.length + seqB.length;
    const workerType = algo === 'global' ? 'ALIGN_GLOBAL' : 'ALIGN_LOCAL';

    try {
      window.withLoading('panel-alignment', (result) => {
        displayAlignment(result, seqA.length, seqB.length, algo);
      }, totalLen, {
        type: workerType,
        payload: { a: seqA, b: seqB, match, mismatch, gap }
      });
    } catch (error) {
      window.handleError(error, 'Alignment');
    }
  }

  function displayAlignment(result, lenA, lenB, algo) {
    const { alnA, alnB, mid, score } = result;
    const len     = alnA.length;
    const matches = (mid.match(/\|/g) || []).length;
    const gaps    = ((alnA + alnB).match(/-/g) || []).length;
    const identity = len ? ((matches / len) * 100).toFixed(1) : '0';
    const coverage = Math.max(lenA, lenB) ? ((len / Math.max(lenA, lenB)) * 100).toFixed(1) : '0';

    window.buildStatCards('alignStatRow', [
      { val: score,          label: 'Alignment Score' },
      { val: identity + '%', label: 'Identity' },
      { val: matches,        label: 'Matches' },
      { val: gaps,           label: 'Gaps' },
      { val: algo === 'global' ? 'Global' : 'Local', label: 'Algorithm' },
    ]);

    const CHUNK = 60;
    const display = document.getElementById('alignDisplay');
    display.textContent = '';

    for (let start = 0; start < len; start += CHUNK) {
      const sliceA = alnA.slice(start, start + CHUNK);
      const sliceM = mid.slice(start, start + CHUNK);
      const sliceB = alnB.slice(start, start + CHUNK);

      let htmlA = `<span style="color:#8899bb;user-select:none">A </span>`;
      let htmlB = `<span style="color:#8899bb;user-select:none">B </span>`;
      let htmlM = `<span style="color:transparent;user-select:none">  </span>`;

      for (let k = 0; k < sliceA.length; k++) {
        const ca = sliceA[k], cb = sliceB[k];
        const escapedCA = window.escapeHTML(ca);
        const escapedCB = window.escapeHTML(cb);
        if (ca === '-' || cb === '-') {
          htmlA += `<span style="color:#4a5c7a">${escapedCA}</span>`;
          htmlB += `<span style="color:#4a5c7a">${escapedCB}</span>`;
          htmlM += `<span style="color:#4a5c7a"> </span>`;
        } else if (ca === cb) {
          htmlA += `<span style="color:#00d4aa">${escapedCA}</span>`;
          htmlB += `<span style="color:#00d4aa">${escapedCB}</span>`;
          htmlM += `<span style="color:#00d4aabb">|</span>`;
        } else {
          htmlA += `<span style="color:#f43f5e">${escapedCA}</span>`;
          htmlB += `<span style="color:#f43f5e">${escapedCB}</span>`;
          htmlM += `<span style="color:#f43f5e88">·</span>`;
        }
      }

      const block = document.createElement('div');
      block.style.marginBottom = '6px';
      const blockHtml = `
        <div>${htmlA}<span style="color:#4a5c7a;font-size:0.75em;margin-left:8px">${start+1}–${Math.min(start+CHUNK,len)}</span></div>
        <div>${htmlM}</div>
        <div>${htmlB}</div>
      `;
      block.innerHTML = window.sanitizeHTML(blockHtml);
      display.appendChild(block);
    }

    document.getElementById('alignResults').classList.remove('hidden');
  }
})();
