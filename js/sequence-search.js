/* ============================================================
   sequence-search.js — BLAST-like Sequence Search (Fixed DOM IDs)
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  function run() {
    const input = document.getElementById('blastInput');
    if (!input) return;

    const raw = input.value;
    const validation = window.validateSequence(raw, 'dna');

    if (!validation.valid) {
      window.showValidationWarning('panel-sequence-search', validation);
      return;
    }

    const query = validation.clean;
    if (query.length < 15) {
      window.showToast('Query must be at least 15bp.');
      return;
    }

    const database = window.EXAMPLE_GENOMES || [];
    if (database.length === 0) {
      window.showToast('No genome database loaded.');
      return;
    }

    try {
      window.withLoading('panel-sequence-search', (results) => {
        displayResults(query, results);
      }, query.length, {
        type: 'SEQ_SEARCH',
        payload: {
          query: query,
          database: database,
          params: { match: 2, mismatch: -1, gap: -2 }
        }
      });
    } catch (error) {
      window.handleError(error, 'Sequence Search');
    }
  }

  function displayResults(query, results) {
    const container = document.getElementById('blastResults');
    if (!container) return;
    container.classList.remove('hidden');

    window.buildStatCards('blastStatRow', [
      { val: results.length, label: 'Hits Found' },
      { val: query.length + ' bp', label: 'Query Length' },
      { val: (window.EXAMPLE_GENOMES || []).length, label: 'DB Entries' },
    ]);

    const tableWrap = document.getElementById('blastTableWrap');
    if (tableWrap) {
      if (results.length === 0) {
        tableWrap.innerHTML = '<p style="color:var(--text-muted);">No significant hits found.</p>';
      } else {
        let html = `<table class="re-table"><thead><tr>
          <th>Target</th><th>Organism</th><th>Score</th><th>E-value</th><th>Identity</th><th>Align Len</th>
        </tr></thead><tbody>`;
        results.forEach(r => {
          const t = r.target;
          const a = r.aln;
          html += `<tr>
            <td><strong style="color:var(--teal)">${window.escapeHTML(t.gene)}</strong><br><small style="color:var(--text-muted)">${window.escapeHTML(t.description)}</small></td>
            <td><em>${window.escapeHTML(t.organism)}</em></td>
            <td>${a.score}</td>
            <td>${a.eValue < 0.001 ? a.eValue.toExponential(1) : a.eValue.toFixed(3)}</td>
            <td>${a.identity}%</td>
            <td>${a.alnLen}</td>
          </tr>`;
        });
        html += '</tbody></table>';
        tableWrap.innerHTML = window.sanitizeHTML(html);
      }
    }

    // Detailed alignments
    const alnWrap = document.getElementById('blastAlignments');
    if (alnWrap && results.length > 0) {
      let html = '';
      results.slice(0, 5).forEach(r => {
        const t = r.target;
        const a = r.aln;
        html += `
          <div style="background:var(--bg-card);border:1px solid var(--border);border-radius:var(--radius);padding:16px;margin-top:16px;">
            <h4 style="color:var(--teal);margin-bottom:8px;">${window.escapeHTML(t.gene)} — ${window.escapeHTML(t.organism)}</h4>
            <div style="font-size:0.8rem;color:var(--text-muted);margin-bottom:12px;">
              Score: ${a.score} | E-value: ${a.eValue < 0.001 ? a.eValue.toExponential(1) : a.eValue.toFixed(3)} | Identity: ${a.identity}% | Query: ${a.qStart}-${a.qEnd} | Subject: ${a.sStart}-${a.sEnd}
            </div>
            <pre style="font-family:var(--mono);font-size:0.8rem;line-height:1.6;overflow-x:auto;">Query  ${window.escapeHTML(a.queryAln)}
       ${window.escapeHTML(a.midAln)}
Sbjct  ${window.escapeHTML(a.subjAln)}</pre>
          </div>`;
      });
      alnWrap.innerHTML = window.sanitizeHTML(html);
    }
  }

  // Events
  document.getElementById('blastAnalyzeBtn')?.addEventListener('click', window.debounce(run, 300));

  document.getElementById('blastSampleBtn')?.addEventListener('click', () => {
    const input = document.getElementById('blastInput');
    if (input) {
      // Short E. coli lacZ fragment
      input.value = 'ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAAC';
      input.dispatchEvent(new Event('input'));
    }
  });

  document.getElementById('blastClearBtn')?.addEventListener('click', () => {
    const input = document.getElementById('blastInput');
    if (input) input.value = '';
    document.getElementById('blastResults')?.classList.add('hidden');
  });

})();
