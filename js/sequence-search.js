/* ============================================================
   sequence-search.js — Custom SW Search Engine
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  let customDatabase = [];

  // ── FASTA parser ──────────────────────────────────────────
  function parseFasta(text) {
    const entries = [];
    let header = null;
    const lines = [];
    for (const raw of text.split(/\r?\n/)) {
      const line = raw.trim();
      if (line.startsWith('>')) {
        if (header !== null) entries.push({ header, sequence: lines.join('') });
        header = line.slice(1).trim();
        lines.length = 0;
      } else if (header !== null && line.length > 0) {
        lines.push(line);
      }
    }
    if (header !== null && lines.length > 0) entries.push({ header, sequence: lines.join('') });
    return entries;
  }

  // ── Custom DB upload ──────────────────────────────────────
  const uploadInput = document.getElementById('customDbUpload');
  const dbStatus = document.getElementById('dbLoadStatus');

  if (uploadInput) {
    uploadInput.addEventListener('change', (e) => {
      const file = e.target.files[0];
      if (!file) return;

      const reader = new FileReader();
      reader.onload = (ev) => {
        customDatabase = parseFasta(ev.target.result);
        if (dbStatus) {
          dbStatus.textContent = customDatabase.length + ' sequence' +
            (customDatabase.length !== 1 ? 's' : '') + ' loaded.';
          dbStatus.style.color = customDatabase.length > 0
            ? 'var(--teal)' : 'var(--text-muted)';
        }
        if (customDatabase.length === 0) {
          window.showToast('No valid sequences found in the uploaded file.');
        }
      };
      reader.onerror = () => window.showToast('Failed to read file.');
      reader.readAsText(file);
    });
  }

  // ── Scoring parameter helpers ─────────────────────────────
  function getScoringParams() {
    return {
      match:    parseInt(document.getElementById('swMatch')?.value   ?? '2',  10),
      mismatch: parseInt(document.getElementById('swMismatch')?.value ?? '-1', 10),
      gap:      parseInt(document.getElementById('swGap')?.value      ?? '-1', 10),
    };
  }

  // ── Main search ───────────────────────────────────────────
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

    if (customDatabase.length === 0) {
      window.showToast('Please upload a target FASTA database first.');
      return;
    }

    const params = getScoringParams();

    try {
      window.withLoading('panel-sequence-search', (results) => {
        displayResults(query, results);
      }, query.length, {
        type: 'SEQ_SEARCH',
        payload: {
          query:         query,
          customDatabase: customDatabase,
          matchScore:    params.match,
          mismatchScore: params.mismatch,
          gapScore:      params.gap
        }
      });
    } catch (error) {
      window.handleError(error, 'Sequence Search');
    }
  }

  // ── Result rendering ──────────────────────────────────────
  function displayResults(query, results) {
    const container = document.getElementById('blastResults');
    if (!container) return;
    container.classList.remove('hidden');

    window.buildStatCards('blastStatRow', [
      { val: results.length,           label: 'Hits Found' },
      { val: query.length + ' bp',     label: 'Query Length' },
      { val: customDatabase.length,    label: 'DB Entries' },
    ]);

    const tableWrap = document.getElementById('blastTableWrap');
    if (tableWrap) {
      if (results.length === 0) {
        tableWrap.innerHTML = '<p style="color:var(--text-muted);">No significant hits found.</p>';
      } else {
        let html = `<table class="re-table"><thead><tr>
          <th>Target</th><th>Organism</th><th>Score</th><th>Identity</th><th>Align Len</th>
        </tr></thead><tbody>`;
        results.forEach(r => {
          const t = r.target;
          const a = r.aln;
          html += `<tr>
            <td><strong style="color:var(--teal)">${window.escapeHTML(t.gene || t.header || '')}</strong><br>
                <small style="color:var(--text-muted)">${window.escapeHTML(t.description || '')}</small></td>
            <td><em>${window.escapeHTML(t.organism || '—')}</em></td>
            <td>${a.score}</td>
            <td>${a.identity}%</td>
            <td>${a.alnLen}</td>
          </tr>`;
        });
        html += '</tbody></table>';
        tableWrap.innerHTML = window.sanitizeHTML(html);
      }
    }

    const alnWrap = document.getElementById('blastAlignments');
    if (alnWrap && results.length > 0) {
      let html = '';
      results.slice(0, 5).forEach(r => {
        const t = r.target;
        const a = r.aln;
        html += `
          <div style="background:var(--bg-card);border:1px solid var(--border);border-radius:var(--radius);padding:16px;margin-top:16px;">
            <h4 style="color:var(--teal);margin-bottom:8px;">${window.escapeHTML(t.gene || t.header || '')} — ${window.escapeHTML(t.organism || '—')}</h4>
            <div style="font-size:0.8rem;color:var(--text-muted);margin-bottom:12px;">
              Score: ${a.score} | Identity: ${a.identity}% | Query: ${a.qStart}-${a.qEnd} | Subject: ${a.sStart}-${a.sEnd}
            </div>
            <pre style="font-family:var(--mono);font-size:0.8rem;line-height:1.6;overflow-x:auto;">Query  ${window.escapeHTML(a.queryAln)}
       ${window.escapeHTML(a.midAln)}
Sbjct  ${window.escapeHTML(a.subjAln)}</pre>
          </div>`;
      });
      alnWrap.innerHTML = window.sanitizeHTML(html);
    }
  }

  // ── Button events ─────────────────────────────────────────
  document.getElementById('blastAnalyzeBtn')?.addEventListener('click', window.debounce(run, 300));

  document.getElementById('blastClearBtn')?.addEventListener('click', () => {
    const input = document.getElementById('blastInput');
    if (input) input.value = '';
    document.getElementById('blastResults')?.classList.add('hidden');
  });

})();
