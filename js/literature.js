/* ============================================================
   literature.js — PubMed Literature Search Module
   Queries /api/pubmed/search (NCBI ESearch + ESummary proxy)
   Returns real article metadata from PubMed: no fabricated data
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── State ──────────────────────────────────────────────────
  let _lastResults = [];
  let _isSearching = false;

  // ── DOM References ─────────────────────────────────────────
  function el(id) { return document.getElementById(id); }

  // ── Search ─────────────────────────────────────────────────
  async function runSearch() {
    const input = el('litSearchInput');
    const btn = el('litSearchBtn');
    const resultsWrap = el('litResults');
    const statusEl = el('litStatus');

    if (!input || _isSearching) return;

    const query = input.value.trim();
    if (!query) {
      if (window.showToast) window.showToast('Enter a search term (e.g. BRCA1 methylation).');
      return;
    }

    _isSearching = true;
    if (btn) { btn.disabled = true; btn.textContent = 'Searching…'; }
    if (statusEl) statusEl.textContent = 'Querying PubMed…';
    if (resultsWrap) resultsWrap.innerHTML = '';

    try {
      const resp = await fetch('/api/pubmed/search', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query, retmax: 20 })
      });

      if (!resp.ok) {
        const errData = await resp.json().catch(() => ({}));
        throw new Error(errData.error || `Server returned HTTP ${resp.status}`);
      }

      const data = await resp.json();
      _lastResults = data.articles || [];

      if (statusEl) {
        statusEl.textContent = _lastResults.length > 0
          ? `${data.total.toLocaleString()} results found · showing top ${_lastResults.length}`
          : 'No results found for this query.';
      }

      renderResults(_lastResults);

    } catch (err) {
      const msg = err.message || 'Search failed.';
      if (statusEl) statusEl.textContent = 'Error — ' + msg;
      if (window.showToast) window.showToast('⚠ ' + msg);
      console.error('[Literature]', err);
    } finally {
      _isSearching = false;
      if (btn) { btn.disabled = false; btn.textContent = 'Search PubMed'; }
    }
  }

  // ── Render ─────────────────────────────────────────────────
  function renderResults(articles) {
    const wrap = el('litResults');
    if (!wrap) return;

    if (articles.length === 0) {
      wrap.innerHTML = '<p style="color:var(--text-muted);padding:24px 0;">No articles matched your query.</p>';
      return;
    }

    const escape = window.escapeHTML || (s => s.replace(/[&<>"']/g, c => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
    }[c])));

    let html = '<div class="lit-article-list">';
    articles.forEach((a, i) => {
      const authorsShort = a.authors.length > 3
        ? escape(a.authors.slice(0, 3).join(', ')) + ' <em>et al.</em>'
        : escape(a.authors.join(', '));

      const doiLink = a.doi
        ? `<a href="https://doi.org/${escape(a.doi)}" target="_blank" rel="noopener" class="lit-doi">DOI ↗</a>`
        : '';

      const pubmedLink = `<a href="https://pubmed.ncbi.nlm.nih.gov/${escape(a.pmid)}/" target="_blank" rel="noopener" class="lit-pmid">PMID: ${escape(a.pmid)}</a>`;

      const journalInfo = [
        a.journal ? `<em>${escape(a.journal)}</em>` : '',
        a.volume ? `<strong>${escape(a.volume)}</strong>` : '',
        a.issue ? `(${escape(a.issue)})` : '',
        a.pages ? `:${escape(a.pages)}` : '',
        a.pubdate ? `(${escape(a.pubdate)})` : ''
      ].filter(Boolean).join(' ');

      html += `
        <div class="lit-article" data-index="${i}">
          <div class="lit-article-num">${i + 1}</div>
          <div class="lit-article-body">
            <div class="lit-article-title">${escape(a.title)}</div>
            <div class="lit-article-authors">${authorsShort || '<span style="color:var(--text-muted)">No authors listed</span>'}</div>
            <div class="lit-article-journal">${journalInfo}</div>
            <div class="lit-article-links">
              ${pubmedLink}
              ${doiLink}
              <button class="lit-cite-btn" data-index="${i}" title="Copy APA citation">📋 Cite</button>
            </div>
          </div>
        </div>`;
    });
    html += '</div>';

    const sanitize = window.sanitizeHTML || (s => s);
    wrap.innerHTML = sanitize(html);

    // Wire cite buttons
    wrap.querySelectorAll('.lit-cite-btn').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const idx = parseInt(e.target.dataset.index);
        copyAPA(_lastResults[idx]);
      });
    });
  }

  // ── APA Citation Generator ─────────────────────────────────
  function copyAPA(article) {
    if (!article) return;

    const authStr = article.authors.length > 0
      ? article.authors.slice(0, 6).join(', ') + (article.authors.length > 6 ? ', ... ' : '')
      : 'Unknown';

    // Extract year from pubdate (e.g. "2023 Jan 15" → "2023")
    const yearMatch = (article.pubdate || '').match(/(\d{4})/);
    const year = yearMatch ? yearMatch[1] : 'n.d.';

    let citation = `${authStr} (${year}). ${article.title} `;
    if (article.journal) citation += `${article.journal}`;
    if (article.volume) citation += `, ${article.volume}`;
    if (article.issue) citation += `(${article.issue})`;
    if (article.pages) citation += `, ${article.pages}`;
    citation += '.';
    if (article.doi) citation += ` https://doi.org/${article.doi}`;

    navigator.clipboard.writeText(citation).then(() => {
      if (window.showToast) window.showToast('📋 APA citation copied to clipboard.');
    }).catch(() => {
      if (window.showToast) window.showToast('⚠ Could not copy — check browser permissions.');
    });
  }

  // ── Event Wiring ──────────────────────────────────────────
  document.addEventListener('DOMContentLoaded', () => {
    el('litSearchBtn')?.addEventListener('click', window.debounce ? window.debounce(runSearch, 300) : runSearch);

    el('litSearchInput')?.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') runSearch();
    });

    el('litClearBtn')?.addEventListener('click', () => {
      const input = el('litSearchInput');
      if (input) input.value = '';
      const wrap = el('litResults');
      if (wrap) wrap.innerHTML = '';
      const status = el('litStatus');
      if (status) status.textContent = '';
      _lastResults = [];
    });
  });

})();
