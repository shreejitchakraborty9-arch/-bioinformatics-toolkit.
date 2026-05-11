/* ============================================================
   literature.js — PubMed Literature Search Module
   Queries /api/pubmed/search (NCBI ESearch + ESummary + EFetch proxy)
   Returns real article metadata from PubMed: no fabricated data
   ============================================================ */

(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── State ──────────────────────────────────────────────────
  let _lastResults = [];
  let _isSearching = false;
  let _currentQuery = "";
  let _retstart = 0;
  let _totalHits = 0;

  // ── DOM References ─────────────────────────────────────────
  function el(id) { return document.getElementById(id); }

  // ── Search ─────────────────────────────────────────────────
  async function runSearch(isLoadMore = false) {
    const input = el('litSearchInput');
    const btn = el('litSearchBtn');
    const resultsWrap = el('litResults');
    const statusEl = el('litStatus');
    const filterDate = el('litFilterDate');
    const filterType = el('litFilterType');

    if (!input || _isSearching) return;

    let baseQuery = input.value.trim();
    if (!baseQuery) {
      if (window.showToast) window.showToast('Enter a search term (e.g. BRCA1 methylation).');
      return;
    }

    if (!isLoadMore) {
      _retstart = 0;
      _lastResults = [];
      if (resultsWrap) resultsWrap.innerHTML = '';
      
      let queryParts = [baseQuery];
      if (filterDate && filterDate.value) {
        queryParts.push(`"last ${filterDate.value} years"[dp]`);
      }
      if (filterType && filterType.value) {
        queryParts.push(`${filterType.value}[ptyp]`);
      }
      _currentQuery = queryParts.join(' AND ');
    }

    _isSearching = true;
    if (btn) { btn.disabled = true; btn.textContent = 'Searching…'; }
    if (statusEl) statusEl.textContent = 'Querying PubMed…';

    const loadMoreBtn = el('litLoadMoreBtn');
    if (loadMoreBtn) {
      loadMoreBtn.disabled = true;
      loadMoreBtn.textContent = 'Loading...';
    }

    try {
      const resp = await fetch('/api/pubmed/search', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query: _currentQuery, retmax: 20, retstart: _retstart })
      });

      if (!resp.ok) {
        const errData = await resp.json().catch(() => ({}));
        throw new Error(errData.error || `Server returned HTTP ${resp.status}`);
      }

      const data = await resp.json();
      const newArticles = data.articles || [];
      _totalHits = data.total || 0;
      
      if (isLoadMore) {
        _lastResults = _lastResults.concat(newArticles);
      } else {
        _lastResults = newArticles;
      }
      
      _retstart += newArticles.length;

      if (statusEl) {
        statusEl.textContent = _lastResults.length > 0
          ? `${_totalHits.toLocaleString()} results found · showing ${_lastResults.length}`
          : 'No results found for this query.';
      }

      renderResults(_lastResults);

    } catch (err) {
      const msg = err.message || 'Search failed.';
      if (statusEl && !isLoadMore) statusEl.textContent = 'Error — ' + msg;
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
      
      const abstractHtml = a.abstract 
        ? `<details style="margin-top:10px; font-size:0.8rem; color:var(--text-secondary);">
             <summary style="cursor:pointer; color:var(--teal); font-weight:600; margin-bottom:6px;">Abstract</summary>
             <div style="line-height:1.5; padding:8px 12px; background:rgba(0,0,0,0.1); border-radius:4px;">${escape(a.abstract).replace(/\n\n/g, '<br><br>')}</div>
           </details>`
        : '';

      html += `
        <div class="lit-article" data-index="${i}">
          <div class="lit-article-num">${i + 1}</div>
          <div class="lit-article-body">
            <div class="lit-article-title">${escape(a.title)}</div>
            <div class="lit-article-authors">${authorsShort || '<span style="color:var(--text-muted)">No authors listed</span>'}</div>
            <div class="lit-article-journal">${journalInfo}</div>
            ${abstractHtml}
            <div class="lit-article-links" style="margin-top:12px;">
              ${pubmedLink}
              ${doiLink}
              <button class="lit-cite-btn" data-index="${i}" title="Copy APA citation">📋 Cite</button>
              <button class="lit-cite-btn lit-ris-btn" data-index="${i}" title="Download RIS format for Zotero/EndNote">💾 RIS</button>
            </div>
          </div>
        </div>`;
    });
    html += '</div>';

    if (_retstart < _totalHits) {
      html += `
        <div style="margin-top:16px; text-align:center;">
          <button id="litLoadMoreBtn" class="btn-secondary">Load More</button>
        </div>
      `;
    }

    const sanitize = window.sanitizeHTML || (s => s);
    wrap.innerHTML = sanitize(html);

    // Wire cite buttons
    wrap.querySelectorAll('.lit-cite-btn:not(.lit-ris-btn)').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const idx = parseInt(e.target.dataset.index);
        copyAPA(_lastResults[idx]);
      });
    });
    
    // Wire RIS buttons
    wrap.querySelectorAll('.lit-ris-btn').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const idx = parseInt(e.target.dataset.index);
        downloadRIS(_lastResults[idx]);
      });
    });

    // Wire Load More button
    const loadMoreBtn = el('litLoadMoreBtn');
    if (loadMoreBtn) {
      loadMoreBtn.addEventListener('click', () => runSearch(true));
    }
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

  // ── RIS Export ─────────────────────────────────────────────
  function downloadRIS(article) {
    if (!article) return;
    
    const yearMatch = (article.pubdate || '').match(/(\d{4})/);
    const year = yearMatch ? yearMatch[1] : '';

    let ris = `TY  - JOUR\n`;
    ris += `TI  - ${article.title}\n`;
    article.authors.forEach(a => {
      ris += `AU  - ${a}\n`;
    });
    if (article.journal) ris += `T2  - ${article.journal}\n`;
    if (year) ris += `PY  - ${year}\n`;
    if (article.volume) ris += `VL  - ${article.volume}\n`;
    if (article.issue) ris += `IS  - ${article.issue}\n`;
    if (article.pages) {
      const p = article.pages.split('-');
      if (p[0]) ris += `SP  - ${p[0]}\n`;
      if (p[1]) ris += `EP  - ${p[1]}\n`;
    }
    if (article.doi) ris += `DO  - ${article.doi}\n`;
    if (article.pmid) ris += `AN  - ${article.pmid}\n`;
    if (article.abstract) ris += `AB  - ${article.abstract.replace(/\n/g, ' ')}\n`;
    ris += `ER  - \n`;

    const blob = new Blob([ris], { type: 'application/x-research-info-systems' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `PMID_${article.pmid}.ris`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
    
    if (window.showToast) window.showToast('💾 RIS citation downloaded.');
  }

  // ── Event Wiring ──────────────────────────────────────────
  // ── Event Wiring ──────────────────────────────────────────
  function wireEvents() {
    el('litSearchBtn')?.addEventListener('click', window.debounce ? window.debounce(() => runSearch(false), 300) : () => runSearch(false));

    el('litSearchInput')?.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') runSearch(false);
    });

    el('litClearBtn')?.addEventListener('click', () => {
      const input = el('litSearchInput');
      if (input) input.value = '';
      const filterDate = el('litFilterDate');
      if (filterDate) filterDate.value = '';
      const filterType = el('litFilterType');
      if (filterType) filterType.value = '';
      const wrap = el('litResults');
      if (wrap) wrap.innerHTML = '';
      const status = el('litStatus');
      if (status) status.textContent = '';
      _lastResults = [];
      _retstart = 0;
      _currentQuery = "";
      _totalHits = 0;
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', wireEvents);
  } else {
    wireEvents();
  }

})();
