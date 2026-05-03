/* ============================================================
 batch.js — High-Throughput Omics Batch Queue
 File upload → async Celery job → short-poll /api/status →
 progress bar → in-app results table → CSV / TSV download
 ============================================================ */

(function () {
  'use strict';

  const POLL_INTERVAL_MS  = 2000;
  const HISTORY_KEY       = 'bioToolkit_batchHistory';
  const HISTORY_MAX       = 20;

  let currentJobId    = null;
  let pollTimer       = null;
  let uploadedContent = null;  // raw FASTA string from file or paste

  // Full dataset fetched from /api/batch/<id>/results after completion
  let _allRows        = [];
  let _columns        = [];
  let _visibleCols    = new Set();
  let _filteredRows   = [];

  // Table state
  let _sortCol        = null;
  let _sortDir        = 'asc';   // 'asc' | 'desc'
  let _currentPage    = 1;
  let _pageSize       = 50;

  // ── Column metadata (tooltips, Phase 6) ─────────────────────
  const COL_META = {
    id:                       { label: 'ID',            tip: 'FASTA record identifier' },
    description:              { label: 'Description',   tip: 'Full FASTA header line' },
    length:                   { label: 'Length',        tip: 'Sequence length in nt (nucleotide) or aa (protein)' },
    type:                     { label: 'Type',          tip: 'Auto-detected sequence type: DNA / RNA / PROTEIN' },
    gc_percent:               { label: 'GC %',          tip: 'G+C content (%) — canonical bases only\nDegenerate IUPAC codes are excluded' },
    at_percent:               { label: 'AT %',          tip: 'A+T/U content (%) — canonical bases only' },
    a_percent:                { label: 'A %',           tip: 'Adenine fraction (%)' },
    t_percent:                { label: 'T/U %',         tip: 'Thymine (DNA) or Uracil (RNA) fraction (%)' },
    g_percent:                { label: 'G %',           tip: 'Guanine fraction (%)' },
    c_percent:                { label: 'C %',           tip: 'Cytosine fraction (%)' },
    cpg_obs_exp:              { label: 'CpG O/E',       tip: 'CpG observed/expected ratio\n(Gardiner-Garden & Frommer, 1987)\n≥ 0.60 indicates a CpG island (promoter-associated)' },
    melting_temp_c:           { label: 'Tm (°C)',       tip: 'Melting temperature — Tm_Wallace formula:\n4·(G+C) + 2·(A+T)\nCalibrated for oligonucleotides; estimate only for > 50 nt' },
    shannon_entropy:          { label: 'Entropy',       tip: 'Normalized Shannon entropy (0–1)\nMeasures sequence complexity\n< 0.70 = low complexity (DUST/SEG region)' },
    molecular_weight_kda:     { label: 'MW (kDa)',      tip: 'Molecular weight in kDa\n(BioPython molecular_weight)' },
    isoelectric_point:        { label: 'pI',            tip: 'Isoelectric point\n(BioPython ProtParam)' },
    gravy:                    { label: 'GRAVY',         tip: 'Grand average of hydropathicity\n(Kyte & Doolittle, 1982)\n> 0 = hydrophobic; < 0 = hydrophilic' },
    instability_index:        { label: 'Instability',   tip: 'Instability index (Guruprasad, 1990)\n< 40 = stable in vitro; ≥ 40 = unstable' },
    instability_label:        { label: 'Stability',     tip: 'Stable (< 40) or Unstable (≥ 40)' },
    aromaticity:              { label: 'Aromaticity',   tip: 'Fraction of aromatic residues (F + W + Y)\n(BioPython ProtParam)' },
    aliphatic_index:          { label: 'Aliphatic',     tip: 'Aliphatic index (Ikai, 1980)\nSurrogate for thermal stability\n(relative volume of aliphatic side chains)' },
    net_charge_ph74:          { label: 'Charge pH7.4',  tip: 'Net charge at physiological pH 7.4\n(BioPython charge_at_pH)' },
    extinction_coeff_reduced: { label: 'ε red',         tip: 'Extinction coefficient (M⁻¹cm⁻¹) at 280 nm\nAll cysteines free (reduced)' },
    extinction_coeff_oxidized:{ label: 'ε ox',          tip: 'Extinction coefficient (M⁻¹cm⁻¹) at 280 nm\nAll cysteines disulfide-bridged (oxidized)' },
    warning:                  { label: 'Warning',       tip: 'Non-fatal analysis notes (duplicate IDs, Tm caveats, removed residues…)' },
    error:                    { label: 'Error',         tip: 'Per-record exception — record failed to compute one or more metrics' },
  };

  // Default column visibility: hide verbose / narrow-use columns on first load
  const COLS_HIDDEN_DEFAULT = new Set([
    'description', 'at_percent', 'a_percent', 't_percent', 'g_percent', 'c_percent',
    'extinction_coeff_reduced', 'extinction_coeff_oxidized',
  ]);

  // ── File / Drop-Zone Wiring ──────────────────────────────────
  const dropZone  = document.getElementById('batchDropZone');
  const fileInput = document.getElementById('batchFileInput');

  if (dropZone) {
    dropZone.addEventListener('dragover', (e) => {
      e.preventDefault();
      dropZone.classList.add('dragover');
    });
    dropZone.addEventListener('dragleave', () => dropZone.classList.remove('dragover'));
    dropZone.addEventListener('drop', (e) => {
      e.preventDefault();
      dropZone.classList.remove('dragover');
      const file = e.dataTransfer?.files[0];
      if (file) handleFileSelect(file);
    });
  }

  fileInput?.addEventListener('change', () => {
    if (fileInput.files[0]) handleFileSelect(fileInput.files[0]);
  });

  function handleFileSelect(file) {
    const MAX_BYTES = 5 * 1024 * 1024; // 5 MB
    if (file.size > MAX_BYTES) {
      window.showToast?.('File too large. Maximum is 5 MB.');
      return;
    }
    const reader = new FileReader();
    reader.onload = (e) => {
      uploadedContent = e.target.result;
      const pasteArea = document.getElementById('batchPasteInput');
      if (pasteArea) pasteArea.value = uploadedContent;
      updateSeqCount(uploadedContent);
      window.showToast?.(`Loaded: ${file.name}`);
    };
    reader.readAsText(file);
  }

  // Live count as the user types in the paste area (Phase 1 client-side count)
  document.getElementById('batchPasteInput')?.addEventListener('input', (e) => {
    uploadedContent = e.target.value;
    updateSeqCount(uploadedContent);
  });

  function updateSeqCount(text) {
    const count = text ? text.split('\n').filter(l => l.trimStart().startsWith('>')).length : 0;
    const el    = document.getElementById('batchSeqCount');
    if (!el) return;
    if (count > 0) {
      el.textContent = `${count.toLocaleString()} sequence${count === 1 ? '' : 's'} detected`;
      el.style.display = '';
    } else {
      el.style.display = 'none';
    }
  }

  // ── Submit ───────────────────────────────────────────────────
  document.getElementById('batchSubmitBtn')?.addEventListener('click', submitBatch);

  async function submitBatch() {
    const pasteArea = document.getElementById('batchPasteInput');
    const content   = (uploadedContent || pasteArea?.value || '').trim();

    if (!content) {
      window.showToast?.('Upload a FASTA file or paste sequences first.');
      return;
    }

    // Client-side pre-flight: require at least one column-0 header line
    const headerLines = content.split('\n').filter(l => l.startsWith('>'));
    if (headerLines.length === 0) {
      window.showToast?.('No FASTA headers found (lines starting with >). Check your input.');
      return;
    }

    stopPolling();
    resetJobUI();
    setBatchStatus('Submitting job…', false);

    const mode = document.getElementById('batchModeSelect')?.value || 'all';

    try {
      const formData = new FormData();
      const blob     = new Blob([content], { type: 'text/plain' });
      formData.append('file', blob, 'batch.fasta');
      formData.append('mode', mode);

      const res = await fetch('/api/analyze/batch', { method: 'POST', body: formData });

      if (res.status === 503) {
        showBatchDegradedPanel(content, headerLines.length);
        return;
      }

      if (!res.ok) {
        const err = await res.json().catch(() => ({ error: `HTTP ${res.status}` }));
        throw new Error(err.error || `Server error ${res.status}`);
      }

      const data    = await res.json();
      currentJobId  = data.job_id;

      saveJobToHistory(currentJobId, { mode, seqs: headerLines.length, status: 'running' });
      showJobArea(currentJobId);
      startPolling(currentJobId);
    } catch (err) {
      setBatchStatus('Submission failed: ' + err.message, true);
      window.showToast?.('Submission failed: ' + err.message);
    }
  }

  // ── Batch 503 Degraded-Mode Panel ───────────────────────────
  const INLINE_THRESHOLD = 20; // max sequences for the inline fallback

  function showBatchDegradedPanel(fastaContent, seqCount) {
    resetJobUI();
    setBatchStatus('Batch worker unavailable', true);

    let panel = document.getElementById('batchDegradedPanel');
    if (!panel) {
      panel = document.createElement('div');
      panel.id = 'batchDegradedPanel';
      panel.style.cssText = [
        'margin-top:16px',
        'padding:16px 20px',
        'border:1px solid var(--warning-border,#f59e0b)',
        'border-radius:8px',
        'background:var(--warning-bg,rgba(245,158,11,0.08))',
        'color:var(--text-primary)',
        'font-size:0.9rem',
        'line-height:1.6',
      ].join(';');

      const jobArea = document.getElementById('batchJobArea');
      if (jobArea && jobArea.parentNode) {
        jobArea.parentNode.insertBefore(panel, jobArea.nextSibling);
      }
    }

    const canRunInline = seqCount > 0 && seqCount <= INLINE_THRESHOLD;

    panel.innerHTML = window.BioKit.utils.sanitizeHTML(`
      <div style="display:flex;align-items:flex-start;gap:12px;">
        <span style="font-size:1.4rem;line-height:1;flex-shrink:0;">⚠</span>
        <div style="flex:1;">
          <strong style="font-size:1rem;">Batch Processing Requires Background Infrastructure</strong>
          <p style="margin:8px 0 0;">The batch queue runs on <strong>Celery</strong> (task worker) + <strong>Redis</strong> (message broker).
          These services are not currently running on this server.</p>
          <details style="margin-top:10px;">
            <summary style="cursor:pointer;font-weight:600;color:var(--teal);">Setup Instructions</summary>
            <div style="margin-top:8px;padding:10px 12px;background:var(--bg-surface,#1e293b);border-radius:6px;font-family:var(--mono);font-size:0.82rem;line-height:1.8;">
              <div><span style="color:var(--text-muted)"># 1. Install dependencies</span></div>
              <div>pip install celery redis</div>
              <div style="margin-top:6px;"><span style="color:var(--text-muted)"># 2. Start Redis (Docker)</span></div>
              <div>docker run -d -p 6379:6379 redis:7-alpine</div>
              <div style="margin-top:6px;"><span style="color:var(--text-muted)"># 3. Start the Celery worker</span></div>
              <div>celery -A worker worker --loglevel=info</div>
              <div style="margin-top:6px;"><span style="color:var(--text-muted)"># 4. Restart the Flask app</span></div>
              <div>python app.py</div>
            </div>
          </details>
          ${canRunInline ? `
          <div style="margin-top:14px;padding-top:12px;border-top:1px solid var(--border-muted,#334155);">
            <strong>Inline Fallback Available</strong>
            <p style="margin:4px 0 8px;color:var(--text-secondary);">
              Your input has ${seqCount} sequence${seqCount !== 1 ? 's' : ''} (≤ ${INLINE_THRESHOLD}).
              Basic metrics (length, GC%) can be computed inline without the worker.
            </p>
            <button id="batchRunInlineBtn" class="btn-primary" style="font-size:0.85rem;padding:6px 16px;">
              Run Inline (basic metrics)
            </button>
          </div>` : `
          <p style="margin-top:10px;color:var(--text-muted);">
            Your input has ${seqCount} sequence${seqCount !== 1 ? 's' : ''}.
            Inline fallback is available for ≤ ${INLINE_THRESHOLD} sequences.
          </p>`}
        </div>
      </div>
    `);

    panel.style.display = '';

    if (canRunInline) {
      document.getElementById('batchRunInlineBtn')?.addEventListener('click', () => {
        runInlineFallback(fastaContent, panel);
      });
    }
  }

  function _parseFasta(text) {
    const records = [];
    let id = null, desc = '', lines = [];
    text.split('\n').forEach(raw => {
      const line = raw.trim();
      if (line.startsWith('>')) {
        if (id !== null) records.push({ id, desc, seq: lines.join('') });
        const header = line.slice(1);
        const space  = header.indexOf(' ');
        id    = space >= 0 ? header.slice(0, space) : header;
        desc  = space >= 0 ? header.slice(space + 1) : '';
        lines = [];
      } else if (line) {
        lines.push(line);
      }
    });
    if (id !== null) records.push({ id, desc, seq: lines.join('') });
    return records;
  }

  async function runInlineFallback(fastaContent, panel) {
    const btn = document.getElementById('batchRunInlineBtn');
    if (btn) { btn.disabled = true; btn.textContent = 'Analyzing…'; }

    const records = _parseFasta(fastaContent);
    const results = [];

    for (const rec of records) {
      try {
        const res = await fetch('/api/analyze/basic', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ sequence: rec.seq.toUpperCase() })
        });
        const data = await res.json().catch(() => ({}));
        results.push({
          id:         rec.id,
          length:     data.length     ?? rec.seq.length,
          gc_percent: data.gc_content != null ? data.gc_content.toFixed(1) : '—',
          error:      data.error || null
        });
      } catch (err) {
        results.push({ id: rec.id, length: rec.seq.length, gc_percent: '—', error: err.message });
      }
    }

    // Render inline results table
    let tableHtml = `
      <div style="margin-top:14px;padding-top:12px;border-top:1px solid var(--border-muted,#334155);">
        <strong>Inline Results (${results.length} sequence${results.length !== 1 ? 's' : ''})</strong>
        <div style="overflow-x:auto;margin-top:8px;">
          <table class="data-table" style="font-size:0.85rem;">
            <thead><tr>
              <th>ID</th><th>Length (nt/aa)</th><th>GC %</th><th>Note</th>
            </tr></thead>
            <tbody>`;
    results.forEach(r => {
      const noteStyle = r.error ? 'color:#f87171;' : 'color:var(--text-muted);';
      tableHtml += `<tr>
        <td style="font-family:var(--mono);">${window.escapeHTML(r.id)}</td>
        <td>${window.escapeHTML(String(r.length))}</td>
        <td>${window.escapeHTML(String(r.gc_percent))}</td>
        <td style="${noteStyle}">${r.error ? window.escapeHTML(r.error) : '✓'}</td>
      </tr>`;
    });
    tableHtml += '</tbody></table></div></div>';

    const existing = panel.querySelector('#batchInlineResults');
    if (existing) existing.remove();

    const div = document.createElement('div');
    div.id = 'batchInlineResults';
    div.innerHTML = window.BioKit.utils.sanitizeHTML(tableHtml);
    panel.appendChild(div);

    if (btn) { btn.disabled = true; btn.textContent = 'Done'; }
    window.showToast?.('Inline analysis complete.');
  }

  // ── Polling ──────────────────────────────────────────────────
  function startPolling(jobId) {
    pollTimer = setInterval(() => pollStatus(jobId), POLL_INTERVAL_MS);
    pollStatus(jobId); // immediate first check
  }

  function stopPolling() {
    if (pollTimer) { clearInterval(pollTimer); pollTimer = null; }
  }

  async function pollStatus(jobId) {
    try {
      const res = await fetch(`/api/status/${encodeURIComponent(jobId)}`);
      if (!res.ok) return; // transient HTTP error — retry next tick
      const data = await res.json();
      updateProgress(data);

      if (data.status === 'complete') {
        stopPolling();
        onJobComplete(jobId, data);
      } else if (data.status === 'failed' || data.status === 'failure') {
        stopPolling();
        onJobFailed(data.error || 'Worker error');
      }
    } catch (_) {
      // Network blip — keep polling
    }
  }

  // ── Progress UI ──────────────────────────────────────────────
  function showJobArea(jobId) {
    document.getElementById('batchJobArea')?.classList.remove('hidden');
    const idEl = document.getElementById('batchJobId');
    if (idEl) idEl.textContent = `JOB: ${jobId.slice(0, 8).toUpperCase()}…`;
    document.getElementById('batchDownloadArea')?.classList.add('hidden');
    document.getElementById('batchResultsArea')?.classList.add('hidden');
    setProgressBar(0);
  }

  function resetJobUI() {
    document.getElementById('batchJobArea')?.classList.add('hidden');
    document.getElementById('batchDownloadArea')?.classList.add('hidden');
    document.getElementById('batchResultsArea')?.classList.add('hidden');
    const degraded = document.getElementById('batchDegradedPanel');
    if (degraded) degraded.style.display = 'none';
    setProgressBar(0);
  }

  function setBatchStatus(msg, isError) {
    const el = document.getElementById('batchJobStatus');
    if (!el) return;
    el.textContent = msg;
    el.style.color = isError ? '#f87171' : '';
  }

  function updateProgress(data) {
    const pct = Math.min(Math.round(data.progress || 0), 100);
    setProgressBar(pct);

    const pctEl = document.getElementById('batchProgressPct');
    if (pctEl) pctEl.textContent = `${pct}%`;

    const labelEl = document.getElementById('batchProgressLabel');

    if (data.status === 'running') {
      const processed = data.processed || 0;
      const total     = data.total     || 0;
      const msg       = total > 0
        ? `Processing ${processed.toLocaleString()} / ${total.toLocaleString()} sequences…`
        : 'Processing…';
      setBatchStatus(msg, false);
      if (labelEl && total > 0) {
        labelEl.textContent = `${processed.toLocaleString()} of ${total.toLocaleString()} sequences processed`;
      }
    } else if (data.status === 'pending') {
      setBatchStatus('Queued — waiting for worker…', false);
    } else if (data.status === 'complete') {
      setBatchStatus('Complete', false);
      if (labelEl) labelEl.textContent = '';
    } else if (data.status === 'failed' || data.status === 'failure') {
      setBatchStatus('Failed', true);
    }

    const spinner = document.getElementById('batchSpinner');
    if (spinner) {
      spinner.style.display =
        (data.status === 'running' || data.status === 'pending') ? '' : 'none';
    }
  }

  function setProgressBar(pct) {
    const bar = document.getElementById('batchProgressBar');
    if (bar) bar.style.width = `${pct}%`;
  }

  // ── Job Complete ─────────────────────────────────────────────
  async function onJobComplete(jobId, data) {
    setProgressBar(100);
    setBatchStatus('Complete', false);
    document.getElementById('batchSpinner') && (document.getElementById('batchSpinner').style.display = 'none');

    const total  = data.total  || 0;
    const failed = data.failed || 0;
    const mode   = data.mode   || 'all';

    const summary = document.getElementById('batchCompleteSummary');
    if (summary) {
      summary.textContent =
        `${total.toLocaleString()} sequences · mode: ${mode}` +
        (failed > 0 ? ` · ${failed} failed` : '');
    }

    document.getElementById('batchDownloadArea')?.classList.remove('hidden');

    // Update history entry with final status
    updateJobInHistory(jobId, { status: 'complete', total, failed, mode });

    // Wire download buttons for this job
    const csvBtn = document.getElementById('batchDownloadBtn');
    if (csvBtn) csvBtn.onclick = () => triggerDownload(jobId, 'csv');

    const tsvBtn = document.getElementById('batchDownloadTsvBtn');
    if (tsvBtn) tsvBtn.onclick = () => triggerDownload(jobId, 'tsv');

    // Fetch and render the in-app results table
    try {
      await fetchAndRenderResults(jobId);
    } catch (err) {
      window.showToast?.('Could not load inline results table: ' + err.message);
    }

    window.showToast?.('Batch analysis complete — results ready.');
  }

  function onJobFailed(errMsg) {
    setBatchStatus('Job failed: ' + errMsg, true);
    const spinner = document.getElementById('batchSpinner');
    if (spinner) spinner.style.display = 'none';
    window.showToast?.('Batch job failed: ' + errMsg);
  }

  // ── Download ─────────────────────────────────────────────────
  function triggerDownload(jobId, fmt) {
    const url = `/api/batch/${encodeURIComponent(jobId)}/download?format=${fmt}&summary=1`;
    window.location.href = url;
  }

  // ── Results Table (Phase 2) ───────────────────────────────────
  async function fetchAndRenderResults(jobId) {
    const res = await fetch(`/api/batch/${encodeURIComponent(jobId)}/results`);
    if (!res.ok) {
      const body = await res.json().catch(() => ({}));
      throw new Error(body.error || `HTTP ${res.status}`);
    }
    const data = await res.json();
    _allRows  = data.rows    || [];
    _columns  = data.columns || [];

    // Set default column visibility
    _visibleCols = new Set(_columns.filter(c => !COLS_HIDDEN_DEFAULT.has(c)));

    buildColPickerGrid();
    applyFiltersAndRender();

    document.getElementById('batchResultsArea')?.classList.remove('hidden');
  }

  // ── Column Picker (Phase 5) ───────────────────────────────────
  function buildColPickerGrid() {
    const grid = document.getElementById('batchColPickerGrid');
    if (!grid) return;
    grid.innerHTML = '';
    _columns.forEach(col => {
      const label = document.createElement('label');
      label.style.cssText = 'display:flex; align-items:center; gap:5px; cursor:pointer; color:var(--text-secondary);';
      const cb = document.createElement('input');
      cb.type    = 'checkbox';
      cb.value   = col;
      cb.checked = _visibleCols.has(col);
      cb.style.accentColor = 'var(--teal)';
      cb.addEventListener('change', () => {
        cb.checked ? _visibleCols.add(col) : _visibleCols.delete(col);
        applyFiltersAndRender();
      });
      const meta = COL_META[col] || {};
      label.appendChild(cb);
      label.appendChild(document.createTextNode(meta.label || col));
      grid.appendChild(label);
    });
  }

  document.getElementById('batchColPickerBtn')?.addEventListener('click', () => {
    document.getElementById('batchColPickerPanel')?.classList.toggle('hidden');
  });

  document.getElementById('batchColSelectAll')?.addEventListener('click', () => {
    _visibleCols = new Set(_columns);
    buildColPickerGrid();
    applyFiltersAndRender();
  });

  document.getElementById('batchColSelectNone')?.addEventListener('click', () => {
    _visibleCols = new Set(['id', 'length', 'type']); // always keep identity
    buildColPickerGrid();
    applyFiltersAndRender();
  });

  // ── Filters ──────────────────────────────────────────────────
  document.getElementById('batchFilterInput')?.addEventListener('input',  () => { _currentPage = 1; applyFiltersAndRender(); });
  document.getElementById('batchFilterType')?.addEventListener('change',  () => { _currentPage = 1; applyFiltersAndRender(); });
  document.getElementById('batchFilterErrors')?.addEventListener('change',() => { _currentPage = 1; applyFiltersAndRender(); });
  document.getElementById('batchPageSize')?.addEventListener('change', (e) => {
    _pageSize    = parseInt(e.target.value, 10);
    _currentPage = 1;
    renderTable();
  });

  function applyFiltersAndRender() {
    const text      = (document.getElementById('batchFilterInput')?.value  || '').toLowerCase();
    const typeFilter= (document.getElementById('batchFilterType')?.value   || '').toUpperCase();
    const errOnly   =  document.getElementById('batchFilterErrors')?.checked || false;

    _filteredRows = _allRows.filter(row => {
      if (typeFilter && (row.type || '').toUpperCase() !== typeFilter) return false;
      if (errOnly   && !row.error) return false;
      if (text) {
        const haystack = ((row.id || '') + ' ' + (row.description || '')).toLowerCase();
        if (!haystack.includes(text)) return false;
      }
      return true;
    });

    // Re-apply sort on filtered set
    if (_sortCol) {
      _filteredRows = sortRows(_filteredRows, _sortCol, _sortDir);
    }

    const countEl = document.getElementById('batchFilterCount');
    if (countEl) {
      countEl.textContent = `${_filteredRows.length.toLocaleString()} of ${_allRows.length.toLocaleString()} sequences`;
    }

    _currentPage = 1;
    renderTable();
  }

  function sortRows(rows, col, dir) {
    return [...rows].sort((a, b) => {
      const av = a[col] ?? '';
      const bv = b[col] ?? '';
      // Numeric sort if both parse as float
      const an = parseFloat(av);
      const bn = parseFloat(bv);
      const cmp = (!isNaN(an) && !isNaN(bn))
        ? an - bn
        : String(av).localeCompare(String(bv), undefined, { numeric: true });
      return dir === 'asc' ? cmp : -cmp;
    });
  }

  // ── Table Render ─────────────────────────────────────────────
  function renderTable() {
    const thead = document.getElementById('batchResultsThead');
    const tbody = document.getElementById('batchResultsTbody');
    if (!thead || !tbody) return;

    const visibleCols = _columns.filter(c => _visibleCols.has(c));

    // ── Header ──
    thead.innerHTML = '';
    const tr = document.createElement('tr');
    visibleCols.forEach(col => {
      const th   = document.createElement('th');
      const meta = COL_META[col] || {};
      th.textContent = meta.label || col;
      if (meta.tip) th.setAttribute('data-tooltip', meta.tip);
      if (col === _sortCol) th.classList.add(_sortDir === 'asc' ? 'sort-asc' : 'sort-desc');
      th.addEventListener('click', () => {
        if (_sortCol === col) {
          _sortDir = _sortDir === 'asc' ? 'desc' : 'asc';
        } else {
          _sortCol = col;
          _sortDir = 'asc';
        }
        _filteredRows = sortRows(_filteredRows, _sortCol, _sortDir);
        renderTable();
      });
      tr.appendChild(th);
    });
    thead.appendChild(tr);

    // ── Body ──
    const total    = _filteredRows.length;
    const totalPgs = Math.max(1, Math.ceil(total / _pageSize));
    _currentPage   = Math.min(_currentPage, totalPgs);

    const start = (_currentPage - 1) * _pageSize;
    const end   = Math.min(start + _pageSize, total);
    const page  = _filteredRows.slice(start, end);

    tbody.innerHTML = '';
    page.forEach(row => {
      const tr = document.createElement('tr');
      if (row.error)   tr.classList.add('row-error');
      else if (row.warning) tr.classList.add('row-warning');

      visibleCols.forEach(col => {
        const td  = document.createElement('td');
        const val = row[col] ?? '';

        if (col === 'instability_label' && val) {
          const span    = document.createElement('span');
          span.textContent = val;
          span.className   = val === 'Stable' ? 'instability-stable' : 'instability-unstable';
          td.appendChild(span);
        } else if (col === 'error' && val) {
          const badge       = document.createElement('span');
          badge.className   = 'badge-error';
          badge.textContent = 'ERR';
          badge.title       = val;
          td.appendChild(badge);
          td.appendChild(document.createTextNode(' ' + val.slice(0, 60) + (val.length > 60 ? '…' : '')));
        } else if (col === 'warning' && val) {
          const badge       = document.createElement('span');
          badge.className   = 'badge-warning';
          badge.textContent = 'WARN';
          badge.title       = val;
          td.appendChild(badge);
        } else {
          td.textContent = val;
        }
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });

    // ── Pagination controls ──
    const pageInfo = document.getElementById('batchPageInfo');
    if (pageInfo) {
      pageInfo.textContent = total > 0
        ? `Page ${_currentPage} of ${totalPgs} (${start + 1}–${end} of ${total.toLocaleString()})`
        : 'No results';
    }
    const prevBtn = document.getElementById('batchPagePrev');
    const nextBtn = document.getElementById('batchPageNext');
    if (prevBtn) prevBtn.disabled = _currentPage <= 1;
    if (nextBtn) nextBtn.disabled = _currentPage >= totalPgs;
  }

  document.getElementById('batchPagePrev')?.addEventListener('click', () => {
    if (_currentPage > 1) { _currentPage--; renderTable(); }
  });
  document.getElementById('batchPageNext')?.addEventListener('click', () => {
    const totalPgs = Math.ceil(_filteredRows.length / _pageSize);
    if (_currentPage < totalPgs) { _currentPage++; renderTable(); }
  });

  // ── Client-side TSV export (with column picker, Phase 5) ─────
  // The TSV button in the toolbar calls triggerDownload() which hits the server.
  // This secondary function is available for column-filtered client-side export.
  function exportFilteredTsv() {
    if (!_filteredRows.length) return;
    const visibleCols = _columns.filter(c => _visibleCols.has(c));
    const lines = [];

    // Summary header block
    lines.push(`# BioToolkit Batch Analysis Results`);
    lines.push(`# Exported: ${new Date().toISOString()}`);
    lines.push(`# Showing ${_filteredRows.length} of ${_allRows.length} sequences (filtered)`);
    lines.push(`# Columns: ${visibleCols.join(', ')}`);
    lines.push('#');

    lines.push(visibleCols.join('\t'));
    _filteredRows.forEach(row => {
      lines.push(visibleCols.map(c => String(row[c] ?? '').replace(/\t/g, ' ')).join('\t'));
    });

    const blob = new Blob([lines.join('\r\n')], { type: 'text/tab-separated-values' });
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement('a');
    a.href     = url;
    a.download = 'batch_filtered.tsv';
    a.click();
    URL.revokeObjectURL(url);
  }

  // Override TSV button to export filtered + column-picked rows client-side
  // (triggered after results are loaded; initial wiring uses server download)
  function rewireTsvBtn() {
    const tsvBtn = document.getElementById('batchDownloadTsvBtn');
    if (tsvBtn) tsvBtn.onclick = exportFilteredTsv;
  }

  // ── Sample Data (Phase 6 — 3 curated datasets) ───────────────
  document.getElementById('batchSampleBtn')?.addEventListener('click', () => {
    const sample = [
      '>P04637|TP53_HUMAN Cellular tumor antigen p53 [Homo sapiens]',
      'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDP',
      'GPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLFRNL',
      'NQTSFQNLSDLQPPPSSQS',
      '>P00533|EGFR_HUMAN Epidermal growth factor receptor [Homo sapiens]',
      'MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCE',
      'VVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYA',
      'LAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMS',
      '>P68871|HBB_HUMAN Hemoglobin subunit beta [Homo sapiens]',
      'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNP',
      'KVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHH',
      'FGKEFTPPVQAAYQKVVAGVANALAHKYH',
    ].join('\n');
    _loadSample(sample, 'Example: 3 human proteins (TP53, EGFR, HBB)');
  });

  document.getElementById('batchSampleDnaBtn')?.addEventListener('click', () => {
    // Curated set with IUPAC degenerate codes — tests the IUPAC classification fix
    const sample = [
      '>SEQ1_BRCA1_promoter BRCA1 promoter region — contains IUPAC degenerate codes from primer consensus',
      'ATGCRYSWKMBDHVNATGCATGCRYSWKMBDHVNATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',
      '>SEQ2_p53_response_element p53 response element with degenerate codes (R=A|G, Y=C|T)',
      'GGACATGCCCGGGCATGTCCGGACATGCCCGGGCATGTCCRRYYRRYYRR',
      '>SEQ3_CpG_island CpG island in MLH1 promoter region',
      'CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG',
      '>SEQ4_primer_consensus PCR primer consensus with multiple degenerate positions',
      'ATGCRYSWKMBDHVNATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',
    ].join('\n');
    _loadSample(sample, 'Example: DNA sequences with IUPAC degenerate codes');
    document.getElementById('batchModeSelect') && (document.getElementById('batchModeSelect').value = 'nucleotide');
  });

  document.getElementById('batchSampleRnaBtn')?.addEventListener('click', () => {
    // mRNA sequences with U instead of T
    const sample = [
      '>NM_000546.6|TP53_mRNA_fragment Human TP53 mRNA CDS fragment [Homo sapiens]',
      'AUGGAGGAGCCGCAGUCAGAUCCUAGCGUUAGUGCCCCAACUGGAGGCCAGAACUUUUUCUGGAGCAGUCAAAAGGUGGAGUUUGAGGUCAGU',
      '>NM_005228.5|EGFR_mRNA_fragment Human EGFR mRNA CDS fragment [Homo sapiens]',
      'AUGCGACCCUCGGGGACGGCCGGGGCAGCGCUCUUCUUCUCCUUGCAGUGCAGCUGUCGGGGAGCGGAGACGCUGAGGAAAAAGGUGGUGCAG',
      '>NM_000518.5|HBB_mRNA_fragment Human HBB mRNA CDS [Homo sapiens]',
      'AUGUGCACUGACUCUGAGGAGAAGUCUGCCGUUACUGCCCUGUGGGGCAAGGUGAACGUGGAUGUGUUGCUGAGUCUCCUAGGGAUGGCGAUGU',
    ].join('\n');
    _loadSample(sample, 'Example: 3 human mRNA fragments (TP53, EGFR, HBB)');
    document.getElementById('batchModeSelect') && (document.getElementById('batchModeSelect').value = 'nucleotide');
  });

  function _loadSample(text, toast) {
    const pasteArea = document.getElementById('batchPasteInput');
    if (pasteArea) pasteArea.value = text;
    uploadedContent = text;
    updateSeqCount(text);
    window.showToast?.(toast);
  }

  // ── Clear ────────────────────────────────────────────────────
  document.getElementById('batchClearBtn')?.addEventListener('click', () => {
    stopPolling();
    uploadedContent = null;
    currentJobId    = null;
    _allRows        = [];
    _columns        = [];
    _filteredRows   = [];
    const pasteArea = document.getElementById('batchPasteInput');
    if (pasteArea)   pasteArea.value = '';
    if (fileInput)   fileInput.value = '';
    document.getElementById('batchSeqCount') && (document.getElementById('batchSeqCount').style.display = 'none');
    resetJobUI();
  });

  // ── localStorage Job History (Phase 7) ───────────────────────
  function saveJobToHistory(jobId, meta) {
    const history = getJobHistory();
    // Avoid duplicates
    const existing = history.findIndex(e => e.jobId === jobId);
    const entry    = { jobId, ...meta, timestamp: Date.now() };
    if (existing >= 0) {
      history[existing] = entry;
    } else {
      history.unshift(entry);
    }
    try {
      localStorage.setItem(HISTORY_KEY, JSON.stringify(history.slice(0, HISTORY_MAX)));
    } catch (_) { /* storage quota — silently skip */ }
    renderJobHistory();
  }

  function updateJobInHistory(jobId, updates) {
    const history = getJobHistory();
    const idx     = history.findIndex(e => e.jobId === jobId);
    if (idx >= 0) {
      Object.assign(history[idx], updates);
      try { localStorage.setItem(HISTORY_KEY, JSON.stringify(history)); } catch (_) {}
    }
    renderJobHistory();
  }

  function getJobHistory() {
    try { return JSON.parse(localStorage.getItem(HISTORY_KEY) || '[]'); }
    catch (_) { return []; }
  }

  function renderJobHistory() {
    const history    = getJobHistory();
    const details    = document.getElementById('batchHistoryDetails');
    const listEl     = document.getElementById('batchHistoryList');
    if (!details || !listEl) return;

    if (history.length === 0) {
      details.classList.add('hidden');
      return;
    }
    details.classList.remove('hidden');
    listEl.innerHTML = '';

    history.forEach(entry => {
      const item = document.createElement('div');
      item.className = 'batch-history-item';

      const ago = _relativeTime(entry.timestamp);
      const seqInfo = entry.total
        ? `${entry.total.toLocaleString()} seqs${entry.failed ? `, ${entry.failed} failed` : ''}`
        : `${(entry.seqs || '?').toLocaleString()} seqs submitted`;

      item.innerHTML = `
        <span class="hist-id">${entry.jobId.slice(0, 8).toUpperCase()}…</span>
        <span class="hist-meta">${seqInfo} · mode: ${entry.mode || 'all'}</span>
        <span class="hist-ts">${ago}</span>
      `;

      item.addEventListener('click', () => restoreJob(entry.jobId, entry));
      listEl.appendChild(item);
    });
  }

  async function restoreJob(jobId, entry) {
    stopPolling();
    currentJobId = jobId;
    resetJobUI();

    try {
      const res  = await fetch(`/api/status/${encodeURIComponent(jobId)}`);
      const data = await res.json();

      if (data.status === 'complete') {
        showJobArea(jobId);
        setProgressBar(100);
        setBatchStatus('Complete (restored)', false);
        document.getElementById('batchSpinner') && (document.getElementById('batchSpinner').style.display = 'none');

        const summaryEl = document.getElementById('batchCompleteSummary');
        if (summaryEl) {
          summaryEl.textContent =
            `${(data.total || entry.total || 0).toLocaleString()} sequences · mode: ${data.mode || entry.mode || 'all'}` +
            (data.failed > 0 ? ` · ${data.failed} failed` : '');
        }

        document.getElementById('batchDownloadArea')?.classList.remove('hidden');
        const csvBtn = document.getElementById('batchDownloadBtn');
        if (csvBtn) csvBtn.onclick = () => triggerDownload(jobId, 'csv');
        const tsvBtn = document.getElementById('batchDownloadTsvBtn');
        if (tsvBtn) tsvBtn.onclick = () => triggerDownload(jobId, 'tsv');

        await fetchAndRenderResults(jobId);
        rewireTsvBtn();
        window.showToast?.('Restored job results.');

      } else if (data.status === 'running' || data.status === 'pending') {
        showJobArea(jobId);
        startPolling(jobId);
        window.showToast?.('Job still running — reconnected.');
      } else {
        window.showToast?.(`Job ${jobId.slice(0, 8)} status: ${data.status}`);
      }
    } catch (err) {
      window.showToast?.('Could not restore job: ' + err.message);
    }
  }

  function _relativeTime(ts) {
    const diffSec = Math.round((Date.now() - ts) / 1000);
    if (diffSec < 60)   return `${diffSec}s ago`;
    if (diffSec < 3600) return `${Math.round(diffSec / 60)}m ago`;
    if (diffSec < 86400)return `${Math.round(diffSec / 3600)}h ago`;
    return `${Math.round(diffSec / 86400)}d ago`;
  }

  // ── Auto-reconnect on load (Phase 7) ────────────────────────
  (function autoReconnect() {
    renderJobHistory();
    const history = getJobHistory();
    if (!history.length) return;
    const last = history[0];
    // Only auto-reconnect if the last job was queued/running recently (< 12 h)
    const age  = (Date.now() - (last.timestamp || 0)) / 1000;
    if (last.status === 'running' && age < 43200) {
      restoreJob(last.jobId, last);
    }
  })();

  // ── ViewManager Cleanup ───────────────────────────────────────
  if (window.ViewManager) {
    window.ViewManager.registerUnmount('batch', () => {
      stopPolling();
    });
  }

})();
