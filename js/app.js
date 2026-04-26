/* ============================================================
   app.js — Tab Routing + Shared Utilities
   ============================================================ */

(function () {
  'use strict';

  // ── Initialize BioKit Namespace ────────────────────────────
  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  // ── Browser Compatibility & Polyfills ───────────────────
  if (typeof CanvasRenderingContext2D !== 'undefined' && !CanvasRenderingContext2D.prototype.roundRect) {
    CanvasRenderingContext2D.prototype.roundRect = function (x, y, w, h, r) {
      if (typeof r === 'number') r = { tl: r, tr: r, br: r, bl: r };
      else if (Array.isArray(r)) {
        const [tl = 0, tr = 0, br = 0, bl = 0] = r;
        r = { tl, tr, br, bl };
      } else {
        r = Object.assign({ tl: 0, tr: 0, br: 0, bl: 0 }, r);
      }
      this.beginPath();
      this.moveTo(x + r.tl, y);
      this.lineTo(x + w - r.tr, y);
      this.quadraticCurveTo(x + w, y, x + w, y + r.tr);
      this.lineTo(x + w, y + h - r.br);
      this.quadraticCurveTo(x + w, y + h, x + w - r.br, y + h);
      this.lineTo(x + r.bl, y + h);
      this.quadraticCurveTo(x, y + h, x, y + h - r.bl);
      this.lineTo(x, y + r.tl);
      this.quadraticCurveTo(x, y, x + r.tl, y);
      this.closePath();
      return this;
    };
  }

  // Initialize Lucide icons
  if (window.lucide) {
    window.lucide.createIcons();
  }

  // ── ViewManager — Lifecycle-aware SPA router ──────────────
  //
  // Each tool module may register an unmount() callback via:
  //   window.ViewManager.registerUnmount('toolId', fn)
  // The ViewManager calls this before activating a new view so that
  // heavy DOM nodes are physically removed and JS references are
  // set to null — preventing the browser OOM documented in the audit.
  //
  const ViewManager = {
    _currentToolId: null,
    _unmountRegistry: {},

    /**
     * Register a teardown callback for a specific tool.
     * Call from within each module's IIFE:
     *   window.ViewManager.registerUnmount('fasta', () => { ... })
     */
    registerUnmount(toolId, fn) {
      this._unmountRegistry[toolId] = fn;
    },

    /**
     * Navigate to a new tool panel.
     * 1. Runs the outgoing tool's unmount() to destroy heavy DOM + null refs.
     * 2. Hides all panels via CSS class (fast).
     * 3. Reveals the target panel and marks its nav item active.
     * 4. Scrolls content area to top.
     */
    navigate(toolId) {
      // ── Teardown outgoing view ──────────────────────────────
      const prev = this._currentToolId;
      if (prev && prev !== toolId && typeof this._unmountRegistry[prev] === 'function') {
        try {
          this._unmountRegistry[prev]();
        } catch (e) {
          console.warn('[ViewManager] unmount error for', prev, e);
        }
      }

      // ── Swap panels ──────────────────────────────────────────
      document.querySelectorAll('.tool-panel').forEach(p => p.classList.remove('active'));
      document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));

      const panel = document.getElementById('panel-' + toolId);
      const nav   = document.getElementById('nav-' + toolId);
      if (panel) {
        panel.classList.add('active');
        // Fix 5: notify tool modules that their panel is now visible
        panel.dispatchEvent(new CustomEvent('panel-shown', { bubbles: false }));
      }
      if (nav)   nav.classList.add('active');

      const main = document.getElementById('mainContent');
      if (main) main.scrollTo({ top: 0, behavior: 'smooth' });

      this._currentToolId = toolId;
    }
  };

  // Expose globally so modules can call ViewManager.registerUnmount()
  window.ViewManager = ViewManager;

  // Wire sidebar nav buttons
  document.querySelectorAll('.nav-item').forEach(btn => {
    btn.addEventListener('click', () => ViewManager.navigate(btn.dataset.tool));
  });

  // Home tool-cards navigate to tool
  document.querySelectorAll('.tool-card').forEach(card => {
    card.addEventListener('click', () => ViewManager.navigate(card.dataset.tool));
  });

  // ── Copy to Clipboard ─────────────────────────────────────
  document.addEventListener('click', (e) => {
    const btn = e.target.closest('.copy-btn');
    if (!btn) return;
    const targetId  = btn.dataset.target;
    const targetEl  = document.getElementById(targetId);
    const text      = targetEl ? targetEl.innerText.trim() : '';
    if (!text) return;
    navigator.clipboard.writeText(text).then(() => showToast('Copied!')).catch(() => {});
  });

  // ── Toast Notification ────────────────────────────────────
  function showToast(msg, duration = 2000) {
    const existing = document.querySelector('.toast');
    if (existing) existing.remove();
    const toast = document.createElement('div');
    toast.className = 'toast';
    toast.textContent = msg;
    document.body.appendChild(toast);
    setTimeout(() => toast.remove(), duration);
  }
  window.showToast = showToast;
  window.BioKit.utils.showToast = showToast;

  // ── Debounce Utility ──────────────────────────────────────
  function debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
      const later = () => {
        clearTimeout(timeout);
        func(...args);
      };
      clearTimeout(timeout);
      timeout = setTimeout(later, wait);
    };
  }
  window.debounce = debounce;
  window.BioKit.utils.debounce = debounce;

  // ── Shared Sequence Cleaner ───────────────────────────────────
  function cleanSeq(seq) {
    return stripFastaHeader(seq).replace(/[^A-Z\*]/g, '');
  }
  window.cleanSeq = cleanSeq;
  window.BioKit.utils.cleanSeq = cleanSeq;

  // ── FASTA Header Stripper ─────────────────────────────────
  /**
   * Strips FASTA header lines (lines starting with '>') and
   * collapses all whitespace from the remaining sequence lines.
   * Safe to call on a plain sequence string with no headers.
   * @param {string} raw - Raw textarea value, may include one or more FASTA headers.
   * @returns {string} Pure sequence string (uppercase, no whitespace, no headers).
   */
  function stripFastaHeader(raw) {
    if (!raw) return '';
    return raw
      .split('\n')
      .filter(line => !line.trimStart().startsWith('>'))
      .join('')
      .replace(/[\s\r\n]/g, '')
      .toUpperCase();
  }
  window.stripFastaHeader = stripFastaHeader;
  window.BioKit.utils.stripFastaHeader = stripFastaHeader;

  // ── Format Auto-Detection ─────────────────────────────────
  function detectFormat(seq) {
    const s = cleanSeq(seq);
    if (!s) return 'dna';
    const dnaCount = (s.match(/[ATGCNU]/g) || []).length;
    return (dnaCount / s.length > 0.8) ? 'dna' : 'protein';
  }
  window.detectFormat = detectFormat;
  window.BioKit.utils.detectFormat = detectFormat;

  // ── Sequence Validation ───────────────────────────────────
  // Minimum length thresholds (DNA only)
  const SEQ_MIN_HARD  = 2;   // absolute floor — reject outright
  const SEQ_MIN_SOFT  = 10;  // soft floor — proceed with inline warning

  function validateSequence(seq, type) {
    const raw = seq || '';
    // Strip FASTA headers first so a pasted FASTA file never triggers a false rejection
    const headerStripped = stripFastaHeader(raw);
    const clean = headerStripped.toUpperCase().replace(/\s/g, '');
    if (!clean) return { valid: false, msg: 'Sequence is empty.' };

    let allowed = /^[ACGTURYSWKMBDHVN]+$/i;
    let label = 'DNA/RNA';
    let limit = 0;
    if (type === 'protein') {
      // IUPAC extended: standard 20 AA + X (unknown) + U (Selenocysteine) + Z (Glx) + stop codon *
      allowed = /^[ACDEFGHIKLMNPQRSTVWYXUZ\*]+$/;
      label = 'Protein';
      limit = 5000000; // 5MB
    } else {
      limit = 10000000; // 10MB
    }

    const invalidChars = [];
    for (let i = 0; i < clean.length; i++) {
      if (!allowed.test(clean[i])) {
        invalidChars.push({ char: clean[i], pos: i + 1 });
        if (invalidChars.length > 5) break;
      }
    }

    if (invalidChars.length > 0 && invalidChars.length / clean.length > 0.05) {
      return {
        valid: false,
        msg: `Sequence rejected: Invalid characters detected (>5%). This does not look like a valid ${label} sequence.`,
        suggestion: `First invalid character '${invalidChars[0].char}' found at position ${invalidChars[0].pos}.`
      };
    }

    if (clean.length > limit) {
      return {
        valid: false,
        msg: `Sequence is too large: ${clean.length.toLocaleString()} ${type === 'protein' ? 'amino acids' : 'base pairs'}.`,
        suggestion: `The limit is ${limit.toLocaleString()} for browser processing. Consider using a backend API.`
      };
    }

    // ── Minimum length checks (DNA/RNA path only) ─────────────
    if (type !== 'protein') {
      if (clean.length < SEQ_MIN_HARD) {
        return {
          valid: false,
          msg: `Sequence too short (${clean.length} bp). Minimum is ${SEQ_MIN_HARD} bp.`,
          suggestion: 'At least 2 bases are required to perform any calculation.'
        };
      }
      if (clean.length < SEQ_MIN_SOFT) {
        // Allow, but tag result with a warning for the UI to display
        return {
          valid: true,
          clean: clean.replace(/[^A-Z\*]/g, ''),
          warning: `Short sequence (${clean.length} bp). Results — especially Tm — may be unreliable. Recommended minimum: ${SEQ_MIN_SOFT} bp.`
        };
      }
    }

    return { valid: true, clean: clean.replace(/[^A-Z\*]/g, '') };
  }
  window.validateSequence = validateSequence;
  window.BioKit.utils.validateSequence = validateSequence;

  function showValidationWarning(panelId, validation) {
    const panel = document.getElementById(panelId);
    if (!panel) return;
    clearValidationAlert(panelId);

    const isSoft = validation.valid === true; // soft advisory vs hard error
    const alert = document.createElement('div');
    alert.className = isSoft ? 'validation-alert validation-advisory' : 'validation-alert';
    alert.innerHTML = `
      <div class="validation-alert-icon">
        <i data-lucide="${isSoft ? 'alert-circle' : 'alert-triangle'}"></i>
      </div>
      <div class="validation-alert-content">
        <div class="validation-alert-title">${isSoft ? '⚠ Advisory' : 'Input Validation Error'}</div>
        <div class="validation-alert-msg">${escapeHTML(validation.msg)}</div>
        ${validation.suggestion ? `<div class="validation-suggestion">${escapeHTML(validation.suggestion)}</div>` : ''}
      </div>
    `;
    const inputBlock = panel.querySelector('.input-block');
    if (inputBlock) inputBlock.parentNode.insertBefore(alert, inputBlock);
    if (window.lucide) window.lucide.createIcons();

    // Only hide results on a hard error — soft warnings let results remain visible
    if (!isSoft) {
      const resultsArea = panel.querySelector('.results-area');
      if (resultsArea) resultsArea.classList.add('hidden');
    }
  }
  window.showValidationWarning = showValidationWarning;
  window.BioKit.utils.showValidationWarning = showValidationWarning;

  function clearValidationAlert(panelId) {
    const panel = document.getElementById(panelId);
    if (!panel) return;
    const existing = panel.querySelector('.validation-alert');
    if (existing) existing.remove();
  }
  window.clearValidationAlert = clearValidationAlert;
  window.BioKit.utils.clearValidationAlert = clearValidationAlert;

  // ── Centralized Error Handling ───────────────────────────
  function handleError(error, context = 'Analysis') {
    console.error(`[${context}]`, error);
    const msg = error.message || 'An unexpected calculation error occurred.';
    showToast(`⚠ ${context} Error: ${msg}`, 4000);
    document.querySelectorAll('.progress-container').forEach(p => p.classList.remove('active'));
  }
  window.handleError = handleError;
  window.BioKit.utils.handleError = handleError;

  function clearAllResults() {
    document.querySelectorAll('.results-area').forEach(a => a.classList.add('hidden'));
    document.querySelectorAll('.progress-container').forEach(p => p.classList.remove('active'));
    showToast('All results cleared.');
  }
  window.clearAllResults = clearAllResults;

  // ── HTML Escaping & Sanitization ──────────────────────────
  function escapeHTML(str) {
    if (!str) return '';
    return String(str).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;').replace(/'/g, '&#39;');
  }
  window.escapeHTML = escapeHTML;
  window.BioKit.utils.escapeHTML = escapeHTML;

  function sanitizeHTML(html) {
    if (typeof DOMPurify !== 'undefined') return DOMPurify.sanitize(html);
    return escapeHTML(html);
  }
  window.sanitizeHTML = sanitizeHTML;
  window.BioKit.utils.sanitizeHTML = sanitizeHTML;

  // ── Stat Card Builder ─────────────────────────────────────
  function buildStatCards(containerId, stats) {
    const c = document.getElementById(containerId);
    if (!c) return;
    const html = stats.map(s => `
      <div class="result-stat">
        <div class="result-stat-val">${escapeHTML(String(s.val))}</div>
        <div class="result-stat-label">${escapeHTML(s.label)}</div>
      </div>`).join('');
    c.innerHTML = sanitizeHTML(html);
  }
  window.buildStatCards = buildStatCards;
  window.BioKit.utils.buildStatCards = buildStatCards;

  // ── Result Tab Switching ──────────────────────────────────
  document.addEventListener('click', (e) => {
    const tab = e.target.closest('.result-tab');
    if (!tab) return;
    const panel = e.target.closest('.results-area, .workspace-results');
    if (!panel) return;
    panel.querySelectorAll('.result-tab').forEach(t => t.classList.remove('active'));
    panel.querySelectorAll('.result-tab-content').forEach(c => c.classList.add('hidden'));
    tab.classList.add('active');
    const rtab = tab.dataset.rtab;
    const parentId = panel.id;
    const prefix   = parentId.replace('Results', '');
    const content  = document.getElementById(prefix + '-' + rtab);
    if (content) content.classList.remove('hidden');
  });

  // ── Download helper ────────────────────────────────────────
  function downloadText(filename, text) {
    const a    = document.createElement('a');
    a.href     = URL.createObjectURL(new Blob([text], { type: 'text/plain' }));
    a.download = filename;
    a.click();
    URL.revokeObjectURL(a.href);
  }
  window.downloadText = downloadText;
  window.BioKit.utils.downloadText = downloadText;

  // ── Expose activateTool globally (delegates to ViewManager) ─
  // Kept for backwards compat — any existing call to window.activateTool() still works.
  window.activateTool = (toolId) => ViewManager.navigate(toolId);

  // ── Fancy Sequence Renderer ───────────────────────────────
  window.renderSequence = function (seqOptions) {
    const seq = seqOptions.sequence || '';
    const containerId = seqOptions.containerId;

    let html = '<div class="seq-viewer-wrap">';
    html += '<div class="seq-index" style="display:flex; width:max-content; position:relative; height: 1.5em; margin-bottom: 0.5em;">';
    for (let i = 0; i < seq.length; i += 10) {
      if (i > 0) {
        html += `<span style="position:absolute; left:${i}ch;">${i}</span>`;
      } else {
        html += `<span style="position:absolute; left:0;">1</span>`;
      }
    }
    html += '</div>';
    html += '<div style="display:inline-block; letter-spacing:0;">';
    for (let c of seq) {
      const cleanChar = escapeHTML(c);
      let cClass = `base-${cleanChar}`;
      if (!['A', 'T', 'G', 'C', 'U'].includes(cleanChar)) cClass = 'base-N';
      html += `<span class="${cClass}">${cleanChar}</span>`;
    }
    html += '</div></div>';

    const container = document.getElementById(containerId);
    if (container) container.innerHTML = sanitizeHTML(html);
  };

  // ── Worker Manager (Web Worker + Server-side Backend) ─────
  window.WorkerManager = {
    worker: null,
    currentTask: null,

    init() {
      if (!window.Worker) {
        showToast('⚠ Web Workers are not supported in this browser.');
        return false;
      }
      if (!this.worker) {
        this.worker = new Worker('js/worker.js');
        this.worker.onmessage = (e) => this.handleMessage(e.data);
        this.worker.onerror = (err) => handleError(err, 'Worker');
      }
      return true;
    },

    /**
     * Runs a task either in the local Web Worker or on the remote Server-side (Python/Biopython).
     * Automatically chooses based on strategy and sequence size if not specified.
     */
    async runTask(type, payload, onProgress, strategy = 'auto') {
      const seqSize = (payload.sequence?.length || 0) + (payload.a?.length || 0) + (payload.b?.length || 0);
      
      // Auto-offload to backend if sequence is large (> 500kb) or it's a scientific analytical tool
      const isComplexTool = ['THERMO_ANALYSIS', 'RESTRICTION_SEARCH', 'SEQ_SEARCH'].includes(type);
      const shouldOffload = strategy === 'backend' || (strategy === 'auto' && (seqSize > 500000 || isComplexTool));
      
      if (shouldOffload) {
        return this.runBackendTask(type, payload, onProgress);
      } else {
        return this.runLocalTask(type, payload, onProgress);
      }
    },

    runLocalTask(type, payload, onProgress) {
      if (!this.init()) return Promise.reject('No Worker Support');
      const taskId = Math.random().toString(36).substr(2, 9);
      return new Promise((resolve, reject) => {
        this.currentTask = { taskId, resolve, reject, onProgress };
        this.worker.postMessage({ type, payload, taskId });
      });
    },

    async runBackendTask(type, payload, onProgress) {
      const taskId = "backend-" + Math.random().toString(36).substr(2, 9);
      
      try {
        // 1. Submit the job
        const response = await fetch('/api/jobs/create', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            tool_type: type.toLowerCase().replace('_', '-'),
            params: payload
          })
        });
        
        if (!response.ok) throw new Error(`Server failed to accept task: ${response.status}`);
        const data = await response.json();
        const serverTaskId = data.task_id;
        
        // 2. Poll for status
        return new Promise((resolve, reject) => {
          const pollInterval = setInterval(async () => {
            try {
              const statusRes = await fetch(`/api/jobs/status/${serverTaskId}`);
              const statusData = await statusRes.json();
              
              if (statusData.state === 'SUCCESS') {
                clearInterval(pollInterval);
                resolve(statusData.result);
              } else if (statusData.state === 'FAILURE') {
                clearInterval(pollInterval);
                reject(new Error(statusData.error || 'Server-side task failed.'));
              } else if (onProgress && statusData.progress) {
                onProgress(statusData.progress);
              }
            } catch (err) {
              clearInterval(pollInterval);
              reject(err);
            }
          }, 2000); // Poll every 2 seconds
        });
      } catch (err) {
        throw err;
      }
    },

    handleMessage(data) {
      if (!this.currentTask || data.taskId !== this.currentTask.taskId) return;
      if (data.type === 'PROGRESS') {
        if (this.currentTask.onProgress) this.currentTask.onProgress(data.progress);
      } else if (data.type === 'RESULT') {
        this.currentTask.resolve(data.result);
        this.currentTask = null;
      } else if (data.type === 'ERROR') {
        this.currentTask.reject(new Error(data.message));
        this.currentTask = null;
      }
    },

    terminate() {
      if (this.worker) this.worker.terminate();
      this.worker = null;
      this.currentTask = null;
    }
  };
  window.BioKit.core.WorkerManager = window.WorkerManager;

  // ── Unified Loading Handler ─────────────────────────────
  window.withLoading = function (panelId, callback, seqLength = 0, workerTask = null) {
    const panel = document.getElementById(panelId);
    if (!panel) { callback(); return; }
    const resultsArea = panel.querySelector('.results-area');
    clearValidationAlert(panelId);

    let progContainer = panel.querySelector('.progress-container');
    if (!progContainer) {
      progContainer = document.createElement('div');
      progContainer.className = 'progress-container active';
      progContainer.innerHTML = `
        <div class="progress-bar-wrap">
          <div class="progress-bar"></div>
        </div>
        <div class="progress-meta">
          <span class="progress-text">Processing...</span>
          <button class="progress-cancel-btn">Cancel</button>
        </div>
      `;
      const btnRow = panel.querySelector('.btn-row');
      if (btnRow) btnRow.parentNode.insertBefore(progContainer, btnRow.nextSibling);
    } else {
      progContainer.classList.add('active');
    }

    const bar = progContainer.querySelector('.progress-bar');
    const text = progContainer.querySelector('.progress-text');
    const cancelBtn = progContainer.querySelector('.progress-cancel-btn');

    bar.style.width = '0%';
    if (resultsArea) resultsArea.classList.add('hidden');

    const resetUI = () => {
      progContainer.classList.remove('active');
      if (resultsArea) resultsArea.classList.remove('hidden');
    };

    cancelBtn.onclick = () => {
      if (workerTask) window.WorkerManager.terminate();
      progContainer.classList.remove('active');
      showToast('Task cancelled.');
    };

    if (workerTask) {
      // Async Worker Task
      text.textContent = 'Running heavy computation...';
      window.WorkerManager.runTask(workerTask.type, workerTask.payload, (p) => {
        bar.style.width = (p * 100) + '%';
      }).then(result => {
        bar.style.width = '100%';
        setTimeout(() => {
          resetUI();
          if (window.lucide) window.lucide.createIcons();
          callback(result);
        }, 200);
      }).catch(err => {
        progContainer.classList.remove('active');
        handleError(err, 'Worker');
      });
    } else {
      // Synchronous/Simulated Task
      const duration = seqLength > 10000 ? 800 : 400;
      bar.style.transition = `width ${duration}ms ease-out`;
      setTimeout(() => { bar.style.width = '100%'; }, 50);
      setTimeout(() => {
        resetUI();
        try {
          if (window.lucide) window.lucide.createIcons();
          callback();
        } catch (err) {
          handleError(err, 'Task');
        }
      }, duration + 100);
    }
  };


  // ── Session & Persistence ──────────────────────────────
  window.BioKit.core.SessionManager = {
    STORAGE_KEY: 'biokit_session',
    hasUnsavedChanges: false,

    init() {
      this.load();
      setInterval(() => this.save(), 30000);
      window.addEventListener('beforeunload', (e) => {
        if (this.hasUnsavedChanges) { e.preventDefault(); e.returnValue = ''; }
      });

      document.getElementById('saveSessionBtn')?.addEventListener('click', () => this.export());
      document.getElementById('loadSessionBtn')?.addEventListener('click', () => document.getElementById('loadSessionInput')?.click());
      document.getElementById('loadSessionInput')?.addEventListener('change', (e) => this.import(e));
      document.getElementById('exportPdfBtn')?.addEventListener('click', () => window.print());
    },

    getState() {
      const state = { version: '1.0', timestamp: new Date(), tools: {} };
      document.querySelectorAll('input, textarea, select').forEach(el => {
        if (el.id && !el.classList.contains('fetch-id-input')) {
          state.tools[el.id] = el.type === 'checkbox' ? el.checked : el.value;
        }
      });
      return state;
    },

    applyState(state) {
      if (!state || !state.tools) return;
      Object.entries(state.tools).forEach(([id, val]) => {
        const el = document.getElementById(id);
        if (el) {
          if (el.type === 'checkbox') el.checked = val;
          else el.value = val;
        }
      });
    },

    save() {
      try {
        const state = JSON.stringify(this.getState());
        if (state.length < 5000000) {
          localStorage.setItem(this.STORAGE_KEY, state);
          this.hasUnsavedChanges = false;
        }
      } catch (e) { console.warn('Auto-save failed:', e); }
    },

    load() {
      const stored = localStorage.getItem(this.STORAGE_KEY);
      if (stored) {
        try {
          const state = JSON.parse(stored);
          this.applyState(state);
        } catch (e) { localStorage.removeItem(this.STORAGE_KEY); }
      }
    },

    export() {
      const s = JSON.stringify(this.getState(), null, 2);
      downloadText(`biokit_session_${Date.now()}.json`, s);
      this.hasUnsavedChanges = false;
    },

    import(e) {
      const file = e.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = (ev) => {
        try {
          const state = JSON.parse(ev.target.result);
          this.applyState(state);
          showToast('Session imported.');
        } catch (err) { showToast('Invalid file.'); }
      };
      reader.readAsText(file);
    },

    markChanged() { this.hasUnsavedChanges = true; }
  };

  // ── History & Shortcuts ────────────────────────────────
  window.BioKit.core.HistoryManager = {
    stack: [],
    index: -1,
    max: 10,

    push(snapshot) {
      this.stack = this.stack.slice(0, this.index + 1);
      this.stack.push(snapshot);
      if (this.stack.length > this.max) this.stack.shift();
      this.index = this.stack.length - 1;
    },

    undo() {
      if (this.index > 0) {
        this.index--;
        this.apply();
      }
    },

    redo() {
      if (this.index < this.stack.length - 1) {
        this.index++;
        this.apply();
      }
    },

    apply() {
      const snap = this.stack[this.index];
      if (snap) {
        window.BioKit.core.SessionManager.applyState(snap);
        showToast(`History: ${this.index + 1}/${this.stack.length}`);
      }
    }
  };

  window.BioKit.core.ShortcutManager = {
    init() {
      window.addEventListener('keydown', (e) => {
        const isCtrl = e.ctrlKey || e.metaKey;
        const isShift = e.shiftKey;
        const target = e.target;
        const isInput = ['INPUT', 'TEXTAREA'].includes(target.tagName) || target.isContentEditable;

        // 1. Run Analysis (Ctrl+Enter)
        if (isCtrl && (e.key === 'Enter' || e.keyCode === 13)) {
          e.preventDefault();
          const activePanel = document.querySelector('.tool-panel.active');
          const runBtn = activePanel?.querySelector('.btn-primary');
          if (runBtn) {
            showToast('🚀 Running Analysis...');
            runBtn.click();
          }
        }

        // 2. Undo (Ctrl+Z)
        if (isCtrl && !isShift && (e.key === 'z' || e.key === 'Z') && !isInput) {
          e.preventDefault();
          window.BioKit.core.HistoryManager.undo();
        }

        // 3. Redo (Ctrl+Y or Ctrl+Shift+Z)
        if (((isCtrl && (e.key === 'y' || e.key === 'Y')) || (isCtrl && isShift && (e.key === 'z' || e.key === 'Z'))) && !isInput) {
          e.preventDefault();
          window.BioKit.core.HistoryManager.redo();
        }

        // 4. Export (Ctrl+D)
        if (isCtrl && (e.key === 'd' || e.key === 'D') && !isInput) {
          const exportBtn = document.querySelector('.tool-panel.active #proExportBtn');
          if (exportBtn) {
            e.preventDefault();
            showToast('📄 Exporting results...');
            exportBtn.click();
          }
        }

        // 5. Close Modals (Escape)
        if (e.key === 'Escape' || e.keyCode === 27) {
          const modals = document.querySelectorAll('.settings-overlay:not(.hidden)');
          if (modals.length > 0) {
            modals.forEach(m => m.classList.add('hidden'));
            showToast('Modals closed');
          }
        }

        // 6. Help Modal (?)
        // We only trigger '?' shortcut if user is NOT typing in an input
        if (e.key === '?' && !isInput) {
          e.preventDefault();
          const modal = document.getElementById('shortcutModal');
          if (modal) {
            modal.classList.toggle('hidden');
          }
        }
      });
    }
  };

  // ── App Boot ───────────────────────────────────────────
  document.addEventListener('DOMContentLoaded', () => {
    if (!window.Worker || !window.Int32Array) {
      const banner = document.getElementById('browserWarning');
      if (banner) banner.classList.remove('hidden');
    }

    document.getElementById('closeBrowserWarning')?.addEventListener('click', () => {
      document.getElementById('browserWarning')?.classList.add('hidden');
    });

    document.getElementById('clearAllBtn')?.addEventListener('click', clearAllResults);
    document.getElementById('shortcutHelpBtn')?.addEventListener('click', () => {
      const modal = document.getElementById('shortcutModal');
      if (modal) modal.classList.toggle('hidden');
    });

    // Initialize Session & Shortcut Managers
    window.BioKit.core.SessionManager.init();
    window.BioKit.core.ShortcutManager.init();

      // ── Theme Management ────────────────────────────────────
    const ThemeManager = {
      STORAGE_KEY: 'biokit-theme',

      init() {
        const toggleBtn = document.getElementById('themeToggleBtn');
        if (!toggleBtn) return;

        toggleBtn.addEventListener('click', () => this.toggle());

        // Initial setup already happened in head script to prevent FOUC,
        // but we ensure consistency here.
        this.apply(this.getStoredTheme() || this.getSystemTheme());
      },

      getStoredTheme() {
        return localStorage.getItem(this.STORAGE_KEY);
      },

      getSystemTheme() {
        return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
      },

      toggle() {
        const current = document.documentElement.classList.contains('dark') ? 'dark' : 'light';
        const next = current === 'dark' ? 'light' : 'dark';
        this.apply(next);
      },

      apply(theme) {
        const isDark = theme === 'dark';
        document.documentElement.classList.toggle('dark', isDark);
        localStorage.setItem(this.STORAGE_KEY, theme);
        
        // Dispatch global event for canvas re-rendering
        window.dispatchEvent(new CustomEvent('biokit-theme-change', { detail: { theme } }));
        
        // Update tooltip if present (accessibility)
        const toggleBtn = document.getElementById('themeToggleBtn');
        if (toggleBtn) {
          toggleBtn.setAttribute('aria-pressed', isDark);
        }
      }
    };
    ThemeManager.init();
    window.BioKit.core.ThemeManager = ThemeManager;

    // Mark changes on input
    document.querySelectorAll('input, textarea, select, [contenteditable]').forEach(el => {
      el.addEventListener('input', () => window.BioKit.core.SessionManager.markChanged());
      if (el.tagName === 'SELECT') {
        el.addEventListener('change', () => window.BioKit.core.SessionManager.markChanged());
      }
    });

    // ── Sidebar Toggle ──────────────────────────────────────
    const SidebarManager = {
      STORAGE_KEY: 'bionised_sidebar_collapsed',

      init() {
        const btn = document.getElementById('sidebarToggleBtn');
        if (!btn) return;

        if (localStorage.getItem(this.STORAGE_KEY) === 'true') {
          document.body.classList.add('sidebar-collapsed');
        }

        btn.addEventListener('click', () => this.toggle());
      },

      toggle() {
        const collapsed = document.body.classList.toggle('sidebar-collapsed');
        localStorage.setItem(this.STORAGE_KEY, collapsed);

        // Fire resize after CSS transition finishes so Mol* and Genome Viewer
        // recalculate their WebGL canvas dimensions at the new viewport width.
        setTimeout(() => window.dispatchEvent(new Event('resize')), 310);
      }
    };
    SidebarManager.init();
    window.BioKit.core.SidebarManager = SidebarManager;
  });

})();
