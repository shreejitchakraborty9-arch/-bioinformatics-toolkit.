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

  // ── Navigation ───────────────────────────────────────────
  const panels   = document.querySelectorAll('.tool-panel');
  const navItems = document.querySelectorAll('.nav-item');

  function activateTool(toolId) {
    panels.forEach(p => p.classList.remove('active'));
    navItems.forEach(n => n.classList.remove('active'));

    const panel = document.getElementById('panel-' + toolId);
    const nav   = document.getElementById('nav-' + toolId);
    if (panel) panel.classList.add('active');
    if (nav)   nav.classList.add('active');
    document.getElementById('mainContent').scrollTo({ top: 0, behavior: 'smooth' });
  }

  navItems.forEach(btn => {
    btn.addEventListener('click', () => activateTool(btn.dataset.tool));
  });

  // Home tool-cards navigate to tool
  document.querySelectorAll('.tool-card').forEach(card => {
    card.addEventListener('click', () => activateTool(card.dataset.tool));
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

  // ── Shared Sequence Cleaner ────────────────────────────────
  function cleanSeq(seq) {
    return (seq || '').toUpperCase().replace(/[^A-Z\*]/g, '');
  }
  window.cleanSeq = cleanSeq;
  window.BioKit.utils.cleanSeq = cleanSeq;

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
  function validateSequence(seq, type) {
    const raw = seq || '';
    const clean = raw.toUpperCase().replace(/\s/g, '');
    if (!clean) return { valid: false, msg: 'Sequence is empty.' };

    let allowed = /^[ATGCNU]+$/;
    let label = 'DNA/RNA';
    if (type === 'protein') {
      allowed = /^[ACDEFGHIKLMNPQRSTVWY\*]+$/;
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

    return { valid: true, clean: clean.replace(/[^A-Z\*]/g, '') };
  }
  window.validateSequence = validateSequence;
  window.BioKit.utils.validateSequence = validateSequence;

  function showValidationWarning(panelId, validation) {
    const panel = document.getElementById(panelId);
    if (!panel) return;
    clearValidationAlert(panelId);
    const alert = document.createElement('div');
    alert.className = 'validation-alert';
    alert.innerHTML = `
      <div class="validation-alert-icon"><i data-lucide="alert-triangle"></i></div>
      <div class="validation-alert-content">
        <div class="validation-alert-title">Input Validation Error</div>
        <div class="validation-alert-msg">${escapeHTML(validation.msg)}</div>
        ${validation.suggestion ? `<div class="validation-suggestion">${escapeHTML(validation.suggestion)}</div>` : ''}
      </div>
    `;
    const inputBlock = panel.querySelector('.input-block');
    if (inputBlock) inputBlock.parentNode.insertBefore(alert, inputBlock);
    if (window.lucide) window.lucide.createIcons();
    const resultsArea = panel.querySelector('.results-area');
    if (resultsArea) resultsArea.classList.add('hidden');
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
    const panel = e.target.closest('.results-area');
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

  // ── Expose activateTool globally ──────────────────────────
  window.activateTool = activateTool;

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

  // ── Worker Manager ────────────────────────────────────────
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

    runTask(type, payload, onProgress) {
      if (!this.init()) return Promise.reject('No Worker Support');
      const taskId = Math.random().toString(36).substr(2, 9);
      return new Promise((resolve, reject) => {
        this.currentTask = { taskId, resolve, reject, onProgress };
        this.worker.postMessage({ type, payload, taskId });
      });
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
          callback();
        } catch (err) {
          handleError(err, 'Task');
        }
      }, duration + 100);
    }
  };

  // ── API Settings Management ────────────────────────────────
  const settingsOverlay = document.getElementById('settingsOverlay');
  const openSettingsBtn = document.getElementById('openSettingsBtn');
  const closeSettingsBtn = document.getElementById('closeSettingsBtn');
  const saveSettingsBtn = document.getElementById('saveSettingsBtn');
  const ncbiKeyInput = document.getElementById('ncbiApiKey');

  if (ncbiKeyInput) {
    ncbiKeyInput.value = localStorage.getItem('biokit_ncbi_key') || '';
  }

  if (openSettingsBtn) {
    openSettingsBtn.addEventListener('click', () => {
      settingsOverlay.classList.add('visible');
    });
  }

  const closeSettings = () => {
    settingsOverlay.classList.remove('visible');
  };

  if (closeSettingsBtn) closeSettingsBtn.addEventListener('click', closeSettings);

  if (saveSettingsBtn) {
    saveSettingsBtn.addEventListener('click', () => {
      const key = ncbiKeyInput.value.trim();
      localStorage.setItem('biokit_ncbi_key', key);
      showToast('Settings saved successfully!');
      closeSettings();
    });
  }

  if (settingsOverlay) {
    settingsOverlay.addEventListener('click', (e) => {
      if (e.target === settingsOverlay) closeSettings();
    });
  }

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

        if (isCtrl && e.key === 'Enter') {
          e.preventDefault();
          document.querySelector('.tool-panel.active .btn-primary')?.click();
        }
        if (isCtrl && e.key === 'z') {
          e.preventDefault();
          window.BioKit.core.HistoryManager.undo();
        }
        if (isCtrl && e.key === 'y') {
          e.preventDefault();
          window.BioKit.core.HistoryManager.redo();
        }
        if (isCtrl && e.key === 'd') {
          e.preventDefault();
          document.querySelector('.tool-panel.active #proExportBtn')?.click();
        }
        if (e.key === 'Escape') {
          document.querySelectorAll('.settings-overlay').forEach(el => el.classList.add('hidden'));
          settingsOverlay?.classList.remove('visible');
        }
        if (e.key === '?' && !['INPUT', 'TEXTAREA'].includes(document.activeElement.tagName)) {
          const modal = document.getElementById('shortcutModal');
          if (modal) modal.classList.toggle('hidden');
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

    // Mark changes on input
    document.querySelectorAll('input, textarea').forEach(el => {
      el.addEventListener('input', () => window.BioKit.core.SessionManager.markChanged());
    });
  });

})();
