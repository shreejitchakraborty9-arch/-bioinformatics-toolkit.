/* ============================================================
   structure-viewer.js — 3D Protein Structure Viewer
   pdbe-molstar (Mol*) + AlphaFold EBI API + RCSB PDB
   Native pLDDT confidence coloring via alphafoldView: true
   PDB IDs (4 alphanumeric chars) routed to RCSB BCIF endpoint
   ============================================================ */

(function () {
  'use strict';

  const ALPHAFOLD_API = 'https://alphafold.ebi.ac.uk/api/prediction/';

  let viewerInstance    = null;
  let activeAccession   = null;
  let molstarInjected   = false;

  // ── Mol* Lazy Loader ──────────────────────────────────────
  function ensureMolstar() {
    return new Promise((resolve, reject) => {
      if (typeof window.PDBeMolstarPlugin !== 'undefined') { resolve(); return; }

      if (molstarInjected) {
        let attempts = 0;
        const poll = setInterval(() => {
          if (typeof window.PDBeMolstarPlugin !== 'undefined') { clearInterval(poll); resolve(); }
          if (++attempts > 150) { clearInterval(poll); reject(new Error('Mol* load timeout.')); }
        }, 80);
        return;
      }

      molstarInjected = true;

      const link    = document.createElement('link');
      link.rel      = 'stylesheet';
      link.href     = 'https://cdn.jsdelivr.net/npm/pdbe-molstar@latest/build/pdbe-molstar-light.css';
      link.crossOrigin = 'anonymous';
      document.head.appendChild(link);

      const s         = document.createElement('script');
      s.src           = 'https://cdn.jsdelivr.net/npm/pdbe-molstar@latest/build/pdbe-molstar-plugin.js';
      s.crossOrigin   = 'anonymous';
      s.onload        = () => resolve();
      s.onerror       = () => reject(new Error('Failed to load Mol* library from CDN. Check your internet connection.'));
      document.head.appendChild(s);
    });
  }

  // ── AlphaFold API Fetch ───────────────────────────────────
  async function fetchAlphaFoldEntry(accession) {
    const res = await fetch(
      ALPHAFOLD_API + encodeURIComponent(accession.toUpperCase()),
      { headers: { Accept: 'application/json' } }
    );

    if (res.status === 404) {
      throw new Error(
        `No AlphaFold prediction found for "${accession}". ` +
        `Verify the UniProt accession (e.g. P04637 for TP53).`
      );
    }
    if (!res.ok) {
      throw new Error(`AlphaFold API returned HTTP ${res.status}. Try again shortly.`);
    }

    const data = await res.json();
    if (!Array.isArray(data) || data.length === 0) {
      throw new Error(`AlphaFold returned no predictions for "${accession}".`);
    }
    return data[0];
  }

  // ── Structure Loading Pipeline ────────────────────────────
  async function loadStructure(accession) {
    // 4 alphanumeric characters = PDB ID; anything else = UniProt/AlphaFold
    const isPDB = /^[A-Za-z0-9]{4}$/.test(accession);

    if (isPDB) {
      setOverlay('loading', `Loading PDB structure ${accession.toUpperCase()}…`);
      setStatusBar('Fetching experimental structure from RCSB…');
    } else {
      setOverlay('loading', `Querying AlphaFold DB for ${accession}…`);
      setStatusBar('Fetching prediction metadata…');
    }

    disableBtn(true);
    disableControls(true);

    try {
      const container = document.getElementById('ngl-viewer-container');
      if (!container) throw new Error('Viewer container element not found.');

      container.style.position = 'relative';

      if (viewerInstance) {
        viewerInstance.plugin?.dispose();
        viewerInstance = null;
        container.innerHTML = '';
      }

      viewerInstance = new window.PDBeMolstarPlugin();

      if (isPDB) {
        const id = accession.toLowerCase();
        await viewerInstance.render(container, {
          customData: {
            url:      'https://models.rcsb.org/' + id + '.bcif',
            format:   'bcif',
            isBinary: true,
          },
          bgColor:      { r: 255, g: 255, b: 255 },
          hideControls: false,
        });

        activeAccession = accession;
        setOverlay('ready', '');
        setStatusBar(
          `PDB: ${accession.toUpperCase()} · Source: RCSB PDB · Experimental Structure`
        );
      } else {
        const entry = await fetchAlphaFoldEntry(accession);

        setOverlay('loading', 'Downloading structure file…');

        await viewerInstance.render(container, {
          customData: {
            url:    'https://alphafold.ebi.ac.uk/files/AF-' + accession.toUpperCase() + '-F1-model_v6.cif',
            format: 'cif',
          },
          alphafoldView: true,
          bgColor:       { r: 255, g: 255, b: 255 },
          hideControls:  false,
        });

        activeAccession = accession;
        renderMetadata(entry);
        setOverlay('ready', '');
        setStatusBar(
          `${entry.uniprotAccession || accession}` +
          (entry.gene                   ? ` · Gene: ${entry.gene}` : '') +
          (entry.organismScientificName ? ` · ${entry.organismScientificName}` : '') +
          (entry.latestVersion          ? ` · Model v${entry.latestVersion}` : '')
        );
      }

      // Only unlock controls on successful render
      disableControls(false);

    } catch (err) {
      const msg = err?.message || (typeof err === 'string' ? err : 'Viewer initialization failed');
      setOverlay('ready', '');
      setStatusBar('Error — ' + msg);
      if (window.showToast) window.showToast('⚠ ' + msg);
      console.error('[StructureViewer]', err);
      // Controls intentionally remain disabled on error
    } finally {
      disableBtn(false);
    }
  }

  // ── Metadata Card ─────────────────────────────────────────
  function renderMetadata(entry) {
    const contentEl = document.getElementById('sv-meta-content');
    if (!contentEl) return;

    const row = (label, value) => {
      if (!value) return '';
      const safe = typeof DOMPurify !== 'undefined' ? DOMPurify.sanitize(String(value)) : String(value);
      return `<div class="sv-meta-row">
        <span class="sv-meta-label">${label}</span>
        <span class="sv-meta-value">${safe}</span>
      </div>`;
    };

    const coverage = (entry.uniprotStart != null && entry.uniprotEnd != null)
      ? `${entry.uniprotStart}–${entry.uniprotEnd}`
      : null;

    contentEl.innerHTML = [
      row('UniProt',  entry.uniprotAccession),
      row('Gene',     entry.gene),
      row('Organism', entry.organismScientificName
        ? `<em>${entry.organismScientificName}</em>` : null),
      row('Coverage', coverage),
      row('Version',  entry.latestVersion ? `v${entry.latestVersion}` : null),
      row('Format',   'mmCIF'),
    ].join('');

    document.getElementById('sv-meta')?.classList.remove('hidden');
    document.getElementById('sv-plddt-legend')?.classList.remove('hidden');
  }

  // ── Overlay & Status ──────────────────────────────────────
  function setOverlay(state, message) {
    const overlay = document.getElementById('sv-overlay');
    if (!overlay) return;

    if (state === 'loading') {
      const safe = typeof DOMPurify !== 'undefined' ? DOMPurify.sanitize(message) : message;
      overlay.innerHTML = `<div class="sv-spinner"></div><p class="sv-overlay-text">${safe}</p>`;
      overlay.classList.remove('hidden');
    } else {
      overlay.classList.add('hidden');
    }
  }

  function setStatusBar(msg) {
    const el = document.getElementById('sv-status-bar');
    if (el) el.textContent = msg;
  }

  function disableBtn(disabled) {
    const btn = document.getElementById('sv-fetch-btn');
    if (!btn) return;
    btn.disabled    = disabled;
    btn.textContent = disabled ? 'Loading…' : 'Load Structure';
  }

  // Prompt 2: lock/unlock all three viewer controls together
  function disableControls(disabled) {
    ['sv-repr-select', 'sv-reset-btn', 'sv-screenshot-btn'].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.disabled = disabled;
    });
  }

  // ── Event Wiring ──────────────────────────────────────────
  (function bindEvents() {
    // Disable controls until a structure is loaded
    disableControls(true);

    document.getElementById('sv-fetch-btn')?.addEventListener('click', async () => {
      const acc = document.getElementById('sv-accession-input')?.value.trim();
      if (!acc) {
        if (window.showToast) window.showToast('Enter a UniProt accession (e.g. P04637) or PDB ID (e.g. 6LU7)');
        return;
      }
      try {
        await ensureMolstar();
        await loadStructure(acc);
      } catch (err) {
        const msg = err?.message || (typeof err === 'string' ? err : 'Viewer initialization failed');
        setStatusBar('Error — ' + msg);
        if (window.showToast) window.showToast('⚠ ' + msg);
      }
    });

    document.getElementById('sv-accession-input')?.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') document.getElementById('sv-fetch-btn')?.click();
    });

    document.getElementById('sv-sample-btn')?.addEventListener('click', () => {
      const inp = document.getElementById('sv-accession-input');
      if (inp) inp.value = 'P04637';
      document.getElementById('sv-fetch-btn')?.click();
    });

    // Representation dropdown — awaits the async Mol* update
    document.getElementById('sv-repr-select')?.addEventListener('change', async (e) => {
      if (!viewerInstance || !viewerInstance.plugin) { console.warn('Viewer not ready'); return; }
      const map = {
        'cartoon':    'cartoon',
        'ribbon':     'cartoon',       // Mol* has no separate ribbon; cartoon is equivalent
        'surface':    'molecular-surface',
        'ball-stick': 'ball-and-stick',
        'spacefill':  'spacefill',
      };
      const styleName = map[e.target.value] ?? e.target.value;
      await viewerInstance.visual?.update({ style: { name: styleName } });
    });

    document.getElementById('sv-reset-btn')?.addEventListener('click', () => {
      if (!viewerInstance || !viewerInstance.plugin) { console.warn('Viewer not ready'); return; }
      viewerInstance.plugin.canvas3d?.requestCameraReset();
    });

    document.getElementById('sv-screenshot-btn')?.addEventListener('click', () => {
      if (!viewerInstance || !viewerInstance.plugin) { console.warn('Viewer not ready'); return; }
      viewerInstance.plugin.helpers?.viewportScreenshot?.download({ filename: 'bionised_structure.png' });
    });
  })();

  // ── Panel Lifecycle ───────────────────────────────────────
  if (window.ViewManager) {
    window.ViewManager.registerUnmount('structure', () => {
      if (viewerInstance) {
        const container = document.getElementById('ngl-viewer-container');
        if (container) container.innerHTML = '';
        viewerInstance  = null;
        activeAccession = null;
      }
    });
  }

})();
