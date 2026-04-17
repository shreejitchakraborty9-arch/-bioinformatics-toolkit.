/* ============================================================
   structure-viewer.js — 3D Protein Structure Viewer
   NGL Viewer (industry-standard, RCSB-used) + AlphaFold EBI API
   pLDDT confidence coloring per official AlphaFold color scale
   ============================================================ */

(function () {
  'use strict';

  const ALPHAFOLD_API = 'https://alphafold.ebi.ac.uk/api/prediction/';

  let stage          = null;   // NGL.Stage instance
  let activeComp     = null;   // currently loaded NGL component
  let activeAccession = null;
  let nglScriptInjected = false;

  // ── NGL Lazy Loader ───────────────────────────────────────
  // NGL (~5 MB) is loaded on first user interaction to keep
  // the initial page load fast.
  function ensureNGL() {
    return new Promise((resolve, reject) => {
      if (typeof NGL !== 'undefined') { resolve(); return; }

      if (nglScriptInjected) {
        // Already injected — wait for it to finish parsing
        let attempts = 0;
        const poll = setInterval(() => {
          if (typeof NGL !== 'undefined') { clearInterval(poll); resolve(); }
          if (++attempts > 100) { clearInterval(poll); reject(new Error('NGL load timeout.')); }
        }, 80);
        return;
      }

      nglScriptInjected = true;
      const s   = document.createElement('script');
      s.src     = 'https://unpkg.com/ngl@2.0.0-dev.37/dist/ngl.js';
      s.onload  = () => resolve();
      s.onerror = () => reject(new Error('Failed to load NGL Viewer library from CDN. Check your internet connection.'));
      document.head.appendChild(s);
    });
  }

  // ── AlphaFold pLDDT Color Scheme Registration ─────────────
  // Official AlphaFold pLDDT confidence scale (stored as B-factor):
  //   ≥90  → #0053D6  Very high confidence
  //   ≥70  → #65CBF3  Confident
  //   ≥50  → #FFDB13  Low confidence
  //   <50  → #FF7D45  Very low confidence
  function registerPLDDTScheme() {
    if (NGL.ColormakerRegistry.hasScheme('plddt')) return;
    NGL.ColormakerRegistry.addScheme(function () {
      this.atomColor = function (atom) {
        const b = atom.bfactor;
        if (b >= 90) return 0x0053D6;
        if (b >= 70) return 0x65CBF3;
        if (b >= 50) return 0xFFDB13;
        return 0xFF7D45;
      };
    }, 'plddt');
  }

  // ── Stage Initialization ──────────────────────────────────
  function initStage() {
    const container = document.getElementById('ngl-viewer-container');
    if (!container || stage) return;

    const dark = document.documentElement.classList.contains('dark');
    stage = new NGL.Stage(container, {
      backgroundColor: dark ? '#101214' : '#f8f9fa',
      quality:         'medium',
      impostor:        true,
      tooltip:         false,
    });

    registerPLDDTScheme();

    // Keep canvas filling container on resize
    if (window.ResizeObserver) {
      new ResizeObserver(() => { if (stage) stage.handleResize(); }).observe(container);
    }
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
    return data[0]; // Most recent model version
  }

  // ── Structure Loading Pipeline ────────────────────────────
  async function loadStructure(accession) {
    setOverlay('loading', `Querying AlphaFold DB for ${accession}…`);
    setStatusBar(`Fetching prediction metadata…`);
    disableBtn(true);

    try {
      const entry = await fetchAlphaFoldEntry(accession);

      setOverlay('loading', 'Downloading structure file…');

      // Remove previous structure
      if (activeComp) {
        stage.removeComponent(activeComp);
        activeComp = null;
      }

      // Prefer mmCIF (canonical AlphaFold format); fall back to PDB
      const url = entry.cifUrl || entry.pdbUrl;
      const ext = entry.cifUrl ? 'cif' : 'pdb';

      if (!url) throw new Error('AlphaFold entry has no structure file URL.');

      const comp = await stage.loadFile(url, { ext, name: accession });
      activeComp      = comp;
      activeAccession = accession;

      // Apply default cartoon + pLDDT representation
      const reprVal = document.getElementById('sv-repr-select')?.value || 'cartoon';
      applyRepresentation(reprVal);

      stage.autoView(600);

      renderMetadata(entry);
      setOverlay('ready', '');
      setStatusBar(
        `${entry.uniprotAccession || accession}` +
        (entry.gene                 ? ` · Gene: ${entry.gene}` : '') +
        (entry.organismScientificName ? ` · ${entry.organismScientificName}` : '') +
        (entry.latestVersion        ? ` · Model v${entry.latestVersion}` : '')
      );
    } catch (err) {
      setOverlay('ready', '');
      setStatusBar('Error — ' + err.message);
      if (window.showToast) window.showToast('⚠ ' + err.message);
      console.error('[StructureViewer]', err);
    } finally {
      disableBtn(false);
    }
  }

  // ── Representation ────────────────────────────────────────
  function applyRepresentation(reprType) {
    if (!activeComp) return;
    activeComp.removeAllRepresentations();

    const plddt = { colorScheme: 'plddt' };

    switch (reprType) {
      case 'cartoon':
        activeComp.addRepresentation('cartoon', {
          ...plddt, smoothSheet: true, aspectRatio: 4.0, quality: 'high',
        });
        break;
      case 'ribbon':
        activeComp.addRepresentation('ribbon', { ...plddt, quality: 'high' });
        break;
      case 'surface':
        activeComp.addRepresentation('surface', { ...plddt, opacity: 0.82, useWorker: true });
        break;
      case 'ball-stick':
        activeComp.addRepresentation('ball+stick', { colorScheme: 'element', multipleBond: true });
        break;
      case 'spacefill':
        activeComp.addRepresentation('spacefill', { ...plddt, radiusScale: 0.6 });
        break;
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
      row('UniProt', entry.uniprotAccession),
      row('Gene',    entry.gene),
      row('Organism', entry.organismScientificName
        ? `<em>${entry.organismScientificName}</em>` : null),
      row('Coverage', coverage),
      row('Version',  entry.latestVersion ? `v${entry.latestVersion}` : null),
      row('Format',   entry.cifUrl ? 'mmCIF' : 'PDB'),
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

  // ── Event Wiring ──────────────────────────────────────────
  (function bindEvents() {
    // Fetch button
    document.getElementById('sv-fetch-btn')?.addEventListener('click', async () => {
      const acc = document.getElementById('sv-accession-input')?.value.trim();
      if (!acc) {
        if (window.showToast) window.showToast('Enter a UniProt accession (e.g. P04637)');
        return;
      }
      try {
        await ensureNGL();
        if (!stage) initStage();
        await loadStructure(acc);
      } catch (err) {
        setStatusBar('Error — ' + err.message);
        if (window.showToast) window.showToast('⚠ ' + err.message);
      }
    });

    // Enter key
    document.getElementById('sv-accession-input')?.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') document.getElementById('sv-fetch-btn')?.click();
    });

    // Representation selector
    document.getElementById('sv-repr-select')?.addEventListener('change', (e) => {
      applyRepresentation(e.target.value);
    });

    // Reset camera
    document.getElementById('sv-reset-btn')?.addEventListener('click', () => {
      stage?.autoView(500);
    });

    // High-res screenshot
    document.getElementById('sv-screenshot-btn')?.addEventListener('click', () => {
      if (!stage || !activeAccession) return;
      stage.makeImage({ factor: 2, antialias: true, trim: false, transparent: false })
        .then(blob => {
          const a   = document.createElement('a');
          a.href    = URL.createObjectURL(blob);
          a.download = `${activeAccession}_alphafold.png`;
          a.click();
          URL.revokeObjectURL(a.href);
        });
    });

    // TP53 example
    document.getElementById('sv-sample-btn')?.addEventListener('click', async () => {
      const inp = document.getElementById('sv-accession-input');
      if (inp) inp.value = 'P04637';
      document.getElementById('sv-fetch-btn')?.click();
    });

    // Theme sync — update NGL background color on theme change
    window.addEventListener('biokit-theme-change', () => {
      if (!stage) return;
      const dark = document.documentElement.classList.contains('dark');
      stage.setParameters({ backgroundColor: dark ? '#101214' : '#f8f9fa' });
    });
  })();

  // ── Panel Lifecycle ───────────────────────────────────────
  // Initialise stage lazily when panel becomes active (avoids wasting
  // resources when user never visits this tool).
  const panel = document.getElementById('panel-structure');
  if (panel) {
    new MutationObserver(() => {
      if (panel.classList.contains('active') && !stage && typeof NGL !== 'undefined') {
        initStage();
      }
    }).observe(panel, { attributes: true, attributeFilter: ['class'] });
  }

  // Tear down NGL stage to free GPU memory when leaving the panel
  if (window.ViewManager) {
    window.ViewManager.registerUnmount('structure', () => {
      if (stage) {
        stage.dispose();
        stage      = null;
        activeComp = null;
      }
    });
  }

})();
