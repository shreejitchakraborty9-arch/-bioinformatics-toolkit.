(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };
  
  // IUPAC code to Regex mapping for DNA
  const IUPAC_REGEX = {
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ATGC]'
  };

  window.BioKit.tools.restriction = {
    run: function() {
      const toolId = 'panel-restriction';
      const raw = document.getElementById('reInput').value;
      const validation = window.validateSequence(raw, 'dna');

      if (!validation.valid) {
        window.showToast(validation.msg);
        return;
      }

      const seq = validation.clean;
      const selectedEnzymes = this.getSelectedEnzymes();
      
      if (selectedEnzymes.length === 0) {
        window.showToast('Please select at least one enzyme.');
        return;
      }

      window.withLoading(toolId, (results) => {
        this.displayResults(seq, results);
      }, seq.length, {
        type: 'RESTRICTION_SEARCH',
        strategy: 'auto',
        payload: { sequence: seq, enzymes: selectedEnzymes }
      });
    },

    getSelectedEnzymes: function() {
      const chips = document.querySelectorAll('.enzyme-chip.selected');
      return Array.from(chips).map(c => c.dataset.enzyme);
    },

    calculateFragments: function(seqLength, results) {
      if (results.length === 0) return [seqLength];
      
      // Get all unique sorted cut positions on the forward strand
      const cutPositions = [...new Set(results.map(r => r.cut))].sort((a, b) => a - b);
      // Filter out cuts outside sequence bounds
      const validCuts = cutPositions.filter(c => c > 0 && c < seqLength);
      
      if (validCuts.length === 0) return [seqLength];

      const fragments = [];
      let lastCut = 0;
      validCuts.forEach(cut => {
        fragments.push(cut - lastCut);
        lastCut = cut;
      });
      fragments.push(seqLength - lastCut);
      return fragments.sort((a, b) => b - a); // Largest to smallest
    },

    displayResults: function(seq, results) {
      const container = document.getElementById('reResults');
      if (!container) return;
      container.classList.remove('hidden');
      
      const fragments = this.calculateFragments(seq.length, results);
      
      const statRow = document.getElementById('reStatRow');
      if (statRow) {
        statRow.innerHTML = window.BioKit.utils.sanitizeHTML(`
          <div class="stat-card"><h3>${results.length}</h3><p>Cut Sites</p></div>
          <div class="stat-card"><h3>${new Set(results.map(r => r.name)).size}</h3><p>Enzymes</p></div>
          <div class="stat-card"><h3>${fragments.length}</h3><p>Fragments</p></div>
        `);
      }

      this.renderSequenceMap(seq, results);
      this.renderTable(results, fragments);
    },

    renderSequenceMap: function(seq, results) {
      const mapEl = document.getElementById('reSeqMap');
      if (!mapEl) return;
      let html = `<div class="re-map-track">`;
      html += `<pre style="font-family:var(--mono); font-size:12px; line-height:1.2; letter-spacing:1px; color:var(--text-secondary);">${seq}</pre>`;
      html += `</div>`;
      mapEl.innerHTML = window.BioKit.utils.sanitizeHTML(html);
    },

    renderTable: function(results, fragments) {
      const tableEl = document.getElementById('reTable');
      if (!tableEl) return;
      
      let html = '';
      
      if (results.length === 0) {
        html += '<p>No sites found.</p>';
      } else {
        html += `
          <table class="data-table">
            <thead><tr><th>Enzyme</th><th>Site</th><th>Strand</th><th>Position</th><th>Cut</th></tr></thead>
            <tbody>
              ${results.map(r => `<tr><td>${r.name}</td><td>${r.site}</td><td>${r.strand || '+'}</td><td>${r.pos}</td><td>${r.cut}</td></tr>`).join('')}
            </tbody>
          </table>
        `;
      }
      
      // Add Fragments list below table
      if (fragments.length > 0) {
        html += `
          <div style="margin-top:20px;">
            <h4 style="color:var(--text-primary); margin-bottom:10px;">Predicted Fragments (bp)</h4>
            <div style="display:flex; flex-wrap:wrap; gap:8px;">
              ${fragments.map(f => `<span class="enzyme-chip" style="background:var(--bg-tertiary); cursor:default;">${f.toLocaleString()} bp</span>`).join('')}
            </div>
          </div>
        `;
      }

      tableEl.innerHTML = window.BioKit.utils.sanitizeHTML(html);
    }
  };

  // Enzyme DB Initialization
  document.addEventListener('DOMContentLoaded', () => {
    const allEnzymes = window.BioKit.data.RESTRICTION_ENZYMES || window.RESTRICTION_ENZYMES || [];
    const chipsContainer = document.getElementById('enzymeChips');
    if (chipsContainer) {
      const fragment = document.createDocumentFragment();
      allEnzymes.forEach(e => {
        const chip = document.createElement('div');
        chip.className = 'enzyme-chip selected';
        chip.dataset.enzyme = e.name;
        chip.textContent = e.name;
        fragment.appendChild(chip);
      });
      chipsContainer.textContent = '';
      chipsContainer.appendChild(fragment);
      
      chipsContainer.addEventListener('click', (e) => {
        const chip = e.target.closest('.enzyme-chip');
        if (chip) chip.classList.toggle('selected');
      });

      // Select All / Select None
      document.getElementById('reSelectAll')?.addEventListener('click', () => {
        chipsContainer.querySelectorAll('.enzyme-chip:not([style*="display: none"])').forEach(c => c.classList.add('selected'));
      });
      document.getElementById('reSelectNone')?.addEventListener('click', () => {
        chipsContainer.querySelectorAll('.enzyme-chip').forEach(c => c.classList.remove('selected'));
      });
      
      // Filter logic
      const filterBtns = document.querySelectorAll('.enzyme-filter-btn');
      filterBtns.forEach(btn => {
        btn.addEventListener('click', (e) => {
          filterBtns.forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          const filter = btn.dataset.filter;
          
          const chips = chipsContainer.querySelectorAll('.enzyme-chip');
          chips.forEach(chip => {
            const enz = allEnzymes.find(x => x.name === chip.dataset.enzyme);
            if (!enz) return;
            
            let show = false;
            if (filter === 'all') show = true;
            else if (filter === 'common') show = enz.common;
            else if (filter === '6-cutter') show = enz.site.length === 6;
            else if (filter === '4-cutter') show = enz.site.length === 4;
            else if (filter === 'rare') show = enz.site.length >= 8;
            else if (filter === 'blunt') show = enz.cut === Math.floor(enz.site.length / 2); // Approximate blunt check
            else if (filter === 'sticky') show = enz.cut !== Math.floor(enz.site.length / 2);
            
            chip.style.display = show ? 'inline-block' : 'none';
          });
        });
      });
    }
  });

  document.getElementById('reAnalyzeBtn')?.addEventListener('click', window.debounce(() => window.BioKit.tools.restriction.run(), 300));
  document.getElementById('reSampleBtn')?.addEventListener('click', () => {
    const input = document.getElementById('reInput');
    if (input) {
      input.value = 'GAATTC' + 'A'.repeat(10) + 'GGTACC';
      input.dispatchEvent(new Event('input'));
    }
  });

})();
