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

      // Read topology toggle — default to linear if toggle element absent
      const isCircular = document.getElementById('reTopologyToggle')?.value === 'circular';

      window.withLoading(toolId, (results) => {
        this.displayResults(seq, results);
      }, seq.length, {
        type: 'RESTRICTION_SEARCH',
        strategy: 'auto',
        payload: { sequence: seq, enzymes: selectedEnzymes, isCircular }
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
      
      // Handle Circular Visualization
      if (document.getElementById('reTopologyToggle')?.value === 'circular') {
        let plasmidContainer = document.getElementById('rePlasmidMap');
        if (!plasmidContainer) {
          plasmidContainer = document.createElement('div');
          plasmidContainer.id = 'rePlasmidMap';
          plasmidContainer.style.marginTop = '20px';
          plasmidContainer.style.display = 'flex';
          plasmidContainer.style.justifyContent = 'center';
          plasmidContainer.style.background = 'var(--bg-secondary)';
          plasmidContainer.style.borderRadius = '8px';
          plasmidContainer.style.padding = '20px';
          
          const seqMap = document.getElementById('reSeqMap');
          if (seqMap && seqMap.parentNode) {
            seqMap.parentNode.insertBefore(plasmidContainer, seqMap.nextSibling);
          }
        }
        this.drawPlasmidMap(seq.length, results, 'rePlasmidMap');
      } else {
        const plasmidContainer = document.getElementById('rePlasmidMap');
        if (plasmidContainer) plasmidContainer.remove();
      }

      this.renderTable(results, fragments);
    },

    drawPlasmidMap: function(sequenceLength, cutResults, containerId) {
      const container = document.getElementById(containerId);
      if (!container) return;
      
      const size = 600;
      const cx = size / 2;
      const cy = size / 2;
      const r = 180; // Radius of backbone
      
      let svg = `<svg width="${size}" height="${size}" viewBox="0 0 ${size} ${size}" xmlns="http://www.w3.org/2000/svg">`;
      
      // Backbone (Dark ring)
      svg += `<circle cx="${cx}" cy="${cy}" r="${r}" fill="none" stroke="var(--border-color)" stroke-width="4"/>`;
      
      // Group by exact cut position logic to mitigate collision
      const groups = {};
      cutResults.forEach(cut => {
        const cutPos = cut.cut;
        if (!groups[cutPos]) groups[cutPos] = [];
        groups[cutPos].push(cut);
      });
      
      for (const pos in groups) {
        const cuts = groups[pos];
        const theta = (parseInt(pos) / sequenceLength) * 2 * Math.PI - Math.PI / 2;
        
        // Base line crossing circle
        const x1 = cx + (r - 10) * Math.cos(theta);
        const y1 = cy + (r - 10) * Math.sin(theta);
        const x2 = cx + (r + 10) * Math.cos(theta);
        const y2 = cy + (r + 10) * Math.sin(theta);
        
        svg += `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="var(--text-secondary)" stroke-width="2"/>`;
        
        // Radially offset labels if multiple enzymes share same coordinates
        cuts.forEach((cut, idx) => {
          // Expand radius safely per duplicate site
          const labelR = r + 20 + (idx * 16);
          const tx = cx + labelR * Math.cos(theta);
          const ty = cy + labelR * Math.sin(theta);
          
          let fill = "var(--text-primary)";
          let decoration = "none";
          let icon = "";
          
          // Methylation Interference visual tracking
          if (cut.isBlocked) {
            fill = "#dc2626"; // Wet lab Warning red
            decoration = "line-through";
            icon = "⛔ ";
          }
          
          let anchor = "middle";
          if (Math.cos(theta) > 0.05) anchor = "start";
          else if (Math.cos(theta) < -0.05) anchor = "end";
          
          // Using strict SVG markup injection (no frameworks payload)
          svg += `<text x="${tx}" y="${ty}" font-family="var(--sans)" font-size="12" fill="${fill}" text-anchor="${anchor}" dominant-baseline="middle" style="text-decoration: ${decoration}; font-weight: 600;">${icon}${cut.name}</text>`;
        });
      }
      
      // Plasmid structural summary in the center
      svg += `<text x="${cx}" y="${cy - 8}" font-family="var(--sans)" font-size="16" fill="var(--text-primary)" text-anchor="middle" font-weight="bold">Plasmid Map</text>`;
      svg += `<text x="${cx}" y="${cy + 15}" font-family="var(--sans)" font-size="14" fill="var(--text-secondary)" text-anchor="middle">${sequenceLength.toLocaleString()} bp</text>`;
      
      svg += `</svg>`;
      container.innerHTML = window.BioKit.utils.sanitizeHTML(svg);
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
            <thead><tr><th>Enzyme</th><th>Site</th><th>Strand</th><th>Position</th><th>5' Cut</th><th>End Type</th><th>Overhang (bp)</th><th>Wet-Lab Warning</th></tr></thead>
            <tbody>
              ${results.map(r => {
                let badge = '';
                if (r.isBlocked) {
                  badge = `<span style="background:#dc2626;color:#fff;padding:2px 7px;border-radius:4px;font-size:0.72rem;font-weight:700;letter-spacing:.04em;" title="${r.warning}">⛔ BLOCKED · ${r.blockType}</span>`;
                }
                
                const rowStyle = r.isBlocked 
                  ? 'background: rgba(220, 38, 38, 0.08); opacity: 0.8;' 
                  : '';
                
                return `<tr style="${rowStyle}">
                  <td>${r.name}</td><td>${r.site}</td><td>${r.strand || '+'}</td>
                  <td>${r.pos}</td><td>${r.cut}</td>
                  <td>${r.endType || '—'}</td>
                  <td>${r.overhang !== undefined ? r.overhang : '—'}</td>
                  <td>${badge || '<span style="color:var(--text-muted);font-size:0.75rem;">—</span>'}</td>
                </tr>`;
              }).join('')}
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
