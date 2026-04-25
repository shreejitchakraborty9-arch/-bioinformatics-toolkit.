(function () {
  'use strict';

  window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };

  window.BioKit.tools.restriction = {
    run: async function () {
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

      if (seq.length > 100000) {
        window.showToast('Error: Sequence exceeds 100,000 bp limit.');
        return;
      }

      const topology = document.getElementById('reTopologyToggle')?.value || 'linear';
      const btn = document.getElementById('reAnalyzeBtn');

      if (btn) { btn.disabled = true; btn.textContent = 'Analyzing…'; }

      try {
        const response = await fetch('/api/analyze/restriction', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ sequence: seq, topology, enzymes: selectedEnzymes })
        });

        const data = await response.json().catch(() => ({}));

        if (!response.ok || data.error) {
          window.showToast(data.error || 'Analysis failed. Please try again.');
          return;
        }

        this.displayResults(seq, topology, data);
      } catch (err) {
        window.showToast('Network error. Could not reach the server.');
        console.error(err);
      } finally {
        if (btn) { btn.disabled = false; btn.textContent = 'Analyze'; }
      }
    },

    getSelectedEnzymes: function () {
      const chips = document.querySelectorAll('.enzyme-chip.selected');
      return Array.from(chips).map(c => c.dataset.enzyme);
    },

    displayResults: function (seq, topology, data) {
      const container = document.getElementById('reResults');
      if (!container) return;
      container.classList.remove('hidden');

      const results = data.results || [];
      const totalCuts = results.reduce((sum, r) => sum + r.cuts.length, 0);
      const totalFragments = results.reduce((sum, r) => sum + r.fragments.length, 0);

      const statRow = document.getElementById('reStatRow');
      if (statRow) {
        statRow.innerHTML = window.BioKit.utils.sanitizeHTML(`
          <div class="stat-card"><h3>${totalCuts}</h3><p>Cut Sites</p></div>
          <div class="stat-card"><h3>${results.length}</h3><p>Enzymes</p></div>
          <div class="stat-card"><h3>${totalFragments}</h3><p>Fragments</p></div>
        `);
      }

      this.renderSequenceMap(seq);

      if (topology === 'circular') {
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

        // Flatten per-enzyme cuts into the shape drawPlasmidMap expects
        const flatCuts = results.flatMap(r => r.cuts.map(pos => ({ cut: pos, name: r.enzyme })));
        this.drawPlasmidMap(seq.length, flatCuts, 'rePlasmidMap');

        let exportBtn = document.getElementById('reExportPdfBtn');
        if (!exportBtn) {
          exportBtn = document.createElement('button');
          exportBtn.id = 'reExportPdfBtn';
          exportBtn.className = 'btn-primary';
          exportBtn.style.marginTop = '15px';
          exportBtn.style.alignSelf = 'center';
          exportBtn.innerHTML = 'Export Publication PDF';
          exportBtn.addEventListener('click', () => {
            const svgNode = plasmidContainer.querySelector('svg');
            if (svgNode) {
              const seqNameMatch = document.querySelector('.session-name')?.textContent || 'plasmid';
              const cleanName = seqNameMatch.replace(/[^a-z0-9]/gi, '_').toLowerCase();
              this.exportPlasmidPDF(svgNode, cleanName);
            }
          });
          plasmidContainer.style.flexDirection = 'column';
          plasmidContainer.appendChild(exportBtn);
        }
      } else {
        const plasmidContainer = document.getElementById('rePlasmidMap');
        if (plasmidContainer) plasmidContainer.remove();
      }

      this.renderTable(results);
    },

    drawPlasmidMap: function (sequenceLength, cutResults, containerId) {
      const container = document.getElementById(containerId);
      if (!container) return;
      container.innerHTML = '';

      const size = 600;
      const cx = size / 2;
      const cy = size / 2;
      const r = 180;
      const xmlns = 'http://www.w3.org/2000/svg';

      const svg = document.createElementNS(xmlns, 'svg');
      svg.setAttribute('width', size);
      svg.setAttribute('height', size);
      svg.setAttribute('viewBox', `0 0 ${size} ${size}`);

      const circle = document.createElementNS(xmlns, 'circle');
      circle.setAttribute('cx', cx);
      circle.setAttribute('cy', cy);
      circle.setAttribute('r', r);
      circle.setAttribute('fill', 'none');
      circle.setAttribute('stroke', 'var(--border-color)');
      circle.setAttribute('stroke-width', '4');
      svg.appendChild(circle);

      const groups = {};
      cutResults.forEach(cut => {
        const key = cut.cut;
        if (!groups[key]) groups[key] = [];
        groups[key].push(cut);
      });

      for (const pos in groups) {
        const cuts = groups[pos];
        const theta = (parseInt(pos) / sequenceLength) * 2 * Math.PI - Math.PI / 2;

        const x1 = cx + (r - 10) * Math.cos(theta);
        const y1 = cy + (r - 10) * Math.sin(theta);
        const x2 = cx + (r + 10) * Math.cos(theta);
        const y2 = cy + (r + 10) * Math.sin(theta);

        const line = document.createElementNS(xmlns, 'line');
        line.setAttribute('x1', x1); line.setAttribute('y1', y1);
        line.setAttribute('x2', x2); line.setAttribute('y2', y2);
        line.setAttribute('stroke', 'var(--text-secondary)');
        line.setAttribute('stroke-width', '2');
        svg.appendChild(line);

        cuts.forEach((cut, idx) => {
          const labelR = r + 20 + (idx * 16);
          const tx = cx + labelR * Math.cos(theta);
          const ty = cy + labelR * Math.sin(theta);

          let anchor = 'middle';
          if (Math.cos(theta) > 0.05) anchor = 'start';
          else if (Math.cos(theta) < -0.05) anchor = 'end';

          const text = document.createElementNS(xmlns, 'text');
          text.setAttribute('x', tx);
          text.setAttribute('y', ty);
          text.setAttribute('font-family', 'var(--sans)');
          text.setAttribute('font-size', '12');
          text.setAttribute('fill', 'var(--text-primary)');
          text.setAttribute('text-anchor', anchor);
          text.setAttribute('dominant-baseline', 'middle');
          text.setAttribute('style', 'font-weight: 600;');
          text.textContent = cut.name;
          svg.appendChild(text);
        });
      }

      const titleText = document.createElementNS(xmlns, 'text');
      titleText.setAttribute('x', cx); titleText.setAttribute('y', cy - 8);
      titleText.setAttribute('font-family', 'var(--sans)');
      titleText.setAttribute('font-size', '16');
      titleText.setAttribute('fill', 'var(--text-primary)');
      titleText.setAttribute('text-anchor', 'middle');
      titleText.setAttribute('font-weight', 'bold');
      titleText.textContent = 'Plasmid Map';
      svg.appendChild(titleText);

      const subtitleText = document.createElementNS(xmlns, 'text');
      subtitleText.setAttribute('x', cx); subtitleText.setAttribute('y', cy + 15);
      subtitleText.setAttribute('font-family', 'var(--sans)');
      subtitleText.setAttribute('font-size', '14');
      subtitleText.setAttribute('fill', 'var(--text-secondary)');
      subtitleText.setAttribute('text-anchor', 'middle');
      subtitleText.textContent = `${sequenceLength.toLocaleString()} bp`;
      svg.appendChild(subtitleText);

      container.appendChild(svg);
    },

    exportPlasmidPDF: async function (svgNode, sequenceName) {
      const btn = document.getElementById('reExportPdfBtn');
      if (btn) btn.textContent = 'Generating PDF… Loading Heavy Modules…';

      const loadScript = (src, checkGlobal) => {
        return new Promise((resolve, reject) => {
          if (window[checkGlobal]) return resolve(window[checkGlobal]);
          const existing = document.querySelector(`script[src="${src}"]`);
          if (existing) {
            const interval = setInterval(() => {
              if (window[checkGlobal]) { clearInterval(interval); resolve(window[checkGlobal]); }
            }, 100);
            return;
          }
          const s = document.createElement('script');
          s.src = src;
          s.onload = () => resolve(window[checkGlobal]);
          s.onerror = reject;
          document.head.appendChild(s);
        });
      };

      try {
        await loadScript('https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js', 'jspdf');
        await loadScript('https://cdnjs.cloudflare.com/ajax/libs/svg2pdf.js/2.2.4/svg2pdf.umd.min.js', 'svg2pdf');

        const { jsPDF } = window.jspdf;
        const clone = svgNode.cloneNode(true);
        const originalElements = svgNode.querySelectorAll('*');
        const cloneElements = clone.querySelectorAll('*');

        for (let i = 0; i < originalElements.length; i++) {
          const computed = window.getComputedStyle(originalElements[i]);
          const cloneEl = cloneElements[i];
          cloneEl.style.fill = computed.fill;
          cloneEl.style.stroke = computed.stroke;
          cloneEl.style.strokeWidth = computed.strokeWidth;
          cloneEl.style.fontFamily = computed.fontFamily;
          cloneEl.style.fontSize = computed.fontSize;
          cloneEl.style.fontWeight = computed.fontWeight;
          cloneEl.style.textAnchor = computed.textAnchor;
          cloneEl.style.dominantBaseline = computed.dominantBaseline;
          cloneEl.style.textDecoration = computed.textDecorationLine || computed.textDecoration;
          cloneEl.style.color = computed.color;
        }

        const doc = new jsPDF({ orientation: 'portrait', unit: 'pt', format: 'letter' });
        const pageWidth = doc.internal.pageSize.getWidth();
        const pageHeight = doc.internal.pageSize.getHeight();
        const svgW = parseInt(clone.getAttribute('width')) || 600;
        const svgH = parseInt(clone.getAttribute('height')) || 600;

        await doc.svg(clone, {
          x: (pageWidth - svgW) / 2,
          y: (pageHeight - svgH) / 2,
          width: svgW,
          height: svgH
        });

        const dateStr = new Date().toISOString().split('T')[0];
        doc.setFontSize(9);
        doc.setTextColor(120);
        doc.text(`Generated by BioToolkit ${dateStr}`, pageWidth - 30, pageHeight - 30, { align: 'right' });
        doc.save(`${sequenceName}_plasmid_map.pdf`);
      } catch (err) {
        console.error('PDF generation failed:', err);
        window.showToast('Failed to export PDF.');
      } finally {
        if (btn) btn.textContent = 'Export Publication PDF';
      }
    },

    renderSequenceMap: function (seq) {
      const mapEl = document.getElementById('reSeqMap');
      if (!mapEl) return;
      mapEl.innerHTML = window.BioKit.utils.sanitizeHTML(
        `<div class="re-map-track"><pre style="font-family:var(--mono); font-size:12px; line-height:1.2; letter-spacing:1px; color:var(--text-secondary);">${seq}</pre></div>`
      );
    },

    renderTable: function (results) {
      const tableEl = document.getElementById('reTable');
      if (!tableEl) return;

      let html = '';

      if (results.length === 0) {
        html = '<p>No cut sites found for the selected enzymes.</p>';
      } else {
        const rows = results.map(r => {
          const cuts = Array.isArray(r.cuts) ? r.cuts : [];
          const fragments = Array.isArray(r.fragments) ? r.fragments : [];
          const cutPositions = cuts.length > 0 ? cuts.join(', ') : 'None';
          const fragmentSizes = fragments.slice().sort((a, b) => b - a).join(', ');
          return `
            <tr>
              <td><strong>${r.enzyme}</strong></td>
              <td>${cuts.length}</td>
              <td>${cutPositions}</td>
              <td>${fragmentSizes}</td>
            </tr>`;
        }).join('');

        html = `
          <table class="data-table">
            <thead>
              <tr>
                <th>Enzyme</th>
                <th>Number of Cuts</th>
                <th>Cut Positions</th>
                <th>Fragment Sizes (bp)</th>
              </tr>
            </thead>
            <tbody>${rows}</tbody>
          </table>`;
      }

      tableEl.innerHTML = window.BioKit.utils.sanitizeHTML(html);
    }
  };

  const initEnzymes = () => {
    const allEnzymes = window.BioKit.data.RESTRICTION_ENZYMES || window.RESTRICTION_ENZYMES || [];
    const chipsContainer = document.getElementById('enzymeChips');
    if (chipsContainer && chipsContainer.children.length === 0) {
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

      document.getElementById('reSelectAll')?.addEventListener('click', () => {
        chipsContainer.querySelectorAll('.enzyme-chip:not([style*="display: none"])').forEach(c => c.classList.add('selected'));
      });
      document.getElementById('reSelectNone')?.addEventListener('click', () => {
        chipsContainer.querySelectorAll('.enzyme-chip').forEach(c => c.classList.remove('selected'));
      });

      const filterBtns = document.querySelectorAll('.enzyme-filter-btn');
      filterBtns.forEach(btn => {
        btn.addEventListener('click', () => {
          filterBtns.forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          const filter = btn.dataset.filter;
          chipsContainer.querySelectorAll('.enzyme-chip').forEach(chip => {
            const enz = allEnzymes.find(x => x.name === chip.dataset.enzyme);
            if (!enz) return;
            let show = false;
            if (filter === 'all') show = true;
            else if (filter === 'common') show = enz.is_common;
            else if (filter === '6-cutter') show = enz.site.length === 6;
            else if (filter === '4-cutter') show = enz.site.length === 4;
            else if (filter === 'rare') show = enz.site.length >= 8;
            else if (filter === 'blunt') show = enz.cut === Math.floor(enz.site.length / 2);
            else if (filter === 'sticky') show = enz.cut !== Math.floor(enz.site.length / 2);
            chip.style.display = show ? 'inline-block' : 'none';
          });
        });
      });
    }
  };

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initEnzymes);
  } else {
    initEnzymes();
  }

  document.getElementById('reAnalyzeBtn')?.addEventListener('click', window.debounce(() => window.BioKit.tools.restriction.run(), 300));
  document.getElementById('reSampleBtn')?.addEventListener('click', () => {
    const input = document.getElementById('reInput');
    if (input) {
      input.value = 'GAATTC' + 'A'.repeat(10) + 'GGTACC';
      input.dispatchEvent(new Event('input'));
    }
  });

})();
