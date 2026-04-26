/* ============================================================
 batch.js — High-Throughput Omics Batch Queue
 File upload → async Celery job → short-poll /api/status →
 progress bar → CSV bulk download
 ============================================================ */

(function () {
  'use strict';

  const POLL_INTERVAL_MS = 2000; // 2-second short-poll cadence

  let currentJobId = null;
  let pollTimer = null;
  let uploadedContent = null; // raw FASTA string from file or paste

  // ── File / Drop-Zone Wiring ───────────────────────────────
  const dropZone = document.getElementById('batchDropZone');
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
      window.showToast?.('File too large. Maximum size is 5 MB.');
      return;
    }
    const reader = new FileReader();
    reader.onload = (e) => {
      uploadedContent = e.target.result;
      const pasteArea = document.getElementById('batchPasteInput');
      if (pasteArea) pasteArea.value = uploadedContent;
      window.showToast?.(`Loaded: ${file.name}`);
    };
    reader.readAsText(file);
  }

  // ── Submit Button ─────────────────────────────────────────
  document.getElementById('batchSubmitBtn')?.addEventListener('click', submitBatch);

  async function submitBatch() {
    const pasteArea = document.getElementById('batchPasteInput');
    const content = (uploadedContent || pasteArea?.value || '').trim();

    if (!content) {
      window.showToast?.('Upload a FASTA file or paste sequences first.');
      return;
    }

    // Stop any prior poll
    stopPolling();
    resetJobUI();

    setBatchStatus('Submitting job…', false);

    try {
      const formData = new FormData();
      const blob = new Blob([content], { type: 'text/plain' });
      formData.append('file', blob, 'batch.fasta');

      const res = await fetch('/api/analyze/batch', {
        method: 'POST',
        body: formData,
      });

      if (!res.ok) {
        const err = await res.json().catch(() => ({ error: `HTTP ${res.status}` }));
        throw new Error(err.error || `Server error ${res.status}`);
      }

      const data = await res.json();
      currentJobId = data.job_id;

      showJobArea(currentJobId);
      startPolling(currentJobId);
    } catch (err) {
      setBatchStatus('Submission failed: ' + err.message, true);
      window.showToast?.('⚠ ' + err.message);
    }
  }

  // ── Polling ───────────────────────────────────────────────
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
      if (!res.ok) return; // transient HTTP error — try again next tick
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

  // ── UI Updates ────────────────────────────────────────────
  function showJobArea(jobId) {
    const area = document.getElementById('batchJobArea');
    if (area) area.classList.remove('hidden');
    const idEl = document.getElementById('batchJobId');
    if (idEl) idEl.textContent = `JOB: ${jobId.slice(0, 8).toUpperCase()}…`;
    document.getElementById('batchDownloadArea')?.classList.add('hidden');
    setProgressBar(0);
  }

  function resetJobUI() {
    document.getElementById('batchJobArea')?.classList.add('hidden');
    document.getElementById('batchDownloadArea')?.classList.add('hidden');
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

    const statusMap = { pending: 'Queued — waiting for worker…', running: null, complete: 'Complete', failed: 'Failed', failure: 'Failed' };
    const statusText = statusMap[data.status] ?? data.status;

    if (data.status === 'running') {
      const processed = data.processed || 0;
      const total = data.total || 0;
      setBatchStatus(
        total > 0
          ? `Processing ${processed.toLocaleString()} / ${total.toLocaleString()} sequences…`
          : 'Processing…',
        false
      );
    } else if (statusText) {
      setBatchStatus(statusText, data.status === 'failed' || data.status === 'failure');
    }

    const pctEl = document.getElementById('batchProgressPct');
    if (pctEl) pctEl.textContent = `${pct}%`;

    const labelEl = document.getElementById('batchProgressLabel');
    if (labelEl && data.total > 0 && data.status === 'running') {
      labelEl.textContent = `${data.processed.toLocaleString()} of ${data.total.toLocaleString()} sequences processed`;
    }

    // Show/hide spinner
    const spinner = document.getElementById('batchSpinner');
    if (spinner) {
      spinner.style.display = (data.status === 'running' || data.status === 'pending') ? '' : 'none';
    }
  }

  function setProgressBar(pct) {
    const bar = document.getElementById('batchProgressBar');
    if (bar) bar.style.width = `${pct}%`;
  }

  function onJobComplete(jobId, data) {
    setProgressBar(100);
    setBatchStatus('Complete', false);
    const spinner = document.getElementById('batchSpinner');
    if (spinner) spinner.style.display = 'none';

    const summary = document.getElementById('batchCompleteSummary');
    if (summary) {
      const total = data.total || 0;
      const failed = data.failed || 0;
      summary.textContent = `${total.toLocaleString()} sequences processed${failed > 0 ? `, ${failed} failed` : ''}.`;
    }

    const downloadArea = document.getElementById('batchDownloadArea');
    if (downloadArea) downloadArea.classList.remove('hidden');

    // Wire download button for this job
    const btn = document.getElementById('batchDownloadBtn');
    if (btn) {
      btn.onclick = () => triggerDownload(jobId);
    }

    window.showToast?.('Batch analysis complete — ready to download.');
  }

  function onJobFailed(errMsg) {
    setBatchStatus('Job failed: ' + errMsg, true);
    const spinner = document.getElementById('batchSpinner');
    if (spinner) spinner.style.display = 'none';
    window.showToast?.('⚠ Batch job failed: ' + errMsg);
  }

  // ── CSV Download ──────────────────────────────────────────
  function triggerDownload(jobId) {
    window.location.href = `/api/batch/${encodeURIComponent(jobId)}/download`;
  }

  // ── Sample Data ───────────────────────────────────────────
  document.getElementById('batchSampleBtn')?.addEventListener('click', () => {
    const sample = [
      '>P04637|TP53_HUMAN Cellular tumor antigen p53 [Homo sapiens]',
      'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPE',
      'AAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLFRNLNQTSFQNLSDLQPPPSSQS',
      '>P00533|EGFR_HUMAN Epidermal growth factor receptor [Homo sapiens]',
      'MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA',
      '>P68871|HBB_HUMAN Hemoglobin subunit beta [Homo sapiens]',
      'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
    ].join('\n');

    const pasteArea = document.getElementById('batchPasteInput');
    if (pasteArea) {
      pasteArea.value = sample;
      uploadedContent = sample;
    }
    window.showToast?.('Example loaded — 3 human proteins (TP53, EGFR, HBB)');
  });

  // ── Clear ─────────────────────────────────────────────────
  document.getElementById('batchClearBtn')?.addEventListener('click', () => {
    stopPolling();
    uploadedContent = null;
    currentJobId = null;
    const pasteArea = document.getElementById('batchPasteInput');
    if (pasteArea) pasteArea.value = '';
    if (fileInput) fileInput.value = '';
    resetJobUI();
  });

  // ── ViewManager Cleanup ───────────────────────────────────
  if (window.ViewManager) {
    window.ViewManager.registerUnmount('batch', () => {
      stopPolling();
    });
  }

})();
