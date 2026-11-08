const $ = id => document.getElementById(id);
let currentJob = null, pollTimer = null, currentResult = null, busy = false;
let watchGeneration = 0, refinementJob = null, pendingResult = null, artifactReady = false, lastSnapshot = '';
let importing = false, exporting = false, editorReady = false, editorTarget = 'query';
let datasets = [], resultRows = [], selectedQuery = null, selectedIndex = null, showAllPairs = false;
const integer = value => Number(value).toLocaleString('en-US');
const lines = id => $(id).value.split(/\r?\n/).map(s => s.trim()).filter(Boolean);
const plural = (n, word) => `${integer(n)} ${word}${n === 1 ? '' : 's'}`;
const demoNames = {
  'DEMO-ALCOHOL': 'Alcohol set (synthetic demo)',
  'DEMO-ACID': 'Carboxylic acid set (synthetic demo)',
  'DEMO-AROMATIC': 'Aromatic set (synthetic demo)',
};
const targetName = target => demoNames[target?.id] || target?.name || 'Ligand set A vs B';
const families = target => target?.families?.length ? target.families : target?.family ? target.family.split(' / ') : [];
const queryKey = row => row.query_index == null ? 'comparison' : String(row.query_index);
const bestTc = row => row.similar_ligands?.length ? Math.max(...row.similar_ligands.map(p => p.tanimoto)) : null;
const methodName = method => ({gumbel: 'Calibrated Gumbel', monte_carlo: 'Monte Carlo', exact_zero: 'Exact zero-score evaluation', pending: 'Statistical evaluation pending'}[method] || method);
const statisticsStatus = row => ({pending: 'Pending', resolution_limited: 'Resolution-limited', uncertain: 'Uncertain', evaluated: 'Evaluated', exact: 'Exact', analytic: 'Analytic model'}[row.statistics?.status] || '');
const scrollTo = element => element.scrollIntoView({behavior: matchMedia('(prefers-reduced-motion: reduce)').matches ? 'instant' : 'smooth', block: 'start'});

function el(tag, className = '', text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}
function message(text = '') { $('message').textContent = text; }
function errorText(detail) {
  if (typeof detail === 'string') return detail;
  if (Array.isArray(detail)) return detail.map(d => `${d.loc?.filter(x => x !== 'body').join(' / ') || 'Input'}: ${d.msg || d.message}`).join('\n');
  if (detail?.records) return [detail.message, ...detail.records.map(r => `Record ${r.record}${r.name ? ` (${r.name})` : ''}: ${r.message}`)].join('\n');
  if (detail?.rows) {
    const prefix = `Query or set ${Number(detail.set_index) + 1}`;
    return typeof detail.rows === 'string' ? `${prefix}: ${detail.rows}` : detail.rows.map(r => `${prefix}, row ${r.row}: ${r.message}`).join('\n');
  }
  if (detail?.code === 'too_many_molecules') return `Too many molecules. The limit is ${integer(detail.limit)}.`;
  return detail?.message || JSON.stringify(detail);
}
async function api(path, options = {}) {
  const res = await fetch('/api/v1' + path, {...options, headers: {'Content-Type': 'application/json', 'X-SEA-Request': '1', ...options.headers}});
  if (!res.ok) {
    const body = await res.json().catch(() => ({detail: res.statusText}));
    const error = new Error(errorText(body.detail)); error.status = res.status; throw error;
  }
  return res.json();
}
function refreshInputControls() {
  ['predict-submit', 'compare-submit'].forEach(id => $(id).disabled = busy || importing || exporting || !$('dataset').value);
  document.querySelectorAll('input[type=file], [data-draw]').forEach(input => input.disabled = importing || exporting);
  ['editor-add', 'editor-clear'].forEach(id => $(id).disabled = !editorReady || importing || exporting);
  $('load-example').disabled = importing || exporting;
  $('dataset').disabled = busy;
  $('statistics-profile').disabled = busy;
}
function setBusy(value) { busy = value; refreshInputControls(); }
const inputPreviews = new Map();
const previewPageSize = 4;
for (const id of ['query', 'left', 'right']) {
  const panel = el('section', 'input-preview'); panel.id = `${id}-input-preview`; panel.hidden = true;
  panel.setAttribute('aria-label', `${id === 'query' ? 'Query' : id === 'left' ? 'Set A' : 'Set B'} input structures`);
  $(id + '-import-status').after(panel);
  inputPreviews.set(id, {panel, text: null, timer: null, page: 0});
}
function scheduleInputPreview(id) {
  const state = inputPreviews.get(id), text = $(id).value;
  if (text === state.text) return;
  state.text = text; state.page = 0; clearTimeout(state.timer);
  state.panel.replaceChildren(); state.panel.hidden = !lines(id).length;
  if (!state.panel.hidden) {
    state.panel.append(el('p', 'preview-note', 'Updating structure preview…'));
    state.timer = setTimeout(() => renderInputPreview(id), 350);
  }
}
function renderInputPreview(id) {
  const state = inputPreviews.get(id), structures = lines(id), start = state.page * previewPageSize;
  state.panel.replaceChildren(); state.panel.hidden = !structures.length;
  if (!structures.length) return;
  const header = el('div', 'preview-heading');
  header.append(el('strong', '', 'Input structures'), el('span', '', `${start + 1}–${Math.min(start + previewPageSize, structures.length)} of ${integer(structures.length)}`));
  const grid = el('div', 'preview-grid');
  for (const [offset, smiles] of structures.slice(start, start + previewPageSize).entries()) {
    const card = el('figure', 'input-molecule'), label = el('figcaption', '', `Molecule ${start + offset + 1}`);
    const image = el('img'); image.width = 240; image.height = 144; image.alt = `Input molecule ${start + offset + 1}: ${smiles}`;
    const status = el('span', 'preview-molecule-status', 'Rendering…'); status.setAttribute('role', 'status');
    image.onload = () => { status.hidden = true; card.dataset.state = 'ready'; };
    image.onerror = () => { image.hidden = true; status.textContent = 'Unable to render. Check the SMILES.'; card.dataset.state = 'error'; };
    const code = el('code', '', smiles); code.title = smiles;
    card.append(label, image, status, code); grid.append(card);
    image.src = '/api/v1/structure?standardize=false&smiles=' + encodeURIComponent(smiles);
  }
  const footer = el('div', 'preview-footer');
  footer.append(el('span', 'preview-note', 'Before search standardization. Salts and charges are shown as entered.'));
  if (structures.length > previewPageSize) {
    const nav = el('div', 'preview-navigation');
    for (const [label, delta] of [['Previous', -1], ['Next', 1]]) {
      const button = el('button', 'quiet', label); button.type = 'button';
      button.disabled = delta < 0 ? state.page === 0 : start + previewPageSize >= structures.length;
      button.setAttribute('aria-label', `${label} input structures`);
      button.onclick = () => { state.page += delta; renderInputPreview(id); };
      nav.append(button);
    }
    footer.append(nav);
  }
  state.panel.append(header, grid, footer);
}
function updateInputCounts() {
  for (const id of ['query', 'left', 'right']) {
    $(id + '-count').textContent = plural(lines(id).length, 'molecule');
    scheduleInputPreview(id);
  }
}
for (const id of ['query', 'left', 'right']) $(id).addEventListener('input', updateInputCounts);
function libraryName(dataset) {
  if (dataset.source === 'synthetic-demo') return 'Synthetic demo · chemical groups';
  const release = String(dataset.release || '').replace(/^chembl[_\s-]*/i, '');
  return `${dataset.source} ${release} · ${plural(dataset.target_count, 'target')}`;
}
async function loadWorkspace() {
  try {
    const previous = $('dataset').value;
    datasets = (await api('/datasets')).sort((a, b) =>
      Number(b.source === 'ChEMBL') - Number(a.source === 'ChEMBL') ||
      String(b.release).localeCompare(String(a.release), 'en', {numeric: true}) || a.id.localeCompare(b.id));
    $('dataset').replaceChildren(...datasets.map(d => {
      const option = el('option', '', libraryName(d)); option.value = d.id; return option;
    }));
    if (datasets.some(d => d.id === previous)) $('dataset').value = previous;
    $('dataset').onchange = () => {
      const d = datasets.find(d => d.id === $('dataset').value);
      $('dataset-info').textContent = d ? `${plural(d.molecule_count, 'reference molecule')} · RDKit ${d.rdkit_version}` : 'No reference libraries. Ask your administrator to build a dataset.';
      $('attribution').textContent = d?.attribution || 'SEA Open';
      refreshInputControls();
    };
    $('dataset').onchange();
    $('login-panel').hidden = true; $('workspace').hidden = false; $('logout').hidden = false;
  } catch (error) {
    $('login-panel').hidden = false; $('workspace').hidden = true; $('logout').hidden = true;
    if (error.status !== 401) message(error.message);
  }
}
$('login').onsubmit = async event => {
  event.preventDefault();
  const button = event.currentTarget.querySelector('button'); button.disabled = true;
  try {
    await api('/session', {method: 'POST', body: JSON.stringify({api_key: $('key').value})});
    $('key').value = ''; message(); await loadWorkspace();
  } catch (error) { message(error.message); }
  finally { button.disabled = false; }
};
$('logout').onclick = async () => {
  try {
    await api('/session', {method: 'DELETE'});
    currentJob = null; currentResult = null; watchGeneration++; refinementJob = null; pendingResult = null; clearTimeout(pollTimer); setBusy(false);
    $('progress-panel').hidden = true; $('result-panel').hidden = true; $('job-list').replaceChildren();
    message(); await loadWorkspace();
  } catch (error) { message(error.message); }
};
document.querySelectorAll('[data-tab]').forEach(button => button.onclick = () => {
  document.querySelectorAll('.panel').forEach(panel => panel.hidden = panel.id !== button.dataset.tab);
  document.querySelectorAll('[data-tab]').forEach(tab => {
    tab.classList.toggle('active', tab === button); tab.setAttribute('aria-pressed', String(tab === button));
  });
  if (button.dataset.tab === 'jobs') loadJobs();
});
$('load-example').onclick = () => {
  const demo = datasets.find(d => d.id === $('dataset').value)?.source === 'synthetic-demo';
  $('query').value = demo ? 'CCO\nCCCO' : 'COc1ccc2[nH]cc(c2c1)C(CC)CN';
  importStatus('query', 'Example loaded. Ready to search.'); updateInputCounts(); $('query').focus();
};
async function submit(path, body) {
  if (busy || importing || exporting) return;
  watchGeneration++; clearTimeout(pollTimer);
  setBusy(true); message();
  try {
    const response = await api(path, {method: 'POST', body: JSON.stringify({dataset_id: $('dataset').value, statistics_profile: $('statistics-profile').value, ...body})});
    await watch(response.job_id);
  } catch (error) { setBusy(false); message(error.message); }
}
$('predict-submit').onclick = () => {
  const query = lines('query');
  if (!query.length) { message('Add at least one molecule to search.'); $('query').focus(); return; }
  submit('/predictions', {queries: $('mode').value === 'ensemble' ? [query] : query.map(s => [s]), top_k: Number($('top-k').value)});
};
$('compare-submit').onclick = () => {
  if (!lines('left').length || !lines('right').length) { message('Add at least one molecule to each ligand set.'); return; }
  submit('/comparisons', {left: lines('left'), right: lines('right')});
};
$('cancel').onclick = async () => {
  if (!currentJob) return;
  try { await api(`/jobs/${refinementJob || currentJob}/cancel`, {method: 'POST'}); }
  catch (error) { message(error.message); }
};
function incomingResult(data, id, ready) {
  const signature = JSON.stringify([id, ready, data.analysis?.phase, data.analysis?.evaluated_comparisons, data.results.map(row => [row.query_index, row.target?.id, row.p_value, row.statistics?.phase, row.statistics?.trials])]);
  if (signature === lastSnapshot) return;
  lastSnapshot = signature;
  if (!currentResult) {
    currentJob = id; artifactReady = ready; render(data); setBusy(false);
  } else {
    pendingResult = {data, id, ready}; $('apply-statistics').hidden = false;
    $('apply-statistics').textContent = data.analysis?.complete ? 'Apply final results' : data.analysis?.phase === 'screening' && ready ? 'Apply screening results' : 'Apply latest statistics';
  }
}
$('apply-statistics').onclick = () => {
  if (!pendingResult) return;
  const update = pendingResult, query = $('query-select').value;
  const selected = resultRows[selectedIndex]?.result;
  const filters = [$('result-search').value, $('family-filter').value, $('result-sort').value];
  const position = {top: window.scrollY, left: window.scrollX, behavior: 'instant'};
  pendingResult = null; currentJob = update.id; artifactReady = update.ready;
  render(update.data);
  if ([...$('query-select').options].some(option => option.value === query)) {
    $('query-select').value = query; changeQuery();
  }
  $('result-search').value = filters[0]; $('result-sort').value = filters[2];
  if ([...$('family-filter').options].some(option => option.value === filters[1])) $('family-filter').value = filters[1];
  selectedIndex = resultRows.find(item => queryKey(item.result) === query && item.result.target?.id === selected?.target?.id)?.index ?? null;
  updateResults(); window.scrollTo(position);
};
$('resume-refinement').onclick = async () => {
  if (!refinementJob) return;
  try { await api(`/jobs/${refinementJob}/resume`, {method: 'POST'}); await watch(refinementJob, true); }
  catch (error) { message(error.message); }
};
async function watch(id, preserve = false) {
  const generation = ++watchGeneration;
  clearTimeout(pollTimer); refinementJob = null; pendingResult = null; lastSnapshot = '';
  $('apply-statistics').hidden = true; $('resume-refinement').hidden = true;
  if (!preserve) { currentJob = id; currentResult = null; artifactReady = false; $('result-panel').hidden = true; }
  setBusy(!preserve); message();
  $('progress-panel').hidden = false;
  $('status').textContent = 'Loading analysis…'; $('job-id').textContent = `Job ${id}`;
  $('progress').value = 0; $('cancel').hidden = false; $('cancel').textContent = 'Cancel job';
  let activeId = id;
  async function tick() {
    if (generation !== watchGeneration) return;
    try {
      const job = await api(`/jobs/${activeId}`);
      if (generation !== watchGeneration) return;
      const refining = job.kind === 'refinement';
      if (refining) { refinementJob = activeId; $('cancel').textContent = 'Pause refinement'; }
      if (refining && job.status !== 'completed' && !currentResult && job.parent_id) {
        const source = await api(`/jobs/${job.parent_id}/result`);
        if (generation !== watchGeneration) return;
        incomingResult(source, job.parent_id, true); setBusy(false);
      }
      const percent = Math.round(job.progress * 100);
      $('status').textContent = job.status === 'queued' ? (refining ? 'Refinement queued · you can continue searching' : 'Queued · waiting for an available worker') : job.status === 'running' ? `${refining ? 'Refining fixed-budget statistics' : 'Preparing evidence and screening statistics'} · ${percent}%` : `${job.status[0].toUpperCase()}${job.status.slice(1)}`;
      $('job-id').textContent = `Job ${activeId}`; $('progress').value = job.progress;
      $('cancel').hidden = !['queued', 'running'].includes(job.status);
      if (job.status === 'completed') {
        const data = await api(`/jobs/${activeId}/result`);
        if (generation !== watchGeneration) return;
        incomingResult(data, activeId, true); setBusy(false);
        if (data.analysis && !data.analysis.complete && job.refinement_job_id) {
          activeId = job.refinement_job_id; refinementJob = activeId;
          pollTimer = setTimeout(tick, 100); return;
        }
        $('progress-panel').hidden = true; return;
      }
      if (job.profile !== 'legacy' || refining) {
        try {
          const data = await api(`/jobs/${activeId}/preview`);
          if (generation !== watchGeneration) return;
          incomingResult(data, activeId, false); setBusy(false);
        } catch (error) { if (error.status !== 409) throw error; }
      }
      if (['failed', 'cancelled'].includes(job.status)) {
        setBusy(false); $('resume-refinement').hidden = !refining;
        $('status').textContent = refining ? 'Refinement incomplete · screening results remain available' : 'Analysis incomplete';
        message(job.error || (refining ? 'Refinement is paused. Resume it to complete the declared statistical budget.' : 'The job was cancelled.')); return;
      }
      pollTimer = setTimeout(tick, 500);
    } catch (error) { if (generation === watchGeneration) { setBusy(false); $('cancel').hidden = true; message(error.message); } }
  }
  await tick();
}

// Use stored natural logarithms when a very small p/E value underflows to zero.
function scientific(value, logValue) {
  if (value == null || !Number.isFinite(value)) return el('span', '', '—');
  if (value === 0 && !Number.isFinite(logValue)) return el('span', '', '0');
  if (value === 0 && logValue === 0) return el('span', '', '0');
  const log10 = value > 0 ? Math.log10(value) : logValue / Math.LN10;
  let exponent = Math.floor(log10), mantissa = 10 ** (log10 - exponent);
  if (Number(mantissa.toPrecision(3)) >= 10) { mantissa = 1; exponent += 1; }
  const span = el('span', 'scientific');
  if (exponent >= -2 && exponent < 3) span.textContent = Number((mantissa * 10 ** exponent).toPrecision(4)).toLocaleString('en-US', {maximumSignificantDigits: 4});
  else {
    span.append(document.createTextNode(`${mantissa.toFixed(2)} × 10`), el('sup', '', exponent));
    span.setAttribute('aria-label', `${mantissa.toFixed(2)} times ten to the power ${exponent}`);
  }
  span.title = `Natural log: ${Number.isFinite(logValue) ? logValue.toPrecision(8) : 'not available'}`;
  return span;
}
function targetMetadata(target, showTaxonomy = false) {
  const fragment = document.createDocumentFragment();
  if (target?.id) fragment.append(el('small', 'target-id', target.id));
  const meta = el('div', 'target-meta');
  families(target).forEach(family => { const badge = el('span', 'family-badge', family); badge.dataset.family = family; meta.append(badge); });
  if (target?.protein_classes?.length) meta.title = target.protein_classes.map(path => path.join(' → ')).join('\n');
  if (target?.organism) {
    const species = el('span', 'species-label', target.organism + (showTaxonomy && target.tax_id ? ` · Taxonomy ${target.tax_id}` : ''));
    if (target.tax_id) species.title = `NCBI taxonomy ID: ${target.tax_id}`;
    meta.append(species);
  }
  fragment.append(meta); return fragment;
}
function render(data) {
  currentResult = data; selectedIndex = null; showAllPairs = false;
  const counts = new Map();
  resultRows = data.results.map((result, index) => {
    const key = queryKey(result), rank = (counts.get(key) || 0) + 1; counts.set(key, rank);
    return {result, index, rank};
  });
  let keys = Array.isArray(data.normalization) ? data.normalization.map((_, i) => String(i)) : ['comparison'];
  if (!keys.length) keys = [...counts.keys()];
  $('query-select').replaceChildren(...keys.map(key => {
    const rows = resultRows.filter(row => queryKey(row.result) === key);
    const size = rows[0]?.result.left_size;
    const label = key === 'comparison' ? 'Ligand set A vs B' : `Query ${Number(key) + 1}${size > 1 ? ` · ${plural(size, 'molecule')}` : ''}`;
    const option = el('option', '', label); option.value = key; return option;
  }));
  $('result-search').value = ''; $('result-sort').value = 'sea'; $('family-filter').value = '';
  $('results-title').textContent = keys[0] === 'comparison' ? 'Set comparison results' : 'Target search results';
  const elapsed = data.execution?.elapsed_seconds;
  const fingerprint = data.dataset.fingerprint;
  $('result-context').textContent = `${libraryName(data.dataset)}${fingerprint ? ` · ${fingerprint.kind} · ${integer(fingerprint.bits)} bits` : ''}${Number.isFinite(elapsed) ? ` · ${elapsed.toFixed(2)} s of computation` : ''}`;
  $('metadata').textContent = JSON.stringify({dataset: data.dataset, threshold: data.threshold, normalization: data.normalization, execution: data.execution, analysis: data.analysis, job_id: currentJob}, null, 2);
  const analysis = data.analysis;
  $('analysis-banner').hidden = !analysis;
  $('apply-statistics').hidden = !pendingResult;
  if (analysis) {
    const labels = {evidence: 'Chemical evidence · statistics pending', screening: 'Screening · preliminary statistics', refining: 'Refinement · partial statistics', final: 'Final · fixed budget complete'};
    $('analysis-stage').textContent = labels[analysis.phase] || analysis.phase;
    $('analysis-note').textContent = `${analysis.profile.name} profile · ${integer(analysis.profile.planned_trials)} planned trials per Monte Carlo comparison. ${analysis.complete ? 'Final refers to computational completion; limited resolution or uncertainty may remain.' : 'You can continue searching. Updates are applied only when you choose; displayed rows stay in place.'}`;
  }
  for (const [id, suffix, extension] of [['json-export', '', 'json'], ['csv-export', '?format=csv', 'csv']]) {
    $(id).hidden = !artifactReady;
    $(id).href = `/api/v1/jobs/${currentJob}/result${suffix}`; $(id).setAttribute('download', `${currentJob}.${extension}`);
  }
  $('result-panel').hidden = false; changeQuery();
}
function changeQuery() {
  selectedQuery = $('query-select').value; selectedIndex = null;
  const rows = resultRows.filter(row => queryKey(row.result) === selectedQuery);
  const isComparison = selectedQuery === 'comparison';
  const normalization = isComparison ? currentResult.normalization?.left : currentResult.normalization?.[Number(selectedQuery)];
  const structures = [...new Set((normalization || []).filter(row => !row.duplicate).map(row => row.smiles))];
  $('query-description').textContent = `${plural(structures.length, 'unique molecule')}${isComparison ? ' in set A' : structures.length > 1 ? ' · ensemble query' : ''}`;
  $('query-smiles').textContent = structures.length === 1 ? structures[0] : structures.join(' · ');
  $('query-smiles').title = structures.join('\n');
  $('query-preview').hidden = !structures.length;
  $('query-preview').alt = structures.length > 1 ? 'First molecule in the query ensemble' : 'Query molecule';
  if (structures.length) $('query-preview').src = '/api/v1/structure?smiles=' + encodeURIComponent(structures[0]);
  else $('query-preview').removeAttribute('src');
  $('summary-count').textContent = integer(rows.length);
  $('summary-count-label').textContent = isComparison ? 'set comparison' : `returned target${rows.length === 1 ? '' : 's'}`;
  const similarities = rows.map(row => bestTc(row.result)).filter(value => value !== null);
  $('summary-tc').textContent = similarities.length ? Math.max(...similarities).toFixed(4) : '—';
  const options = [...new Set(rows.flatMap(row => families(row.result.target)))].sort((a, b) => a.localeCompare(b, 'en'));
  const previousFamily = $('family-filter').value;
  $('family-filter').replaceChildren(el('option', '', 'All families'), ...options.map(family => { const option = el('option', '', family); option.value = family; return option; }));
  $('family-filter').firstChild.value = '';
  if (options.includes(previousFamily)) $('family-filter').value = previousFamily;
  $('family-filter').parentElement.hidden = isComparison;
  $('results-layout').querySelector('thead th:nth-child(2)').textContent = isComparison ? 'Ligand sets' : 'Target / family / species';
  $('result-search').placeholder = isComparison ? 'Filter set comparison…' : 'Filter by target, ID or species…';
  $('list-title').textContent = isComparison ? 'Ligand set comparison' : currentResult.analysis?.phase === 'evidence' ? 'Candidates by raw score · statistics pending' : 'Ranked targets';
  updateResults();
}
function updateResults() {
  const text = $('result-search').value.trim().toLocaleLowerCase('en');
  const family = $('family-filter').value;
  const all = resultRows.filter(row => queryKey(row.result) === selectedQuery);
  const filtered = all.filter(({result: row}) => {
    const target = row.target;
    const searchable = [targetName(target), target?.id, target?.organism, target?.tax_id, ...families(target)].filter(Boolean).join(' ').toLocaleLowerCase('en');
    return (!text || searchable.includes(text)) && (!family || families(target).includes(family));
  });
  if ($('result-sort').value === 'similarity') filtered.sort((a, b) => (bestTc(b.result) ?? -1) - (bestTc(a.result) ?? -1) || a.rank - b.rank);
  if ($('result-sort').value === 'name') filtered.sort((a, b) => targetName(a.result.target).localeCompare(targetName(b.result.target), 'en') || a.rank - b.rank);
  if (!filtered.some(row => row.index === selectedIndex)) selectedIndex = filtered[0]?.index ?? null;
  $('visible-count').textContent = `${integer(filtered.length)} of ${integer(all.length)} shown`;
  $('result-rows').replaceChildren();
  for (const item of filtered) {
    const {result: row, index, rank} = item;
    const tr = el('tr'); tr.dataset.index = index; tr.dataset.selected = String(index === selectedIndex);
    tr.append(el('td', '', String(rank).padStart(2, '0')));
    const target = el('td'), button = el('button', 'target-button', targetName(row.target));
    button.type = 'button'; button.setAttribute('aria-pressed', String(index === selectedIndex));
    button.setAttribute('aria-controls', 'target-detail');
    target.append(button, targetMetadata(row.target)); tr.append(target);
    const e = el('td', 'numeric-cell e-cell'); e.dataset.label = 'E-value'; e.append(scientific(row.e_value, row.log_e_value));
    if (row.statistics) e.append(el('small', 'statistical-state', `${row.statistics.phase === 'screening' ? 'Preliminary · ' : ''}${statisticsStatus(row)}`));
    tr.append(e);
    const tc = el('td', 'numeric-cell'), maximum = bestTc(row); tc.dataset.label = 'Best Tanimoto';
    tc.append(el('span', 'tc-value', maximum === null ? '—' : maximum.toFixed(4)));
    if (maximum !== null) { const track = el('span', 'tc-track'), fill = el('span', 'tc-fill'); fill.style.width = `${Math.max(0, Math.min(1, maximum))*100}%`; track.append(fill); tc.append(track); }
    tr.append(tc);
    tr.onclick = () => {
      selectedIndex = index; showAllPairs = false;
      document.querySelectorAll('#result-rows tr[data-index]').forEach(r => {
        r.dataset.selected = String(Number(r.dataset.index) === index);
        r.querySelector('button').setAttribute('aria-pressed', r.dataset.selected);
      });
      evidence(row);
      if (matchMedia('(max-width: 1000px)').matches) scrollTo($('target-detail'));
    };
    $('result-rows').append(tr);
  }
  $('target-detail').hidden = selectedIndex === null;
  $('results-layout').classList.toggle('no-selection', selectedIndex === null);
  if (selectedIndex === null) {
    const row = el('tr'), cell = el('td', 'empty-state', all.length ? 'No returned results match these filters. Try another name or choose All families.' : 'No results were returned for this query.');
    cell.colSpan = 4; row.append(cell); $('result-rows').append(row); $('evidence').replaceChildren();
  } else { showAllPairs = false; evidence(resultRows[selectedIndex].result); }
}
function stat(label, value, description) {
  const node = el('div', 'stat'); node.append(el('label', '', label));
  const strong = el('strong'); strong.append(value instanceof Node ? value : document.createTextNode(value)); node.append(strong);
  if (description) node.title = description; return node;
}
function evidence(row) {
  $('target-detail').hidden = false;
  $('selected-title').textContent = targetName(row.target);
  $('selected-identity').replaceChildren(targetMetadata(row.target, true));
  $('selected-stats').replaceChildren(
    stat('SEA E-value', scientific(row.e_value, row.log_e_value), 'Expected random hits across all searched targets.'),
    stat('p-value', scientific(row.p_value, row.log_p_value)),
    stat('Raw score', Number(row.raw_score).toFixed(3), 'Sum of Tanimoto values at or above the calibrated threshold.'),
    stat('Z-score', row.z_score == null ? '—' : row.z_score.toFixed(3)),
    stat('Contributing pairs', integer(row.contributing_pairs)),
    stat(row.target ? 'Reference ligands' : 'Set sizes (A / B)', row.target ? integer(row.right_size) : `${integer(row.left_size)} / ${integer(row.right_size)}`),
  );
  let method = methodName(row.method);
  if (row.p_resolution) method += ` · p-value resolution ${row.p_resolution}. Equal E-values are ordered by raw score.`;
  else if (row.method === 'exact_zero') method += ' · No ligand pairs reached the similarity threshold.';
  $('method-note').textContent = method;
  let simulation = $('simulation-note');
  if (!simulation) { simulation = el('p', 'simulation-note'); simulation.id = 'simulation-note'; $('method-note').after(simulation); }
  simulation.hidden = !row.statistics;
  if (row.statistics) {
    const stats = row.statistics;
    const interval = stats.e_tail_interval;
    simulation.textContent = `${statisticsStatus(row)} · ${stats.phase}. ${integer(stats.trials)} completed Monte Carlo trials${stats.exceedances == null ? '' : ` · ${integer(stats.exceedances)} null scores at least as large`}. Testing scope: ${integer(row.comparisons)} target comparisons for this query.`;
    if (interval) simulation.textContent += ` Pointwise 95% interval for the limiting null-tail E-value: ${interval.map(value => Number(value.toPrecision(4))).join('–')}. This measures simulation uncertainty, not binding probability.`;
    if (stats.status === 'resolution_limited') simulation.textContent += ' Zero exceedances: the reported value is the finite-sampling floor, not zero or an upper bound.';
    if (stats.status === 'uncertain') simulation.textContent += ' The interval crosses the E = 0.05 reference threshold.';
  }
  const pairs = row.similar_ligands || [], displayed = showAllPairs ? pairs : pairs.slice(0, 3);
  $('evidence-count').textContent = pairs.length ? `${displayed.length} of ${pairs.length} saved pairs · ${plural(row.contributing_pairs, 'contributing pair')} in total` : 'No pairs reached the calibrated similarity threshold.';
  $('evidence').replaceChildren(...displayed.map((pair, i) => {
    const card = el('article', 'pair'), header = el('div', 'pair-heading');
    header.append(el('span', '', `Match ${String(i + 1).padStart(2, '0')}`), el('strong', '', `Tanimoto ${pair.tanimoto.toFixed(4)}`));
    const molecules = el('div', 'pair-molecules');
    for (const [label, smiles] of [[row.target ? 'Query molecule' : 'Set A molecule', pair.query_smiles], [row.target ? 'Reference ligand' : 'Set B molecule', pair.reference_smiles]]) {
      const molecule = el('div', 'molecule-card'), image = el('img');
      image.src = '/api/v1/structure?smiles=' + encodeURIComponent(smiles); image.alt = `${label}: ${smiles}`;
      image.loading = 'lazy'; image.width = 180; image.height = 110;
      const details = el('details', 'smiles-disclosure'); details.append(el('summary', '', 'SMILES'), el('code', '', smiles));
      const copy = el('button', 'text-button', 'Copy SMILES'); copy.type = 'button';
      copy.onclick = async () => { try { await navigator.clipboard.writeText(smiles); copy.textContent = 'Copied'; } catch { copy.textContent = 'Select the SMILES to copy'; } };
      details.append(copy); molecule.append(el('span', 'molecule-label', label), image, details); molecules.append(molecule);
    }
    card.append(header, molecules); return card;
  }));
  $('show-pairs').hidden = pairs.length <= 3;
  $('show-pairs').textContent = showAllPairs ? 'Show fewer pairs' : `Show all ${pairs.length} saved pairs`;
  $('show-pairs').onclick = () => { showAllPairs = !showAllPairs; evidence(row); };
}
$('query-select').onchange = changeQuery;
$('result-search').oninput = updateResults;
$('family-filter').onchange = updateResults;
$('result-sort').onchange = updateResults;
$('back-results').onclick = () => scrollTo($('results-layout'));
$('result-panel').querySelector('a[href="#calculation-details"]').onclick = () => { $('calculation-details').open = true; };

async function loadJobs() {
  try {
    const jobs = await api('/jobs'); $('job-list').replaceChildren();
    if (!jobs.length) $('job-list').append(el('p', 'empty-state', 'No analyses yet. Start a target search or compare two ligand sets.'));
    for (const job of jobs) {
      const row = el('div', 'job'), info = el('div');
      info.append(el('strong', '', `${job.kind === 'prediction' ? 'Target search' : job.kind === 'refinement' ? 'Statistical refinement' : 'Set comparison'} · ${job.profile}`), el('small', '', job.id));
      if (job.parent_id) info.append(el('small', '', `Screening source: ${job.parent_id}`));
      const status = el('span', 'status-badge', job.status[0].toUpperCase() + job.status.slice(1)); status.dataset.status = job.status;
      const date = el('time', '', new Date(job.created * 1000).toLocaleString('en-US', {dateStyle: 'medium', timeStyle: 'short'})); date.dateTime = new Date(job.created * 1000).toISOString();
      const button = el('button', 'quiet', job.status === 'completed' ? 'View results →' : 'Open job →'); button.onclick = () => watch(job.id);
      row.append(info, date, status, button); $('job-list').append(row);
    }
  } catch (error) { message(error.message); }
}
$('refresh').onclick = loadJobs;
function importStatus(id, text, error = false) { $(id + '-import-status').textContent = text; $(id + '-import-status').classList.toggle('error', error); }
async function importInto(id, format, content, append = false) {
  if (importing) throw new Error('Wait for the current import to finish.');
  importing = true; refreshInputControls(); $(id).disabled = true; importStatus(id, 'Reading structures…'); message();
  try {
    const data = await api('/structures/import', {method: 'POST', body: JSON.stringify({format, content})});
    const smiles = data.molecules.map(m => m.smiles).join('\n');
    $(id).value = append && $(id).value.trim() ? $(id).value.trimEnd() + '\n' + smiles : smiles;
    importStatus(id, `${plural(data.count, 'molecule')} ${append ? 'added' : 'imported'}.`); updateInputCounts(); return data.count;
  } catch (error) { importStatus(id, error.message, true); throw error; }
  finally { importing = false; $(id).disabled = false; refreshInputControls(); }
}
const fileFormats = {sdf: 'sdf', sd: 'sdf', mol: 'mol', smi: 'smiles', smiles: 'smiles', txt: 'smiles', csv: 'csv'};
document.querySelectorAll('input[type=file]').forEach(input => input.onchange = async () => {
  try {
    const file = input.files[0]; if (!file) return;
    if (file.size > 10485760) throw new Error('Choose a structure file no larger than 10 MB.');
    const format = fileFormats[file.name.split('.').pop().toLowerCase()];
    if (!format) throw new Error('Choose an SDF, MOL, SMI or CSV file.');
    let content;
    try { content = new TextDecoder('utf-8', {fatal: true}).decode(await file.arrayBuffer()); }
    catch { throw new Error('Choose a UTF-8 text file.'); }
    await importInto(input.dataset.for, format, content);
  } catch (error) { importStatus(input.dataset.for, error.message, true); }
  finally { input.value = ''; }
});
function editorStatus(text, error = false) { $('editor-status').textContent = text; $('editor-status').classList.toggle('error', error); }
function editor() { return $('structure-editor').contentWindow.ketcher; }
const editorPrompt = 'Draw a molecule, then choose Add to input.';
document.querySelectorAll('[data-draw]').forEach(button => button.onclick = () => {
  editorTarget = button.dataset.draw;
  $('editor-title').textContent = 'Draw a structure · ' + ({query: 'Target search', left: 'Ligand set A', right: 'Ligand set B'}[editorTarget]);
  $('structure-dialog').showModal();
  if (editorReady) { editor().editor.clear(); editorStatus(editorPrompt); }
  else if (!$('structure-editor').hasAttribute('src')) $('structure-editor').src = '/static/vendor/ketcher/index.html?disableMacromoleculesEditor=true';
});
window.addEventListener('message', event => {
  if (event.origin !== location.origin || event.source !== $('structure-editor').contentWindow || event.data?.eventType !== 'init') return;
  editorReady = true; refreshInputControls(); editorStatus(editorPrompt);
});
$('editor-close').onclick = () => $('structure-dialog').close();
$('editor-clear').onclick = () => { editor().editor.clear(); editorStatus('Drawing cleared.'); };
$('editor-add').onclick = async () => {
  if (!editorReady || importing || exporting) return;
  exporting = true; refreshInputControls();
  try {
    if (editor().editor.struct().atoms.size === 0) throw new Error('Draw a molecule before adding it.');
    editorStatus('Reading the structure…');
    const target = editorTarget; await importInto(target, 'mol', await editor().getMolfile(), true); $('structure-dialog').close();
  } catch (error) { editorStatus(error.message, true); }
  finally { exporting = false; refreshInputControls(); }
};
updateInputCounts(); refreshInputControls(); loadWorkspace();
