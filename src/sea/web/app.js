const $=id=>document.getElementById(id); let currentJob=null,pollTimer=null,currentResult=null,busy=false,importing=false,exporting=false,editorReady=false,editorTarget='query';
function refreshInputControls(){['predict-submit','compare-submit'].forEach(id=>$(id).disabled=busy||importing||exporting);document.querySelectorAll('input[type=file],[data-draw]').forEach(input=>input.disabled=importing||exporting);['editor-add','editor-clear'].forEach(id=>$(id).disabled=!editorReady||importing||exporting);}
function setBusy(value){busy=value;refreshInputControls();}
function message(text=''){ $('message').textContent=text; }
function errorText(detail){if(typeof detail==='string')return detail;if(detail?.records)return [detail.message,...detail.records.map(r=>`レコード ${r.record}${r.name?` (${r.name})`:''}: ${r.message}`)].join('\n');return detail?.message||JSON.stringify(detail);}
async function api(path,options={}){const res=await fetch('/api/v1'+path,{...options,headers:{'Content-Type':'application/json','X-SEA-Request':'1',...options.headers}});if(!res.ok){const error=await res.json().catch(()=>({detail:res.statusText}));throw new Error(errorText(error.detail));}return res.json();}
const lines=id=>$(id).value.split(/\r?\n/).map(s=>s.trim()).filter(Boolean);
async function loadWorkspace(){try{const datasets=await api('/datasets');$('dataset').replaceChildren(...datasets.map(d=>{const o=document.createElement('option');o.value=d.id;o.textContent=`${d.source} ${d.release} · ${d.target_count}標的 · ${d.id}`;return o;}));$('dataset').onchange=()=>{const d=datasets.find(d=>d.id===$('dataset').value);$('dataset-info').textContent=d?`${d.molecule_count}分子 / RDKit ${d.rdkit_version} / ${d.model_id}`:'管理者がデータセットを構築してください。';$('attribution').textContent=d?.attribution||'SEA Open';};$('dataset').onchange();$('login-panel').hidden=true;$('workspace').hidden=false;}catch(e){$('login-panel').hidden=false;$('workspace').hidden=true;}}
$('login').onsubmit=async e=>{e.preventDefault();try{await api('/session',{method:'POST',body:JSON.stringify({api_key:$('key').value})});$('key').value='';message();await loadWorkspace();}catch(e){message(e.message);}};
$('logout').onclick=async()=>{currentJob=null;setBusy(false);await api('/session',{method:'DELETE'});clearTimeout(pollTimer);await loadWorkspace();};
document.querySelectorAll('[data-tab]').forEach(b=>b.onclick=()=>{document.querySelectorAll('.panel').forEach(p=>p.hidden=p.id!==b.dataset.tab);document.querySelectorAll('[data-tab]').forEach(t=>t.classList.toggle('active',t===b));if(b.dataset.tab==='jobs')loadJobs();});
async function submit(path,body){if(busy||importing||exporting)return;setBusy(true);message();try{const r=await api(path,{method:'POST',body:JSON.stringify({dataset_id:$('dataset').value,...body})});watch(r.job_id);}catch(e){setBusy(false);message(e.message);}}
$('predict-submit').onclick=()=>{const q=lines('query');submit('/predictions',{queries:$('mode').value==='ensemble'?[q]:q.map(s=>[s])});};
$('compare-submit').onclick=()=>submit('/comparisons',{left:lines('left'),right:lines('right')});
$('cancel').onclick=async()=>{if(currentJob)await api(`/jobs/${currentJob}/cancel`,{method:'POST'});};
async function watch(id){clearTimeout(pollTimer);currentJob=id;setBusy(true);$('progress-panel').hidden=false;$('result-panel').hidden=true;async function tick(){if(currentJob!==id)return;try{const j=await api(`/jobs/${id}`);if(currentJob!==id)return;const label=j.status==='queued'?'待機中（先行ジョブの完了を待っています）':j.status==='running'?`計算中 ${Math.round(j.progress*100)}%`:j.status;$('status').textContent=`${label} · ${id}`;$('progress').value=j.progress;$('cancel').hidden=!['queued','running'].includes(j.status);if(j.status==='completed'){currentResult=await api(`/jobs/${id}/result`);render(currentResult);setBusy(false);return;}if(['failed','cancelled'].includes(j.status)){setBusy(false);message(j.error||'ジョブを取り消しました');return;}pollTimer=setTimeout(tick,1000);}catch(e){setBusy(false);message(e.message);}}await tick();}
const num=x=>x===null?'—':Number(x).toExponential(3);
function render(data){$('result-panel').hidden=false;$('result-rows').replaceChildren();data.results.forEach((r,i)=>{const tr=document.createElement('tr');tr.dataset.index=i;const values=[r.query_index===null?'A ↔ B':r.query_index+1,r.target?.name||'集合比較',r.target?.family||'—',r.target?.organism?(r.target.organism+(r.target.tax_id?` (${r.target.tax_id})`:'')):'—',r.raw_score.toFixed(3),r.z_score?.toFixed(3)??'—',num(r.p_value),num(r.e_value),r.method+(r.p_resolution?` (分解能 ${r.p_resolution})`:'')];values.forEach(v=>{const td=document.createElement('td');td.textContent=v;if(tr.children.length===1&&r.target){const id=document.createElement('small');id.className='target-id';id.textContent=r.target.id;td.append(id);}if(tr.children.length===2&&r.target?.protein_classes){td.title=r.target.protein_classes.map(p=>p.join(' → ')).join('\n');}tr.append(td);});tr.onclick=()=>evidence(r);$('result-rows').append(tr);});$('metadata').textContent=JSON.stringify({dataset:data.dataset,threshold:data.threshold,normalization:data.normalization},null,2);$('json-export').href=`/api/v1/jobs/${currentJob}/result`;$('json-export').setAttribute('download',`${currentJob}.json`);$('csv-export').href=`/api/v1/jobs/${currentJob}/result?format=csv`;$('evidence').replaceChildren();if(data.results.length)evidence(data.results[0]);}
function evidence(r){$('evidence').replaceChildren();const h=document.createElement('h3');h.textContent='スコアに寄与した類似分子ペア';$('evidence').append(h);r.similar_ligands.slice(0,8).forEach(p=>{const row=document.createElement('div');row.className='pair';[p.query_smiles,p.reference_smiles].forEach(s=>{const img=document.createElement('img');img.src='/api/v1/structure?smiles='+encodeURIComponent(s);img.alt=s;img.title=s;row.append(img);});const span=document.createElement('span');span.textContent='Tc = '+p.tanimoto.toFixed(3);row.append(span);$('evidence').append(row);});}
async function loadJobs(){try{const jobs=await api('/jobs');$('job-list').replaceChildren();jobs.forEach(j=>{const row=document.createElement('div');row.className='job';const text=document.createElement('span');text.textContent=`${new Date(j.created*1000).toLocaleString()} · ${j.kind} · ${j.status}`;const button=document.createElement('button');button.textContent='開く';button.onclick=()=>watch(j.id);row.append(text,button);$('job-list').append(row);});}catch(e){message(e.message);}}
$('refresh').onclick=loadJobs;
function importStatus(id,text,error=false){const status=$(id+'-import-status');status.textContent=text;status.classList.toggle('error',error);}
async function importInto(id,format,content,append=false){
  if(importing)throw new Error('読み込み完了を待ってください');
  importing=true;refreshInputControls();$(id).disabled=true;importStatus(id,'構造を読み込んでいます…');message();
  try{
    const data=await api('/structures/import',{method:'POST',body:JSON.stringify({format,content})});
    const smiles=data.molecules.map(m=>m.smiles).join('\n');
    $(id).value=append&&$(id).value.trim()?$(id).value.trimEnd()+'\n'+smiles:smiles;
    importStatus(id,`${data.count}分子を${append?'追加':'読み込み'}しました。`);
    return data.count;
  }catch(e){importStatus(id,e.message,true);throw e;}
  finally{importing=false;$(id).disabled=false;refreshInputControls();}
}
const fileFormats={sdf:'sdf',sd:'sdf',mol:'mol',smi:'smiles',smiles:'smiles',txt:'smiles',csv:'csv'};
document.querySelectorAll('input[type=file]').forEach(input=>input.onchange=async()=>{
  try{
    const file=input.files[0];if(!file)return;
    if(file.size>10485760)throw new Error('構造ファイルは10 MB以下にしてください');
    const format=fileFormats[file.name.split('.').pop().toLowerCase()];
    if(!format)throw new Error('SDF・MOL・SMI・CSVファイルを選択してください');
    let content;try{content=new TextDecoder('utf-8',{fatal:true}).decode(await file.arrayBuffer());}catch{throw new Error('UTF-8のテキストファイルを選択してください');}
    await importInto(input.dataset.for,format,content);
  }catch(e){importStatus(input.dataset.for,e.message,true);}
  finally{input.value='';}
});
function editorStatus(text,error=false){$('editor-status').textContent=text;$('editor-status').classList.toggle('error',error);}
function editor(){return $('structure-editor').contentWindow.ketcher;}
document.querySelectorAll('[data-draw]').forEach(button=>button.onclick=()=>{
  editorTarget=button.dataset.draw;
  $('editor-title').textContent='構造を描く · '+({query:'標的予測',left:'集合A',right:'集合B'}[editorTarget]);
  $('structure-dialog').showModal();
  if(editorReady){editor().editor.clear();editorStatus('構造を描画し、「入力に追加」を押してください。');}
  else if(!$('structure-editor').hasAttribute('src'))$('structure-editor').src='/static/vendor/ketcher/index.html?disableMacromoleculesEditor=true';
});
window.addEventListener('message',event=>{
  if(event.origin!==location.origin||event.source!==$('structure-editor').contentWindow||event.data?.eventType!=='init')return;
  editorReady=true;refreshInputControls();editorStatus('構造を描画し、「入力に追加」を押してください。');
});
$('editor-close').onclick=()=>$('structure-dialog').close();
$('editor-clear').onclick=()=>{editor().editor.clear();editorStatus('構造を消去しました。');};
$('editor-add').onclick=async()=>{
  if(!editorReady||importing||exporting)return;
  exporting=true;refreshInputControls();
  try{
    if(editor().editor.struct().atoms.size===0)throw new Error('構造を描いてから追加してください');
    editorStatus('構造を読み込んでいます…');
    const target=editorTarget;
    await importInto(target,'mol',await editor().getMolfile(),true);
    $('structure-dialog').close();
  }catch(e){editorStatus(e.message,true);}
  finally{exporting=false;refreshInputControls();}
};
loadWorkspace();
