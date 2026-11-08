// Optional development smoke check: Node >=22 and Chrome. Services started by smoke.py.
import fs from 'node:fs';
const base='http://127.0.0.1:19021';
const targets=await(await fetch('http://127.0.0.1:19022/json')).json();
const target=targets.find(t=>t.type==='page');
const ws=new WebSocket(target.webSocketDebuggerUrl);await new Promise((r,j)=>{ws.onopen=r;ws.onerror=j;});
let sequence=0;const callbacks=new Map(),externalRequests=[],runtimeErrors=[];
ws.onmessage=e=>{const data=JSON.parse(e.data);if(data.method==='Runtime.exceptionThrown')runtimeErrors.push(data.params.exceptionDetails);if(data.method==='Network.requestWillBeSent'&&/^https?:/.test(data.params.request.url)&&!data.params.request.url.startsWith(base+'/'))externalRequests.push(data.params.request.url);if(data.id){const cb=callbacks.get(data.id);callbacks.delete(data.id);data.error?cb.reject(Error(JSON.stringify(data.error))):cb.resolve(data.result);}};
function call(method,params={}){return new Promise((resolve,reject)=>{const id=++sequence;callbacks.set(id,{resolve,reject});ws.send(JSON.stringify({id,method,params}));});}
async function evaluate(expression){const r=await call('Runtime.evaluate',{expression,awaitPromise:true,returnByValue:true});if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));return r.result.value;}
const delay=ms=>new Promise(r=>setTimeout(r,ms));
await call('Network.enable');await call('Network.setCacheDisabled',{cacheDisabled:true});
await call('Runtime.enable');
await call('Page.enable');await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
await call('Page.navigate',{url:base});await delay(1000);
if(!await evaluate("document.documentElement.lang==='en'&&!/[ぁ-んァ-ン一-龯]/.test(document.body.innerText)"))throw Error('Login is not English');
const loginScreenshot=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync('var/smoke-login.png',Buffer.from(loginScreenshot.data,'base64'));
await evaluate("document.getElementById('key').value='invalid';document.getElementById('login').requestSubmit()");
await delay(250);
if(!await evaluate("document.getElementById('message').textContent.includes('Invalid API key')"))throw Error('Authentication error is not English');
await evaluate(`document.getElementById('key').value=${JSON.stringify(process.env.SEA_SMOKE_KEY)};document.getElementById('login').requestSubmit();`);
for(let i=0;i<30;i++){if(await evaluate("!document.getElementById('workspace').hidden"))break;await delay(200);}
if(!await evaluate("!document.getElementById('workspace').hidden"))throw Error('Login failed');
async function waitFor(expression){for(let i=0;i<200;i++){if(await evaluate(expression))return;await delay(100);}throw Error('Timed out: '+expression);}
async function upload(id,name,content){
  await evaluate(`(()=>{document.getElementById('${id}-import-status').textContent='';const input=document.querySelector('input[data-for="${id}"]');const files=new DataTransfer();files.items.add(new File([${JSON.stringify(content)}],${JSON.stringify(name)},{type:'text/plain'}));input.files=files.files;input.dispatchEvent(new Event('change',{bubbles:true}));})()`);
  await waitFor(`!importing&&document.getElementById('${id}-import-status').textContent!==''`);
}
await upload('left','compounds.smi','CCO ethanol\nCCCO propanol\n');
await upload('right','compounds.csv','smiles,name\r\nCCO,"quoted, name"\r\nCCCO,propanol\r\n');
await upload('query','ethanol.mol',process.env.SEA_SMOKE_MOL);
if(await evaluate("document.getElementById('query').value")!=='CCO')throw Error('MOL import failed');
await upload('query','compounds.sdf',process.env.SEA_SMOKE_SDF);
if(await evaluate("document.getElementById('query').value")!=='CCO\nCCCO')throw Error('SDF import failed');
await upload('query','invalid.sdf',process.env.SEA_SMOKE_MOL+'$$$$\ninvalid\n$$$$\n');
if(!await evaluate("document.getElementById('query-import-status').classList.contains('error')&&document.getElementById('query-import-status').textContent.includes('Record 2')&&document.getElementById('query').value==='CCO\\nCCCO'"))throw Error('Invalid SDF changed the input or lost its record number');
for(const id of ['left','right'])if(await evaluate(`document.getElementById('${id}').value`)!=='CCO\nCCCO')throw Error('Set import failed: '+id);
await evaluate("document.querySelector('[data-draw=query]').click()");
await waitFor('editorReady');
await evaluate("document.getElementById('editor-add').click()");
if(!await evaluate("document.getElementById('editor-status').classList.contains('error')"))throw Error('Empty drawing was accepted');
// Draw a carbon-carbon bond using the real toolbar and browser pointer events.
await evaluate("document.getElementById('structure-editor').contentDocument.querySelector('[data-testid=bonds]').click()");
const point=await evaluate("(()=>{const r=document.getElementById('structure-editor').getBoundingClientRect();return {x:r.x+r.width/2,y:r.y+r.height/2};})()");
await call('Input.dispatchMouseEvent',{type:'mousePressed',...point,button:'left',clickCount:1});
await call('Input.dispatchMouseEvent',{type:'mouseReleased',...point,button:'left',clickCount:1});
await waitFor('editor().editor.struct().atoms.size===2');
if(await evaluate('editor().getSmiles()')!=='CC')throw Error('Drawing tool failed');
await evaluate("editorStatus(editorPrompt)");
const editorScreenshot=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync('var/smoke-editor-desktop.png',Buffer.from(editorScreenshot.data,'base64'));
await evaluate("document.getElementById('editor-add').click();document.getElementById('editor-add').click()");
await waitFor("!document.getElementById('structure-dialog').open");
if(await evaluate("document.getElementById('query').value")!=='CCO\nCCCO\nCC')throw Error('Drawing was not appended exactly once');
for(const id of ['left','right']){
  await evaluate(`document.querySelector('[data-draw=${id}]').click()`);
  await delay(200);
  await evaluate(`editor().setMolecule(${JSON.stringify(process.env.SEA_SMOKE_MOL)})`);
  await evaluate("document.getElementById('editor-add').click()");
  await waitFor("!document.getElementById('structure-dialog').open");
  if(await evaluate(`document.getElementById('${id}').value`)!=='CCO\nCCCO\nCCO')throw Error('Drawing went to the wrong input: '+id);
}
// Restore the two-molecule file input before exercising the actual SEA job.
await upload('query','compounds.sdf',process.env.SEA_SMOKE_SDF);
await evaluate("document.getElementById('predict-submit').click();document.getElementById('predict-submit').click()");
for(let i=0;i<60;i++){if(await evaluate("!document.getElementById('result-panel').hidden"))break;await delay(500);}

await waitFor('!busy');
const rows=await evaluate("document.querySelectorAll('#result-rows tr[data-index]').length");
const jobCount=await evaluate("fetch('/api/v1/jobs').then(r=>r.json()).then(jobs=>jobs.length)");
if(jobCount!==1)throw Error('Duplicate submission created extra jobs: '+jobCount);
if(rows!==3)throw Error('Expected three demo targets for the selected query: '+rows);
if(!await evaluate("document.getElementById('query-select').options.length===2"))throw Error('Queries were not grouped');
await evaluate("document.getElementById('query-select').value='1';document.getElementById('query-select').dispatchEvent(new Event('change'))");
if(!await evaluate("Array.from(document.querySelectorAll('#result-rows tr[data-index]')).every(r=>currentResult.results[Number(r.dataset.index)].query_index===1)&&document.getElementById('query-smiles').textContent==='CCCO'"))throw Error('Switching query mixed its results or structures');
await evaluate("document.getElementById('query-select').value='0';document.getElementById('query-select').dispatchEvent(new Event('change'))");
// Transient rendering fixtures exercise metadata, ties, saved-pair limits and
// underflow. They are never written to a dataset or a stored analysis.
await evaluate(`(()=>{
  window.smokeOriginal=currentResult;
  const families=['GPCR','Kinase','Ion channel'];
  window.smokeFixture={...currentResult,results:families.map((family,i)=>({
    ...currentResult.results[0],query_index:0,
    target:{id:'RENDER-'+i,name:['Serotonin 2A receptor','Epidermal growth factor receptor','Potassium channel subfamily H member 2'][i],family,families:[family],protein_classes:[[family]],organism:'Homo sapiens',tax_id:9606},
    e_value:i===0?0:i*.15,log_e_value:i===0?-1000:Math.log(i*.15),
    p_value:i===0?0:i*.05,log_p_value:i===0?-1000-Math.log(3):Math.log(i*.05),
    contributing_pairs:8,similar_ligands:Array.from({length:8},(_,j)=>({...currentResult.results[0].similar_ligands[0],tanimoto:.6+i*.1-j*.01}))
  }))};
  render(window.smokeFixture);
})()`);
if(!await evaluate("document.querySelectorAll('#result-rows tr[data-index]').length===3&&Array.from(document.querySelectorAll('#result-rows .family-badge')).map(n=>n.textContent).join(',')==='GPCR,Kinase,Ion channel'&&document.getElementById('selected-identity').textContent.includes('Taxonomy 9606')&&document.querySelector('#result-rows tr td:nth-child(3) sup').textContent==='-435'"))throw Error('Metadata or underflow display failed');
await evaluate("document.getElementById('family-filter').value='Kinase';document.getElementById('family-filter').dispatchEvent(new Event('change'))");
if(!await evaluate("document.querySelectorAll('#result-rows tr[data-index]').length===1&&document.getElementById('selected-title').textContent==='Epidermal growth factor receptor'&&document.getElementById('visible-count').textContent==='1 of 3 shown'"))throw Error('Family filter or selection did not update');
await evaluate("document.getElementById('family-filter').value='';document.getElementById('family-filter').dispatchEvent(new Event('change'));document.getElementById('result-search').value='RENDER-2';document.getElementById('result-search').dispatchEvent(new Event('input'))");
if(!await evaluate("document.querySelectorAll('#result-rows tr[data-index]').length===1&&document.getElementById('selected-identity').textContent.includes('RENDER-2')"))throw Error('ID search failed');
await evaluate("document.getElementById('result-search').value='nothing matches';document.getElementById('result-search').dispatchEvent(new Event('input'))");
if(!await evaluate("document.querySelectorAll('#result-rows tr[data-index]').length===0&&document.getElementById('target-detail').hidden&&document.getElementById('visible-count').textContent==='0 of 3 shown'"))throw Error('Empty filter did not clear the selected evidence');
await evaluate("document.getElementById('result-search').value='';document.getElementById('result-search').dispatchEvent(new Event('input'));document.getElementById('result-sort').value='similarity';document.getElementById('result-sort').dispatchEvent(new Event('change'))");
if(!await evaluate("document.querySelector('#result-rows tr').dataset.index==='2'&&document.querySelector('#result-rows tr td').textContent==='03'"))throw Error('Similarity sort changed the original SEA rank');
await evaluate("document.querySelector('#result-rows .target-button').focus()");
await call('Input.dispatchKeyEvent',{type:'keyDown',key:'Enter',code:'Enter',text:'\r',unmodifiedText:'\r',windowsVirtualKeyCode:13});
await call('Input.dispatchKeyEvent',{type:'keyUp',key:'Enter',code:'Enter',windowsVirtualKeyCode:13});
if(!await evaluate("document.getElementById('selected-identity').textContent.includes('RENDER-2')&&document.querySelectorAll('#evidence .pair').length===3"))throw Error('Keyboard target selection failed: '+await evaluate("JSON.stringify({active:document.activeElement.outerHTML,selected:document.getElementById('selected-identity').textContent,pairs:document.querySelectorAll('#evidence .pair').length})"));
await evaluate("document.getElementById('show-pairs').click()");
if(!await evaluate("document.querySelectorAll('#evidence .pair').length===8&&document.getElementById('evidence-count').textContent.includes('8 of 8 saved pairs')"))throw Error('Saved pair expansion failed');
await evaluate("document.getElementById('show-pairs').click()");
if(!await evaluate("document.querySelectorAll('#evidence .pair').length===3"))throw Error('Saved pair collapse failed');
await evaluate("render({...smokeFixture,results:[{...smokeFixture.results[0],raw_score:0,contributing_pairs:0,similar_ligands:[],method:'exact_zero',e_value:3,log_e_value:Math.log(3),p_value:1,log_p_value:0,p_resolution:null}]})");
if(!await evaluate("document.getElementById('summary-tc').textContent==='—'&&document.querySelectorAll('#evidence img').length===0&&document.getElementById('method-note').textContent.includes('No ligand pairs')"))throw Error('Zero score fabricated a maximum similarity or evidence');
await evaluate("render({...smokeOriginal,results:[]})");
if(!await evaluate("document.getElementById('summary-count').textContent==='0'&&document.getElementById('target-detail').hidden"))throw Error('Empty result bundle failed');
await evaluate("render(smokeOriginal)");
if(!await evaluate("Promise.all(['json-export','csv-export'].map(id=>fetch(document.getElementById(id).href).then(r=>r.ok))).then(r=>r.every(Boolean))"))throw Error('Export links failed');
if(!await evaluate("!/[ぁ-んァ-ン一-龯]/.test(document.body.innerText)"))throw Error('Workspace is not English');

// Exercise the set comparison and ensemble paths through the actual worker.
await evaluate("document.querySelector('[data-tab=compare]').click();document.getElementById('compare-submit').click()");
await waitFor("!busy&&document.getElementById('results-title').textContent==='Set comparison results'");
if(!await evaluate("document.getElementById('query-select').options.length===1&&document.getElementById('selected-title').textContent==='Ligand set A vs B'&&document.getElementById('selected-stats').textContent.includes('Set sizes (A / B)')"))throw Error('Set comparison presentation failed');
await evaluate("document.querySelector('[data-tab=predict]').click();document.getElementById('mode').value='ensemble';document.getElementById('predict-submit').click()");
await waitFor("!busy&&document.getElementById('results-title').textContent==='Target search results'");
if(!await evaluate("document.getElementById('query-select').options.length===1&&document.getElementById('query-description').textContent.includes('2 unique molecules')"))throw Error('Ensemble query presentation failed');
await evaluate("document.querySelector('[data-tab=jobs]').click()");
await waitFor("document.querySelectorAll('#job-list .job').length===3");
await evaluate("document.querySelector('#job-list .job button').click()");
await waitFor("!busy&&!document.getElementById('result-panel').hidden");
await evaluate("document.querySelector('[data-tab=predict]').click();document.getElementById('result-panel').scrollIntoView({behavior:'instant'})");
await delay(300);
const desktop=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:true});fs.writeFileSync('var/smoke-desktop.png',Buffer.from(desktop.data,'base64'));
await call('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
await evaluate("document.getElementById('result-panel').scrollIntoView({behavior:'instant'})");
await delay(300);
if(!await evaluate("document.documentElement.scrollWidth<=window.innerWidth"))throw Error('Results overflow on mobile');
const mobile=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:true});fs.writeFileSync('var/smoke-mobile.png',Buffer.from(mobile.data,'base64'));
await evaluate("document.getElementById('results-layout').scrollIntoView({behavior:'instant'})");
await delay(250);
const mobileList=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync('var/smoke-results-mobile.png',Buffer.from(mobileList.data,'base64'));
await evaluate("document.querySelector('#result-rows .target-button').click()");
await delay(400);
const mobileDetail=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync('var/smoke-detail-mobile.png',Buffer.from(mobileDetail.data,'base64'));
if(!await evaluate("document.documentElement.scrollWidth<=window.innerWidth"))throw Error('Selected target overflows on mobile');
// Optional read-only visual review of an existing ChEMBL result bundle.
if(process.env.SEA_SMOKE_PREVIEW){
  const preview=JSON.parse(fs.readFileSync(process.env.SEA_SMOKE_PREVIEW,'utf8'));
  await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
  await evaluate(`render(${JSON.stringify(preview)});document.getElementById('result-panel').scrollIntoView({behavior:'instant'})`);
  await waitFor("Array.from(document.querySelectorAll('#evidence img')).every(img=>img.complete&&img.naturalWidth>0)");
  if(!await evaluate("document.querySelector('.table-wrap').scrollWidth<=document.querySelector('.table-wrap').clientWidth"))throw Error('ChEMBL results need horizontal scrolling on desktop');
  const previewDesktop=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync('var/smoke-chembl-desktop.png',Buffer.from(previewDesktop.data,'base64'));
  await call('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
  await evaluate("document.getElementById('results-layout').scrollIntoView({behavior:'instant'})");
  await delay(300);
  if(!await evaluate("document.documentElement.scrollWidth<=window.innerWidth"))throw Error('ChEMBL results overflow on mobile');
  const previewMobile=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync('var/smoke-chembl-mobile.png',Buffer.from(previewMobile.data,'base64'));
  await evaluate("render(smokeOriginal)");
}
await evaluate("document.querySelector('[data-draw=query]').click()");
await delay(300);
await evaluate(`editor().setMolecule(${JSON.stringify(process.env.SEA_SMOKE_MOL)})`);
await delay(300);const mobileEditor=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync('var/smoke-editor-mobile.png',Buffer.from(mobileEditor.data,'base64'));
if(!await evaluate("document.getElementById('structure-dialog').scrollWidth<=document.getElementById('structure-dialog').clientWidth"))throw Error('Editor dialog overflows on mobile');
if(externalRequests.length)throw Error('Unexpected external resources: '+externalRequests.join(', '));
if(runtimeErrors.length)throw Error('Browser exceptions: '+JSON.stringify(runtimeErrors));
console.log(JSON.stringify({browser_prediction_rows:rows,queries:2,real_jobs:3,checks:['English UI and errors','query switching','family and text filters','retained SEA ranks','keyboard selection','pair expansion','underflow','zero/empty results','exports','comparison','ensemble','history','mobile overflow'],structure_formats:['sdf','mol','smiles','csv'],editor:'local Ketcher 3.18.0',desktop:'var/smoke-desktop.png',mobile:'var/smoke-mobile.png'}));ws.close();
