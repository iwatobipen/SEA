// Optional development smoke check: Node >=22 and Chrome. Services started by smoke.py.
import fs from 'node:fs';
const base='http://127.0.0.1:19021';
const targets=await(await fetch('http://127.0.0.1:19022/json')).json();
const target=targets.find(t=>t.type==='page');
const ws=new WebSocket(target.webSocketDebuggerUrl);await new Promise((r,j)=>{ws.onopen=r;ws.onerror=j;});
let sequence=0;const callbacks=new Map(),externalRequests=[];
ws.onmessage=e=>{const data=JSON.parse(e.data);if(data.method==='Network.requestWillBeSent'&&/^https?:/.test(data.params.request.url)&&!data.params.request.url.startsWith(base+'/'))externalRequests.push(data.params.request.url);if(data.id){const cb=callbacks.get(data.id);callbacks.delete(data.id);data.error?cb.reject(Error(JSON.stringify(data.error))):cb.resolve(data.result);}};
function call(method,params={}){return new Promise((resolve,reject)=>{const id=++sequence;callbacks.set(id,{resolve,reject});ws.send(JSON.stringify({id,method,params}));});}
async function evaluate(expression){const r=await call('Runtime.evaluate',{expression,awaitPromise:true,returnByValue:true});if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));return r.result.value;}
const delay=ms=>new Promise(r=>setTimeout(r,ms));
await call('Network.enable');await call('Network.setCacheDisabled',{cacheDisabled:true});
await call('Page.enable');await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
await call('Page.navigate',{url:base});await delay(1000);
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
if(!await evaluate("document.getElementById('query-import-status').classList.contains('error')&&document.getElementById('query-import-status').textContent.includes('レコード 2')&&document.getElementById('query').value==='CCO\\nCCCO'"))throw Error('Invalid SDF changed the input or lost its record number');
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
await evaluate("editorStatus('構造を描画し、「入力に追加」を押してください。')");
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
const rows=await evaluate("document.querySelectorAll('#result-rows tr').length");
const jobCount=await evaluate("fetch('/api/v1/jobs').then(r=>r.json()).then(jobs=>jobs.length)");
if(jobCount!==1)throw Error('Duplicate submission created extra jobs: '+jobCount);
if(rows<1)throw Error('No results: '+await evaluate("document.getElementById('message').textContent"));
// Exercise metadata cells with a transient rendering fixture; never store fake annotations.
const metadataRendered=await evaluate(`(()=>{const original=currentResult;const families=['GPCR','Kinase','Ion channel'];const fixture={...original,results:families.map((family,i)=>({...original.results[0],target:{id:'RENDER-'+i,name:'Metadata rendering fixture',family,protein_classes:[[family]],organism:'Homo sapiens',tax_id:9606}}))};render(fixture);const rows=Array.from(document.querySelectorAll('#result-rows tr'));const headers=Array.from(document.querySelectorAll('thead th')).map(th=>th.textContent);const valid=headers.length===9&&headers[2]==='ファミリー'&&headers[3]==='生物種'&&rows.every((row,i)=>row.children.length===9&&row.children[2].textContent===families[i]&&row.children[3].textContent==='Homo sapiens (9606)'&&row.children[1].textContent.includes('RENDER-'+i));render(original);return valid;})()`);
if(!metadataRendered)throw Error('Target family/species columns failed to render');

const desktop=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:true});fs.writeFileSync('var/smoke-desktop.png',Buffer.from(desktop.data,'base64'));
await call('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
await delay(300);const mobile=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:true});fs.writeFileSync('var/smoke-mobile.png',Buffer.from(mobile.data,'base64'));
await evaluate("document.querySelector('[data-draw=query]').click()");
await delay(300);
await evaluate(`editor().setMolecule(${JSON.stringify(process.env.SEA_SMOKE_MOL)})`);
await delay(300);const mobileEditor=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync('var/smoke-editor-mobile.png',Buffer.from(mobileEditor.data,'base64'));
if(!await evaluate("document.getElementById('structure-dialog').scrollWidth<=document.getElementById('structure-dialog').clientWidth"))throw Error('Editor dialog overflows on mobile');
if(externalRequests.length)throw Error('Structure editor requested external resources: '+externalRequests.join(', '));
console.log(JSON.stringify({browser_prediction_rows:rows,structure_formats:['sdf','mol','smiles','csv'],editor:'local Ketcher 3.18.0',desktop:'var/smoke-desktop.png',mobile:'var/smoke-mobile.png'}));ws.close();
