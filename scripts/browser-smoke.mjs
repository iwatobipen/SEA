// Optional development smoke check: Node >=22 and Chrome. Services started by smoke.py.
import fs from 'node:fs';
const base='http://127.0.0.1:19021';
const targets=await(await fetch('http://127.0.0.1:19022/json')).json();
const target=targets.find(t=>t.type==='page');
const ws=new WebSocket(target.webSocketDebuggerUrl);await new Promise((r,j)=>{ws.onopen=r;ws.onerror=j;});
let sequence=0;const callbacks=new Map();
ws.onmessage=e=>{const data=JSON.parse(e.data);if(data.id){const cb=callbacks.get(data.id);callbacks.delete(data.id);data.error?cb.reject(Error(JSON.stringify(data.error))):cb.resolve(data.result);}};
function call(method,params={}){return new Promise((resolve,reject)=>{const id=++sequence;callbacks.set(id,{resolve,reject});ws.send(JSON.stringify({id,method,params}));});}
async function evaluate(expression){const r=await call('Runtime.evaluate',{expression,awaitPromise:true,returnByValue:true});if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));return r.result.value;}
const delay=ms=>new Promise(r=>setTimeout(r,ms));
await call('Network.enable');await call('Network.setCacheDisabled',{cacheDisabled:true});
await call('Page.enable');await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
await call('Page.navigate',{url:base});await delay(1000);
await evaluate(`document.getElementById('key').value=${JSON.stringify(process.env.SEA_SMOKE_KEY)};document.getElementById('login').requestSubmit();`);
for(let i=0;i<30;i++){if(await evaluate("!document.getElementById('workspace').hidden"))break;await delay(200);}
if(!await evaluate("!document.getElementById('workspace').hidden"))throw Error('Login failed');
await evaluate("document.getElementById('query').value='CCO\\nCCCO';document.getElementById('predict-submit').click();document.getElementById('predict-submit').click()");
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
console.log(JSON.stringify({browser_prediction_rows:rows,desktop:'var/smoke-desktop.png',mobile:'var/smoke-mobile.png'}));ws.close();
