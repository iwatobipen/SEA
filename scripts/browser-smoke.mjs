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
await call('Page.enable');await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
await call('Page.navigate',{url:base});await delay(1000);
await evaluate(`document.getElementById('key').value=${JSON.stringify(process.env.SEA_SMOKE_KEY)};document.getElementById('login').requestSubmit();`);
for(let i=0;i<30;i++){if(await evaluate("!document.getElementById('workspace').hidden"))break;await delay(200);}
if(!await evaluate("!document.getElementById('workspace').hidden"))throw Error('Login failed');
await evaluate("document.getElementById('query').value='CCO\\nCCCO';document.getElementById('predict-submit').click()");
for(let i=0;i<60;i++){if(await evaluate("!document.getElementById('result-panel').hidden"))break;await delay(500);}
const rows=await evaluate("document.querySelectorAll('#result-rows tr').length");
if(rows<1)throw Error('No results: '+await evaluate("document.getElementById('message').textContent"));
const desktop=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:true});fs.writeFileSync('var/smoke-desktop.png',Buffer.from(desktop.data,'base64'));
await call('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
await delay(300);const mobile=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:true});fs.writeFileSync('var/smoke-mobile.png',Buffer.from(mobile.data,'base64'));
console.log(JSON.stringify({browser_prediction_rows:rows,desktop:'var/smoke-desktop.png',mobile:'var/smoke-mobile.png'}));ws.close();
