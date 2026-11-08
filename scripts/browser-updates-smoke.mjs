// Deterministic display/polling regressions; fixtures never reach the job API.
export async function checkUpdates({evaluate, waitFor, delay}) {
  await evaluate(`(async()=>{
    window.updateTest={api,setTimeout,original:structuredClone(smokeOriginal),job:currentJob,stage:0,calls:0,failures:{},errors:0,aborts:0};
    const state=updateTest;
    state.screen=structuredClone(state.original);
    state.screen.analysis={...state.screen.analysis,phase:'screening',complete:false,source_job_id:'ui-screen',refinement_job_id:'ui-final'};
    for(const row of state.screen.results)if(row.method==='monte_carlo')row.statistics.phase='screening';
    state.screen.results.reverse();
    state.evidence=structuredClone(state.screen);state.evidence.analysis.phase='evidence';
    for(const row of state.evidence.results){
      Object.assign(row,{method:'pending',p_value:null,log_p_value:null,e_value:null,log_e_value:null,z_score:null});
      Object.assign(row.statistics,{phase:'pending',status:'pending',trials:0,exceedances:null});
    }
    state.final=structuredClone(state.original);state.final.analysis.source_job_id='ui-screen';state.final.analysis.refinement_job_id='ui-final';
    state.refining=structuredClone(state.screen);state.refining.analysis.phase='refining';
    api=async(path,options)=>{
      if(!path.startsWith('/jobs/ui-'))return state.api(path,options);
      state.calls++;
      const route=path.endsWith('/preview')?'preview':path.endsWith('/result')?'result':'status';
      if(state.hangStatus&&route==='status'){
        state.hangStatus=false;
        return new Promise((resolve,reject)=>options.signal.addEventListener('abort',()=>{state.aborts++;reject(new DOMException('Simulated timeout','AbortError'));},{once:true}));
      }
      const failure=state.failures[route];
      if(failure?.remaining>0){failure.remaining--;state.errors++;const error=new TypeError('Simulated interrupted update');if(failure.status)error.status=failure.status;throw error;}
      const final=path.startsWith('/jobs/ui-final');
      if(route==='status')return {id:final?'ui-final':'ui-screen',kind:final?'refinement':'prediction',profile:'balanced',status:(final?state.stage===2:state.stage>=1)?'completed':'running',progress:0.5,parent_id:final?'ui-screen':null,refinement_job_id:final?null:'ui-final'};
      return structuredClone(final?(state.stage===2?state.final:state.refining):route==='preview'&&state.stage===0?state.evidence:state.screen);
    };
    await watch('ui-screen');
  })()`);
  try {
    if (!await evaluate("currentResult.results.some(r=>r.p_value===null)&&currentResult.analysis.phase==='evidence'")) throw Error('Evidence must identify statistics as pending');
    await evaluate(`document.getElementById('query-select').value='1';changeQuery();
      document.getElementById('result-sort').value='name';document.getElementById('result-search').value='Alcohol';updateResults();
      updateTest.selected=resultRows[selectedIndex].result.target.id;updateTest.stage=1;`);
    await waitFor("currentResult?.analysis?.phase==='screening'&&currentResult.results.every(r=>r.p_value!==null&&r.e_value!==null)");
    if (!await evaluate("document.getElementById('query-select').value==='1'&&document.getElementById('result-search').value==='Alcohol'&&document.getElementById('result-sort').value==='name'&&resultRows[selectedIndex].result.target.id===updateTest.selected&&!document.getElementById('json-export').hidden")) throw Error('First statistics lost the selected query, target, filters or export');
    await evaluate("updateTest.rows=document.getElementById('result-rows').innerHTML;updateTest.stage=2;");
    await waitFor('Boolean(pendingResult?.data.analysis?.complete)');
    if (!await evaluate("currentResult.analysis.phase==='screening'&&document.getElementById('result-rows').innerHTML===updateTest.rows")) throw Error('Later refinement changed the displayed numeric result without applying it');
    await evaluate("document.getElementById('apply-statistics').click()");
    if (!await evaluate("currentResult.analysis.complete&&resultRows[selectedIndex].result.target.id===updateTest.selected")) throw Error('Applying final statistics lost selection');

    // Lose status, result and preview reads at different transitions.
    await evaluate("(async()=>{updateTest.stage=0;await watch('ui-screen');updateTest.failures={status:{remaining:1},result:{remaining:1},preview:{remaining:1,status:503}};})()");
    await waitFor('updateTest.errors===1');
    await evaluate('updateTest.stage=1');
    await waitFor("updateTest.errors===3&&currentResult.results.every(r=>r.p_value!==null)");
    await evaluate('updateTest.stage=2');
    await waitFor('Boolean(pendingResult?.data.analysis?.complete)');
    if (!await evaluate("!document.getElementById('message').textContent&&currentResult.results.every(r=>r.e_value!==null)")) throw Error('Update reads did not recover automatically');

    // A new watch must invalidate a queued retry from the previous job.
    await evaluate("(async()=>{updateTest.failures={status:{remaining:1}};await watch('ui-screen');updateTest.failures={};await watch('ui-final');updateTest.afterNewWatch=updateTest.calls;})()");
    await delay(1200);
    if (!await evaluate("currentJob==='ui-final'&&currentResult.analysis.complete&&updateTest.calls===updateTest.afterNewWatch")) throw Error('An obsolete retry replaced the newer analysis');

    // Authentication failures must stop, rather than retrying indefinitely.
    await evaluate("(async()=>{updateTest.failures={status:{remaining:5,status:401}};await watch('ui-screen');updateTest.afterUnauthorized=updateTest.calls;})()");
    await delay(1200);
    if (!await evaluate("!busy&&updateTest.calls===updateTest.afterUnauthorized&&document.getElementById('message').textContent")) throw Error('Authentication failure was silently retried');
    // Exercise the real request AbortController without waiting 15 wall seconds.
    await evaluate("(async()=>{updateTest.failures={};updateTest.hangStatus=true;window.setTimeout=(fn,ms,...args)=>updateTest.setTimeout.call(window,fn,ms===15000?30:ms,...args);await watch('ui-final');})()");
    await waitFor("updateTest.aborts===1&&currentJob==='ui-final'&&currentResult?.analysis?.complete");
    return ['automatic first statistics', 'preserved query/target/filters', 'explicit later refinement', 'status/result/preview retry', 'obsolete retry isolation', 'terminal authentication', 'request timeout recovery'];
  } finally {
    await evaluate("watchGeneration++;clearTimeout(pollTimer);api=updateTest.api;window.setTimeout=updateTest.setTimeout;pendingResult=null;refinementJob=null;currentJob=updateTest.job;artifactReady=true;render(updateTest.original);setBusy(false);message();document.getElementById('progress-panel').hidden=true;");
  }
}
