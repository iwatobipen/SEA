# Pending statistics and automatic UI updates

Investigated on November 8, 2026. A short-lived **Pending** state is expected:
the evidence preview contains molecular matches before statistical evaluation
finishes. A completed screening or final result contains numerical P/E values.
Pending should not require a page reload once the first statistics are ready.

## Findings

Two UI paths could leave an old evidence preview visible:

1. Every update after the initial preview required an explicit **Apply…** click.
   This included the first completed screening result, so the screen could keep
   showing Pending after both screening and refinement had finished. Opening a
   completed job again loaded its numerical values directly.
2. A failed status, preview or result request stopped polling after one error.
   Worker calculations continued, but the browser did not retry. A reload or
   reopening the job started a new watch.

All 38 completed lab results inspected during this investigation contained
P/E values with no pending rows. The display path was
reproduced independently using a deterministic browser fixture: the unmodified
UI never replaced its evidence preview after a completed screening response.
The original browser request log was not available, so it is not possible to
identify which path the reported session encountered.

## Corrected behavior

- A committed result with P/E values for all returned rows automatically
  replaces a preview containing pending values. In ordinary Balanced/Research
  searches this is the completed 999-trial screening artifact. Partial previews
  do not silently become completed statistical results.
- The update preserves the query, filters, display sort, scroll position,
  the saved-pair expansion setting and selected target when it remains in the returned set.
  The first statistical ranking can differ from the initial raw-score list.
- After numerical values are visible, further refinement remains explicit via
  **Apply latest statistics** or **Apply final results**. The banner identifies
  when final results are available. No statistical estimate is fabricated while
  sampling is still pending or cancelled.
- Read-only job status/result/preview requests time out after 15 seconds.
  Network failures, timeouts, HTTP 408/409/429 and server errors retry with
  delays of 1, 2, 4, 8 and at most 10 seconds. The UI indicates the interruption
  and clears that message after recovery. Submissions are not retried.
- Authentication and other permanent client errors stop with an error message.
  Starting another watch or logging out invalidates earlier timers and results.

Scoring, statistical budgets, reference data and stored job artifacts are
unchanged. Frontend asset version 8 loads the correction; an already-open tab
needs one reload to load the new JavaScript. Subsequent searches do not need a
reload to move from Pending to the first numerical statistics.

## Regression checks

The [browser regression helper](../scripts/browser-updates-smoke.mjs), called
by the existing isolated smoke check, exercises evidence → screening → final,
selection/filter retention, manual later updates, failed status/result/preview
reads, retry cancellation when switching jobs, terminal authentication errors
and an aborted request followed by recovery. It uses in-memory response fixtures;
they are never submitted or saved as scientific jobs.

```bash
pixi run test
pixi run python scripts/smoke.py
```

The optional smoke check requires local Node and Chrome and uses temporary
API/worker services, keys, jobs and browser profiles. Lab records are not used
for these regression submissions.
