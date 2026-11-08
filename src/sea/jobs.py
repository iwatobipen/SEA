import hashlib
import hmac
import json
import secrets
import sqlite3
import time
from functools import lru_cache
from .data import home, load_dataset
from .chemistry import normalize, raw_score
from .statistics import significance
from .annotations import target_info, load_annotations
from .background import NullCache


class Cancelled(Exception):
    pass


def connection():
    con = sqlite3.connect(home()/"service.sqlite", timeout=30)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA journal_mode=WAL")
    con.executescript("""
    CREATE TABLE IF NOT EXISTS keys (id TEXT PRIMARY KEY, name TEXT, hash TEXT UNIQUE, created REAL);
    CREATE TABLE IF NOT EXISTS sessions (hash TEXT PRIMARY KEY, owner TEXT, expires REAL);
    CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY, owner TEXT, kind TEXT, payload TEXT,
      status TEXT, progress REAL DEFAULT 0, result TEXT, error TEXT, cancel INTEGER DEFAULT 0,
      created REAL, updated REAL, worker TEXT);
    """)
    columns = {row[1] for row in con.execute("PRAGMA table_info(jobs)")}
    additions = {"parent_id": "TEXT", "refinement_job_id": "TEXT", "phase": "TEXT",
                 "preview": "TEXT", "checkpoint": "TEXT", "profile": "TEXT DEFAULT 'legacy'"}
    if set(additions) - columns:
        con.execute("BEGIN IMMEDIATE")
        columns = {row[1] for row in con.execute("PRAGMA table_info(jobs)")}
        for name, declaration in additions.items():
            if name not in columns:
                con.execute(f"ALTER TABLE jobs ADD COLUMN {name} {declaration}")
        con.execute("CREATE UNIQUE INDEX IF NOT EXISTS refinement_per_profile ON jobs(parent_id,profile) WHERE kind='refinement'")
        con.commit()
    return con


def hashed(token):
    return hashlib.sha256(token.encode()).hexdigest()


def create_key(name):
    token = secrets.token_urlsafe(32)
    with connection() as con:
        con.execute("INSERT INTO keys VALUES (?,?,?,?)", (secrets.token_hex(12), name, hashed(token), time.time()))
    return token


def authenticate(token):
    if not token:
        return None
    with connection() as con:
        row = con.execute("SELECT id FROM keys WHERE hash=?", (hashed(token),)).fetchone()
    return row["id"] if row else None


def submit(owner, kind, payload):
    jid, now = secrets.token_hex(16), time.time()
    with connection() as con:
        con.execute("INSERT INTO jobs (id,owner,kind,payload,status,created,updated,profile,phase) VALUES (?,?,?,?,?,?,?,?,?)",
                    (jid,owner,kind,json.dumps(payload),"queued",now,now,payload.get("statistics_profile","legacy"),"queued"))
    return jid


def claim(worker, role="any"):
    with connection() as con:
        con.execute("BEGIN IMMEDIATE")
        restriction = " AND kind!='refinement'" if role == "foreground" else ""
        row = con.execute("SELECT * FROM jobs WHERE status='queued' AND cancel=0" + restriction +
                          " ORDER BY CASE WHEN kind='refinement' THEN 1 ELSE 0 END,created LIMIT 1").fetchone()
        if row:
            con.execute("UPDATE jobs SET status='running',worker=?,updated=? WHERE id=?", (worker,time.time(),row["id"]))
    return dict(row) if row else None


@lru_cache(maxsize=2)
def background_cache(root, dataset_id):
    """Reuse null samples in memory and on disk across workers and restarts.

    The dataset ID includes calibration settings, seed and the background pool.
    The root isolates deployments with separate SEA_HOME directories.
    """
    return NullCache(root, dataset_id)


def process(job):
    if job["kind"] == "refinement" or json.loads(job["payload"]).get("statistics_profile", "legacy") != "legacy":
        from .phased_jobs import process as phased_process
        return phased_process(job)
    return _legacy_process(job)


def _legacy_process(job):
    jid = job["id"]
    last_check = [0.]
    last_progress = [0.]
    def check():
        if time.monotonic()-last_check[0] < .2:
            return
        last_check[0] = time.monotonic()
        with connection() as con:
            row = con.execute("SELECT cancel FROM jobs WHERE id=?", (jid,)).fetchone()
            con.execute("UPDATE jobs SET updated=? WHERE id=?", (time.time(),jid))
        if row["cancel"]:
            raise Cancelled()
    started = time.perf_counter()
    try:
        payload = json.loads(job["payload"])
        data, structures = load_dataset(payload["dataset_id"])
        pool = list(structures.values())
        model = data["model"]
        annotations = load_annotations(payload["dataset_id"])
        output = {"dataset": data["manifest"], "threshold": model["threshold"], "results": []}
        if job["kind"] == "comparison":
            a, ar = normalize(payload["left"])
            b, br = normalize(payload["right"])
            comparisons = [(None,a,b)]
            comparison_count = 1
            output["normalization"] = {"left": ar, "right": br}
            ndb = 1
        else:
            queries = payload["queries"]
            normalized = [normalize(q) for q in queries]
            targets = [(t, [structures[s] for s in t["smiles"]]) for t in data["targets"]]
            comparisons = ((i,a,b,t) for i,(a,_) in enumerate(normalized) for t,b in targets)
            comparison_count = len(normalized) * len(targets)
            output["normalization"] = [rows for _, rows in normalized]
            ndb = len(data["targets"])
        null_cache = background_cache(str(home().resolve()), payload["dataset_id"])
        for i, entry in enumerate(comparisons):
            check()
            index,a,b,*target = entry
            score,count,evidence = raw_score(a,b,model["threshold"],check=check)
            stat = significance(score,len(a),len(b),pool,model,ndb,check,null_cache)
            result = {"query_index": index, "raw_score": score, "contributing_pairs": count,
                      "left_size": len(a), "right_size": len(b), "similar_ligands": evidence, **stat}
            if target:
                result["target"] = target_info(target[0], annotations)
            output["results"].append(result)
            if i + 1 == comparison_count or time.monotonic() - last_progress[0] >= .2:
                with connection() as con:
                    con.execute("UPDATE jobs SET progress=?, updated=? WHERE id=?", ((i+1)/comparison_count,time.time(),jid))
                last_progress[0] = time.monotonic()
        if annotations:
            output["target_annotations"] = {k:v for k,v in annotations.items() if k != "targets"}
        output["results"].sort(key=lambda r: (r["query_index"] or 0, r["log_e_value"], -r["raw_score"]))
        if job["kind"] != "comparison":
            top = payload.get("top_k",20)
            output["results"] = [r for q in range(len(payload["queries"])) for r in [r for r in output["results"] if r["query_index"]==q][:top]]
        check()
        import resource
        output["execution"] = {"elapsed_seconds":time.perf_counter()-started,"worker_peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        with connection() as con:
            con.execute("UPDATE jobs SET status=CASE WHEN cancel=1 THEN 'cancelled' ELSE 'completed' END,result=?,updated=? WHERE id=?",
                        (json.dumps(output,ensure_ascii=False,allow_nan=False),time.time(),jid))
    except Cancelled:
        with connection() as con:
            con.execute("UPDATE jobs SET status='cancelled',updated=? WHERE id=?", (time.time(),jid))
    except Exception as exc:
        with connection() as con:
            con.execute("UPDATE jobs SET status='failed',error=?,updated=? WHERE id=?", (str(exc),time.time(),jid))


def work(stop, name, role="any"):
    import signal
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while not stop.is_set():
        job = claim(name, role)
        if job:
            process(job)
        else:
            stop.wait(.5)


def run_workers(count=2):
    import threading
    import signal
    # A service-wide lock prevents a second supervisor from reclaiming live jobs.
    import fcntl
    lock = open(home()/"worker.lock", "w")
    fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
    with connection() as con:
        con.execute("UPDATE jobs SET status=CASE WHEN cancel=1 THEN 'cancelled' ELSE 'queued' END,worker=NULL WHERE status='running'")
    from .data import preload_datasets
    preload_datasets()
    # Processes isolate RDKit workloads and avoid sharing the Python interpreter.
    import multiprocessing
    stop = multiprocessing.Event()
    shutdown = [False]
    def request_shutdown(*_):
        # Signal handlers must not acquire multiprocessing synchronization locks.
        shutdown[0] = True
    signal.signal(signal.SIGTERM, request_shutdown)
    signal.signal(signal.SIGINT, request_shutdown)
    workers = [multiprocessing.Process(target=work,args=(stop,f"worker-{i}",
               "foreground" if count > 1 and i == 0 else "any")) for i in range(count)]
    for worker in workers:
        worker.start()
    try:
        while not shutdown[0]:
            time.sleep(.2)
            if any(not p.is_alive() for p in workers):
                raise RuntimeError("A worker stopped unexpectedly. Restart the worker supervisor.")
    finally:
        stop.set()
        for worker in workers:
            worker.join(timeout=5)
            if worker.is_alive():
                worker.terminate()
                worker.join()
        lock.close()
