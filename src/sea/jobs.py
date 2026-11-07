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
        con.execute("INSERT INTO jobs (id,owner,kind,payload,status,created,updated) VALUES (?,?,?,?,?,?,?)",
                    (jid,owner,kind,json.dumps(payload),"queued",now,now))
    return jid


def claim(worker):
    with connection() as con:
        con.execute("BEGIN IMMEDIATE")
        row = con.execute("SELECT * FROM jobs WHERE status='queued' AND cancel=0 ORDER BY created LIMIT 1").fetchone()
        if row:
            con.execute("UPDATE jobs SET status='running',worker=?,updated=? WHERE id=?", (worker,time.time(),row["id"]))
    return dict(row) if row else None


@lru_cache(maxsize=2)
def background_cache(root, dataset_id):
    """Reuse deterministic null samples across jobs in this worker process.

    The dataset ID includes calibration settings, seed and the background pool.
    The root isolates deployments with separate SEA_HOME directories.
    """
    return {}


def process(job):
    jid = job["id"]
    last_check = [0.]
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
            output["normalization"] = {"left": ar, "right": br}
            ndb = 1
        else:
            queries = payload["queries"]
            normalized = [normalize(q) for q in queries]
            comparisons = [(i,a,[structures[s] for s in t["smiles"]],t)
                           for i,(a,_) in enumerate(normalized) for t in data["targets"]]
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
            with connection() as con:
                con.execute("UPDATE jobs SET progress=?, updated=? WHERE id=?", ((i+1)/len(comparisons),time.time(),jid))
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


def work(stop, name):
    import signal
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while not stop.is_set():
        job = claim(name)
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
    # Processes isolate RDKit workloads and avoid sharing the Python interpreter.
    import multiprocessing
    stop = multiprocessing.Event()
    shutdown = [False]
    def request_shutdown(*_):
        # Signal handlers must not acquire multiprocessing synchronization locks.
        shutdown[0] = True
    signal.signal(signal.SIGTERM, request_shutdown)
    signal.signal(signal.SIGINT, request_shutdown)
    workers = [multiprocessing.Process(target=work,args=(stop,f"worker-{i}")) for i in range(count)]
    for worker in workers:
        worker.start()
    try:
        while not shutdown[0]:
            time.sleep(.2)
            if any(not p.is_alive() for p in workers):
                raise RuntimeError("workerが停止しました。supervisorを再起動してください")
    finally:
        stop.set()
        for worker in workers:
            worker.join(timeout=5)
            if worker.is_alive():
                worker.terminate()
                worker.join()
        lock.close()
