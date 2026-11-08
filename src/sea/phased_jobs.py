"""Evidence previews, immutable screening results, and linked precision jobs."""
from functools import lru_cache
import json
import resource
import secrets
import time

from .annotations import load_annotations, target_info
from .chemistry import normalize, raw_score
from .data import home, load_dataset
from .jobs import Cancelled, background_cache, connection
from .precision import details, pending, profile
from .precision_cache import PrecisionCache
from .statistics import significance


def encode(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False)


def displayed(bundle, top_k):
    output = {**bundle, "results": sorted(bundle["results"], key=lambda r: (
        r["query_index"] or 0, r["log_e_value"] if r["log_e_value"] is not None else float("inf"), -r["raw_score"]))}
    if top_k is not None:
        counts, rows = {}, []
        for row in output["results"]:
            key = row["query_index"]
            counts[key] = counts.get(key, 0) + 1
            if counts[key] <= top_k:
                rows.append(row)
        output["results"] = rows
    return output


def create_refinement(source_id, owner, name, con):
    source = con.execute("SELECT * FROM jobs WHERE id=? AND owner=?", (source_id, owner)).fetchone()
    if not source or source["status"] != "completed" or not source["checkpoint"] or source["kind"] == "refinement":
        raise ValueError("Refinement requires a completed balanced/research screening job.")
    existing = con.execute("SELECT id FROM jobs WHERE parent_id=? AND profile=? AND kind='refinement'", (source_id, name)).fetchone()
    if existing:
        return existing["id"]
    payload = json.loads(source["payload"])
    payload.update(source_job_id=source_id, base_kind=source["kind"], statistics_profile=name)
    jid, now = secrets.token_hex(16), time.time()
    con.execute("INSERT INTO jobs(id,owner,kind,payload,status,created,updated,parent_id,profile,phase) VALUES(?,?,?,?,?,?,?,?,?,?)",
                (jid, owner, "refinement", encode(payload), "queued", now, now, source_id, name, "queued"))
    con.execute("UPDATE jobs SET refinement_job_id=COALESCE(refinement_job_id,?) WHERE id=?", (jid, source_id))
    return jid


@lru_cache(maxsize=2)
def precision_cache(root, dataset_id):
    return PrecisionCache(root, dataset_id, background_cache(root, dataset_id))


def process(job):
    jid, started = job["id"], time.perf_counter()
    last_check, last_publish = [0.], [0.]
    payload = json.loads(job["payload"])
    top = payload.get("top_k", 20) if payload.get("base_kind", job["kind"]) == "prediction" else None
    bundle = None

    def check():
        if time.monotonic() - last_check[0] < .2:
            return
        last_check[0] = time.monotonic()
        with connection() as con:
            row = con.execute("SELECT cancel FROM jobs WHERE id=?", (jid,)).fetchone()
            con.execute("UPDATE jobs SET updated=? WHERE id=?", (time.time(), jid))
        if row["cancel"]:
            raise Cancelled()

    def publish(fraction, force=False):
        if not force and time.monotonic() - last_publish[0] < .5:
            return
        bundle["execution"] = {"elapsed_seconds": time.perf_counter() - started,
                               "worker_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        bundle["analysis"]["evaluated_comparisons"] = sum(r.get("statistics", {}).get("phase") == "final" for r in bundle["results"])
        with connection() as con:
            con.execute("UPDATE jobs SET preview=?,checkpoint=?,phase=?,progress=?,updated=? WHERE id=?",
                        (encode(displayed(bundle, top)), encode(bundle), bundle["analysis"]["phase"], fraction, time.time(), jid))
        last_publish[0] = time.monotonic()

    try:
        check()
        data, structures = load_dataset(payload["dataset_id"])
        pool, model = list(structures.values()), data["model"]
        effective = profile(data, payload["statistics_profile"])
        annotations = load_annotations(payload["dataset_id"])
        if job["kind"] == "refinement":
            with connection() as con:
                current = con.execute("SELECT checkpoint FROM jobs WHERE id=?", (jid,)).fetchone()[0]
                source = con.execute("SELECT checkpoint FROM jobs WHERE id=?", (job["parent_id"],)).fetchone()[0]
            bundle = json.loads(current or source)
            bundle["analysis"].update(phase="refining", complete=False, profile=effective,
                                      source_job_id=job["parent_id"], refinement_job_id=jid)
            remaining = [r for r in bundle["results"] if r["method"] == "monte_carlo" and
                         not (r.get("statistics", {}).get("phase") == "final" and
                              r["statistics"].get("profile_id") == effective["id"])]
            remaining.sort(key=lambda r: (r["log_e_value"], -r["raw_score"]))
            cache = precision_cache(str(home().resolve()), payload["dataset_id"])
            final_model = {**model, "empirical_trials": effective["planned_trials"]}
            publish(0., True)
            try:
                for index, row in enumerate(remaining):
                    check()
                    cache.on_progress = lambda done, total: publish((index + done / total) / max(1, len(remaining)))
                    stat = significance(row["raw_score"], row["left_size"], row["right_size"], pool,
                                        final_model, row["comparisons"], check, cache)
                    row.update(stat)
                    row["statistics"] = details(row, effective["planned_trials"], "final", effective)
                    publish((index + 1) / len(remaining), True)
            finally:
                cache.on_progress = lambda *_: None
            for row in bundle["results"]:
                if row["method"] != "monte_carlo":
                    row["statistics"] = details(row, 0, "final", effective)
            check()
            with connection() as con:
                con.execute("BEGIN IMMEDIATE")
                if con.execute("SELECT cancel FROM jobs WHERE id=?", (jid,)).fetchone()[0]:
                    raise Cancelled()
                bundle["analysis"].update(phase="final", complete=True)
                result = encode(displayed(bundle, top))
                con.execute("UPDATE jobs SET status='completed',result=?,preview=?,checkpoint=?,progress=1,phase='final',updated=? WHERE id=?",
                            (result, result, encode(bundle), time.time(), jid))
            return

        bundle = {"dataset": data["manifest"], "threshold": model["threshold"], "results": [],
                  "execution": {}, "analysis": {"phase": "evidence", "complete": False,
                  "profile": effective, "source_job_id": jid, "refinement_job_id": None,
                  "total_comparisons": 0, "evaluated_comparisons": 0}}
        if job["kind"] == "comparison":
            left, ar = normalize(payload["left"], data['manifest']['fingerprint'])
            right, br = normalize(payload["right"], data['manifest']['fingerprint'])
            entries = [(None, left, right, None)]
            bundle["normalization"] = {"left": ar, "right": br}
            ndb = 1
        else:
            normalized = [normalize(q, data['manifest']['fingerprint']) for q in payload["queries"]]
            targets = [(t, [structures[s] for s in t["smiles"]]) for t in data["targets"]]
            entries = ((i, a, b, t) for i, (a, _) in enumerate(normalized) for t, b in targets)
            bundle["normalization"] = [records for _, records in normalized]
            ndb = len(targets)
        for index, left, right, target in entries:
            check()
            score, count, evidence = raw_score(left, right, model["threshold"], check=check)
            row = {"query_index": index, "raw_score": score, "contributing_pairs": count,
                   "left_size": len(left), "right_size": len(right), "similar_ligands": evidence, **pending(ndb)}
            if target:
                row["target"] = target_info(target, annotations)
            row["statistics"] = details(row, 0, "pending", effective)
            bundle["results"].append(row)
        if annotations:
            bundle["target_annotations"] = {k: v for k, v in annotations.items() if k != "targets"}
        bundle["analysis"]["total_comparisons"] = len(bundle["results"])
        publish(0., True)
        bundle["analysis"]["phase"] = "screening"
        # The screening budget is fixed at 999 regardless of the dataset's
        # original default. Its cache is isolated by that effective budget.
        screening_model = {**model, "empirical_trials": 999}
        screening_cache = precision_cache(str(home().resolve()), payload["dataset_id"])
        try:
            for index, row in enumerate(bundle["results"]):
                check()
                screening_cache.on_progress = lambda *_: publish(index / max(1, len(bundle["results"])))
                row.update(significance(row["raw_score"], row["left_size"], row["right_size"], pool,
                                        screening_model, ndb, check, screening_cache))
                row["statistics"] = details(row, 999, "screening", effective)
                publish((index + 1) / len(bundle["results"]))
        finally:
            screening_cache.on_progress = lambda *_: None
        check()
        publish(1., True)
        with connection() as con:
            con.execute("BEGIN IMMEDIATE")
            if con.execute("SELECT cancel FROM jobs WHERE id=?", (jid,)).fetchone()[0]:
                raise Cancelled()
            con.execute("UPDATE jobs SET status='completed',checkpoint=? WHERE id=?", (encode(bundle), jid))
            child_id = create_refinement(jid, job["owner"], effective["name"], con)
            bundle["analysis"]["refinement_job_id"] = child_id
            con.execute("UPDATE jobs SET result=?,preview=?,checkpoint=?,progress=1,phase='screening',updated=? WHERE id=?",
                        (encode(displayed(bundle, top)), encode(displayed(bundle, top)), encode(bundle), time.time(), jid))
    except Cancelled:
        with connection() as con:
            con.execute("UPDATE jobs SET status='cancelled',updated=? WHERE id=?", (time.time(), jid))
    except Exception as exc:
        with connection() as con:
            con.execute("UPDATE jobs SET status='failed',error=?,updated=? WHERE id=?", (str(exc), time.time(), jid))
