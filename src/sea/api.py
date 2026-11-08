import csv
import io
import json
import os
from pathlib import Path
import secrets
import time
from contextlib import asynccontextmanager
from typing import Literal
from fastapi import FastAPI, Depends, HTTPException, Request, Response
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.security import APIKeyHeader
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from .chemistry import normalize, molecule, InputError
from .data import list_datasets, load_dataset, preload_datasets
from .jobs import connection, authenticate, hashed, submit
from .annotations import enrich_results, target_info, load_annotations
from .structure_inputs import StructureFormat, StructureImportError, import_structures

@asynccontextmanager
async def lifespan(app):
    preload_datasets()
    yield


app = FastAPI(title="SEA Open", version="0.1.0", lifespan=lifespan, description="Ligand ensemble similarity. E-values are expected random hits, not binding probabilities.")
WEB = Path(__file__).parent/"web"
app.mount("/static", StaticFiles(directory=WEB), name="static")


@app.middleware("http")
async def request_limits(request: Request, call_next):
    maximum = int(os.getenv("SEA_MAX_UPLOAD_BYTES",10485760))
    if request.method in ("POST", "PUT"):
        # Buffer with a hard bound even for chunked requests.
        body = bytearray()
        async for chunk in request.stream():
            body.extend(chunk)
            if len(body) > maximum:
                return Response('{"detail":"Request body too large"}',413,media_type="application/json")
        request._body = bytes(body)
    response = await call_next(request)
    if request.url.path == "/" or request.url.path.startswith("/static/"):
        response.headers["Cache-Control"] = "no-cache"
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["Referrer-Policy"] = "same-origin"
    return response


api_key_header = APIKeyHeader(name="X-API-Key", auto_error=False)


def owner(request: Request, token: str | None = Depends(api_key_header)):
    if token:
        user = authenticate(token)
    else:
        session = request.cookies.get("sea_session")
        with connection() as con:
            row = con.execute("SELECT owner FROM sessions WHERE hash=? AND expires>?", (hashed(session or ""),time.time())).fetchone()
        user = row["owner"] if row else None
        if user and request.method not in ("GET","HEAD"):
            origin = request.headers.get("origin")
            if origin and origin != str(request.base_url).rstrip("/"):
                raise HTTPException(403,"Origin mismatch")
            if request.headers.get("X-SEA-Request") != "1":
                raise HTTPException(403,"Missing request header")
    if not user:
        raise HTTPException(401,"An API key or an active session is required.")
    return user


class Login(BaseModel):
    api_key: str


@app.post("/api/v1/session", operation_id="create_session")
def login(body: Login, response: Response):
    user = authenticate(body.api_key)
    if not user:
        raise HTTPException(401,"Invalid API key.")
    token = secrets.token_urlsafe(32)
    with connection() as con:
        con.execute("DELETE FROM sessions WHERE expires<?",(time.time(),))
        con.execute("INSERT INTO sessions VALUES (?,?,?)",(hashed(token),user,time.time()+28800))
    response.set_cookie("sea_session",token,httponly=True,samesite="strict",secure=os.getenv("SEA_SECURE_COOKIE","0")=="1",max_age=28800)
    return {"authenticated":True}


@app.delete("/api/v1/session",operation_id="delete_session")
def logout(request:Request,response:Response,user=Depends(owner)):
    with connection() as con:
        con.execute("DELETE FROM sessions WHERE hash=?",(hashed(request.cookies.get("sea_session","")),))
    response.delete_cookie("sea_session")
    return {"authenticated":False}


class Prediction(BaseModel):
    dataset_id: str
    queries: list[list[str]] = Field(min_length=1, description="Each inner list is one ligand ensemble. Use singleton lists for separate compounds.")
    top_k: int = Field(default=20,ge=1,le=1000)
    statistics_profile: Literal["legacy", "balanced", "research"] = "legacy"


class Comparison(BaseModel):
    dataset_id: str
    left: list[str] = Field(min_length=1)
    right: list[str] = Field(min_length=1)
    statistics_profile: Literal["legacy", "balanced", "research"] = "legacy"


class StructureImport(BaseModel):
    format: StructureFormat
    content: str = Field(min_length=1, description="UTF-8 SDF, MOL (V2000/V3000), SMILES lines, or CSV text with a smiles column.")


class ImportedMolecule(BaseModel):
    record: int
    name: str
    smiles: str


class StructureImportResult(BaseModel):
    format: StructureFormat
    count: int
    molecules: list[ImportedMolecule]


@app.post("/api/v1/structures/import", response_model=StructureImportResult, operation_id="import_structures")
def structures_import(body: StructureImport, user=Depends(owner)):
    """Convert all records to isomeric SMILES; invalid records reject the whole import."""
    try:
        return import_structures(body.content, body.format, int(os.getenv("SEA_MAX_SET_MOLECULES", 1000)))
    except StructureImportError as exc:
        raise HTTPException(422, detail=exc.args[0]) from exc


class JobAccepted(BaseModel):
    job_id: str


class JobStatus(BaseModel):
    id: str
    kind: str
    status: Literal["queued", "running", "completed", "failed", "cancelled"]
    progress: float
    error: str | None
    created: float
    updated: float
    parent_id: str | None = None
    refinement_job_id: str | None = None
    phase: str | None = None
    profile: str = "legacy"


class PairEvidence(BaseModel):
    query_smiles: str
    reference_smiles: str
    tanimoto: float


class TargetInfo(BaseModel):
    id: str
    name: str
    family: str | None = None
    families: list[str] = Field(default_factory=list)
    protein_classes: list[list[str]] = Field(default_factory=list)
    organism: str | None = None
    tax_id: int | None = None


class TargetSummary(TargetInfo):
    ligand_count: int


class PrecisionProfile(BaseModel):
    id: str
    name: Literal["balanced", "research"]
    version: str
    planned_trials: int
    dataset_id: str
    base_model_id: str
    null_definition: str
    sampler_version: str
    stream_policy: str
    seed: int
    testing_family: Literal["per_query"]
    confidence_level: float
    interval_method: str


class StatisticalDetails(BaseModel):
    phase: Literal["pending", "screening", "final"]
    status: Literal["pending", "resolution_limited", "uncertain", "evaluated", "exact", "analytic"]
    profile_id: str
    planned_trials: int
    trials: int
    exceedances: int | None
    tail_probability_interval: list[float] | None
    e_tail_interval: list[float] | None
    confidence_level: float
    interval_method: str | None
    decision_threshold_e: float


class AnalysisState(BaseModel):
    phase: Literal["evidence", "screening", "refining", "final"]
    complete: bool
    profile: PrecisionProfile
    source_job_id: str
    refinement_job_id: str | None
    total_comparisons: int
    evaluated_comparisons: int


class ScoreResult(BaseModel):
    query_index: int | None
    target: TargetInfo | None = None
    raw_score: float
    contributing_pairs: int
    left_size: int
    right_size: int
    similar_ligands: list[PairEvidence]
    z_score: float | None
    p_value: float | None
    log_p_value: float | None
    e_value: float | None
    log_e_value: float | None
    comparisons: int
    method: Literal["gumbel", "monte_carlo", "exact_zero", "pending"]
    p_resolution: float | None
    statistics: StatisticalDetails | None = None


class ResultBundle(BaseModel):
    dataset: dict
    threshold: float
    normalization: dict | list
    results: list[ScoreResult]
    execution: dict
    target_annotations: dict | None = None
    analysis: AnalysisState | None = None


def validate_sets(dataset_id, sets, maximum):
    if sum(map(len,sets)) > maximum:
        raise HTTPException(422,detail={"code":"too_many_molecules","limit":maximum})
    try:
        data, structures = load_dataset(dataset_id)
    except (ValueError,FileNotFoundError) as exc:
        raise HTTPException(422,str(exc))
    for i,s in enumerate(sets):
        try:
            mols,_ = normalize(s, data['manifest']['fingerprint'])
        except InputError as exc:
            raise HTTPException(422,detail={"code":"invalid_smiles","set_index":i,"rows":exc.args[0]})
        if len(mols)>len(structures):
            raise HTTPException(422,"A ligand set exceeds the background molecule count.")


@app.post("/api/v1/predictions",status_code=202,response_model=JobAccepted,operation_id="submit_prediction")
def prediction(body:Prediction,user=Depends(owner)):
    validate_sets(body.dataset_id,body.queries,int(os.getenv("SEA_MAX_QUERY_MOLECULES",100)))
    return {"job_id":submit(user,"prediction",body.model_dump())}


@app.post("/api/v1/comparisons",status_code=202,response_model=JobAccepted,operation_id="submit_comparison")
def comparison(body:Comparison,user=Depends(owner)):
    limit = int(os.getenv("SEA_MAX_SET_MOLECULES",1000))
    if max(len(body.left),len(body.right))>limit:
        raise HTTPException(422,"A ligand set exceeds the molecule limit.")
    validate_sets(body.dataset_id,[body.left,body.right],2*limit)
    return {"job_id":submit(user,"comparison",body.model_dump())}


def get_job(jid,user):
    with connection() as con:
        row = con.execute("SELECT * FROM jobs WHERE id=? AND owner=?",(jid,user)).fetchone()
    if not row:
        raise HTTPException(404,"Job not found.")
    return dict(row)


@app.get("/api/v1/jobs",operation_id="list_jobs")
def jobs(user=Depends(owner)):
    with connection() as con:
        return [dict(r) for r in con.execute("SELECT id,kind,status,progress,error,created,updated,parent_id,refinement_job_id,phase,profile FROM jobs WHERE owner=? ORDER BY created DESC LIMIT 100",(user,))]


@app.get("/api/v1/jobs/{job_id}",response_model=JobStatus,operation_id="get_job_status")
def job_status(job_id:str,user=Depends(owner)):
    row = get_job(job_id,user)
    return {k:row[k] for k in ("id","kind","status","progress","error","created","updated","parent_id","refinement_job_id","phase","profile")}


@app.get("/api/v1/jobs/{job_id}/preview", response_model=ResultBundle, operation_id="get_job_preview")
def preview(job_id:str,user=Depends(owner)):
    row = get_job(job_id,user)
    snapshot = row["preview"] or row["result"]
    if not snapshot:
        raise HTTPException(409, {"status": row["status"], "error": row["error"]})
    return enrich_results(json.loads(snapshot))


class RefinementRequest(BaseModel):
    statistics_profile: Literal["balanced", "research"] = "balanced"


@app.post("/api/v1/jobs/{job_id}/refine", status_code=202, response_model=JobAccepted, operation_id="refine_job")
def refine(job_id:str,body:RefinementRequest,user=Depends(owner)):
    row = get_job(job_id,user)
    source_id = row["parent_id"] or job_id
    from .phased_jobs import create_refinement
    try:
        with connection() as con:
            con.execute("BEGIN IMMEDIATE")
            child = create_refinement(source_id, user, body.statistics_profile, con)
    except ValueError as exc:
        raise HTTPException(409,str(exc)) from exc
    return {"job_id": child}


@app.post("/api/v1/jobs/{job_id}/resume", response_model=JobStatus, operation_id="resume_refinement")
def resume(job_id:str,user=Depends(owner)):
    get_job(job_id,user)
    with connection() as con:
        con.execute("BEGIN IMMEDIATE")
        row = con.execute("SELECT kind,status FROM jobs WHERE id=?",(job_id,)).fetchone()
        if row["kind"] != "refinement" or row["status"] not in ("failed","cancelled"):
            raise HTTPException(409,"Only failed or cancelled refinement jobs can be resumed.")
        con.execute("UPDATE jobs SET status='queued',cancel=0,error=NULL,worker=NULL,updated=? WHERE id=?",(time.time(),job_id))
    return job_status(job_id,user)


@app.post("/api/v1/jobs/{job_id}/cancel",operation_id="cancel_job")
def cancel(job_id:str,user=Depends(owner)):
    get_job(job_id,user)
    with connection() as con:
        con.execute("UPDATE jobs SET cancel=1,status=CASE WHEN status='queued' THEN 'cancelled' ELSE status END WHERE id=?",(job_id,))
    return job_status(job_id,user)


@app.get("/api/v1/jobs/{job_id}/result",response_model=ResultBundle,operation_id="get_job_result")
def result(job_id:str,format:Literal["json","csv"]="json",user=Depends(owner)):
    row = get_job(job_id,user)
    if row["status"]!="completed":
        raise HTTPException(409,{"status":row["status"],"error":row["error"]})
    data = enrich_results(json.loads(row["result"]))
    if format=="json":
        return data
    out = io.StringIO()
    fields = ["query_index","target_id","target_name","target_family","target_families","target_organism","target_tax_id","target_protein_classes","raw_score","z_score","p_value","log_p_value","e_value","log_e_value","method","p_resolution","comparisons","dataset_id","model_id",
              "analysis_phase","analysis_complete","statistics_profile_id","statistics_profile","statistics_version","sampler_version","stream_policy","null_definition","seed","planned_trials","statistics_phase","trials","exceedances","statistics_status","p_tail_lower","p_tail_upper","e_tail_lower","e_tail_upper","confidence_level","interval_method","source_job_id","refinement_job_id","testing_family"]
    writer = csv.DictWriter(out,fieldnames=fields,extrasaction="ignore")
    writer.writeheader()
    for r in data["results"]:
        target = r.get("target") or {}
        record = {**r,"target_id":target.get("id"),"target_name":target.get("name"),"dataset_id":data["dataset"]["id"],"model_id":data["dataset"]["model_id"]}
        record.update({"target_family":target.get("family"),"target_families":"; ".join(target.get("families",[])),
                       "target_organism":target.get("organism"),"target_tax_id":target.get("tax_id"),
                       "target_protein_classes":json.dumps(target.get("protein_classes",[]),ensure_ascii=False)})
        analysis = data.get("analysis") or {}
        effective = analysis.get("profile") or {}
        stat = r.get("statistics") or {}
        p_interval = stat.get("tail_probability_interval") or [None, None]
        e_interval = stat.get("e_tail_interval") or [None, None]
        record.update(analysis_phase=analysis.get("phase"), analysis_complete=analysis.get("complete"),
                      statistics_profile_id=effective.get("id"), statistics_profile=effective.get("name"),
                      statistics_version=effective.get("version"), sampler_version=effective.get("sampler_version"),
                      stream_policy=effective.get("stream_policy"), null_definition=effective.get("null_definition"), seed=effective.get("seed"),
                      statistics_phase=stat.get("phase"),
                      planned_trials=stat.get("planned_trials"), trials=stat.get("trials"), exceedances=stat.get("exceedances"),
                      statistics_status=stat.get("status"), p_tail_lower=p_interval[0], p_tail_upper=p_interval[1],
                      e_tail_lower=e_interval[0], e_tail_upper=e_interval[1], confidence_level=stat.get("confidence_level"),
                      interval_method=stat.get("interval_method"), source_job_id=analysis.get("source_job_id"),
                      refinement_job_id=analysis.get("refinement_job_id"), testing_family=effective.get("testing_family"))
        # Prevent spreadsheet formulas in user/data-supplied strings.
        writer.writerow({k:("'"+v if isinstance(v,str) and v.startswith(("=","+","-","@")) else v) for k,v in record.items()})
    return Response(out.getvalue(),media_type="text/csv",headers={"Content-Disposition":f'attachment; filename="{job_id}.csv"'})


@app.get("/api/v1/datasets",operation_id="list_datasets")
def datasets(user=Depends(owner)):
    return list_datasets()


@app.get("/api/v1/targets",response_model=list[TargetSummary],operation_id="list_targets")
def targets(dataset_id:str,q:str="",offset:int=0,limit:int=50,user=Depends(owner)):
    if offset<0 or not 1<=limit<=1000:
        raise HTTPException(422,"Invalid pagination")
    try:
        data,_ = load_dataset(dataset_id)
    except (ValueError,FileNotFoundError) as exc:
        raise HTTPException(404,str(exc))
    annotations = load_annotations(dataset_id)
    enriched = [{**target_info(t, annotations),"ligand_count":len(t["smiles"])} for t in data["targets"]]
    return [t for t in enriched if q.casefold() in " ".join(str(t.get(k) or "") for k in ("id","name","family","organism","tax_id")).casefold()][offset:offset+limit]


@app.get("/api/v1/structure",operation_id="render_structure")
def structure(smiles:str,standardize:bool=True,user=Depends(owner)):
    try:
        if standardize:
            mol = Chem.MolFromSmiles(molecule(smiles).smiles)
        else:
            if not smiles.strip() or len(smiles) > 10000:
                raise ValueError("Invalid SMILES")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or mol.GetNumAtoms() == 0:
                raise ValueError("Invalid SMILES")
    except Exception:
        raise HTTPException(422,"Invalid SMILES")
    drawer = rdMolDraw2D.MolDraw2DSVG(300,180)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return Response(drawer.GetDrawingText(),media_type="image/svg+xml")


@app.get("/health",operation_id="health")
def health():
    return {"status":"ok"}


@app.get("/",include_in_schema=False)
def index():
    return FileResponse(WEB/"index.html")
