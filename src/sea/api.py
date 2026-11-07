import csv
import io
import json
import os
from pathlib import Path
import secrets
import time
from typing import Literal
from fastapi import FastAPI, Depends, HTTPException, Request, Response
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.security import APIKeyHeader
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from .chemistry import normalize, molecule, InputError
from .data import list_datasets, load_dataset
from .jobs import connection, authenticate, hashed, submit

app = FastAPI(title="SEA Open", version="0.1.0", description="Ligand ensemble similarity. E-values are expected random hits, not binding probabilities.")
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
        raise HTTPException(401,"APIキーまたはログインが必要です")
    return user


class Login(BaseModel):
    api_key: str


@app.post("/api/v1/session", operation_id="create_session")
def login(body: Login, response: Response):
    user = authenticate(body.api_key)
    if not user:
        raise HTTPException(401,"APIキーが不正です")
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


class Comparison(BaseModel):
    dataset_id: str
    left: list[str] = Field(min_length=1)
    right: list[str] = Field(min_length=1)


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


class PairEvidence(BaseModel):
    query_smiles: str
    reference_smiles: str
    tanimoto: float


class ScoreResult(BaseModel):
    query_index: int | None
    target: dict[str, str] | None = None
    raw_score: float
    contributing_pairs: int
    left_size: int
    right_size: int
    similar_ligands: list[PairEvidence]
    z_score: float | None
    p_value: float
    log_p_value: float
    e_value: float
    log_e_value: float
    comparisons: int
    method: Literal["gumbel", "monte_carlo", "exact_zero"]
    p_resolution: float | None


class ResultBundle(BaseModel):
    dataset: dict
    threshold: float
    normalization: dict | list
    results: list[ScoreResult]
    execution: dict


def validate_sets(dataset_id, sets, maximum):
    if sum(map(len,sets)) > maximum:
        raise HTTPException(422,detail={"code":"too_many_molecules","limit":maximum})
    try:
        data, structures = load_dataset(dataset_id)
    except (ValueError,FileNotFoundError) as exc:
        raise HTTPException(422,str(exc))
    for i,s in enumerate(sets):
        try:
            mols,_ = normalize(s)
        except InputError as exc:
            raise HTTPException(422,detail={"code":"invalid_smiles","set_index":i,"rows":exc.args[0]})
        if len(mols)>len(structures):
            raise HTTPException(422,"背景分子数を超える集合です")


@app.post("/api/v1/predictions",status_code=202,response_model=JobAccepted,operation_id="submit_prediction")
def prediction(body:Prediction,user=Depends(owner)):
    validate_sets(body.dataset_id,body.queries,int(os.getenv("SEA_MAX_QUERY_MOLECULES",100)))
    return {"job_id":submit(user,"prediction",body.model_dump())}


@app.post("/api/v1/comparisons",status_code=202,response_model=JobAccepted,operation_id="submit_comparison")
def comparison(body:Comparison,user=Depends(owner)):
    limit = int(os.getenv("SEA_MAX_SET_MOLECULES",1000))
    if max(len(body.left),len(body.right))>limit:
        raise HTTPException(422,"集合サイズ上限を超えています")
    validate_sets(body.dataset_id,[body.left,body.right],2*limit)
    return {"job_id":submit(user,"comparison",body.model_dump())}


def get_job(jid,user):
    with connection() as con:
        row = con.execute("SELECT * FROM jobs WHERE id=? AND owner=?",(jid,user)).fetchone()
    if not row:
        raise HTTPException(404,"ジョブがありません")
    return dict(row)


@app.get("/api/v1/jobs",operation_id="list_jobs")
def jobs(user=Depends(owner)):
    with connection() as con:
        return [dict(r) for r in con.execute("SELECT id,kind,status,progress,error,created,updated FROM jobs WHERE owner=? ORDER BY created DESC LIMIT 100",(user,))]


@app.get("/api/v1/jobs/{job_id}",response_model=JobStatus,operation_id="get_job_status")
def job_status(job_id:str,user=Depends(owner)):
    row = get_job(job_id,user)
    return {k:row[k] for k in ("id","kind","status","progress","error","created","updated")}


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
    data = json.loads(row["result"])
    if format=="json":
        return data
    out = io.StringIO()
    fields = ["query_index","target_id","target_name","raw_score","z_score","p_value","log_p_value","e_value","log_e_value","method","p_resolution","comparisons","dataset_id","model_id"]
    writer = csv.DictWriter(out,fieldnames=fields,extrasaction="ignore")
    writer.writeheader()
    for r in data["results"]:
        record = {**r,"target_id":r.get("target",{}).get("id"),"target_name":r.get("target",{}).get("name"),"dataset_id":data["dataset"]["id"],"model_id":data["dataset"]["model_id"]}
        # Prevent spreadsheet formulas in user/data-supplied strings.
        writer.writerow({k:("'"+v if isinstance(v,str) and v.startswith(("=","+","-","@")) else v) for k,v in record.items()})
    return Response(out.getvalue(),media_type="text/csv",headers={"Content-Disposition":f'attachment; filename="{job_id}.csv"'})


@app.get("/api/v1/datasets",operation_id="list_datasets")
def datasets(user=Depends(owner)):
    return list_datasets()


@app.get("/api/v1/targets",operation_id="list_targets")
def targets(dataset_id:str,q:str="",offset:int=0,limit:int=50,user=Depends(owner)):
    if offset<0 or not 1<=limit<=1000:
        raise HTTPException(422,"Invalid pagination")
    try:
        data,_ = load_dataset(dataset_id)
    except (ValueError,FileNotFoundError) as exc:
        raise HTTPException(404,str(exc))
    return [{"id":t["id"],"name":t["name"],"ligand_count":len(t["smiles"])} for t in data["targets"] if q.lower() in (t["id"]+" "+t["name"]).lower()][offset:offset+limit]


@app.get("/api/v1/structure",operation_id="render_structure")
def structure(smiles:str,user=Depends(owner)):
    try:
        canonical = molecule(smiles).smiles
    except Exception:
        raise HTTPException(422,"Invalid SMILES")
    drawer = rdMolDraw2D.MolDraw2DSVG(300,180)
    drawer.DrawMolecule(Chem.MolFromSmiles(canonical))
    drawer.FinishDrawing()
    return Response(drawer.GetDrawingText(),media_type="image/svg+xml")


@app.get("/health",operation_id="health")
def health():
    return {"status":"ok"}


@app.get("/",include_in_schema=False)
def index():
    return FileResponse(WEB/"index.html")
