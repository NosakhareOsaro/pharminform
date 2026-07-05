"""FastAPI endpoint: upload a VCF, get back a JSON PGx report.

Run with: uvicorn pharminform.api:app --reload
"""

from __future__ import annotations

import os
import tempfile

from fastapi import FastAPI, File, HTTPException, Query, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse

from . import __version__
from .pipeline import run_pipeline
from .report import render_html

app = FastAPI(
    title="PharmInform API",
    version=__version__,
    description="Upload a VCF, get back star allele calls, phenotypes, "
    "CPIC recommendations, and experimental composite multi-gene risk scores.",
)


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "version": __version__}


@app.post("/report")
async def create_report(
    vcf: UploadFile = File(...),
    min_depth: int = Query(20, ge=1),
    patient_id: str | None = Query(None),
) -> JSONResponse:
    if not vcf.filename or not vcf.filename.endswith((".vcf", ".vcf.gz")):
        raise HTTPException(status_code=400, detail="Expected a .vcf or .vcf.gz file")

    suffix = ".vcf.gz" if vcf.filename.endswith(".gz") else ".vcf"
    tmp_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            tmp.write(await vcf.read())
            tmp_path = tmp.name

        context = run_pipeline(tmp_path, min_depth=min_depth, patient_id=patient_id)
    except Exception as exc:  # surfaced as a 422 rather than a bare 500 traceback
        raise HTTPException(status_code=422, detail=f"Could not process VCF: {exc}") from exc
    finally:
        if tmp_path is not None:
            os.unlink(tmp_path)

    return JSONResponse(content=context.to_dict())


@app.post("/report/html", response_class=HTMLResponse)
async def create_report_html(
    vcf: UploadFile = File(...),
    min_depth: int = Query(20, ge=1),
    patient_id: str | None = Query(None),
) -> HTMLResponse:
    if not vcf.filename or not vcf.filename.endswith((".vcf", ".vcf.gz")):
        raise HTTPException(status_code=400, detail="Expected a .vcf or .vcf.gz file")

    suffix = ".vcf.gz" if vcf.filename.endswith(".gz") else ".vcf"
    tmp_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            tmp.write(await vcf.read())
            tmp_path = tmp.name

        context = run_pipeline(tmp_path, min_depth=min_depth, patient_id=patient_id)
    except Exception as exc:
        raise HTTPException(status_code=422, detail=f"Could not process VCF: {exc}") from exc
    finally:
        if tmp_path is not None:
            os.unlink(tmp_path)

    return HTMLResponse(content=render_html(context))
