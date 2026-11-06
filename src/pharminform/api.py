"""FastAPI endpoint: upload a VCF, get back a JSON PGx report.

Run with: uvicorn pharminform.api:app --reload
"""

from __future__ import annotations

import os
import tempfile
from typing import Annotated

from fastapi import FastAPI, File, HTTPException, Query, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse

from . import __version__
from .fhir_export import report_to_bundle
from .pipeline import run_pipeline
from .report import ReportContext, render_html

app = FastAPI(
    title="PharmInform API",
    version=__version__,
    description="Upload a VCF, get back star allele calls, phenotypes, "
    "CPIC recommendations, and experimental composite multi-gene risk scores.",
)


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "version": __version__}


async def _context_from_upload(
    vcf: UploadFile, min_depth: int, patient_id: str | None, sample: str | None
) -> ReportContext:
    if not vcf.filename or not vcf.filename.endswith((".vcf", ".vcf.gz")):
        raise HTTPException(status_code=400, detail="Expected a .vcf or .vcf.gz file")

    suffix = ".vcf.gz" if vcf.filename.endswith(".gz") else ".vcf"
    tmp_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            tmp.write(await vcf.read())
            tmp_path = tmp.name

        return run_pipeline(tmp_path, min_depth=min_depth, patient_id=patient_id, sample=sample)
    except HTTPException:
        raise
    except Exception as exc:  # surfaced as a 422 rather than a bare 500 traceback
        raise HTTPException(status_code=422, detail=f"Could not process VCF: {exc}") from exc
    finally:
        if tmp_path is not None:
            os.unlink(tmp_path)


@app.post("/report")
async def create_report(
    vcf: Annotated[UploadFile, File()],
    min_depth: int = Query(20, ge=1),
    patient_id: str | None = Query(None),
    sample: str | None = Query(None, description="Sample column; required for multi-sample VCFs"),
) -> JSONResponse:
    context = await _context_from_upload(vcf, min_depth, patient_id, sample)
    return JSONResponse(content=context.to_dict())


@app.post("/report/fhir")
async def create_report_fhir(
    vcf: Annotated[UploadFile, File()],
    min_depth: int = Query(20, ge=1),
    sample: str | None = Query(None, description="Sample column; required for multi-sample VCFs"),
    patient_reference: str | None = Query(
        None,
        description=(
            "FHIR reference for the subject, e.g. Patient/abc. Omitted, the bundle "
            "carries no patient data."
        ),
    ),
) -> JSONResponse:
    """The report as a FHIR R4 Bundle following the HL7 Genomics Reporting IG."""
    context = await _context_from_upload(vcf, min_depth, None, sample)
    return JSONResponse(content=report_to_bundle(context, patient_reference))


@app.post("/report/html", response_class=HTMLResponse)
async def create_report_html(
    vcf: Annotated[UploadFile, File()],
    min_depth: int = Query(20, ge=1),
    patient_id: str | None = Query(None),
    sample: str | None = Query(None, description="Sample column; required for multi-sample VCFs"),
) -> HTMLResponse:
    context = await _context_from_upload(vcf, min_depth, patient_id, sample)
    return HTMLResponse(content=render_html(context))
