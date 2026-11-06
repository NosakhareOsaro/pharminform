"""PharmInform CLI.

    pharminform report --vcf sample.vcf --out report.html [--pdf report.pdf]
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from . import __version__
from .caller import VCFIngestError
from .pipeline import run_pipeline
from .report import PDFExportUnavailable, write_report

app = typer.Typer(help="PharmInform: clinical pharmacogenomics CLI")


@app.command("version")
def show_version() -> None:
    """Print the installed PharmInform version.

    Also has the side effect of keeping Typer from collapsing this into a
    single-command CLI without subcommands (Typer's default behavior with
    only one @app.command() registered is to drop the subcommand name
    entirely), so `pharminform report ...` stays the invocation shape.
    """
    typer.echo(__version__)


@app.command("report")
def make_report(
    vcf: Annotated[Path, typer.Option("--vcf", exists=True, readable=True, help="Input VCF file")],
    out: Annotated[Path, typer.Option("--out", help="Output HTML report path")],
    pdf: Annotated[
        Path | None, typer.Option("--pdf", help="Also write a PDF report to this path")
    ] = None,
    min_depth: Annotated[
        int, typer.Option("--min-depth", help="Minimum read depth to trust a defining variant")
    ] = 20,
    patient_id: Annotated[
        str | None,
        typer.Option("--patient-id", help="Optional patient/sample label for the report"),
    ] = None,
    sample: Annotated[
        str | None,
        typer.Option("--sample", help="Sample column to report on (required for multi-sample VCFs)"),
    ] = None,
) -> None:
    """Call star alleles from a VCF and generate a PGx clinical report."""
    typer.echo(f"Loading and calling variants from {vcf} (min depth {min_depth})...")
    try:
        context = run_pipeline(vcf, min_depth=min_depth, patient_id=patient_id, sample=sample)
    except VCFIngestError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=2) from exc

    for row in context.gene_results:
        flag = " [LOW CONFIDENCE]" if row.confidence_level == "low" else ""
        typer.echo(f"  {row.gene}: {row.diplotype} -> {row.phenotype}{flag}")

    try:
        write_report(context, out, pdf)
    except PDFExportUnavailable as exc:
        typer.echo(f"HTML report written; PDF export failed: {exc}", err=True)
        write_report(context, out, None)
        raise typer.Exit(code=0)

    typer.echo(f"Report written to {out}" + (f" and {pdf}" if pdf else ""))


if __name__ == "__main__":
    app()
