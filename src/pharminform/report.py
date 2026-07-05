"""Jinja2 HTML clinical report + PDF export.

PDF export via weasyprint requires native Pango/GObject/cairo libraries
that are not bundled with the Python package. On this project's
development machine (macOS + Homebrew), plain `pip install weasyprint` and
even `brew install pango` were not enough on their own -- weasyprint
raised ``OSError: cannot load library 'libgobject-2.0-0'`` until
``DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib`` was set, because Homebrew
on Apple Silicon installs to ``/opt/homebrew`` rather than a default dyld
search path. ``render_pdf`` below catches that failure mode and raises a
clear, actionable error rather than a confusing traceback; it never
silently skips PDF generation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from . import __version__

# templates/ lives at the repo root (per the project layout), not inside
# the installed package -- fine for editable/dev installs and for running
# from a source checkout. A real PyPI wheel would need this bundled as
# package data (e.g. copied under src/pharminform/templates and referenced
# via importlib.resources) instead; tracked as a packaging TODO rather than
# silently broken.
TEMPLATE_DIR = Path(__file__).parent.parent.parent / "templates"


class PDFExportUnavailable(RuntimeError):
    """Raised when weasyprint (or its native dependencies) can't be loaded."""


@dataclass
class GeneReportRow:
    gene: str
    diplotype: str
    phenotype: str | None
    calling_method: str
    confidence_level: str
    confidence_reasons: list[str] = field(default_factory=list)


@dataclass
class ReportContext:
    gene_results: list[GeneReportRow]
    cpic_recommendations: list[dict]
    composite_results: list[dict]
    patient_id: str | None = None
    generated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    version: str = __version__

    def to_dict(self) -> dict:
        return {
            "patient_id": self.patient_id,
            "generated_at": self.generated_at,
            "version": self.version,
            "gene_results": [
                {
                    "gene": g.gene,
                    "diplotype": g.diplotype,
                    "phenotype": g.phenotype,
                    "calling_method": g.calling_method,
                    "confidence_level": g.confidence_level,
                    "confidence_reasons": g.confidence_reasons,
                }
                for g in self.gene_results
            ],
            "cpic_recommendations": self.cpic_recommendations,
            "composite_risk": self.composite_results,
        }


def _env() -> Environment:
    return Environment(
        loader=FileSystemLoader(str(TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
    )


def render_html(context: ReportContext, template_name: str = "report_template.html") -> str:
    template = _env().get_template(template_name)
    return template.render(
        patient_id=context.patient_id,
        generated_at=context.generated_at,
        version=context.version,
        gene_results=context.gene_results,
        cpic_recommendations=context.cpic_recommendations,
        composite_results=context.composite_results,
    )


def render_pdf(html: str) -> bytes:
    """Render an HTML report string to PDF bytes via weasyprint.

    Raises PDFExportUnavailable with actionable instructions if weasyprint
    or its native libraries can't be loaded, instead of letting a raw
    OSError/ImportError propagate.
    """
    try:
        import weasyprint
    except (ImportError, OSError) as exc:
        raise PDFExportUnavailable(
            "weasyprint could not be loaded (native Pango/GObject/cairo "
            "libraries missing or not on the dynamic linker search path). "
            "On macOS with Homebrew: `brew install pango` and then run with "
            "`DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib` set (Apple "
            "Silicon) or `/usr/local/lib` (Intel). On Linux: install "
            "libpango-1.0-0, libpangoft2-1.0-0, and libgdk-pixbuf2.0-0. See "
            "https://doc.courtbouillon.org/weasyprint/stable/first_steps.html#installation"
        ) from exc
    return weasyprint.HTML(string=html).write_pdf()


def write_report(
    context: ReportContext,
    html_out: str | Path,
    pdf_out: str | Path | None = None,
) -> None:
    html = render_html(context)
    Path(html_out).write_text(html)
    if pdf_out is not None:
        pdf_bytes = render_pdf(html)
        Path(pdf_out).write_bytes(pdf_bytes)
