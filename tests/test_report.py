import os

import pytest

from pharminform.report import (
    GeneReportRow,
    PDFExportUnavailable,
    ReportContext,
    render_html,
    render_pdf,
    write_report,
)


def make_context():
    return ReportContext(
        patient_id="SAMPLE-001",
        gene_results=[
            GeneReportRow(
                gene="CYP2D6",
                diplotype="*4/*4",
                phenotype="poor metabolizer",
                calling_method="rule_based_vcf",
                confidence_level="low",
                confidence_reasons=["VCF-only rule-based call: no copy-number or phasing information available"],
            ),
        ],
        cpic_recommendations=[
            {
                "drug": "codeine",
                "phenotypes": {"CYP2D6": "Poor Metabolizer"},
                "recommendation": "Avoid codeine; use non-CYP2D6-dependent analgesic.",
                "classification": "Strong",
                "citations": ["22205192"],
            }
        ],
        composite_results=[
            {
                "drug": "amitriptyline",
                "composite_function_score": 0.25,
                "risk_tier": "High",
                "directional_flag": None,
                "weighting_method": "equal-weight fallback (no published fraction for at least one contributing gene)",
                "contributing_genes": [
                    {"gene": "CYP2D6", "phenotype": "poor metabolizer", "weight": 0.5},
                    {"gene": "CYP2C19", "phenotype": "intermediate metabolizer", "weight": 0.5},
                ],
            }
        ],
    )


def test_render_html_contains_key_sections():
    html = render_html(make_context())
    assert "PharmInform Pharmacogenomic Report" in html
    assert "CYP2D6" in html
    assert "*4/*4" in html
    assert "codeine" in html
    assert "amitriptyline" in html
    assert "Exploratory" in html or "experimental" in html.lower()
    assert "SAMPLE-001" in html


def test_render_html_flags_low_confidence_row():
    html = render_html(make_context())
    assert "confidence-low" in html


def test_write_report_html_only(tmp_path):
    out = tmp_path / "report.html"
    write_report(make_context(), out)
    assert out.exists()
    assert "PharmInform" in out.read_text()


def test_render_pdf_works_or_fails_clearly():
    html = render_html(make_context())
    # This environment needs DYLD_FALLBACK_LIBRARY_PATH set for weasyprint's
    # native deps to load (see report.py docstring) -- assert either a real
    # PDF comes back, or a clear, actionable PDFExportUnavailable is raised.
    try:
        pdf_bytes = render_pdf(html)
        assert pdf_bytes[:5] == b"%PDF-"
    except PDFExportUnavailable as exc:
        assert "brew install pango" in str(exc) or "libpango" in str(exc)
