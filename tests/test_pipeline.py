from pathlib import Path

from pharminform.pipeline import run_pipeline

FIXTURE_VCF = Path(__file__).parent / "fixtures" / "synthetic_sample.vcf"


def test_run_pipeline_end_to_end():
    ctx = run_pipeline(FIXTURE_VCF, min_depth=10, patient_id="TEST-1")
    genes_called = {row.gene for row in ctx.gene_results}
    assert genes_called == {"CYP2D6", "CYP2C19", "CYP2C9", "DPYD", "UGT1A1", "TPMT"}

    tpmt_row = next(r for r in ctx.gene_results if r.gene == "TPMT")
    assert tpmt_row.diplotype == "*3C/*3C"
    assert tpmt_row.phenotype == "poor metabolizer"

    # low-depth CYP2D6*4 variant in the fixture should be filtered, and every
    # row should be flagged low-confidence since this is the VCF-only path
    assert all(r.confidence_level == "low" for r in ctx.gene_results)

    # azathioprine/mercaptopurine should show up in composite results since
    # TPMT was called poor metabolizer
    composite_drugs = {c["drug"] for c in ctx.composite_results}
    assert "azathioprine" in composite_drugs
    aza = next(c for c in ctx.composite_results if c["drug"] == "azathioprine")
    assert aza["risk_tier"] == "High"


def test_run_pipeline_includes_cpic_recommendations_when_matched():
    ctx = run_pipeline(FIXTURE_VCF, min_depth=10)
    drugs_with_recs = {r["drug"] for r in ctx.cpic_recommendations}
    assert "azathioprine" in drugs_with_recs
