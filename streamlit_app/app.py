"""PharmInform Streamlit demo.

Run with: streamlit run streamlit_app/app.py
"""

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from pharminform.caller import VCFIngestError
from pharminform.pipeline import run_pipeline
from pharminform.report import PDFExportUnavailable, render_pdf

st.set_page_config(page_title="PharmInform", page_icon="\U0001f9ec", layout="wide")

st.title("PharmInform")
st.caption(
    "Clinical pharmacogenomics demo: star allele calling, CPIC dosing "
    "recommendations, and an experimental multi-gene composite risk score."
)

with st.sidebar:
    st.header("Input")
    uploaded = st.file_uploader("VCF file", type=["vcf"])
    patient_id = st.text_input("Patient/sample label (optional)")
    sample = st.text_input("Sample column (needed for multi-sample VCFs)")
    min_depth = st.slider("Minimum read depth", min_value=1, max_value=100, value=20)
    run_clicked = st.button("Run PharmInform", type="primary", disabled=uploaded is None)

st.info(
    "Research/educational demo only. Star allele calls from a VCF alone use "
    "a simplified rule-based matcher with no copy-number detection -- see "
    "the Confidence column below and caller.py's docstring. Composite risk "
    "scores are experimental, not a CPIC/DPWG guideline output.",
    icon="⚠️",
)

if run_clicked and uploaded is not None:
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmp:
        tmp.write(uploaded.getvalue())
        tmp_path = tmp.name

    try:
        with st.spinner("Calling star alleles and scoring..."):
            context = run_pipeline(
                tmp_path, min_depth=min_depth, patient_id=patient_id or None, sample=sample or None
            )
    except VCFIngestError as exc:
        st.error(str(exc))
        st.stop()
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    st.header("Star allele calls and phenotypes")
    gene_df = pd.DataFrame(
        [
            {
                "Gene": g.gene,
                "Diplotype": g.diplotype,
                "Phenotype": g.phenotype,
                "Method": g.calling_method,
                "Confidence": g.confidence_level,
                "Reasons": "; ".join(g.confidence_reasons) if g.confidence_reasons else "",
            }
            for g in context.gene_results
        ]
    )
    st.dataframe(gene_df, use_container_width=True)

    st.header("CPIC dosing recommendations")
    if context.cpic_recommendations:
        cpic_df = pd.DataFrame(
            [
                {
                    "Drug": r["drug"],
                    "Phenotype(s)": r["phenotypes"],
                    "Recommendation": r["recommendation"],
                    "Classification": r["classification"],
                    "Citations": ", ".join(f"PMID:{c}" for c in r["citations"]),
                }
                for r in context.cpic_recommendations
            ]
        )
        st.dataframe(cpic_df, use_container_width=True)
    else:
        st.write("No matching CPIC recommendations found for the called phenotypes.")

    st.header("Composite multi-gene risk (experimental)")
    st.warning(
        "Exploratory composite score. Not a CPIC/DPWG guideline output. For "
        "research/educational use only -- clinical decisions should rely on "
        "individual gene-level CPIC/DPWG recommendations above.",
        icon="⚠️",
    )
    if context.composite_results:
        composite_df = pd.DataFrame(
            [
                {
                    "Drug": c["drug"],
                    "Composite score": round(c["composite_function_score"], 3),
                    "Risk tier": c["risk_tier"],
                    "Directional flag": c["directional_flag"] or "--",
                    "Weighting method": c["weighting_method"],
                }
                for c in context.composite_results
            ]
        )
        st.dataframe(composite_df, use_container_width=True)
    else:
        st.write("No composite risk drugs matched the called genes.")

    st.header("Export")
    html = None
    from pharminform.report import render_html

    html = render_html(context)
    st.download_button("Download HTML report", data=html, file_name="pharminform_report.html", mime="text/html")
    try:
        pdf_bytes = render_pdf(html)
        st.download_button("Download PDF report", data=pdf_bytes, file_name="pharminform_report.pdf", mime="application/pdf")
    except PDFExportUnavailable as exc:
        st.caption(f"PDF export unavailable in this environment: {exc}")
elif uploaded is None:
    st.write("Upload a VCF file in the sidebar to get started.")
