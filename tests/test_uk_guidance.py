"""The UK rules layer: UK Chemotherapy Board DPYD and CERSI-PGx UGT1A1.

These tests assert the dose numbers the published UK documents state, so
a later edit to data/uk_guidance.json that changes a dose, or quietly
extends a rule to a genotype the source does not table, fails here.
"""

from __future__ import annotations

import pytest

from pharminform import regimens, uk_guidance

# --- DPYD (UK Chemotherapy Board 2020) ---------------------------------


def test_no_actionable_variant_is_full_dose_with_the_documents_caveat():
    recs = uk_guidance.dpyd_recommendations("*1/*1")
    assert {r.drug for r in recs} == {"fluorouracil", "capecitabine", "tegafur"}
    rec = next(r for r in recs if r.drug == "capecitabine")
    assert rec.start_dose_percent == 100
    assert rec.titration is None
    assert any("does not exclude DPD deficiency" in c for c in rec.caveats)


@pytest.mark.parametrize("diplotype", ["*2A/*1", "*13/*1"])
def test_heterozygous_no_function_starts_at_50_and_escalates_to_75(diplotype):
    rec = uk_guidance.dpyd_recommendations(diplotype)[0]
    assert rec.dpd_activity == "50%"
    assert rec.start_dose_percent == 50
    assert rec.max_dose_percent == 75
    # The 85% ceiling belongs to the decreased-function group only.
    assert rec.max_dose_percent_with_caution is None
    assert "12.5%" in rec.titration and "two increments" in rec.titration


@pytest.mark.parametrize("diplotype", ["c.2846A>T/*1", "c.1236G>A(HapB3)/*1"])
def test_heterozygous_decreased_function_may_reach_85_percent_with_caution(diplotype):
    rec = uk_guidance.dpyd_recommendations(diplotype)[0]
    assert rec.dpd_activity == "50-75%"
    assert rec.start_dose_percent == 50
    assert rec.max_dose_percent == 75
    assert rec.max_dose_percent_with_caution == 85


@pytest.mark.parametrize(
    "diplotype",
    ["*2A/*2A", "*2A/*13", "*2A/c.2846A>T", "c.2846A>T/c.1236G>A(HapB3)"],
)
def test_two_variant_alleles_give_no_dose_and_name_alternatives(diplotype):
    rec = uk_guidance.dpyd_recommendations(diplotype)[0]
    # Additive loss of activity: a compound genotype must not fall through
    # to a heterozygous rule and be dosed at 50%.
    assert rec.start_dose_percent is None
    assert rec.alternatives == ["raltitrexed", "trifluridine/tipiracil"]
    assert "avoided" in rec.advice


def test_allele_outside_the_uk_tables_gets_no_uk_advice():
    # *9A is normal function in CPIC and is not tabled by the UK document;
    # inventing a dose for it would be an extrapolation.
    assert uk_guidance.dpyd_recommendations("*9A/*1") == []
    assert uk_guidance.classify_dpyd("*9A/*1") is None


def test_dpyd_recommendation_carries_its_source_and_version():
    rec = uk_guidance.dpyd_recommendations("*2A/*1")[0]
    assert rec.source.issuer == "UK Chemotherapy Board"
    assert rec.source.version == "2020"
    assert rec.source.label == "UK Chemotherapy Board 2020"


# --- UGT1A1 (CERSI-PGx 2026) -------------------------------------------


def test_poor_metaboliser_gets_a_30_percent_cycle_one_reduction():
    rec = uk_guidance.ugt1a1_recommendation("poor metabolizer", "*28/*28")
    assert rec.start_dose_percent == 70
    assert "30%" in rec.advice and "cycle 1" in rec.advice
    assert rec.snomed_code == "738536002"
    assert rec.source.doi == "10.1002/bcp.70659"
    assert any("epithelial" in c or "non-pegylated" in c for c in rec.caveats)


def test_intermediate_metaboliser_gets_no_dose_reduction():
    rec = uk_guidance.ugt1a1_recommendation("intermediate metabolizer", "*28/*1")
    assert rec.start_dose_percent == 100
    assert "No dose reduction" in rec.advice


def test_star36_is_flagged_as_dosed_like_wild_type():
    rec = uk_guidance.ugt1a1_recommendation("normal metabolizer", "*36/*1")
    assert any("*36" in c for c in rec.caveats)


def test_indeterminate_phenotype_gets_no_uk_advice():
    assert uk_guidance.ugt1a1_recommendation("indeterminate", "*80/*1") is None


def test_sources_are_listed_with_dates():
    sources = {s.id: s for s in uk_guidance.available_sources()}
    assert set(sources) == {"ukcb_dpyd_2020", "cersi_ugt1a1_2026"}
    assert sources["cersi_ugt1a1_2026"].date == "2026-07-02"


# --- Regimen view ------------------------------------------------------


def test_folfiri_gives_one_decision_per_drug_from_the_right_gene():
    assessment = regimens.assess_regimen(
        "FOLFIRI",
        {"DPYD": "intermediate metabolizer", "UGT1A1": "poor metabolizer"},
        {"DPYD": "*2A/*1", "UGT1A1": "*28/*28"},
    )
    assert assessment.complete
    by_drug = {d.drug: d for d in assessment.decisions}
    assert set(by_drug) == {"fluorouracil", "irinotecan"}
    assert by_drug["fluorouracil"].gene == "DPYD"
    assert by_drug["fluorouracil"].recommendations[0].start_dose_percent == 50
    assert by_drug["irinotecan"].recommendations[0].start_dose_percent == 70
    # Two separate dose decisions, not one blended one.
    assert all(d.status == "advice" for d in assessment.decisions)
    assert "folinic acid" in assessment.unmodified_drugs


def test_combined_regimen_says_to_wait_for_both_results():
    assessment = regimens.assess_regimen(
        "FOLFIRI",
        {"DPYD": "normal metabolizer", "UGT1A1": "normal metabolizer"},
        {"DPYD": "*1/*1", "UGT1A1": "*1/*1"},
    )
    assert any("wait for both" in n.lower() for n in assessment.notes)


def test_a_missing_gene_is_reported_not_assumed_normal():
    assessment = regimens.assess_regimen(
        "FOLFIRI", {"DPYD": "normal metabolizer"}, {"DPYD": "*1/*1"}
    )
    assert not assessment.complete
    assert assessment.missing_genes == ["UGT1A1"]
    irinotecan = next(d for d in assessment.decisions if d.drug == "irinotecan")
    assert irinotecan.status == "missing result"
    assert irinotecan.recommendations == []
    assert "not a normal-metaboliser result" in irinotecan.note
    assert any(n.startswith("Incomplete") for n in assessment.notes)


def test_fluoropyrimidine_only_regimen_has_no_irinotecan_decision():
    assessment = regimens.assess_regimen("CAPOX", {"DPYD": "normal metabolizer"}, {"DPYD": "*1/*1"})
    assert [d.drug for d in assessment.decisions] == ["capecitabine"]
    assert assessment.complete
    assert "oxaliplatin" in assessment.unmodified_drugs
    assert not any("wait for both" in n.lower() for n in assessment.notes)


def test_assess_all_skips_regimens_no_called_gene_informs():
    assessments = regimens.assess_all_regimens(
        {"UGT1A1": "poor metabolizer"}, {"UGT1A1": "*28/*28"}
    )
    names = {a.regimen for a in assessments}
    assert "IRINOTECAN" in names
    assert "CAPOX" not in names  # DPYD-only regimen, DPYD not called
    assert "FOLFIRI" in names and not next(a for a in assessments if a.regimen == "FOLFIRI").complete


def test_unknown_regimen_raises():
    with pytest.raises(KeyError):
        regimens.assess_regimen("FOLFWHAT", {}, {})


# --- Rendering ---------------------------------------------------------


def test_html_report_shows_uk_guidance_and_the_regimen_view():
    from pharminform.report import GeneReportRow, ReportContext, render_html

    phenotypes = {"DPYD": "intermediate metabolizer", "UGT1A1": "poor metabolizer"}
    diplotypes = {"DPYD": "*2A/*1", "UGT1A1": "*28/*28"}
    context = ReportContext(
        gene_results=[
            GeneReportRow(g, diplotypes[g], phenotypes[g], "vcf", "high") for g in diplotypes
        ],
        cpic_recommendations=[],
        uk_recommendations=[r.to_dict() for r in uk_guidance.uk_recommendations(phenotypes, diplotypes)],
        regimens=[a.to_dict() for a in regimens.assess_all_regimens(phenotypes, diplotypes)],
        composite_results=[],
    )
    html = render_html(context)
    assert "UK Chemotherapy Board" in html
    assert "CERSI-PGx" in html
    assert "50%" in html and "70%" in html
    assert "wait for both" in html
    assert "FOLFIRI" in html
    assert "raltitrexed" not in html  # heterozygous: no alternative agent is named
