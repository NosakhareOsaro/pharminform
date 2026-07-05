from datetime import date

from pharminform.omop_mapping import (
    records_to_rows,
    to_drug_exposure,
    to_pgx_annotation,
)


def test_to_drug_exposure_basic_fields():
    rec = to_drug_exposure(
        person_id=1,
        drug_name="codeine",
        drug_rxnorm_id="RxNorm:2670",
        recommendation_text="Avoid codeine use because of potential for serious toxicity.",
        exposure_start=date(2026, 1, 1),
        exposure_end=date(2026, 1, 30),
    )
    assert rec.person_id == 1
    assert rec.drug_concept_id == 0  # unresolved without OMOP vocabulary tables
    assert rec.drug_source_concept_code == "RxNorm:2670"
    assert "Avoid" in rec.sig


def test_to_pgx_annotation_carries_citations():
    ann = to_pgx_annotation(
        person_id=1,
        gene="CYP2D6",
        diplotype="*1x2/*1x2",
        phenotype="ultrarapid metabolizer",
        calling_method="aldy",
        confidence="high",
        cpic_recommendation="Avoid codeine",
        cpic_classification="Strong",
        linked_drug="codeine",
        citations=["22205192"],
    )
    assert ann.gene_symbol == "CYP2D6"
    assert ann.citations == ["22205192"]
    assert ann.linked_drug_exposure == "codeine"


def test_records_to_rows_flattens_to_plain_dicts():
    rec = to_drug_exposure(
        person_id=2,
        drug_name="warfarin",
        drug_rxnorm_id="RxNorm:11289",
        recommendation_text="Use a validated pharmacogenetic dosing algorithm.",
        exposure_start=date(2026, 2, 1),
        exposure_end=date(2026, 2, 28),
    )
    rows = records_to_rows([rec])
    assert rows[0]["drug_source_value"] == "warfarin"
    assert rows[0]["drug_exposure_start_date"] == "2026-02-01"
