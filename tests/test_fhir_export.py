"""FHIR R4 output against the HL7 Genomics Reporting IG profiles.

These tests pin the profile URLs and the codes the IG requires, so a
refactor cannot quietly emit a bundle a conformant consumer would
reject, and cannot invent a SNOMED code for a phenotype this project has
not checked against a SNOMED release.
"""

from __future__ import annotations

from pharminform import fhir_export
from pharminform.report import GeneReportRow, ReportContext
from pharminform.uk_guidance import dpyd_recommendations, ugt1a1_recommendation


def _context() -> ReportContext:
    return ReportContext(
        patient_id="LOCAL-123",
        gene_results=[
            GeneReportRow(
                gene="DPYD", diplotype="*2A/*1", phenotype="intermediate metabolizer",
                calling_method="vcf", confidence_level="high",
            ),
            GeneReportRow(
                gene="UGT1A1", diplotype="*28/*28", phenotype="poor metabolizer",
                calling_method="vcf", confidence_level="high",
            ),
            GeneReportRow(
                gene="CYP2D6", diplotype="not called", phenotype=None,
                calling_method="none", confidence_level="low",
                confidence_reasons=["CYP2D6 cannot be called from a VCF"],
            ),
        ],
        cpic_recommendations=[
            {
                "drug": "fluorouracil",
                "phenotypes": {"DPYD": "Intermediate Metabolizer"},
                "recommendation": "Reduce starting dose by 50%.",
                "classification": "Strong",
                "citations": ["PMID:29152729"],
                "assumed_phenotypes": {},
            }
        ],
        uk_recommendations=[
            r.to_dict() for r in dpyd_recommendations("*2A/*1", drugs=["fluorouracil"])
        ]
        + [ugt1a1_recommendation("poor metabolizer", "*28/*28").to_dict()],
        composite_results=[],
    )


def _resources(bundle: dict, resource_type: str) -> list[dict]:
    return [e["resource"] for e in bundle["entry"] if e["resource"]["resourceType"] == resource_type]


def test_bundle_shape_and_report_profile():
    bundle = fhir_export.report_to_bundle(_context())
    assert bundle["resourceType"] == "Bundle"
    assert bundle["type"] == "collection"
    report = _resources(bundle, "DiagnosticReport")[0]
    assert report["meta"]["profile"] == [
        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomic-report"
    ]
    # Required pattern on the genomic-report profile.
    assert report["code"]["coding"][0] == {
        "system": "http://loinc.org",
        "code": "51969-4",
        "display": "Genetic analysis report",
    }
    # Every observation in the bundle is referenced from the report.
    observation_ids = {o["id"] for o in _resources(bundle, "Observation")}
    assert {r["reference"].split("/", 1)[1] for r in report["result"]} == observation_ids


def test_genotype_observations_use_the_ig_profile_and_gene_studied_component():
    bundle = fhir_export.report_to_bundle(_context())
    genotypes = [
        o
        for o in _resources(bundle, "Observation")
        if o["code"]["coding"][0]["code"] == "84413-4"
    ]
    assert {o["valueCodeableConcept"]["text"] for o in genotypes} == {
        "DPYD *2A/*1",
        "UGT1A1 *28/*28",
    }
    dpyd = next(o for o in genotypes if o["valueCodeableConcept"]["text"].startswith("DPYD"))
    assert dpyd["meta"]["profile"] == [
        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genotype"
    ]
    gene = next(c for c in dpyd["component"] if c["code"]["coding"][0]["code"] == "48018-6")
    assert gene["valueCodeableConcept"]["coding"][0]["code"] == "HGNC:3012"


def test_a_gene_that_was_not_called_produces_no_observation():
    bundle = fhir_export.report_to_bundle(_context())
    assert not any(
        "CYP2D6" in o.get("valueCodeableConcept", {}).get("text", "")
        for o in _resources(bundle, "Observation")
    )


def test_ugt1a1_phenotype_carries_the_cersi_snomed_code_first():
    bundle = fhir_export.report_to_bundle(_context())
    ugt = next(
        o
        for o in _resources(bundle, "Observation")
        if o.get("valueCodeableConcept", {}).get("text", "").startswith("UGT1A1")
    )
    effect = next(c for c in ugt["component"] if c["code"]["coding"][0]["code"] == "functional-effect")
    codings = effect["valueCodeableConcept"]["coding"]
    assert codings[0] == {
        "system": "http://snomed.info/sct",
        "code": "738536002",
        "display": "poor metabolizer",
    }


def test_no_snomed_code_is_invented_for_a_gene_without_a_checked_one():
    bundle = fhir_export.report_to_bundle(_context())
    dpyd = next(
        o
        for o in _resources(bundle, "Observation")
        if o.get("valueCodeableConcept", {}).get("text", "").startswith("DPYD")
    )
    effect = next(c for c in dpyd["component"] if c["code"]["coding"][0]["code"] == "functional-effect")
    systems = {c["system"] for c in effect["valueCodeableConcept"]["coding"]}
    assert systems == {"http://loinc.org"}
    assert fhir_export._phenotype_snomed("DPYD", "intermediate metabolizer") is None


def test_implications_are_typed_and_linked_to_the_genotype_they_rest_on():
    bundle = fhir_export.report_to_bundle(_context())
    implications = [
        o
        for o in _resources(bundle, "Observation")
        if o["code"]["coding"][0]["code"] == "therapeutic-implication"
    ]
    assert implications
    for obs in implications:
        assert obs["code"]["coding"][0]["system"] == (
            "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs"
        )
        assert obs["meta"]["profile"] == [
            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/therapeutic-implication"
        ]
        drug = next(c for c in obs["component"] if c["code"]["coding"][0]["code"] == "51963-7")
        assert drug["valueCodeableConcept"]["text"]
        assert obs["derivedFrom"]


def test_both_cpic_and_uk_advice_appear_each_labelled_with_its_source():
    bundle = fhir_export.report_to_bundle(_context())
    notes = [
        n["text"]
        for o in _resources(bundle, "Observation")
        for n in o.get("note", [])
        if n["text"].startswith("Source: ")
    ]
    assert "Source: CPIC" in notes
    assert "Source: UK Chemotherapy Board 2020" in notes
    assert any(n.startswith("Source: CERSI-PGx") for n in notes)


def test_uk_implication_states_the_dose_percentages():
    bundle = fhir_export.report_to_bundle(_context())
    texts = [n["text"] for o in _resources(bundle, "Observation") for n in o.get("note", [])]
    assert any("Starting dose: 50%" in t and "at most 75%" in t for t in texts)
    assert any("Starting dose: 70%" in t for t in texts)


def test_a_bundle_without_a_patient_reference_carries_no_patient_data():
    bundle = fhir_export.report_to_bundle(_context())
    serialized = repr(bundle)
    assert "LOCAL-123" not in serialized
    assert "subject" not in serialized
    assert not _resources(bundle, "Patient")


def test_patient_reference_is_used_verbatim_as_the_subject():
    bundle = fhir_export.report_to_bundle(_context(), "Patient/abc")
    for entry in bundle["entry"]:
        assert entry["resource"]["subject"] == {"reference": "Patient/abc"}


def test_research_use_only_is_stated_on_the_report():
    report = _resources(fhir_export.report_to_bundle(_context()), "DiagnosticReport")[0]
    assert "Research use only" in report["conclusion"]
    assert "not been validated" in report["conclusion"]


def test_resource_ids_are_stable_for_the_same_content():
    first = fhir_export.report_to_bundle(_context())
    context = _context()
    context.generated_at = first["timestamp"]
    second = fhir_export.report_to_bundle(context)
    assert [e["resource"]["id"] for e in first["entry"]] == [
        e["resource"]["id"] for e in second["entry"]
    ]
