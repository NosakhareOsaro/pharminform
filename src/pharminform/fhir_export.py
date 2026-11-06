"""FHIR R4 output following the HL7 Genomics Reporting IG (STU3).

A PDF is the end of the road for a result: nobody can act on it
automatically, and an oncology prescribing system cannot read it. The
CERSI-PGx UGT1A1 guideline says as much -- results should sit alongside
other pharmacogenomic results as structured data, with SNOMED-CT codes
(its Table 6). This module emits the report as a FHIR Bundle so a
Genomic Laboratory Hub or an EPR can ingest it.

Profiles used (canonical URLs from
https://hl7.org/fhir/uv/genomics-reporting/STU3/artifacts.html):

* ``genomic-report`` -- DiagnosticReport, with the profile's required
  ``DiagnosticReport.code`` pattern of LOINC 51969-4.
* ``genotype`` -- one Observation per gene, Observation.code LOINC
  84413-4, with the gene in the LOINC 48018-6 ``gene-studied`` component.
* ``therapeutic-implication`` -- one Observation per dosing
  recommendation, Observation.code ``therapeutic-implication`` from the
  IG's ``tbd-codes-cs`` CodeSystem, ``derivedFrom`` the genotype
  Observation(s) it rests on, with the drug in the LOINC 51963-7
  ``medication-assessed`` component and the advice text as the
  ``conclusion-string`` component.

Two honest limitations, surfaced rather than hidden:

* The metabolizer phenotype is emitted as an extra coding on the
  genotype Observation's value using LOINC 53040-2 ("Genetic variation's
  effect on drug metabolism"). Some LOINC mirrors list that code as
  retained for historical use. Where this project has a SNOMED-CT code
  for the phenotype -- currently the three UGT1A1 findings named in
  CERSI-PGx Table 6 -- that coding is emitted too, and it is the one a
  UK consumer should prefer.
* This output has not been validated against the IG's own validator, and
  the project is research use only. The bundle carries that statement as
  a ``conclusion-string`` on the DiagnosticReport.

Nothing here invents a code. A phenotype with no SNOMED code gets no
SNOMED coding, and a recommendation with no source gets no citation.
"""

from __future__ import annotations

import hashlib
from typing import Any

from . import __version__
from .report import ReportContext

IG_BASE = "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition"
TBD_CODES = "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs"
LOINC = "http://loinc.org"
SNOMED = "http://snomed.info/sct"

GENOMIC_REPORT_PROFILE = f"{IG_BASE}/genomic-report"
GENOTYPE_PROFILE = f"{IG_BASE}/genotype"
THERAPEUTIC_IMPLICATION_PROFILE = f"{IG_BASE}/therapeutic-implication"

REPORT_CODE = "51969-4"  # Required pattern on the genomic-report profile.
GENOTYPE_CODE = "84413-4"  # Genotype display name.
GENE_STUDIED_CODE = "48018-6"  # Gene studied [ID].
MEDICATION_ASSESSED_CODE = "51963-7"  # Medication assessed [ID].
DRUG_METABOLISM_CODE = "53040-2"  # Genetic variation's effect on drug metabolism.

# HGNC gene ids, so a consumer does not have to parse the symbol.
HGNC_IDS = {
    "DPYD": "HGNC:3012",
    "TPMT": "HGNC:12014",
    "UGT1A1": "HGNC:12530",
    "CYP2C9": "HGNC:2623",
    "CYP2C19": "HGNC:2621",
    "CYP2D6": "HGNC:2625",
    "NUDT15": "HGNC:23063",
}

RESEARCH_USE_ONLY = (
    "Research use only. PharmInform is not a medical device and this output has "
    "not been validated against the HL7 Genomics Reporting IG validator. Do not "
    "use it to make a prescribing decision without an accredited laboratory result."
)


def _id(*parts: str) -> str:
    """A stable, non-identifying resource id derived from its content."""
    digest = hashlib.sha256("|".join(parts).encode()).hexdigest()[:12]
    return f"{parts[0].lower().replace('_', '-')}-{digest}"


def _codeable(system: str, code: str, display: str | None = None) -> dict:
    coding: dict[str, Any] = {"system": system, "code": code}
    if display:
        coding["display"] = display
    return {"coding": [coding]}


def _component(code: dict, **value: Any) -> dict:
    return {"code": code, **value}


def _genotype_observation(row: dict, subject: dict | None) -> dict:
    gene = row["gene"]
    obs: dict[str, Any] = {
        "resourceType": "Observation",
        "id": _id("genotype", gene, row["diplotype"]),
        "meta": {"profile": [GENOTYPE_PROFILE]},
        "status": "final",
        "category": [_codeable("http://terminology.hl7.org/CodeSystem/observation-category",
                               "laboratory")],
        "code": _codeable(LOINC, GENOTYPE_CODE, "Genotype display name"),
        "valueCodeableConcept": {"text": f"{gene} {row['diplotype']}"},
        "component": [
            _component(
                _codeable(LOINC, GENE_STUDIED_CODE, "Gene studied [ID]"),
                valueCodeableConcept=(
                    _codeable("http://www.genenames.org/geneId", HGNC_IDS[gene], gene)
                    if gene in HGNC_IDS
                    else {"text": gene}
                ),
            )
        ],
    }
    if subject:
        obs["subject"] = subject

    phenotype = row.get("phenotype")
    if phenotype:
        value = {
            "coding": [{"system": LOINC, "code": DRUG_METABOLISM_CODE,
                        "display": "Genetic variation's effect on drug metabolism"}],
            "text": phenotype,
        }
        snomed = _phenotype_snomed(gene, phenotype)
        if snomed:
            value["coding"].insert(0, {"system": SNOMED, "code": snomed, "display": phenotype})
        obs["component"].append(
            _component(_codeable(TBD_CODES, "functional-effect", "Functional Effect"),
                       valueCodeableConcept=value)
        )

    notes = list(row.get("confidence_reasons") or [])
    notes.insert(0, f"Confidence: {row.get('confidence_level')}; "
                    f"calling method: {row.get('calling_method')}")
    obs["note"] = [{"text": n} for n in notes]
    return obs


_UGT1A1_SNOMED = {
    "normal metabolizer": "738538001",
    "intermediate metabolizer": "738537006",
    "poor metabolizer": "738536002",
}


def _phenotype_snomed(gene: str, phenotype: str) -> str | None:
    """SNOMED-CT finding code for a phenotype, where this project has one.

    Only the three UGT1A1 findings listed in CERSI-PGx Table 6 are
    encoded. Other genes have SNOMED findings, but they have not been
    checked against a SNOMED release here, so none is emitted.
    """
    if gene != "UGT1A1":
        return None
    return _UGT1A1_SNOMED.get(phenotype.lower())


def _implication_observation(
    drug: str,
    advice: str,
    derived_from: list[str],
    subject: dict | None,
    source_label: str,
    extra_notes: list[str] | None = None,
) -> dict:
    obs: dict[str, Any] = {
        "resourceType": "Observation",
        "id": _id("implication", source_label, drug, advice),
        "meta": {"profile": [THERAPEUTIC_IMPLICATION_PROFILE]},
        "status": "final",
        "code": _codeable(TBD_CODES, "therapeutic-implication", "Therapeutic Implication"),
        "component": [
            _component(
                _codeable(LOINC, MEDICATION_ASSESSED_CODE, "Medication assessed [ID]"),
                valueCodeableConcept={"text": drug},
            ),
            _component(
                _codeable(TBD_CODES, "conclusion-string", "Conclusion Text"),
                valueString=advice,
            ),
        ],
        "note": [{"text": f"Source: {source_label}"}]
        + [{"text": n} for n in (extra_notes or []) if n],
    }
    if subject:
        obs["subject"] = subject
    if derived_from:
        obs["derivedFrom"] = [{"reference": f"Observation/{ref}"} for ref in derived_from]
    return obs


def report_to_bundle(context: ReportContext, patient_reference: str | None = None) -> dict:
    """The report as a FHIR R4 Bundle of type ``collection``.

    ``patient_reference`` is used verbatim as a FHIR reference (e.g.
    ``"Patient/abc"``). No patient resource is created and no identifier
    from the report is written into the bundle, so a bundle built without
    it carries no patient data at all.
    """
    subject = {"reference": patient_reference} if patient_reference else None
    entries: list[dict] = []
    genotype_ids: dict[str, str] = {}

    for row in context.to_dict()["gene_results"]:
        if row["diplotype"] == "not called":
            continue
        obs = _genotype_observation(row, subject)
        genotype_ids[row["gene"]] = obs["id"]
        entries.append(obs)

    for rec in context.cpic_recommendations:
        genes = sorted(rec.get("phenotypes", {}))
        derived = [genotype_ids[g] for g in genes if g in genotype_ids]
        notes = []
        if rec.get("assumed_phenotypes"):
            notes.append(
                "Assumes, not tested: "
                + ", ".join(f"{g} {p}" for g, p in sorted(rec["assumed_phenotypes"].items()))
            )
        notes.extend(rec.get("citations") or [])
        entries.append(
            _implication_observation(
                rec["drug"], rec["recommendation"], derived, subject, "CPIC", notes
            )
        )

    for rec in context.uk_recommendations:
        src = rec["source"]
        label = f"{src['issuer']} {src['version']}"
        notes = []
        if rec.get("start_dose_percent") is not None:
            dose = f"Starting dose: {rec['start_dose_percent']}% of the full dose"
            if rec.get("max_dose_percent"):
                dose += f"; escalate to at most {rec['max_dose_percent']}%"
            if rec.get("max_dose_percent_with_caution"):
                dose += f" ({rec['max_dose_percent_with_caution']}% with caution)"
            notes.append(dose + ".")
        if rec.get("dpd_activity"):
            notes.append(f"Estimated DPD activity: {rec['dpd_activity']}.")
        notes.extend(rec.get("caveats") or [])
        if rec.get("titration"):
            notes.append(rec["titration"])
        if rec.get("alternatives"):
            notes.append("Alternatives named by the guidance: " + ", ".join(rec["alternatives"]))
        for key in ("doi", "pmid", "url"):
            if src.get(key):
                notes.append(f"{key.upper()}: {src[key]}")
        derived = [genotype_ids[rec["gene"]]] if rec["gene"] in genotype_ids else []
        entries.append(
            _implication_observation(
                rec["drug"], rec["advice"], derived, subject, label, notes
            )
        )

    report: dict[str, Any] = {
        "resourceType": "DiagnosticReport",
        "id": _id("report", context.generated_at),
        "meta": {"profile": [GENOMIC_REPORT_PROFILE]},
        "status": "final",
        "category": [_codeable("http://terminology.hl7.org/CodeSystem/v2-0074", "GE",
                               "Genetics")],
        "code": _codeable(LOINC, REPORT_CODE, "Genetic analysis report"),
        "effectiveDateTime": context.generated_at,
        "issued": context.generated_at,
        "result": [{"reference": f"Observation/{e['id']}"} for e in entries],
        "conclusion": RESEARCH_USE_ONLY,
        "extension": [
            {
                "url": f"{IG_BASE}/genomic-report-note",
                "valueString": f"Generated by PharmInform {__version__}.",
            }
        ],
    }
    if subject:
        report["subject"] = subject

    return {
        "resourceType": "Bundle",
        "type": "collection",
        "timestamp": context.generated_at,
        "entry": [{"resource": r} for r in [report, *entries]],
    }
