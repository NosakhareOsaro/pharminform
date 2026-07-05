"""OMOP CDM drug_exposure mapping.

Maps a PGx report (star allele calls, phenotypes, CPIC recommendation) onto
an OMOP Common Data Model v5.4-shaped ``drug_exposure`` record, so PGx
findings can be joined into an existing OMOP-based clinical data
warehouse.

Two honest limitations, documented rather than papered over:

1. OMOP's standard vocabulary uses an internal surrogate ``concept_id``
   (an integer) for every drug, resolved via the OHDSI vocabulary tables
   (distributed through Athena, https://athena.ohdsi.org -- a free but
   registration-gated, multi-GB download, not fetchable at runtime from
   this package). What this module CAN populate directly is the RxNorm
   identifier itself (``drug_source_value`` / ``drug_source_concept_code``),
   which CPIC's own data already carries for every drug in
   ``cpic_dpwg.py``'s cache (e.g. "RxNorm:704" for amitriptyline). Real
   ``drug_concept_id`` resolution is left as an explicit TODO for sites
   that have the OMOP vocabulary tables loaded locally.

2. Standard OMOP CDM has no native fields for star alleles, diplotypes, or
   metabolizer phenotype. Rather than inventing non-standard columns on
   ``drug_exposure`` itself, PGx context is attached as a structured
   ``pgx_annotation`` side-record, deliberately shaped like a minimal
   version of OHDSI's Genomic CDM extension (GENOMIC_TEST /
   VARIANT_OCCURRENCE-style annotation tables), without claiming to
   implement that extension in full.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date


@dataclass
class DrugExposureRecord:
    """Shaped after OMOP CDM v5.4's drug_exposure table. Fields not
    resolvable without a loaded OMOP vocabulary are left as None/0 rather
    than guessed, per OMOP convention (concept_id 0 = "no matching concept")."""

    person_id: int
    drug_concept_id: int  # 0 = unmapped; see module docstring limitation (1)
    drug_exposure_start_date: date
    drug_exposure_end_date: date
    drug_type_concept_id: int  # e.g. 32838 = "Drug recommendation" custom type; site-specific
    drug_source_value: str  # human-readable drug name
    drug_source_concept_code: str | None = None  # e.g. "RxNorm:704"
    quantity: float | None = None
    days_supply: int | None = None
    sig: str | None = None  # free-text dosing instructions (CPIC recommendation text)


@dataclass
class PGxAnnotation:
    """Minimal genomic-context side-record for a drug_exposure, not part of
    standard OMOP CDM. Shaped loosely after OHDSI's Genomic CDM extension
    but not a claimed implementation of it."""

    person_id: int
    gene_symbol: str
    diplotype: str
    phenotype: str
    calling_method: str  # "aldy" or "rule_based_vcf"
    confidence: str  # "high" or "low", see confidence.py
    cpic_recommendation: str | None
    cpic_classification: str | None
    linked_drug_exposure: str | None = None  # drug_source_value this annotation informs
    citations: list[str] = field(default_factory=list)


def to_drug_exposure(
    person_id: int,
    drug_name: str,
    drug_rxnorm_id: str | None,
    recommendation_text: str,
    exposure_start: date,
    exposure_end: date,
    drug_type_concept_id: int = 32838,
) -> DrugExposureRecord:
    """Build a drug_exposure record for a PGx-informed prescribing decision."""
    return DrugExposureRecord(
        person_id=person_id,
        drug_concept_id=0,  # not resolvable without local OMOP vocabulary tables
        drug_exposure_start_date=exposure_start,
        drug_exposure_end_date=exposure_end,
        drug_type_concept_id=drug_type_concept_id,
        drug_source_value=drug_name,
        drug_source_concept_code=drug_rxnorm_id,
        sig=recommendation_text,
    )


def to_pgx_annotation(
    person_id: int,
    gene: str,
    diplotype: str,
    phenotype: str,
    calling_method: str,
    confidence: str,
    cpic_recommendation: str | None = None,
    cpic_classification: str | None = None,
    linked_drug: str | None = None,
    citations: list[str] | None = None,
) -> PGxAnnotation:
    return PGxAnnotation(
        person_id=person_id,
        gene_symbol=gene,
        diplotype=diplotype,
        phenotype=phenotype,
        calling_method=calling_method,
        confidence=confidence,
        cpic_recommendation=cpic_recommendation,
        cpic_classification=cpic_classification,
        linked_drug_exposure=linked_drug,
        citations=citations or [],
    )


def records_to_rows(records: list[DrugExposureRecord]) -> list[dict]:
    """Flatten to plain dicts suitable for writing to CSV/a drug_exposure
    table via pandas.DataFrame(...) or a DB bulk insert."""
    return [
        {
            "person_id": r.person_id,
            "drug_concept_id": r.drug_concept_id,
            "drug_exposure_start_date": r.drug_exposure_start_date.isoformat(),
            "drug_exposure_end_date": r.drug_exposure_end_date.isoformat(),
            "drug_type_concept_id": r.drug_type_concept_id,
            "drug_source_value": r.drug_source_value,
            "drug_source_concept_code": r.drug_source_concept_code,
            "quantity": r.quantity,
            "days_supply": r.days_supply,
            "sig": r.sig,
        }
        for r in records
    ]
