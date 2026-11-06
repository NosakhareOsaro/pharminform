"""Regimen-level view: one combination chemotherapy regimen, several genes.

Gene-by-gene reports leave the last step to the reader. An oncology
pharmacist checking a FOLFIRI prescription needs two separate dose
decisions -- the fluoropyrimidine dose from DPYD and the irinotecan dose
from UGT1A1 -- and needs to know that CERSI-PGx says to wait for *both*
results before prescribing, because acting on the first one alone can
mean a second dose change a cycle later.

This module does no new pharmacology. It groups the per-drug decisions
the UK guidance layer already produces (``uk_guidance``) under the
regimen a patient is actually being prescribed, names the drugs in that
regimen that no gene in the panel informs, and states plainly which
required results are still missing. A regimen whose genes were not all
called is reported as incomplete; it is never scored as if the missing
gene were normal.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .uk_guidance import UKRecommendation, load_uk_guidance, uk_recommendations

# Drug -> the panel gene whose result changes its dose.
DRUG_GENE = {
    "fluorouracil": "DPYD",
    "capecitabine": "DPYD",
    "tegafur": "DPYD",
    "irinotecan": "UGT1A1",
}

REGIMENS: dict[str, dict] = {
    "FOLFIRI": {
        "label": "FOLFIRI (folinic acid, 5-fluorouracil, irinotecan)",
        "drugs": ["fluorouracil", "irinotecan"],
        "other_drugs": ["folinic acid"],
    },
    "CAPIRI": {
        "label": "CAPIRI / XELIRI (capecitabine, irinotecan)",
        "drugs": ["capecitabine", "irinotecan"],
        "other_drugs": [],
    },
    "FOLFOXIRI": {
        "label": "FOLFOXIRI (folinic acid, 5-fluorouracil, oxaliplatin, irinotecan)",
        "drugs": ["fluorouracil", "irinotecan"],
        "other_drugs": ["folinic acid", "oxaliplatin"],
    },
    "FOLFIRINOX": {
        "label": "FOLFIRINOX (folinic acid, 5-fluorouracil, irinotecan, oxaliplatin)",
        "drugs": ["fluorouracil", "irinotecan"],
        "other_drugs": ["folinic acid", "oxaliplatin"],
    },
    "FOLFOX": {
        "label": "FOLFOX (folinic acid, 5-fluorouracil, oxaliplatin)",
        "drugs": ["fluorouracil"],
        "other_drugs": ["folinic acid", "oxaliplatin"],
    },
    "CAPOX": {
        "label": "CAPOX / XELOX (capecitabine, oxaliplatin)",
        "drugs": ["capecitabine"],
        "other_drugs": ["oxaliplatin"],
    },
    "IRINOTECAN": {
        "label": "Single-agent irinotecan",
        "drugs": ["irinotecan"],
        "other_drugs": [],
    },
}


@dataclass
class DrugDecision:
    drug: str
    gene: str
    status: str  # "advice", "missing result", or "no panel gene"
    recommendations: list[UKRecommendation] = field(default_factory=list)
    note: str = ""

    def to_dict(self) -> dict:
        return {
            "drug": self.drug,
            "gene": self.gene,
            "status": self.status,
            "note": self.note,
            "recommendations": [r.to_dict() for r in self.recommendations],
        }


@dataclass
class RegimenAssessment:
    regimen: str
    label: str
    decisions: list[DrugDecision]
    unmodified_drugs: list[str]
    missing_genes: list[str]
    notes: list[str] = field(default_factory=list)

    @property
    def complete(self) -> bool:
        return not self.missing_genes

    def to_dict(self) -> dict:
        return {
            "regimen": self.regimen,
            "label": self.label,
            "complete": self.complete,
            "missing_genes": list(self.missing_genes),
            "unmodified_drugs": list(self.unmodified_drugs),
            "notes": list(self.notes),
            "decisions": [d.to_dict() for d in self.decisions],
        }


def available_regimens() -> list[str]:
    return sorted(REGIMENS)


def assess_regimen(
    regimen: str,
    gene_phenotypes: dict[str, str],
    gene_diplotypes: dict[str, str] | None = None,
    data: dict | None = None,
) -> RegimenAssessment:
    """Per-drug UK dosing decisions for one regimen.

    ``gene_phenotypes`` and ``gene_diplotypes`` hold only genes that were
    actually called with a usable result; a gene absent from both is
    reported as a missing result, never assumed normal.
    """
    key = regimen.upper()
    if key not in REGIMENS:
        raise KeyError(f"Unknown regimen {regimen!r}; known: {', '.join(available_regimens())}")
    spec = REGIMENS[key]
    data = data if data is not None else load_uk_guidance()
    gene_diplotypes = gene_diplotypes or {}

    all_recs = uk_recommendations(gene_phenotypes, gene_diplotypes, data)
    by_drug: dict[str, list[UKRecommendation]] = {}
    for rec in all_recs:
        by_drug.setdefault(rec.drug, []).append(rec)

    decisions: list[DrugDecision] = []
    missing: list[str] = []
    for drug in spec["drugs"]:
        gene = DRUG_GENE[drug]
        recs = by_drug.get(drug, [])
        if recs:
            decisions.append(DrugDecision(drug, gene, "advice", recs))
        else:
            missing.append(gene)
            decisions.append(
                DrugDecision(
                    drug,
                    gene,
                    "missing result",
                    note=(
                        f"No usable {gene} result, so no dose advice is given for {drug}. "
                        f"This is not a normal-metaboliser result."
                    ),
                )
            )

    notes: list[str] = []
    genes_in_regimen = {DRUG_GENE[d] for d in spec["drugs"]}
    if {"DPYD", "UGT1A1"} <= genes_in_regimen:
        notes.append(data["ugt1a1"]["combined_regimen_note"])
    if missing:
        notes.append(
            "Incomplete: "
            + ", ".join(sorted(set(missing)))
            + " not available. Do not prescribe on a partial result."
        )
    if "irinotecan" in spec["drugs"]:
        notes.append(data["sources"]["cersi_ugt1a1_2026"]["notes"])

    return RegimenAssessment(
        regimen=key,
        label=spec["label"],
        decisions=decisions,
        unmodified_drugs=list(spec["other_drugs"]),
        missing_genes=sorted(set(missing)),
        notes=notes,
    )


def assess_all_regimens(
    gene_phenotypes: dict[str, str],
    gene_diplotypes: dict[str, str] | None = None,
    data: dict | None = None,
) -> list[RegimenAssessment]:
    """Every regimen for which at least one gene in the panel was called."""
    data = data if data is not None else load_uk_guidance()
    called = {g for g in gene_diplotypes or {}} | set(gene_phenotypes)
    out = []
    for name, spec in REGIMENS.items():
        if not {DRUG_GENE[d] for d in spec["drugs"]} & called:
            continue
        out.append(assess_regimen(name, gene_phenotypes, gene_diplotypes, data))
    return sorted(out, key=lambda a: a.regimen)
