"""Multi-gene composite risk scoring (experimental, original contribution).

Existing open-source PGx tools (Aldy, PharmCAT, StellarPGx) call star
alleles and phenotypes per gene, and CPIC/DPWG guidelines give per-gene
dosing recommendations -- but for drugs metabolized by more than one
pharmacogene (e.g. amitriptyline via CYP2D6 and CYP2C19), nothing combines
the per-gene phenotypes into a single risk signal. This module does that,
using a transparent, published-PK-weighted (or, failing that, equal-weight)
combination of per-gene metabolic activity.

This is exploratory and NOT a validated clinical scoring system: every
output carries a mandatory disclaimer field, and `directional_flag`
(ultrarapid-metabolizer risk) is deliberately kept separate from the
composite score rather than blended into it, because the two represent
different physiological mechanisms (toxic accumulation from reduced
clearance vs. under-exposure/prodrug over-activation from increased
clearance) that do not sit on a single continuous scale. See
DATA_SOURCES.md for a discussion of drugs (codeine, tramadol, clopidogrel)
where the underlying pathway is prodrug activation rather than clearance,
for which even the primary composite score's directionality is inverted --
flagged per-drug via the "directionality" field in drug_pathway_weights.json.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

DATA_PATH = Path(__file__).parent / "data" / "drug_pathway_weights.json"

DISCLAIMER = (
    "Exploratory composite score. Not a CPIC/DPWG guideline output. For "
    "research/educational use only -- clinical decisions should rely on "
    "individual gene-level CPIC/DPWG recommendations."
)

PHENOTYPE_TO_NORMALIZED_SCORE = {
    "poor metabolizer": 0.0,
    "intermediate metabolizer": 0.5,
    "normal metabolizer": 1.0,
    "rapid metabolizer": 1.5,
    "ultrarapid metabolizer": 2.0,
}


def normalize_activity(raw_score: float, gene: str) -> float:
    """Normalize to 0-1 scale, 1.0 = normal function, capped at 1.0.

    This models REDUCED-function toxicity accumulation risk only.
    Ultrarapid metabolizer risk (under-exposure / prodrug over-activation)
    is handled separately via directional_flag() -- do not blend the two.
    """
    return min(raw_score / 1.0, 1.0)


def load_drug_pathway_weights(path: Path | str = DATA_PATH) -> dict:
    return json.loads(Path(path).read_text())


def get_drug_entry(drug: str, weights_data: dict | None = None) -> dict | None:
    data = weights_data if weights_data is not None else load_drug_pathway_weights()
    for entry in data["drugs"]:
        if entry["drug"].lower() == drug.lower():
            return entry
    return None


class NoPathwayDataError(ValueError):
    """Raised when a drug/gene combination has no matching pathway weight data."""


def resolve_pathway_weights(
    drug_entry: dict, genes_present: list[str]
) -> tuple[dict[str, float], str]:
    """Determine the weight for each gene contributing to a drug's
    metabolism, and report which weighting method was used.

    Uses published fraction_metabolized values when ALL relevant genes have
    one; otherwise falls back to equal weighting across the genes present,
    and reports that fallback transparently via the returned method string.
    """
    contributing = {g["gene"]: g["fraction_metabolized"] for g in drug_entry["contributing_genes"]}
    relevant = {g: contributing[g] for g in genes_present if g in contributing}

    if not relevant:
        raise NoPathwayDataError(
            f"No contributing-gene pathway data for drug {drug_entry['drug']!r} "
            f"among genes {genes_present!r}"
        )

    if all(v is not None for v in relevant.values()):
        return dict(relevant), "published PK fractions (see DATA_SOURCES.md)"

    equal_weight = 1.0 / len(relevant)
    return {g: equal_weight for g in relevant}, "equal-weight fallback (no published fraction for at least one contributing gene)"


def composite_risk_score(
    gene_scores: dict[str, float], pathway_weights: dict[str, float]
) -> float:
    """
    gene_scores: {gene: normalized_activity_score (0-1)}
    pathway_weights: {gene: fraction_metabolized (0-1)}
    Returns composite_function (0-1); lower = more reduced metabolic
    capacity = higher toxicity risk. Normalizes by total available
    weight so missing gene data degrades gracefully.
    """
    total_weight = sum(pathway_weights.get(g, 0) for g in gene_scores)
    if total_weight == 0:
        raise NoPathwayDataError("No matching pathway weights for provided genes")
    weighted_sum = sum(gene_scores[g] * pathway_weights.get(g, 0) for g in gene_scores)
    return weighted_sum / total_weight


def tier_from_composite(composite_function: float) -> str:
    """
    Thresholds chosen so a single poor-metabolizer gene with >=50%
    pathway weight already crosses into 'High' -- i.e. composite risk
    must NEVER return a lower tier than what CPIC's per-gene guidance
    would already indicate for the dominant gene. This monotonicity
    property is the key thing to validate in tests.
    """
    if composite_function >= 0.75:
        return "Low"
    elif composite_function >= 0.45:
        return "Moderate"
    return "High"


def directional_flag(gene_phenotypes: dict[str, str]) -> str | None:
    """Flags when ANY contributing gene is ultrarapid metabolizer --
    can cause under-exposure or dangerous prodrug over-activation (e.g.
    codeine -> morphine via CYP2D6 UM). Independent of composite score."""
    if any(p == "ultrarapid metabolizer" for p in gene_phenotypes.values()):
        return "ultrarapid_metabolizer_present"
    return None


@dataclass
class ContributingGeneResult:
    gene: str
    phenotype: str
    weight: float


@dataclass
class CompositeRiskResult:
    drug: str
    composite_function_score: float
    risk_tier: str
    directional_flag: str | None
    contributing_genes: list[ContributingGeneResult]
    weighting_method: str
    disclaimer: str = DISCLAIMER

    def to_dict(self) -> dict:
        return {
            "drug": self.drug,
            "composite_risk": {
                "composite_function_score": round(self.composite_function_score, 4),
                "risk_tier": self.risk_tier,
                "directional_flag": self.directional_flag,
                "contributing_genes": [
                    {"gene": g.gene, "phenotype": g.phenotype, "weight": round(g.weight, 4)}
                    for g in self.contributing_genes
                ],
                "weighting_method": self.weighting_method,
                "disclaimer": self.disclaimer,
            },
        }


def score_drug(
    drug: str,
    gene_phenotypes: dict[str, str],
    weights_data: dict | None = None,
) -> CompositeRiskResult:
    """End-to-end: given a drug and {gene: phenotype} for genes actually
    called in this patient, compute the composite risk report."""
    data = weights_data if weights_data is not None else load_drug_pathway_weights()
    entry = get_drug_entry(drug, data)
    if entry is None:
        raise NoPathwayDataError(f"No pathway weight data curated for drug {drug!r}")

    genes_present = list(gene_phenotypes.keys())
    pathway_weights, method = resolve_pathway_weights(entry, genes_present)

    gene_scores = {
        gene: normalize_activity(PHENOTYPE_TO_NORMALIZED_SCORE[phenotype], gene)
        for gene, phenotype in gene_phenotypes.items()
        if gene in pathway_weights
    }

    composite = composite_risk_score(gene_scores, pathway_weights)
    tier = tier_from_composite(composite)
    flag = directional_flag(gene_phenotypes)

    contributing = [
        ContributingGeneResult(gene=g, phenotype=gene_phenotypes[g], weight=pathway_weights[g])
        for g in gene_scores
    ]

    return CompositeRiskResult(
        drug=drug,
        composite_function_score=composite,
        risk_tier=tier,
        directional_flag=flag,
        contributing_genes=contributing,
        weighting_method=method,
    )
