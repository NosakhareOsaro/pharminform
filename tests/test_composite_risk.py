import pytest

from pharminform.composite_risk import (
    NoPathwayDataError,
    composite_risk_score,
    directional_flag,
    load_drug_pathway_weights,
    normalize_activity,
    resolve_pathway_weights,
    score_drug,
    tier_from_composite,
)

WEIGHTS_DATA = load_drug_pathway_weights()
ALL_DRUGS = [entry["drug"] for entry in WEIGHTS_DATA["drugs"]]


# --- Monotonicity: a poor metabolizer gene with >=50% weight must never
# --- resolve to "Low" risk, for every drug in the reference dataset.

@pytest.mark.parametrize("drug", ALL_DRUGS)
def test_monotonicity_poor_metabolizer_dominant_gene_never_low(drug):
    entry = next(e for e in WEIGHTS_DATA["drugs"] if e["drug"] == drug)
    genes = [g["gene"] for g in entry["contributing_genes"]]
    weights, _ = resolve_pathway_weights(entry, genes)

    # composite_risk_score normalizes by the total weight of genes actually
    # scored, so what matters for the monotonicity guarantee is each gene's
    # EFFECTIVE share of that total -- not its raw fraction_metabolized
    # (e.g. codeine's CYP2D6 is its only modeled gene, at a raw fraction of
    # just 0.10, but is 100% of what's being scored).
    total = sum(weights.values())
    dominant_genes = [g for g, w in weights.items() if (w / total) >= 0.5]
    assert dominant_genes, f"{drug} has no dominant (>=50% effective weight) gene to test"

    for dominant in dominant_genes:
        phenotypes = {g: "normal metabolizer" for g in genes}
        phenotypes[dominant] = "poor metabolizer"
        result = score_drug(drug, phenotypes, WEIGHTS_DATA)
        assert result.risk_tier != "Low", (
            f"{drug}: poor metabolizer {dominant} (weight={weights[dominant]}) "
            f"produced tier={result.risk_tier}, composite={result.composite_function_score}"
        )


# --- Missing pathway data handled gracefully.

def test_composite_score_ignores_genes_without_weights_instead_of_crashing():
    gene_scores = {"CYP2D6": 0.0, "CYP2C19": 1.0}
    pathway_weights = {"CYP2D6": 0.6}  # CYP2C19 has no weight entry
    score = composite_risk_score(gene_scores, pathway_weights)
    # only CYP2D6 contributes: weighted_sum = 0.0*0.6 = 0, total_weight = 0.6
    assert score == 0.0


def test_composite_score_raises_clear_error_when_no_overlap():
    with pytest.raises(NoPathwayDataError):
        composite_risk_score({"TPMT": 1.0}, {"CYP2D6": 0.6})


def test_score_drug_unknown_drug_raises_clear_error():
    with pytest.raises(NoPathwayDataError):
        score_drug("not-a-real-drug", {"CYP2D6": "normal metabolizer"}, WEIGHTS_DATA)


def test_score_drug_no_genes_present_for_drug_raises():
    with pytest.raises(NoPathwayDataError):
        score_drug("codeine", {"TPMT": "normal metabolizer"}, WEIGHTS_DATA)


# --- Equal-weight fallback.

def test_equal_weight_fallback_triggers_and_is_flagged():
    entry = next(e for e in WEIGHTS_DATA["drugs"] if e["drug"] == "amitriptyline")
    genes = [g["gene"] for g in entry["contributing_genes"]]
    weights, method = resolve_pathway_weights(entry, genes)
    assert weights == {"CYP2D6": 0.5, "CYP2C19": 0.5}
    assert "equal-weight fallback" in method

    result = score_drug(
        "amitriptyline",
        {"CYP2D6": "normal metabolizer", "CYP2C19": "normal metabolizer"},
        WEIGHTS_DATA,
    )
    assert "equal-weight fallback" in result.weighting_method


def test_published_fraction_used_when_available():
    result = score_drug("warfarin", {"CYP2C9": "normal metabolizer"}, WEIGHTS_DATA)
    assert result.weighting_method == "published PK fractions (see DATA_SOURCES.md)"
    assert result.contributing_genes[0].weight == 0.80


# --- Directional flag independence.

def test_directional_flag_fires_regardless_of_composite_tier():
    # Ultrarapid metabolizer normalizes to score 1.0 (capped), same as normal,
    # so composite tier reads "Low"/no elevated risk, but the flag must still fire.
    result = score_drug("codeine", {"CYP2D6": "ultrarapid metabolizer"}, WEIGHTS_DATA)
    assert result.directional_flag == "ultrarapid_metabolizer_present"
    assert result.risk_tier == "Low"


def test_directional_flag_absent_when_no_um_gene():
    assert directional_flag({"CYP2D6": "poor metabolizer", "CYP2C19": "normal metabolizer"}) is None


def test_directional_flag_present_with_one_um_among_several():
    flag = directional_flag({
        "CYP2D6": "poor metabolizer",
        "CYP2C19": "ultrarapid metabolizer",
    })
    assert flag == "ultrarapid_metabolizer_present"


# --- normalize_activity / tier_from_composite unit checks.

def test_normalize_activity_caps_at_one():
    assert normalize_activity(2.0, "CYP2D6") == 1.0


def test_tier_boundaries():
    assert tier_from_composite(0.75) == "Low"
    assert tier_from_composite(0.74) == "Moderate"
    assert tier_from_composite(0.45) == "Moderate"
    assert tier_from_composite(0.44) == "High"
