import responses

from pharminform.cpic_dpwg import CPICClient, GuidelineLookup, load_cached_guidelines


def test_load_cached_guidelines_has_expected_drugs():
    cache = load_cached_guidelines()
    for drug in ["codeine", "clopidogrel", "fluorouracil", "irinotecan", "azathioprine", "warfarin"]:
        assert drug in cache["drugs"]


def test_warfarin_has_no_recommendation_rows():
    # Real finding from the live CPIC API: warfarin dosing is algorithm-based,
    # not a phenotype -> recommendation lookup table.
    cache = load_cached_guidelines()
    assert cache["drugs"]["warfarin"]["recommendations"] == []


def test_lookup_codeine_ultrarapid():
    lookup = GuidelineLookup()
    result = lookup.lookup("codeine", {"CYP2D6": "ultrarapid metabolizer"})
    assert result is not None
    assert "avoid" in result.recommendation.lower()
    assert result.source == "CPIC"
    assert len(result.citations) > 0


def test_lookup_unknown_drug_returns_none():
    lookup = GuidelineLookup()
    assert lookup.lookup("not-a-real-drug", {"CYP2D6": "normal metabolizer"}) is None


def test_lookup_no_matching_phenotype_returns_none():
    lookup = GuidelineLookup()
    assert lookup.lookup("codeine", {"CYP2D6": "not-a-real-phenotype"}) is None


def test_azathioprine_tpmt_poor_metabolizer_recommends_reduction():
    lookup = GuidelineLookup()
    result = lookup.lookup("azathioprine", {"TPMT": "poor metabolizer"})
    assert result is not None
    assert result.classification.lower() in ("strong", "moderate")


@responses.activate
def test_cpic_client_get_recommendations_mocked():
    responses.add(
        responses.GET,
        "https://api.cpicpgx.org/v1/recommendation",
        json=[{"phenotypes": {"CYP2D6": "Poor Metabolizer"}, "drugrecommendation": "avoid"}],
        status=200,
    )
    client = CPICClient()
    rows = client.get_recommendations("RxNorm:2670")
    assert rows[0]["drugrecommendation"] == "avoid"
