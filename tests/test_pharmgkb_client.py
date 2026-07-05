import responses

from pharminform.pharmgkb_client import (
    PharmGKBClient,
    diplotype_activity_score,
    diplotype_to_phenotype,
)


@responses.activate
def test_get_gene_mocked():
    responses.add(
        responses.GET,
        "https://api.pharmgkb.org/v1/data/gene",
        json={"data": [{"symbol": "CYP2D6", "id": "PA128"}]},
        status=200,
    )
    client = PharmGKBClient()
    gene = client.get_gene("CYP2D6")
    assert gene["id"] == "PA128"


@responses.activate
def test_get_gene_no_match_returns_none():
    responses.add(
        responses.GET,
        "https://api.pharmgkb.org/v1/data/gene",
        json={"data": []},
        status=200,
    )
    client = PharmGKBClient()
    assert client.get_gene("NOTAGENE") is None


def test_activity_score_wildtype_normal():
    score = diplotype_activity_score("*1/*1", {})
    assert score == 2.0


def test_activity_score_one_no_function_allele():
    score = diplotype_activity_score("*4/*1", {"*4": "no function"})
    assert score == 1.0


def test_diplotype_to_phenotype_poor_metabolizer():
    phenotype, score = diplotype_to_phenotype(
        "CYP2D6", "*4/*4", {"*4": "no function"}
    )
    assert phenotype == "poor metabolizer"
    assert score == 0.0


def test_diplotype_to_phenotype_ultrarapid():
    phenotype, score = diplotype_to_phenotype(
        "CYP2D6", "*1x2/*1x2", {"*1x2": "increased function"}
    )
    assert phenotype == "ultrarapid metabolizer"
    assert score == 4.0


def test_diplotype_to_phenotype_tpmt_three_tier():
    phenotype, _ = diplotype_to_phenotype(
        "TPMT", "*3C/*1", {"*3C": "no function"}
    )
    assert phenotype == "intermediate metabolizer"
