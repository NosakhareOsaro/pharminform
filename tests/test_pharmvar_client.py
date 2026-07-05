import pytest
import responses

from pharminform.pharmvar_client import PharmVarAuthError, PharmVarClient


def test_no_api_key_raises_clear_error(monkeypatch):
    monkeypatch.delenv("PHARMVAR_API_KEY", raising=False)
    client = PharmVarClient(api_key=None)
    with pytest.raises(PharmVarAuthError):
        client.list_genes()


@responses.activate
def test_rejected_key_raises_auth_error():
    responses.add(
        responses.GET,
        "https://www.pharmvar.org/api-service/genes",
        json={"errorMessage": "API Key is invalid or missing.", "errorCode": 401},
        status=401,
    )
    client = PharmVarClient(api_key="bad-key")
    with pytest.raises(PharmVarAuthError):
        client.list_genes()


@responses.activate
def test_valid_key_returns_genes():
    responses.add(
        responses.GET,
        "https://www.pharmvar.org/api-service/genes",
        json=[{"geneSymbol": "CYP2D6"}, {"geneSymbol": "CYP2C19"}],
        status=200,
    )
    client = PharmVarClient(api_key="good-key")
    genes = client.list_genes()
    assert {g["geneSymbol"] for g in genes} == {"CYP2D6", "CYP2C19"}


@responses.activate
def test_get_gene_alleles():
    responses.add(
        responses.GET,
        "https://www.pharmvar.org/api-service/genes/CYP2D6",
        json={"alleles": [{"name": "*4"}]},
        status=200,
    )
    client = PharmVarClient(api_key="good-key")
    result = client.get_gene_alleles("CYP2D6")
    assert result["alleles"][0]["name"] == "*4"
