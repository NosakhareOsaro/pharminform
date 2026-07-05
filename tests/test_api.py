from pathlib import Path

from fastapi.testclient import TestClient

from pharminform.api import app

FIXTURE_VCF = Path(__file__).parent / "fixtures" / "synthetic_sample.vcf"

client = TestClient(app)


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_create_report_json():
    with open(FIXTURE_VCF, "rb") as f:
        resp = client.post(
            "/report",
            files={"vcf": ("synthetic_sample.vcf", f, "text/plain")},
            params={"min_depth": 10, "patient_id": "API-TEST"},
        )
    assert resp.status_code == 200
    body = resp.json()
    assert body["patient_id"] == "API-TEST"
    genes = {g["gene"] for g in body["gene_results"]}
    assert genes == {"CYP2D6", "CYP2C19", "CYP2C9", "DPYD", "UGT1A1", "TPMT"}
    assert any(c["drug"] == "azathioprine" for c in body["composite_risk"])


def test_create_report_html():
    with open(FIXTURE_VCF, "rb") as f:
        resp = client.post(
            "/report/html",
            files={"vcf": ("synthetic_sample.vcf", f, "text/plain")},
            params={"min_depth": 10},
        )
    assert resp.status_code == 200
    assert "PharmInform" in resp.text


def test_rejects_non_vcf_file():
    resp = client.post(
        "/report",
        files={"vcf": ("notavcf.txt", b"hello", "text/plain")},
    )
    assert resp.status_code == 400
