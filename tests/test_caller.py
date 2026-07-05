from pathlib import Path

import pytest

from pharminform import caller

FIXTURE_VCF = Path(__file__).parent / "fixtures" / "synthetic_sample.vcf"


def test_load_variants_filters_low_depth():
    variants = caller.load_variants(FIXTURE_VCF, min_depth=10)
    cyp2d6_positions = variants["CYP2D6"]
    # the synthetic CYP2D6*4 record has DP=4 and must be filtered out
    assert cyp2d6_positions == []


def test_load_variants_keeps_passing_records():
    variants = caller.load_variants(FIXTURE_VCF, min_depth=10)
    cyp2c19 = variants["CYP2C19"]
    assert len(cyp2c19) == 1
    assert cyp2c19[0].rsid == "rs4244285"
    assert cyp2c19[0].depth == 45


def test_call_diplotype_from_vcf_heterozygous():
    variants = caller.load_variants(FIXTURE_VCF, min_depth=10)
    call = caller.call_diplotype_from_vcf("CYP2C19", variants["CYP2C19"])
    assert call.diplotype == "*2/*1"
    assert call.method == "rule_based_vcf"


def test_call_diplotype_from_vcf_homozygous():
    variants = caller.load_variants(FIXTURE_VCF, min_depth=10)
    call = caller.call_diplotype_from_vcf("TPMT", variants["TPMT"])
    assert call.diplotype == "*3C/*3C"


def test_call_diplotype_from_vcf_defaults_to_wildtype():
    call = caller.call_diplotype_from_vcf("UGT1A1", [])
    assert call.diplotype == "*1/*1"
    assert call.ambiguous_candidates == []


def test_vcf_not_found_raises():
    with pytest.raises(caller.VCFIngestError):
        caller.load_variants("does_not_exist.vcf")


def test_aldy_not_available_raises_clear_error(tmp_path, monkeypatch):
    monkeypatch.setattr(caller, "aldy_available", lambda: False)
    fake_bam = tmp_path / "sample.bam"
    fake_bam.write_bytes(b"not a real bam")
    with pytest.raises(caller.AldyNotAvailable):
        caller.run_aldy(fake_bam, "CYP2D6")
