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


def test_load_variants_excludes_homozygous_reference_calls():
    # Real bug this guards against: a gene-region VCF slice lists every
    # polymorphic site in the population, so a homozygous-reference (0/0)
    # call at a defining position must NOT be treated as "this sample
    # carries the variant" -- caught by the real GeT-RM validation run,
    # invisible against a hand-authored fixture that never included a 0/0
    # row before.
    variants = caller.load_variants(FIXTURE_VCF, min_depth=10)
    cyp2c9 = variants["CYP2C9"]
    assert all(v.rsid != "rs1799853" for v in cyp2c9)


def test_call_diplotype_from_vcf_ignores_homozygous_reference_variant():
    hom_ref = caller.VariantCall(
        gene="CYP2C9", chrom="chr10", pos=1, ref="C", alt="T",
        genotype="0/0", depth=40, qual=99.0, filter="PASS", rsid="rs1799853",
    )
    call = caller.call_diplotype_from_vcf("CYP2C9", [hom_ref])
    assert call.diplotype == "*1/*1"


def test_vcf_not_found_raises():
    with pytest.raises(caller.VCFIngestError):
        caller.load_variants("does_not_exist.vcf")


def test_aldy_not_available_raises_clear_error(tmp_path, monkeypatch):
    monkeypatch.setattr(caller, "aldy_available", lambda: False)
    fake_bam = tmp_path / "sample.bam"
    fake_bam.write_bytes(b"not a real bam")
    with pytest.raises(caller.AldyNotAvailable):
        caller.run_aldy(fake_bam, "CYP2D6")
