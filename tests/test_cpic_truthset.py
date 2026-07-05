"""Pytest suite against a CPIC/PharmVar-consistent star allele truth set.

Each row below is a (gene, rsID, star allele, function) truth-set entry
consistent with CPIC gene-specific allele function tables and PharmVar
allele definitions -- the same rsIDs already embedded in
``caller.BOOTSTRAP_ALLELE_DEFINITIONS`` (see that module's docstring for
the "bootstrap fixture, not a clinical source" caveat). This suite exists
to pin the calling pipeline's behavior against >=10 named, individually
identifiable variants rather than only exercising a couple of them
incidentally through other tests.
"""

import pytest

from pharminform.caller import VariantCall, call_diplotype_from_vcf

TRUTH_SET = [
    # (gene, rsid, ref, alt, star_allele, zygosity, expected_diplotype)
    ("CYP2D6", "rs3892097", "G", "A", "*4", "het", "*4/*1"),
    ("CYP2D6", "rs1065852", "G", "A", "*10", "het", "*10/*1"),
    ("CYP2D6", "rs28371706", "C", "T", "*17", "het", "*17/*1"),
    ("CYP2D6", "rs28371725", "C", "T", "*41", "het", "*41/*1"),
    ("CYP2D6", "rs35742686", "A", "AGT", "*3", "hom", "*3/*3"),
    ("CYP2C19", "rs4244285", "G", "A", "*2", "het", "*2/*1"),
    ("CYP2C19", "rs4986893", "G", "A", "*3", "hom", "*3/*3"),
    ("CYP2C19", "rs12248560", "C", "T", "*17", "het", "*17/*1"),
    ("CYP2C9", "rs1799853", "C", "T", "*2", "het", "*2/*1"),
    ("CYP2C9", "rs1057910", "A", "C", "*3", "het", "*3/*1"),
    ("DPYD", "rs3918290", "C", "T", "*2A", "het", "*2A/*1"),
    ("DPYD", "rs67376798", "A", "T", "c.2846A>T", "het", "c.2846A>T/*1"),
    ("TPMT", "rs1800462", "C", "G", "*2", "hom", "*2/*2"),
    ("TPMT", "rs1142345", "A", "G", "*3C", "het", "*3C/*1"),
    ("UGT1A1", "rs4148323", "G", "A", "*6", "het", "*6/*1"),
]


@pytest.mark.parametrize(
    "gene,rsid,ref,alt,star_allele,zygosity,expected", TRUTH_SET,
    ids=[f"{g}-{a}" for g, _, _, _, a, _, _ in TRUTH_SET],
)
def test_truthset_variant_calls_expected_star_allele(gene, rsid, ref, alt, star_allele, zygosity, expected):
    genotype = "1/1" if zygosity == "hom" else "0/1"
    variant = VariantCall(
        gene=gene, chrom="chrN", pos=1, ref=ref, alt=alt,
        genotype=genotype, depth=50, qual=99.0, filter="PASS", rsid=rsid,
    )
    call = call_diplotype_from_vcf(gene, [variant])
    assert call.diplotype == expected


def test_truthset_covers_at_least_ten_variants():
    assert len(TRUTH_SET) >= 10


def test_truthset_covers_all_six_target_genes():
    assert {row[0] for row in TRUTH_SET} == {"CYP2D6", "CYP2C19", "CYP2C9", "DPYD", "TPMT", "UGT1A1"}
