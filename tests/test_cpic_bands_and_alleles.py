"""Regression tests for CPIC phenotype bands, allele coverage and gene
regions in the oncology genes (DPYD, UGT1A1, TPMT).

Each case below was wrong before the fix on this branch: DPYD and UGT1A1
phenotypes came out one tier too mild, DPYD *13 and UGT1A1 *28 were not
detected, TPMT *3B + *3C was always called *3B/*3C (poor metabolizer), and
the GRCh38 DPYD and CYP2C19 regions missed *2A, c.2846A>T and *17.
"""

import pytest
from tests.getrm_harness import GRCH37_GENE_REGIONS, compare_phenotypes

from pharminform import caller
from pharminform.caller import VariantCall, call_diplotype_from_vcf, normalize_diplotype
from pharminform.pharmgkb_client import diplotype_to_phenotype


def _functions(gene):
    return {name: info["function"] for name, info in caller.BOOTSTRAP_ALLELE_DEFINITIONS[gene].items()}


def _variant(gene, rsid, genotype="0/1", ref="A", alt="G", alt_index=1):
    return VariantCall(
        gene=gene, chrom="chrN", pos=1, ref=ref, alt=alt, genotype=genotype,
        depth=50, qual=99.0, filter="PASS", rsid=rsid, alt_index=alt_index,
    )


# (gene, diplotype, expected phenotype) -- CPIC DPYD (PMID 29152729),
# TPMT/NUDT15 (PMID 30447069) and UGT1A1 (PMID 26417955) guidance, plus
# the CYP genes below.
PHENOTYPE_CASES = [
    ("DPYD", "*1/*1", "normal metabolizer"),
    ("DPYD", "c.2846A>T/*1", "intermediate metabolizer"),
    ("DPYD", "c.1236G>A(HapB3)/*1", "intermediate metabolizer"),
    ("DPYD", "*2A/*1", "intermediate metabolizer"),
    ("DPYD", "*13/*1", "intermediate metabolizer"),
    ("DPYD", "c.2846A>T/c.2846A>T", "intermediate metabolizer"),
    ("DPYD", "*2A/c.2846A>T", "poor metabolizer"),
    ("DPYD", "*2A/*2A", "poor metabolizer"),
    ("UGT1A1", "*1/*1", "normal metabolizer"),
    ("UGT1A1", "*36/*1", "normal metabolizer"),
    ("UGT1A1", "*6/*1", "intermediate metabolizer"),
    ("UGT1A1", "*28/*1", "intermediate metabolizer"),
    ("UGT1A1", "*6/*6", "poor metabolizer"),
    ("UGT1A1", "*28/*28", "poor metabolizer"),
    ("UGT1A1", "*28/*6", "poor metabolizer"),
    ("TPMT", "*1/*1", "normal metabolizer"),
    ("TPMT", "*3A/*1", "intermediate metabolizer"),
    ("TPMT", "*3B/*3C", "poor metabolizer"),
    ("TPMT", "*3A/*3A", "poor metabolizer"),
    # CYP2C9 (PMID 32189324), CYP2C19 (PMID 35034351), CYP2D6 (PMID 31647186).
    ("CYP2C9", "*1/*1", "normal metabolizer"),
    ("CYP2C9", "*1/*2", "intermediate metabolizer"),
    ("CYP2C9", "*1/*3", "intermediate metabolizer"),
    ("CYP2C9", "*2/*2", "intermediate metabolizer"),
    ("CYP2C9", "*2/*3", "poor metabolizer"),
    ("CYP2C9", "*3/*3", "poor metabolizer"),
    ("CYP2C19", "*1/*1", "normal metabolizer"),
    ("CYP2C19", "*1/*2", "intermediate metabolizer"),
    ("CYP2C19", "*2/*17", "intermediate metabolizer"),
    ("CYP2C19", "*3/*17", "intermediate metabolizer"),
    ("CYP2C19", "*2/*2", "poor metabolizer"),
    ("CYP2C19", "*1/*17", "rapid metabolizer"),
    ("CYP2C19", "*17/*17", "ultrarapid metabolizer"),
    ("CYP2D6", "*1/*1", "normal metabolizer"),
    ("CYP2D6", "*4/*1", "intermediate metabolizer"),
    ("CYP2D6", "*4/*4", "poor metabolizer"),
]


@pytest.mark.parametrize("gene,diplotype,expected", PHENOTYPE_CASES, ids=[f"{g}-{d}" for g, d, _ in PHENOTYPE_CASES])
def test_cpic_phenotype_bands(gene, diplotype, expected):
    phenotype, _ = diplotype_to_phenotype(gene, diplotype, _functions(gene))
    assert phenotype == expected


def test_dpyd_star13_detected():
    call = call_diplotype_from_vcf("DPYD", [_variant("DPYD", "rs55886062", ref="A", alt="C")])
    assert call.diplotype == "*13/*1"


def test_dpyd_hapb3_detected_from_causal_intronic_variant():
    call = call_diplotype_from_vcf("DPYD", [_variant("DPYD", "rs75017182", ref="G", alt="C")])
    assert call.diplotype == "c.1236G>A(HapB3)/*1"


def test_dpyd_hapb3_two_markers_count_as_one_allele():
    variants = [_variant("DPYD", "rs56038477", ref="C", alt="T"), _variant("DPYD", "rs75017182", ref="G", alt="C")]
    call = call_diplotype_from_vcf("DPYD", variants)
    assert call.diplotype == "c.1236G>A(HapB3)/*1"


def test_ugt1a1_star28_as_recorded_in_1000_genomes():
    # 1000 Genomes phase 3 records the TA7 allele as rs34983651 C>CAT.
    call = call_diplotype_from_vcf("UGT1A1", [_variant("UGT1A1", "rs34983651", ref="C", alt="CAT")])
    assert call.diplotype == "*28/*1"
    assert call.ambiguous_candidates == []


def test_ugt1a1_star28_homozygous():
    call = call_diplotype_from_vcf("UGT1A1", [_variant("UGT1A1", "rs3064744", "1/1", ref="CAT", alt="CATAT")])
    assert call.diplotype == "*28/*28"


def test_ugt1a1_repeat_length_distinguishes_star36_from_star28():
    call = call_diplotype_from_vcf("UGT1A1", [_variant("UGT1A1", "rs3064744", ref="CAT", alt="C")])
    assert call.diplotype == "*36/*1"


def test_ugt1a1_star80_proxy_reports_star28_with_warning():
    call = call_diplotype_from_vcf("UGT1A1", [_variant("UGT1A1", "rs887829", ref="C", alt="T")])
    assert call.diplotype == "*28/*1"
    assert any("proxy" in note for note in call.ambiguous_candidates)


def test_tpmt_3b_and_3c_unphased_called_3a_with_warning():
    variants = [_variant("TPMT", "rs1800460", ref="C", alt="T"), _variant("TPMT", "rs1142345", ref="T", alt="C")]
    call = call_diplotype_from_vcf("TPMT", variants)
    assert call.diplotype == "*3A/*1"
    assert any("*3B/*3C in trans" in note for note in call.ambiguous_candidates)


def test_tpmt_3b_and_3c_phased_in_cis_is_3a():
    variants = [_variant("TPMT", "rs1800460", "1|0"), _variant("TPMT", "rs1142345", "1|0")]
    call = call_diplotype_from_vcf("TPMT", variants)
    assert call.diplotype == "*3A/*1"
    assert call.ambiguous_candidates == []


def test_tpmt_3b_and_3c_phased_in_trans_is_3b_3c():
    variants = [_variant("TPMT", "rs1800460", "1|0"), _variant("TPMT", "rs1142345", "0|1")]
    call = call_diplotype_from_vcf("TPMT", variants)
    assert call.diplotype == "*3B/*3C"
    assert call.ambiguous_candidates == []


def test_tpmt_3c_homozygous_with_3b_het_is_3a_3c():
    variants = [_variant("TPMT", "rs1800460", "0/1"), _variant("TPMT", "rs1142345", "1/1")]
    call = call_diplotype_from_vcf("TPMT", variants)
    assert call.diplotype == "*3A/*3C"


def test_multiallelic_site_only_counts_carried_alt():
    # GT 0/2: the sample carries ALT 2, so a call for ALT 1 is not evidence.
    call = call_diplotype_from_vcf("CYP2C19", [_variant("CYP2C19", "rs4244285", "0/2", alt_index=1)])
    assert call.diplotype == "*1/*1"


def test_semicolon_joined_ids_match():
    call = call_diplotype_from_vcf("DPYD", [_variant("DPYD", "rs3918290;rs199469548", ref="C", alt="T")])
    assert call.diplotype == "*2A/*1"


def test_normalize_diplotype():
    assert normalize_diplotype("*2A/*1") == "*1/*2A"
    assert normalize_diplotype("*41/*4") == "*4/*41"
    assert normalize_diplotype("c.2846A>T/*1") == "*1/c.2846A>T"


# Ensembl positions of each defining variant (GRCh38, GRCh37).
DEFINING_VARIANT_POSITIONS = {
    "CYP2D6": {"rs35742686": (42128242, None), "rs3892097": (42128945, None), "rs1065852": (42130692, None),
               "rs28371706": (42129770, None), "rs28371725": (42127803, None)},
    "CYP2C19": {"rs4244285": (94781859, 96541616), "rs4986893": (94780653, 96540410),
                "rs12248560": (94761900, 96521657)},
    "CYP2C9": {"rs1799853": (94942290, 96702047), "rs1057910": (94981296, 96741053)},
    "DPYD": {"rs3918290": (97450058, 97915614), "rs55886062": (97515787, 97981343),
             "rs67376798": (97082391, 97547947), "rs56038477": (97573863, 98039419),
             "rs75017182": (97579893, 98045449)},
    "UGT1A1": {"rs4148323": (233760498, 234669144), "rs3064744": (233760234, 234668880),
               "rs887829": (233759924, 234668570)},
    "TPMT": {"rs1800462": (18143724, 18143955), "rs1800460": (18138997, 18139228),
             "rs1142345": (18130687, 18130918)},
}

REGION_CASES = [(g, r, p) for g, rs in DEFINING_VARIANT_POSITIONS.items() for r, p in rs.items()]


@pytest.mark.parametrize("gene,rsid,positions", REGION_CASES, ids=[f"{g}-{r}" for g, r, _ in REGION_CASES])
def test_defining_variants_inside_gene_regions(gene, rsid, positions):
    grch38, grch37 = positions
    _, span = caller.GENE_REGIONS_GRCH38[gene].split(":")
    start, end = (int(x) for x in span.split("-"))
    assert start <= grch38 <= end
    if grch37 is not None:
        _, start37, end37 = GRCH37_GENE_REGIONS[gene]
        assert start37 <= grch37 <= end37


def test_load_variants_finds_grch38_dpyd_and_cyp2c19_star17(tmp_path):
    vcf = tmp_path / "grch38.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=248956422>\n"
        "##contig=<ID=chr10,length=133797422>\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "chr1\t97082391\trs67376798\tT\tA\t99\tPASS\tDP=40\tGT\t0/1\n"
        "chr1\t97450058\trs3918290\tC\tT\t99\tPASS\tDP=40\tGT\t0/1\n"
        "chr10\t94761900\trs12248560\tC\tT\t99\tPASS\tDP=40\tGT\t0/1\n"
    )
    variants = caller.load_variants(vcf, genes=["DPYD", "CYP2C19"])
    assert {v.rsid for v in variants["DPYD"]} == {"rs67376798", "rs3918290"}
    assert call_diplotype_from_vcf("CYP2C19", variants["CYP2C19"]).diplotype == "*17/*1"
    dpyd = call_diplotype_from_vcf("DPYD", variants["DPYD"])
    assert dpyd.diplotype == "*2A/c.2846A>T"
    phenotype, _ = diplotype_to_phenotype("DPYD", dpyd.diplotype, _functions("DPYD"))
    assert phenotype == "poor metabolizer"


def test_harness_phenotype_comparison():
    assert compare_phenotypes("DPYD", "*1/*2A", "*2A/*1") == (
        "intermediate metabolizer", "intermediate metabolizer", "match",
    )
    # CYP2D6 *5 (whole-gene deletion) is not in the bootstrap table.
    assert compare_phenotypes("CYP2D6", "*1/*5", "*1/*1")[2] == "not_assessable"


def test_low_depth_defining_variant_is_reported_not_silent():
    # The synthetic fixture's CYP2D6 *4 record has DP=4: it must not be used
    # for the call, but the confidence result must say it was there.
    from pathlib import Path

    from pharminform.confidence import assess_confidence

    fixture = Path(__file__).parent / "fixtures" / "synthetic_sample.vcf"
    dropped = {}
    variants = caller.load_variants(fixture, min_depth=10, dropped=dropped)
    assert variants["CYP2D6"] == []
    assert [v.rsid for v in dropped["CYP2D6"]] == ["rs3892097"]
    call = call_diplotype_from_vcf("CYP2D6", variants["CYP2D6"])
    assert call.diplotype == "*1/*1"
    conf = assess_confidence(call, all_gene_variants=variants["CYP2D6"], min_depth=10,
                             dropped_variants=dropped["CYP2D6"])
    assert any("excluded for low depth" in r and "rs3892097" in r for r in conf.reasons)


def test_pipeline_surfaces_low_depth_defining_variant():
    # The fixture's only low-depth defining variant is CYP2D6 *4, and the
    # pipeline no longer calls CYP2D6 from a VCF, so check the confidence
    # step directly on what load_variants dropped.
    from pathlib import Path

    from pharminform.confidence import assess_confidence

    vcf = Path(__file__).parent / "fixtures" / "synthetic_sample.vcf"
    dropped: dict = {}
    variants = caller.load_variants(vcf, min_depth=10, dropped=dropped)
    call = caller.call_diplotype_from_vcf("CYP2D6", variants.get("CYP2D6", []))
    conf = assess_confidence(call, variants.get("CYP2D6", []), min_depth=10, dropped_variants=dropped.get("CYP2D6"))
    assert any("rs3892097" in r for r in conf.reasons)


def test_pipeline_does_not_call_cyp2d6_from_vcf():
    from pathlib import Path

    from pharminform.pipeline import run_pipeline

    ctx = run_pipeline(Path(__file__).parent / "fixtures" / "synthetic_sample.vcf", min_depth=10)
    row = next(r for r in ctx.gene_results if r.gene == "CYP2D6")
    assert (row.diplotype, row.phenotype, row.confidence_level) == ("not called", None, "low")
    assert "Aldy" in row.confidence_reasons[0]
    # Amitriptyline (CYP2D6 + CYP2C19) can still match on CYP2C19, but must
    # say it assumes a normal CYP2D6 rather than silently picking one.
    for rec in ctx.cpic_recommendations:
        if "CYP2D6" in rec["phenotypes"]:
            assert rec["assumed_phenotypes"] == {"CYP2D6": "Normal Metabolizer"}


def test_cyp2d6_activity_values_follow_cpic():
    # CPIC gives *10 and *41 an activity value of 0.25, not the generic 0.5.
    functions = {"*4": "no function", "*10": "decreased function", "*41": "decreased function"}
    assert diplotype_to_phenotype("CYP2D6", "*1/*10", functions) == ("normal metabolizer", 1.25)
    assert diplotype_to_phenotype("CYP2D6", "*10/*41", functions) == ("intermediate metabolizer", 0.5)
    assert diplotype_to_phenotype("CYP2D6", "*4/*10", functions) == ("intermediate metabolizer", 0.25)


def test_ugt1a1_star36_does_not_offset_a_decreased_allele():
    functions = {n: i["function"] for n, i in caller.BOOTSTRAP_ALLELE_DEFINITIONS["UGT1A1"].items()}
    assert functions["*36"] == "increased function"
    assert diplotype_to_phenotype("UGT1A1", "*28/*36", functions)[0] == "intermediate metabolizer"
    assert diplotype_to_phenotype("UGT1A1", "*1/*36", functions)[0] == "normal metabolizer"
    assert diplotype_to_phenotype("UGT1A1", "*36/*36", functions)[0] == "normal metabolizer"


def test_ugt1a1_proxy_only_star28_is_indeterminate(tmp_path):
    # CPIC lists *80 without the TA repeat itself as indeterminate.
    from pharminform.pipeline import run_pipeline

    v = caller.VariantCall(
        gene="UGT1A1", chrom="chr2", pos=233759924, ref="C", alt="T",
        genotype="0/1", depth=40, qual=50.0, filter="PASS", rsid="rs887829",
    )
    call = caller.call_diplotype_from_vcf("UGT1A1", [v])
    assert call.proxy_only_alleles == ["*28"]

    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr2>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "chr2\t233759924\trs887829\tC\tT\t50\tPASS\t.\tGT:DP\t0/1:40",
    ]
    vcf = tmp_path / "proxy.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    ctx = run_pipeline(vcf, min_depth=10, genes=["UGT1A1"])
    row = ctx.gene_results[0]
    assert (row.diplotype, row.phenotype) == ("*28/*1", "indeterminate")
    assert not ctx.cpic_recommendations


def test_cyp2d6_duplication_is_ultrarapid_not_rapid():
    # CYP2D6 has no "rapid" tier in CPIC; AS 3 (e.g. *1/*1x2) is ultrarapid.
    phenotype, score = diplotype_to_phenotype("CYP2D6", "*1/*1x2", {"*1x2": "increased function"})
    assert (phenotype, score) == ("ultrarapid metabolizer", 3.0)
