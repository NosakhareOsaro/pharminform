from pharminform.caller import DiplotypeCall, VariantCall
from pharminform.confidence import assess_confidence


def make_variant(gene="CYP2C19", chrom="chr10", pos=94781859, rsid="rs4244285", depth=45):
    return VariantCall(
        gene=gene, chrom=chrom, pos=pos, ref="G", alt="A", genotype="0/1",
        depth=depth, qual=99.0, filter="PASS", rsid=rsid,
    )


def test_high_confidence_aldy_call_no_issues():
    call = DiplotypeCall(gene="CYP2C19", diplotype="*2/*1", method="aldy")
    result = assess_confidence(call)
    assert result.level == "high"
    assert result.reasons == []


def test_rule_based_vcf_call_on_cyp2d6_is_low_confidence():
    call = DiplotypeCall(gene="CYP2D6", diplotype="*4/*1", method="rule_based_vcf")
    result = assess_confidence(call)
    assert result.is_low
    assert any("structural variation" in r for r in result.reasons)


def test_rule_based_vcf_call_on_non_structural_gene_still_flagged_but_no_cn_reason():
    call = DiplotypeCall(gene="CYP2C19", diplotype="*2/*1", method="rule_based_vcf")
    result = assess_confidence(call)
    assert result.is_low
    assert not any("structural variation" in r for r in result.reasons)


def test_low_depth_variant_flagged():
    v = make_variant(depth=5)
    call = DiplotypeCall(gene="CYP2C19", diplotype="*2/*1", method="aldy", matched_variants=[v])
    result = assess_confidence(call, min_depth=20)
    assert result.is_low
    assert any("low coverage" in r for r in result.reasons)


def test_phasing_ambiguity_flagged():
    call = DiplotypeCall(
        gene="CYP2D6", diplotype="*4/*10", method="rule_based_vcf",
        ambiguous_candidates=["phase unresolved among: *10, *4"],
    )
    result = assess_confidence(call)
    assert result.is_low
    assert any("phasing ambiguity" in r for r in result.reasons)


def test_novel_variant_not_in_definitions_flagged():
    known = make_variant()
    novel = make_variant(rsid="rs999999999", pos=94781900)
    call = DiplotypeCall(gene="CYP2C19", diplotype="*2/*1", method="aldy", matched_variants=[known])
    result = assess_confidence(call, all_gene_variants=[known, novel])
    assert result.is_low
    assert any("not in the allele definition table" in r for r in result.reasons)


def test_unnamed_variant_flagged_as_possible_novel():
    unnamed = make_variant(rsid=None)
    call = DiplotypeCall(gene="CYP2C19", diplotype="*1/*1", method="aldy")
    result = assess_confidence(call, all_gene_variants=[unnamed])
    assert result.is_low
    assert any("no rsID" in r for r in result.reasons)
