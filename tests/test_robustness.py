"""Inputs that used to give a confident-looking wrong answer: VCFs with no
rsIDs, on GRCh37, with several samples or none, and defining sites that
were not genotyped or were dropped for depth."""

import pytest
from typer.testing import CliRunner

from pharminform import caller
from pharminform.cli import app
from pharminform.pipeline import run_pipeline

GRCH38_CONTIGS = [
    "##contig=<ID=chr1,length=248956422>",
    "##contig=<ID=chr2,length=242193529>",
    "##contig=<ID=chr6,length=170805979>",
    "##contig=<ID=chr10,length=133797422>",
]
GRCH37_CONTIGS = ["##contig=<ID=1,length=249250621>", "##contig=<ID=10,length=135534747>"]


def write_vcf(tmp_path, records, samples=("S1",), contigs=GRCH38_CONTIGS, name="in.vcf"):
    lines = [
        "##fileformat=VCFv4.2",
        *contigs,
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", *(["FORMAT"] if samples else []), *samples]),
        *("\t".join(r) for r in records),
    ]
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


def rec(chrom, pos, ref, alt, *genotypes, rsid=".", fmt="GT:DP"):
    return [chrom, str(pos), rsid, ref, alt, "50", "PASS", ".", *([fmt] if genotypes else []), *genotypes]


def row(ctx, gene):
    return next(r for r in ctx.gene_results if r.gene == gene)


def test_defining_variants_found_by_position_without_rsids(tmp_path):
    vcf = write_vcf(tmp_path, [
        rec("chr1", 97450058, "C", "T", "0/1:40"),         # DPYD *2A
        rec("chr6", 18130687, "T", "C", "1/1:40"),         # TPMT *3C
        rec("chr2", 233760233, "C", "CAT", "0/1:40"),      # UGT1A1 *28 (TA7)
        rec("chr10", 94781859, "G", "A", "0/1:40"),        # CYP2C19 *2
    ])
    ctx = run_pipeline(vcf, min_depth=10)
    assert row(ctx, "DPYD").diplotype == "*2A/*1"
    assert row(ctx, "TPMT").diplotype == "*3C/*3C"
    assert row(ctx, "UGT1A1").diplotype == "*28/*1"
    assert row(ctx, "CYP2C19").diplotype == "*2/*1"


def test_grch37_vcf_is_refused(tmp_path):
    vcf = write_vcf(tmp_path, [rec("1", 97915614, "C", "T", "0/1:40")], contigs=GRCH37_CONTIGS)
    with pytest.raises(caller.VCFIngestError, match="GRCh37"):
        run_pipeline(vcf)


def test_multi_sample_vcf_needs_a_sample(tmp_path):
    vcf = write_vcf(tmp_path, [rec("chr1", 97450058, "C", "T", "0/0:40", "0/1:40")], samples=("A", "B"))
    with pytest.raises(caller.VCFIngestError, match="2 samples"):
        run_pipeline(vcf)
    assert row(run_pipeline(vcf, min_depth=10, sample="B"), "DPYD").diplotype == "*2A/*1"
    assert row(run_pipeline(vcf, min_depth=10, sample="A"), "DPYD").diplotype == "*1/*1"
    with pytest.raises(caller.VCFIngestError, match="not in"):
        run_pipeline(vcf, sample="C")


def test_sites_only_vcf_is_refused(tmp_path):
    vcf = write_vcf(tmp_path, [rec("chr1", 97450058, "C", "T")], samples=())
    with pytest.raises(caller.VCFIngestError, match="no sample columns"):
        run_pipeline(vcf)


def test_no_call_at_a_defining_site_is_indeterminate_not_normal(tmp_path):
    # A missing genotype at DPYD *2A must not read as *1/*1 normal, which
    # would mean a full fluoropyrimidine dose.
    vcf = write_vcf(tmp_path, [rec("chr1", 97450058, "C", "T", "./.:40")])
    dpyd = row(run_pipeline(vcf, min_depth=10), "DPYD")
    assert dpyd.phenotype == "indeterminate"
    assert any("no genotype call" in r and "rs3918290" in r for r in dpyd.confidence_reasons)


def test_low_depth_defining_variant_is_indeterminate_not_normal(tmp_path):
    vcf = write_vcf(tmp_path, [rec("chr1", 97450058, "C", "T", "0/1:4")])
    ctx = run_pipeline(vcf, min_depth=10)
    dpyd = row(ctx, "DPYD")
    assert (dpyd.diplotype, dpyd.phenotype) == ("*1/*1", "indeterminate")
    assert not any("DPYD" in r["phenotypes"] for r in ctx.cpic_recommendations)


def test_absent_sites_and_unstated_build_are_reported(tmp_path):
    vcf = write_vcf(tmp_path, [rec("chr6", 18130687, "T", "C", "0/1:40")], contigs=[])
    ctx = run_pipeline(vcf, min_depth=10)
    tpmt = row(ctx, "TPMT")
    assert tpmt.phenotype == "intermediate metabolizer"
    assert any("no VCF record at 2 of 3" in r for r in tpmt.confidence_reasons)
    assert any("reference build not stated" in r for r in tpmt.confidence_reasons)


def test_cli_reports_ingest_errors_cleanly(tmp_path):
    vcf = write_vcf(tmp_path, [rec("1", 97915614, "C", "T", "0/1:40")], contigs=GRCH37_CONTIGS)
    result = CliRunner().invoke(app, ["report", "--vcf", str(vcf), "--out", str(tmp_path / "r.html")])
    assert result.exit_code == 2
    assert "GRCh37" in result.output
    assert "Traceback" not in result.output
