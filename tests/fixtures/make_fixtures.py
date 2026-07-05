"""Generate small synthetic VCF fixtures used by the test suite.

These are NOT real patient data -- genotypes are synthetically constructed
to exercise specific star-allele-defining positions listed in
``caller.BOOTSTRAP_ALLELE_DEFINITIONS``. Run this script to regenerate the
fixture VCFs if the bootstrap allele table changes.
"""

from __future__ import annotations

from pathlib import Path

FIXTURES_DIR = Path(__file__).parent

HEADER = """##fileformat=VCFv4.2
##source=pharminform-synthetic-fixture
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr6,length=170805979>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""

# (chrom, pos, rsid, ref, alt, genotype, depth) -- positions are placeholders
# within the correct gene region, not the real genomic coordinates of each
# rsID (those require a PharmVar/dbSNP lookup at import time). This is
# sufficient for exercising the rule-based matcher, which keys off rsID.
VARIANTS = [
    # CYP2C19*2 het -> intermediate metabolizer candidate
    ("chr10", 94781859, "rs4244285", "G", "A", "0/1", 45),
    # CYP2C9*3 het
    ("chr10", 94981296, "rs1057910", "A", "C", "0/1", 38),
    # TPMT*3C hom -> poor metabolizer
    ("chr6", 18143955, "rs1142345", "A", "G", "1/1", 30),
    # A low-depth CYP2D6*4 call -> should be filtered out by min_depth
    ("chr22", 42126611, "rs3892097", "G", "A", "0/1", 4),
]


def write_vcf(path: Path) -> None:
    lines = [HEADER]
    for chrom, pos, rsid, ref, alt, gt, dp in VARIANTS:
        lines.append(
            f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t99\tPASS\tDP={dp}\tGT:DP\t{gt}:{dp}\n"
        )
    path.write_text("".join(lines))


if __name__ == "__main__":
    write_vcf(FIXTURES_DIR / "synthetic_sample.vcf")
    print("wrote", FIXTURES_DIR / "synthetic_sample.vcf")
