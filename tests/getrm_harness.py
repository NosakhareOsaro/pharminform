"""GeT-RM validation harness.

Compares PharmInform's star allele calls against published GeT-RM
(Genetic Testing Reference Materials Coordination Program) consensus
genotypes for public Coriell reference samples, using real sequencing
data from the 1000 Genomes Project where a GeT-RM Coriell sample is
also a 1000 Genomes participant.

Two real, live findings from getting this running end to end:

1. Pulling 1000 Genomes-hosted sequencing data directly into cyvcf2 over
   HTTPS fails in this environment -- ``cyvcf2.VCF(url)`` raises
   ``OSError: Failed to open file ... No such file or directory``, which
   traces to the bundled htslib not being built with remote (libcurl) URL
   support. ``bcftools`` (installed separately, confirmed working via
   ``bcftools --version``) is used instead for the remote fetch step,
   with the result then parsed locally by cyvcf2 -- both are exercised
   for real here.

2. The AWS Open Data mirror's phase 3 VCF filenames use a "v5a" suffix
   for every chromosome (confirmed via the bucket's own listing API),
   not "v5b" as an earlier version of this file assumed by analogy with
   the EBI FTP mirror's naming for some chromosomes -- that assumption
   was never actually verified against this specific mirror and was
   wrong; fixed here.

See validation/getrm_concordance_report.md for the full writeup,
including which of the 11 bundled CYP2D6 consensus samples are and are
not present in the 1000 Genomes panel, and the real concordance results.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import cyvcf2

from pharminform import caller

CONSENSUS_PATH = Path(__file__).parent.parent / "src" / "pharminform" / "data" / "getrm_consensus_genotypes.json"

# 1000 Genomes phase 3 GRCh37 VCFs, per-chromosome, publicly readable.
# Confirmed live via the AWS Open Data mirror's own bucket listing
# (https://1000genomes.s3.amazonaws.com/?list-type=2&prefix=release/20130502/ALL.chr<N>) --
# every chromosome uses a "v5a" filename suffix on this specific mirror.
ONE_KG_BASE_URL = "https://1000genomes.s3.amazonaws.com/release/20130502"
ONE_KG_CHROM_FILES = {
    str(n): f"{ONE_KG_BASE_URL}/ALL.chr{n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    for n in list(range(1, 23))
}

# The official 1000 Genomes phase 3 sample panel (sample ID, population,
# super-population, sex) -- 2504 samples, confirmed matching the actual
# VCF header's sample columns exactly (cross-checked both sources).
ONE_KG_PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)

# Approximate GRCh37 region for CYP2D6/CYP2D7 on chr22, wide enough to
# catch the defining SNVs in caller.BOOTSTRAP_ALLELE_DEFINITIONS. GRCh37,
# NOT GRCh38 -- 1000 Genomes phase 3 is GRCh37, unlike
# caller.GENE_REGIONS_GRCH38. This only matters for the initial region
# slice; the actual diplotype match is by rsID, which is build-independent.
GRCH37_GENE_REGIONS = {
    "CYP2D6": ("22", 42_500_000, 42_550_000),
}


class SequencingDataUnavailable(RuntimeError):
    """Raised when the sequencing data needed for a real concordance check
    cannot be fetched in the current environment -- never caught silently."""


@dataclass
class ConcordanceOutcome:
    coriell_id: str
    gene: str
    consensus_diplotype: str
    called_diplotype: str | None
    status: str  # "match", "mismatch", "skipped_no_data"
    detail: str = ""
    n_variants_considered: int = 0


def load_consensus_table(path: Path | str = CONSENSUS_PATH) -> list[dict]:
    data = json.loads(Path(path).read_text())
    return data["samples"]


def bcftools_available() -> bool:
    return shutil.which("bcftools") is not None


def fetch_1000genomes_panel_samples() -> set[str]:
    """Fetch the real 1000 Genomes phase 3 sample panel and return the set
    of sample IDs. Raises SequencingDataUnavailable on any network failure
    -- never returns a guessed/cached list presented as live data."""
    proc = subprocess.run(["curl", "-s", "-f", ONE_KG_PANEL_URL], capture_output=True, text=True, timeout=60)
    if proc.returncode != 0 or not proc.stdout:
        raise SequencingDataUnavailable(
            f"Could not fetch the 1000 Genomes sample panel from {ONE_KG_PANEL_URL}"
        )
    lines = proc.stdout.strip().splitlines()
    samples = {line.split("\t")[0] for line in lines[1:] if line.strip()}
    return samples


def fetch_1000genomes_region(chrom: str, start: int, end: int, sample: str) -> str:
    """Fetch a small region of a 1000 Genomes phase 3 VCF for one sample via
    bcftools (remote HTTPS + tabix range requests). Requires bcftools on
    PATH. Raises SequencingDataUnavailable with actionable instructions
    otherwise -- this function must never return fabricated variant data.
    """
    if not bcftools_available():
        raise SequencingDataUnavailable(
            "bcftools is not installed/on PATH. Install it (e.g. `brew install "
            "bcftools` or `conda install -c bioconda bcftools`) to run live "
            "GeT-RM concordance checks against 1000 Genomes-hosted sequencing "
            "data, or download the relevant region locally and adapt this "
            "function's input source."
        )
    url = ONE_KG_CHROM_FILES.get(chrom)
    if url is None:
        raise SequencingDataUnavailable(
            f"No indexed 1000 Genomes URL configured for chr{chrom} in "
            "ONE_KG_CHROM_FILES; add it to extend coverage."
        )
    cmd = ["bcftools", "view", "-r", f"{chrom}:{start}-{end}", "-s", sample, url]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        raise SequencingDataUnavailable(
            f"bcftools failed fetching {chrom}:{start}-{end} for {sample}: {proc.stderr.strip()}"
        )
    return proc.stdout


def call_diplotype_from_1000genomes(gene: str, sample: str) -> tuple[caller.DiplotypeCall, int]:
    """Fetch a sample's real 1000 Genomes genotypes for a gene region and
    run them through PharmInform's actual rule-based VCF caller.

    Returns (DiplotypeCall, n_variants_considered). Raises
    SequencingDataUnavailable if the region isn't configured or the fetch
    fails. Genotype build is GRCh37; matching is by rsID (build-independent),
    per GRCH37_GENE_REGIONS' docstring note.
    """
    region = GRCH37_GENE_REGIONS.get(gene)
    if region is None:
        raise SequencingDataUnavailable(
            f"No GRCh37 region configured for {gene} in GRCH37_GENE_REGIONS; "
            "add it to extend real concordance testing beyond CYP2D6."
        )
    chrom, start, end = region
    vcf_text = fetch_1000genomes_region(chrom, start, end, sample)

    tmp_path = Path(f"/tmp/getrm_{sample}_{gene}.vcf")
    tmp_path.write_text(vcf_text)
    try:
        vcf = cyvcf2.VCF(str(tmp_path))
        variants = []
        for rec in vcf:
            if rec.ID is None:
                continue
            gt_str = caller._genotype_string(rec)
            if not caller._has_alt_allele(gt_str):
                continue
            for alt in rec.ALT:
                variants.append(
                    caller.VariantCall(
                        gene=gene, chrom=rec.CHROM, pos=rec.POS, ref=rec.REF, alt=alt,
                        # 1000 Genomes phase 3's integrated call set has no
                        # per-sample DP in FORMAT (population-level
                        # integrated/imputed callset) -- depth genuinely
                        # unavailable here, not a bug.
                        genotype=gt_str, depth=None, qual=rec.QUAL, filter=rec.FILTER, rsid=rec.ID,
                    )
                )
    finally:
        tmp_path.unlink(missing_ok=True)

    call = caller.call_diplotype_from_vcf(gene, variants)
    return call, len(variants)


def run_concordance_check(samples: list[dict] | None = None) -> list[ConcordanceOutcome]:
    """Run a real concordance check for each (sample, gene) pair in the
    consensus table against actual 1000 Genomes sequencing data.

    Samples not present in the 1000 Genomes panel, or blocked by a real
    tooling limitation (no bcftools, no region configured), are reported
    as "skipped_no_data" with the specific reason -- never silently
    marked as a match, and never given a fabricated diplotype.
    """
    samples = samples if samples is not None else load_consensus_table()
    outcomes: list[ConcordanceOutcome] = []

    try:
        panel_samples = fetch_1000genomes_panel_samples()
    except SequencingDataUnavailable as exc:
        return [
            ConcordanceOutcome(
                coriell_id=entry["coriell_id"], gene=entry["gene"],
                consensus_diplotype=entry["consensus_diplotype"], called_diplotype=None,
                status="skipped_no_data", detail=str(exc),
            )
            for entry in samples
        ]

    for entry in samples:
        coriell_id, gene, consensus = entry["coriell_id"], entry["gene"], entry["consensus_diplotype"]

        if coriell_id not in panel_samples:
            outcomes.append(
                ConcordanceOutcome(
                    coriell_id=coriell_id, gene=gene, consensus_diplotype=consensus, called_diplotype=None,
                    status="skipped_no_data",
                    detail="not present in the 1000 Genomes phase 3 sample panel",
                )
            )
            continue

        try:
            call, n_variants = call_diplotype_from_1000genomes(gene, coriell_id)
        except SequencingDataUnavailable as exc:
            outcomes.append(
                ConcordanceOutcome(
                    coriell_id=coriell_id, gene=gene, consensus_diplotype=consensus, called_diplotype=None,
                    status="skipped_no_data", detail=str(exc),
                )
            )
            continue

        status = "match" if call.diplotype == consensus else "mismatch"
        outcomes.append(
            ConcordanceOutcome(
                coriell_id=coriell_id, gene=gene, consensus_diplotype=consensus,
                called_diplotype=call.diplotype, status=status,
                detail="" if status == "match" else "; ".join(call.ambiguous_candidates),
                n_variants_considered=n_variants,
            )
        )

    return outcomes


def summarize(outcomes: list[ConcordanceOutcome]) -> dict:
    total = len(outcomes)
    matches = sum(1 for o in outcomes if o.status == "match")
    mismatches = sum(1 for o in outcomes if o.status == "mismatch")
    skipped = sum(1 for o in outcomes if o.status == "skipped_no_data")
    return {
        "total": total,
        "matches": matches,
        "mismatches": mismatches,
        "skipped_no_data": skipped,
        "concordance_pct_of_attempted": (
            round(100 * matches / (matches + mismatches), 1) if (matches + mismatches) else None
        ),
    }
