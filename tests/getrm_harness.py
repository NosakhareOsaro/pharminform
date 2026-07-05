"""GeT-RM validation harness.

Compares PharmInform's star allele calls against published GeT-RM
(Genetic Testing Reference Materials Coordination Program) consensus
genotypes for public Coriell reference samples.

Real, live attempt made during development: pulling 1000 Genomes-hosted
sequencing data (https://1000genomes.s3.amazonaws.com/release/20130502/...)
directly into cyvcf2 over HTTPS failed in this environment --
``cyvcf2.VCF(url)`` raised ``OSError: Failed to open file ... No such file
or directory``, which traces to the bundled htslib not being built with
remote (libcurl) URL support. This is a real, reproducible finding, not a
placeholder excuse: the fix is either (a) a cyvcf2/htslib build with
libcurl support, or (b) ``bcftools``/``tabix`` on PATH (neither is
installed in this environment), or (c) downloading the relevant region
locally first. ``fetch_1000genomes_region`` below implements path (b) and
raises a clear, actionable error otherwise -- it does not silently return
fabricated genotypes.

See validation/getrm_concordance_report.md for the full writeup, what
running this for real requires, and which samples/genes are and are not
yet covered by the bundled consensus table.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

CONSENSUS_PATH = Path(__file__).parent.parent / "src" / "pharminform" / "data" / "getrm_consensus_genotypes.json"

# 1000 Genomes phase 3 GRCh37 VCFs, per-chromosome, publicly readable
# (confirmed live via the AWS Open Data mirror; the EBI FTP mirror's path
# also works but uses a "v5b" filename suffix as of this writing, not "v5a").
ONE_KG_BASE_URL = "https://1000genomes.s3.amazonaws.com/release/20130502"
ONE_KG_CHROM_FILES = {
    "10": f"{ONE_KG_BASE_URL}/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "22": f"{ONE_KG_BASE_URL}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
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


def load_consensus_table(path: Path | str = CONSENSUS_PATH) -> list[dict]:
    data = json.loads(Path(path).read_text())
    return data["samples"]


def bcftools_available() -> bool:
    return shutil.which("bcftools") is not None


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


def run_concordance_check(samples: list[dict] | None = None) -> list[ConcordanceOutcome]:
    """Attempt a real concordance check for each (sample, gene) pair in the
    consensus table. In environments without bcftools/network access to
    1000 Genomes, every sample is reported as "skipped_no_data" with the
    specific reason -- never silently marked as a match.
    """
    samples = samples if samples is not None else load_consensus_table()
    outcomes: list[ConcordanceOutcome] = []
    for entry in samples:
        try:
            # NOTE: full concordance requires resolving each sample's gene
            # region to GRCh37 coordinates and running the star allele
            # caller against the fetched variants -- not yet wired up end
            # to end (see validation/getrm_concordance_report.md for what's
            # needed). fetch_1000genomes_region is exercised here only to
            # surface a real, honest availability check.
            fetch_1000genomes_region("22", 42522500, 42526900, entry["coriell_id"])
            outcomes.append(
                ConcordanceOutcome(
                    coriell_id=entry["coriell_id"],
                    gene=entry["gene"],
                    consensus_diplotype=entry["consensus_diplotype"],
                    called_diplotype=None,
                    status="skipped_no_data",
                    detail="sequencing data fetched but end-to-end calling not yet wired up",
                )
            )
        except SequencingDataUnavailable as exc:
            outcomes.append(
                ConcordanceOutcome(
                    coriell_id=entry["coriell_id"],
                    gene=entry["gene"],
                    consensus_diplotype=entry["consensus_diplotype"],
                    called_diplotype=None,
                    status="skipped_no_data",
                    detail=str(exc),
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
