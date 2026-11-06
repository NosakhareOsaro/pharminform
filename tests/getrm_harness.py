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
import tempfile
from dataclasses import dataclass
from pathlib import Path

import cyvcf2

from pharminform import caller
from pharminform.pharmgkb_client import diplotype_to_phenotype

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

# Approximate GRCh37 regions, wide enough to catch the defining variants in
# caller.BOOTSTRAP_ALLELE_DEFINITIONS (checked against Ensembl GRCh37
# positions, including CYP2C19 *17 upstream of the gene and UGT1A1 *80 in
# the promoter). GRCh37, NOT GRCh38 -- 1000 Genomes phase 3 is GRCh37,
# unlike caller.GENE_REGIONS_GRCH38. This only matters for the initial
# region slice; the actual diplotype match is by rsID, which is
# build-independent.
GRCH37_GENE_REGIONS = {
    "CYP2D6": ("22", 42_500_000, 42_550_000),
    "CYP2C19": ("10", 96_515_000, 96_615_000),
    "CYP2C9": ("10", 96_695_000, 96_750_000),
    "DPYD": ("1", 97_540_000, 98_390_000),
    "UGT1A1": ("2", 234_665_000, 234_682_000),
    "TPMT": ("6", 18_125_000, 18_158_000),
}

# GRCh37 positions of the defining variants in DPYD, TPMT, UGT1A1, CYP2C19
# and CYP2C9
# (Ensembl GRCh37). DPYD spans ~850 kb, and each remote fetch pulls every
# sample's data for the requested range, so for these genes only small
# windows around the defining sites are fetched. The rule-based caller
# only uses those sites, so the call is the same.
GRCH37_DEFINING_SITES = {
    "DPYD": [97_915_614, 97_981_343, 97_547_947, 98_039_419, 98_045_449],
    "TPMT": [18_143_955, 18_139_228, 18_130_918],
    "UGT1A1": [234_669_144, 234_668_879, 234_668_880, 234_668_570],
    "CYP2C19": [96_541_616, 96_540_410, 96_521_657],
    "CYP2C9": [96_702_047, 96_741_053],
}
SITE_WINDOW = 25

# 1000 Genomes 30x high-coverage phased panel (NYGC; Byrska-Bishop 2022,
# PMID 36055201), GRCh38, 3202 samples: the 2504 phase 3 samples plus 698
# relatives, so more GeT-RM samples are covered, and each sample is
# sequenced directly rather than at low coverage with imputation. Hosted on
# the EBI FTP mirror; the AWS bucket has no copy of this directory.
ONE_KG_30X_BASE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
    "working/20220422_3202_phased_SNV_INDEL_SV"
)
ONE_KG_30X_CHROM_FILES = {
    str(n): f"{ONE_KG_30X_BASE_URL}/1kGP_high_coverage_Illumina.chr{n}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    for n in range(1, 23)
}

# The 30x VCF has no rsIDs (IDs read "10:94781859:G:A") and the caller
# matches by rsID, so the defining sites are labelled from this table:
# (GRCh38 position, REF, ALT, rsID), checked against Ensembl GRCh38 and
# against the records in the 30x VCF itself. REF/ALT of None matches any
# indel at that position: the UGT1A1 TA repeat appears there as C>CAT
# (TA7, *28), C>CATAT (TA8, *37) and CAT>C (TA5, *36), which the caller
# tells apart by length.
GRCH38_DEFINING_SITES = {
    "DPYD": [
        (97_450_058, "C", "T", "rs3918290"),
        (97_515_787, "A", "C", "rs55886062"),
        (97_082_391, "T", "A", "rs67376798"),
        (97_573_863, "C", "T", "rs56038477"),
        (97_579_893, "G", "C", "rs75017182"),
    ],
    "TPMT": [
        (18_143_724, "C", "G", "rs1800462"),
        (18_138_997, "C", "T", "rs1800460"),
        (18_130_687, "T", "C", "rs1142345"),
    ],
    "UGT1A1": [
        (233_760_498, "G", "A", "rs4148323"),
        (233_760_233, None, None, "rs3064744"),
        (233_759_924, "C", "T", "rs887829"),
    ],
    "CYP2C19": [
        (94_781_859, "G", "A", "rs4244285"),
        (94_780_653, "G", "A", "rs4986893"),
        (94_761_900, "C", "T", "rs12248560"),
    ],
    "CYP2C9": [
        (94_942_290, "C", "T", "rs1799853"),
        (94_981_296, "A", "C", "rs1057910"),
    ],
}
GRCH38_CHROMS = {"DPYD": "1", "TPMT": "6", "UGT1A1": "2", "CYP2C19": "10", "CYP2C9": "10"}


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
    consensus_phenotype: str | None = None
    called_phenotype: str | None = None
    # "match", "mismatch", or "not_assessable" when an allele in either
    # diplotype has no function in the bootstrap table (e.g. CYP2D6 *5).
    phenotype_status: str | None = None


def load_consensus_table(path: Path | str = CONSENSUS_PATH) -> list[dict]:
    data = json.loads(Path(path).read_text())
    return data["samples"]


def bcftools_available() -> bool:
    return shutil.which("bcftools") is not None


def fetch_1000genomes_panel_samples() -> set[str]:
    """Fetch the real 1000 Genomes phase 3 sample panel and return the set
    of sample IDs. Raises SequencingDataUnavailable on any network failure
    -- never returns a guessed/cached list presented as live data."""
    proc = subprocess.run(
        ["curl", "-s", "-f", ONE_KG_PANEL_URL], capture_output=True, text=True, timeout=60, check=False
    )
    if proc.returncode != 0 or not proc.stdout:
        raise SequencingDataUnavailable(
            f"Could not fetch the 1000 Genomes sample panel from {ONE_KG_PANEL_URL}"
        )
    lines = proc.stdout.strip().splitlines()
    samples = {line.split("\t")[0] for line in lines[1:] if line.strip()}
    return samples


def fetch_1000genomes_region(chrom: str, start: int, end: int, sample: str) -> str:
    return fetch_1000genomes_regions(chrom, [(start, end)], sample)


def fetch_1000genomes_regions(chrom: str, spans: list[tuple[int, int]], sample: str) -> str:
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
    merged: list[list[int]] = []
    for start, end in sorted(spans):
        # Overlapping regions would make bcftools emit a record twice.
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    regions = ",".join(f"{chrom}:{start}-{end}" for start, end in merged)
    cmd = ["bcftools", "view", "-r", regions, "-s", sample, url]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=False)
    if proc.returncode != 0:
        raise SequencingDataUnavailable(
            f"bcftools failed fetching {regions} for {sample}: {proc.stderr.strip()}"
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
    sites = GRCH37_DEFINING_SITES.get(gene)
    if sites:
        spans = [(pos - SITE_WINDOW, pos + SITE_WINDOW) for pos in sites]
        vcf_text = fetch_1000genomes_regions(chrom, spans, sample)
    else:
        vcf_text = fetch_1000genomes_region(chrom, start, end, sample)

    tmp_path = Path(f"/tmp/getrm_{sample}_{gene}.vcf")
    tmp_path.write_text(vcf_text)
    try:
        variants = _variants_from_vcf(gene, tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)

    call = caller.call_diplotype_from_vcf(gene, variants)
    return call, len(variants)


def _site_rsid(gene: str, pos: int, ref: str, alt: str) -> str | None:
    for site_pos, site_ref, site_alt, rsid in GRCH38_DEFINING_SITES.get(gene, []):
        if pos != site_pos:
            continue
        if site_ref is None and len(ref) != len(alt):
            return rsid
        if (ref, alt) == (site_ref, site_alt):
            return rsid
    return None


def fetch_30x_sites(gene: str, samples: list[str]) -> str:
    """Fetch a gene's defining sites for many samples from the 30x panel in
    one remote call; the per-sample calls are then made locally."""
    if not bcftools_available():
        raise SequencingDataUnavailable("bcftools is not installed/on PATH.")
    chrom = GRCH38_CHROMS.get(gene)
    if chrom is None:
        raise SequencingDataUnavailable(f"No GRCh38 defining sites configured for {gene}.")
    regions = ",".join(f"chr{chrom}:{pos}" for pos, *_ in sorted(GRCH38_DEFINING_SITES[gene]))
    cmd = ["bcftools", "view", "-r", regions, "-s", ",".join(samples), ONE_KG_30X_CHROM_FILES[chrom]]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600, check=False)
    if proc.returncode != 0:
        raise SequencingDataUnavailable(f"bcftools failed fetching {regions}: {proc.stderr.strip()}")
    return proc.stdout


def call_diplotypes_from_30x(gene: str, samples: list[str]) -> dict[str, tuple[caller.DiplotypeCall, int]]:
    """Call each sample from the 30x panel, labelling the defining sites with
    their rsIDs (GRCH38_DEFINING_SITES) before the caller sees them."""
    tmp_dir = Path(tempfile.mkdtemp(prefix=f"getrm30x_{gene}_"))
    try:
        multi = tmp_dir / "sites.vcf"
        multi.write_text(fetch_30x_sites(gene, samples))
        calls = {}
        for sample in samples:
            single = tmp_dir / f"{sample}.vcf"
            subprocess.run(
                ["bcftools", "view", "-s", sample, "-o", str(single), str(multi)],
                check=True, capture_output=True, text=True,
            )
            variants = _variants_from_vcf(gene, single, label_sites=True)
            calls[sample] = (caller.call_diplotype_from_vcf(gene, variants), len(variants))
        return calls
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _variants_from_vcf(gene: str, path: Path, label_sites: bool = False) -> list[caller.VariantCall]:
    variants = []
    for rec in cyvcf2.VCF(str(path)):
        if label_sites:
            rsid = _site_rsid(gene, rec.POS, rec.REF, rec.ALT[0]) if len(rec.ALT) == 1 else None
        else:
            rsid = rec.ID
        if rsid is None:
            continue
        gt_str = caller._genotype_string(rec)
        if not caller._has_alt_allele(gt_str):
            continue
        for alt_index, alt in enumerate(rec.ALT, start=1):
            if not caller._alt_copies(gt_str, alt_index):
                continue
            variants.append(
                caller.VariantCall(
                    gene=gene, chrom=rec.CHROM, pos=rec.POS, ref=rec.REF, alt=alt,
                    # Neither 1000 Genomes call set has per-sample DP in
                    # FORMAT (phase 3 is an integrated/imputed callset; the
                    # 30x phased panel carries GT only) -- depth genuinely
                    # unavailable here, not a bug.
                    genotype=gt_str, depth=None, qual=rec.QUAL, filter=rec.FILTER, rsid=rsid,
                    alt_index=alt_index,
                )
            )
    return variants


def _phenotype_if_assessable(gene: str, diplotype: str) -> str | None:
    """Phenotype for a diplotype, or None if any allele's function is not
    in the bootstrap table -- an unknown allele would otherwise silently
    default to normal function and make the comparison meaningless."""
    defs = caller.BOOTSTRAP_ALLELE_DEFINITIONS.get(gene, {})
    alleles = diplotype.split("/")
    if len(alleles) != 2 or any(a != "*1" and a not in defs for a in alleles):
        return None
    functions = {name: info["function"] for name, info in defs.items()}
    phenotype, _ = diplotype_to_phenotype(gene, diplotype, functions)
    return phenotype


def compare_phenotypes(gene: str, consensus: str, called: str) -> tuple[str | None, str | None, str]:
    consensus_pheno = _phenotype_if_assessable(gene, consensus)
    called_pheno = _phenotype_if_assessable(gene, called)
    if consensus_pheno is None or called_pheno is None:
        return consensus_pheno, called_pheno, "not_assessable"
    return consensus_pheno, called_pheno, "match" if consensus_pheno == called_pheno else "mismatch"


def fetch_30x_panel_samples() -> set[str]:
    proc = subprocess.run(
        ["bcftools", "query", "-l", ONE_KG_30X_CHROM_FILES["10"]],
        capture_output=True, text=True, timeout=120, check=False,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        raise SequencingDataUnavailable(f"Could not list 30x panel samples: {proc.stderr.strip()}")
    return set(proc.stdout.split())


def mask_untested(diplotype: str, untested: list[str] | tuple[str, ...]) -> str:
    """Read alleles the consensus source never tested as *1, so the
    comparison covers only what both sides looked at."""
    return "/".join("*1" if a in untested else a for a in diplotype.split("/"))


def run_concordance_check(samples: list[dict] | None = None, source: str = "phase3") -> list[ConcordanceOutcome]:
    """Run a real concordance check for each (sample, gene) pair in the
    consensus table against actual 1000 Genomes sequencing data.

    Samples not present in the 1000 Genomes panel, or blocked by a real
    tooling limitation (no bcftools, no region configured), are reported
    as "skipped_no_data" with the specific reason -- never silently
    marked as a match, and never given a fabricated diplotype.

    ``source`` is "phase3" (GRCh37, per-sample fetches) or "30x" (the
    GRCh38 high-coverage panel, fetched once per gene).
    """
    samples = samples if samples is not None else load_consensus_table()
    outcomes: list[ConcordanceOutcome] = []

    try:
        panel_samples = fetch_30x_panel_samples() if source == "30x" else fetch_1000genomes_panel_samples()
    except SequencingDataUnavailable as exc:
        return [
            ConcordanceOutcome(
                coriell_id=entry["coriell_id"], gene=entry["gene"],
                consensus_diplotype=entry["consensus_diplotype"], called_diplotype=None,
                status="skipped_no_data", detail=str(exc),
            )
            for entry in samples
        ]

    batch_calls: dict[tuple[str, str], tuple[caller.DiplotypeCall, int]] = {}
    batch_errors: dict[str, str] = {}
    if source == "30x":
        for gene in sorted({e["gene"] for e in samples}):
            ids = sorted({e["coriell_id"] for e in samples if e["gene"] == gene} & panel_samples)
            if not ids:
                continue
            try:
                for sample, result in call_diplotypes_from_30x(gene, ids).items():
                    batch_calls[(sample, gene)] = result
            except SequencingDataUnavailable as exc:
                batch_errors[gene] = str(exc)
    panel_name = "1000 Genomes 30x" if source == "30x" else "1000 Genomes phase 3"

    for entry in samples:
        coriell_id, gene, consensus = entry["coriell_id"], entry["gene"], entry["consensus_diplotype"]

        if coriell_id not in panel_samples:
            outcomes.append(
                ConcordanceOutcome(
                    coriell_id=coriell_id, gene=gene, consensus_diplotype=consensus, called_diplotype=None,
                    status="skipped_no_data",
                    detail=f"not present in the {panel_name} sample panel",
                )
            )
            continue

        try:
            if source == "30x":
                if gene in batch_errors:
                    raise SequencingDataUnavailable(batch_errors[gene])
                call, n_variants = batch_calls[(coriell_id, gene)]
            else:
                call, n_variants = call_diplotype_from_1000genomes(gene, coriell_id)
        except SequencingDataUnavailable as exc:
            outcomes.append(
                ConcordanceOutcome(
                    coriell_id=coriell_id, gene=gene, consensus_diplotype=consensus, called_diplotype=None,
                    status="skipped_no_data", detail=str(exc),
                )
            )
            continue

        untested = entry.get("untested_alleles") or ()
        compared = mask_untested(call.diplotype, untested)
        # Compare in canonical allele order: "*4/*1" and "*1/*4" are the
        # same diplotype.
        same = caller.normalize_diplotype(compared) == caller.normalize_diplotype(consensus)
        status = "match" if same else "mismatch"
        consensus_pheno, called_pheno, pheno_status = compare_phenotypes(gene, consensus, compared)
        notes = [] if status == "match" else list(call.ambiguous_candidates)
        if compared != call.diplotype:
            notes.append(
                f"called {call.diplotype}; {', '.join(a for a in untested if a in call.diplotype)} "
                "not tested by the consensus source, so compared as *1"
            )
        outcomes.append(
            ConcordanceOutcome(
                coriell_id=coriell_id, gene=gene, consensus_diplotype=consensus,
                called_diplotype=call.diplotype, status=status,
                detail="; ".join(notes),
                n_variants_considered=n_variants,
                consensus_phenotype=consensus_pheno, called_phenotype=called_pheno,
                phenotype_status=pheno_status,
            )
        )

    return outcomes


def summarize(outcomes: list[ConcordanceOutcome]) -> dict:
    total = len(outcomes)
    matches = sum(1 for o in outcomes if o.status == "match")
    mismatches = sum(1 for o in outcomes if o.status == "mismatch")
    skipped = sum(1 for o in outcomes if o.status == "skipped_no_data")
    pheno_matches = sum(1 for o in outcomes if o.phenotype_status == "match")
    pheno_mismatches = sum(1 for o in outcomes if o.phenotype_status == "mismatch")
    return {
        "total": total,
        "matches": matches,
        "mismatches": mismatches,
        "skipped_no_data": skipped,
        "concordance_pct_of_attempted": (
            round(100 * matches / (matches + mismatches), 1) if (matches + mismatches) else None
        ),
        "phenotype_matches": pheno_matches,
        "phenotype_mismatches": pheno_mismatches,
        "phenotype_not_assessable": sum(1 for o in outcomes if o.phenotype_status == "not_assessable"),
    }
