"""Star allele calling.

Two calling paths are implemented, because Aldy and a plain VCF solve
different problems:

* ``run_aldy`` shells out to the real Aldy CLI (https://github.com/0xTCG/aldy)
  against aligned reads (BAM/CRAM). Aldy performs copy-number-aware calling,
  which is required to correctly resolve structural variation in genes like
  CYP2D6 (whole-gene deletions/duplications, CYP2D6/CYP2D7 hybrids). This is
  the clinically-preferred path whenever aligned reads are available.

* ``call_diplotype_from_vcf`` is a simplified, rule-based matcher that looks
  up called variants against a small embedded table of core star-allele
  defining variants. A VCF alone cannot represent copy number or phasing, so
  this path is a best-effort fallback for VCF-only input and is always
  reported at reduced confidence (see ``confidence.py``). The embedded table
  in this module is a bootstrap/test fixture, NOT an authoritative allele
  definition source -- ``pharmvar_client.py`` is the intended source of truth
  and should be preferred once network access to PharmVar is available.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import cyvcf2

TARGET_GENES = ["CYP2D6", "CYP2C19", "DPYD", "UGT1A1", "TPMT", "CYP2C9"]

# Approximate gene regions (GRCh38), used only to subset the VCF to the
# relevant loci before matching. Widened slightly beyond the coding sequence
# to catch upstream/intronic defining variants (e.g. splice-site alleles).
GENE_REGIONS_GRCH38 = {
    "CYP2D6": "chr22:42125000-42132000",
    "CYP2C19": "chr10:94762000-94856000",
    "DPYD": "chr1:97543000-98387000",
    "UGT1A1": "chr2:233757000-233774000",
    "TPMT": "chr6:18128000-18156000",
    "CYP2C9": "chr10:94938000-94991000",
}

# Minimal, illustrative set of core star-allele defining variants per gene,
# keyed by rsID. This is NOT an exhaustive or guaranteed-current allele
# definition table -- it exists so the rule-based VCF path and the test
# suite have something concrete to run against offline. Production use
# should source definitions from PharmVar (see pharmvar_client.py).
#
# rsIDs below are widely cited in the PGx literature for these alleles;
# verify against PharmVar/CPIC before relying on them clinically.
BOOTSTRAP_ALLELE_DEFINITIONS: dict[str, dict[str, dict]] = {
    "CYP2D6": {
        "*3": {"rsid": "rs35742686", "ref": "A", "alt": "AGT" , "function": "no function"},
        "*4": {"rsid": "rs3892097", "ref": "G", "alt": "A", "function": "no function"},
        "*10": {"rsid": "rs1065852", "ref": "G", "alt": "A", "function": "decreased function"},
        "*17": {"rsid": "rs28371706", "ref": "C", "alt": "T", "function": "decreased function"},
        "*41": {"rsid": "rs28371725", "ref": "C", "alt": "T", "function": "decreased function"},
    },
    "CYP2C19": {
        "*2": {"rsid": "rs4244285", "ref": "G", "alt": "A", "function": "no function"},
        "*3": {"rsid": "rs4986893", "ref": "G", "alt": "A", "function": "no function"},
        "*17": {"rsid": "rs12248560", "ref": "C", "alt": "T", "function": "increased function"},
    },
    "CYP2C9": {
        "*2": {"rsid": "rs1799853", "ref": "C", "alt": "T", "function": "decreased function"},
        "*3": {"rsid": "rs1057910", "ref": "A", "alt": "C", "function": "no function"},
    },
    "DPYD": {
        "*2A": {"rsid": "rs3918290", "ref": "C", "alt": "T", "function": "no function"},
        "c.2846A>T": {"rsid": "rs67376798", "ref": "A", "alt": "T", "function": "decreased function"},
        "c.1236G>A(HapB3)": {"rsid": "rs56038477", "ref": "G", "alt": "A", "function": "decreased function"},
    },
    "UGT1A1": {
        "*6": {"rsid": "rs4148323", "ref": "G", "alt": "A", "function": "decreased function"},
        # *28 is a TA-repeat expansion in the promoter (TA6>TA7), a repeat
        # variant rather than a clean SNV -- flagged here rather than matched
        # by position, since naive REF/ALT matching does not represent it well.
    },
    "TPMT": {
        "*2": {"rsid": "rs1800462", "ref": "C", "alt": "G", "function": "no function"},
        "*3B": {"rsid": "rs1800460", "ref": "G", "alt": "A", "function": "no function"},
        "*3C": {"rsid": "rs1142345", "ref": "A", "alt": "G", "function": "no function"},
    },
}


class VCFIngestError(RuntimeError):
    """Raised when the input VCF cannot be parsed or has no usable records."""


class AldyNotAvailable(RuntimeError):
    """Raised when Aldy is requested but not installed/importable."""


@dataclass
class VariantCall:
    gene: str
    chrom: str
    pos: int
    ref: str
    alt: str
    genotype: str  # e.g. "0/1", "1/1"
    depth: int | None
    qual: float | None
    filter: str | None
    rsid: str | None = None


@dataclass
class DiplotypeCall:
    gene: str
    diplotype: str
    method: str  # "aldy" or "rule_based_vcf"
    matched_variants: list[VariantCall] = field(default_factory=list)
    ambiguous_candidates: list[str] = field(default_factory=list)
    raw_aldy_output: str | None = None


def load_variants(
    vcf_path: str | Path,
    genes: list[str] | None = None,
    min_depth: int = 10,
    require_pass: bool = True,
) -> dict[str, list[VariantCall]]:
    """Parse a VCF with cyvcf2, apply quality filters, and return variant
    calls grouped by target gene.

    Records are dropped (not silently kept) if they fail the PASS filter
    (when ``require_pass``) or fall below ``min_depth``; both are common
    causes of an unreliable star-allele call and are surfaced later via
    confidence.py rather than hidden here.
    """
    genes = genes or TARGET_GENES
    path = Path(vcf_path)
    if not path.exists():
        raise VCFIngestError(f"VCF not found: {path}")

    try:
        vcf = cyvcf2.VCF(str(path))
    except Exception as exc:  # pragma: no cover - cyvcf2 raises various C-level errors
        raise VCFIngestError(f"Could not parse VCF {path}: {exc}") from exc

    results: dict[str, list[VariantCall]] = {gene: [] for gene in genes}

    for gene in genes:
        region = GENE_REGIONS_GRCH38.get(gene)
        if region is None:
            continue
        try:
            records = list(vcf(region))
        except Exception:
            # No tabix/CSI index, or region string didn't match this VCF's
            # contig naming (e.g. "22" vs "chr22") -- fall back to a full scan.
            records = _fallback_scan(path, region)

        for rec in records:
            if require_pass and rec.FILTER not in (None, "PASS", "."):
                continue
            depth = _extract_depth(rec)
            if depth is not None and depth < min_depth:
                continue
            gt_str = _genotype_string(rec)  # per-sample; assume single-sample VCF
            for alt in rec.ALT:
                results[gene].append(
                    VariantCall(
                        gene=gene,
                        chrom=rec.CHROM,
                        pos=rec.POS,
                        ref=rec.REF,
                        alt=alt,
                        genotype=gt_str,
                        depth=depth,
                        qual=rec.QUAL,
                        filter=rec.FILTER,
                        rsid=rec.ID,
                    )
                )
    return results


def _fallback_scan(path: Path, region: str) -> list:
    chrom, rest = region.split(":")
    start, end = (int(x) for x in rest.split("-"))
    alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    out = []
    for candidate in (chrom, alt_chrom):
        try:
            vcf = cyvcf2.VCF(str(path))
            for rec in vcf:
                if rec.CHROM == candidate and start <= rec.POS <= end:
                    out.append(rec)
        except Exception:
            continue
        if out:
            break
    return out


def _extract_depth(rec) -> int | None:
    try:
        dp = rec.INFO.get("DP")
        if dp is not None:
            return int(dp)
    except Exception:
        pass
    try:
        fmt_dp = rec.format("DP")
        if fmt_dp is not None and len(fmt_dp):
            return int(fmt_dp[0][0])
    except Exception:
        pass
    return None


def _genotype_string(rec) -> str:
    try:
        gt = rec.genotypes[0]
        alleles = "/".join(str(a) for a in gt[:-1])
        phased = "|" if len(gt) > 2 and gt[-1] else "/"
        return alleles.replace("/", phased) if phased == "|" else alleles
    except Exception:
        return "./."


def aldy_available() -> bool:
    return shutil.which("aldy") is not None


def run_aldy(
    bam_path: str | Path,
    gene: str,
    genome: str = "hg38",
    profile: str = "illumina",
    solver: str = "any",
    output_dir: str | Path | None = None,
) -> DiplotypeCall:
    """Run the real Aldy CLI against aligned reads for one gene.

    Requires ``aldy`` to be installed (``pip install pharminform[aldy]``) and
    on PATH; raises ``AldyNotAvailable`` otherwise rather than silently
    returning a fake result.
    """
    if not aldy_available():
        raise AldyNotAvailable(
            "aldy is not installed or not on PATH. Install with "
            "`pip install pharminform[aldy]` (requires an ILP solver -- "
            "Gurobi or Google OR-Tools/CBC -- see Aldy's documentation)."
        )
    bam_path = Path(bam_path)
    if not bam_path.exists():
        raise VCFIngestError(f"BAM/CRAM not found: {bam_path}")

    out_dir = Path(output_dir) if output_dir else bam_path.parent
    out_file = out_dir / f"{bam_path.stem}.{gene}.aldy"

    cmd = [
        "aldy", "genotype",
        "-g", gene,
        "--genome", genome,
        "--profile", profile,
        "--solver", solver,
        "-o", str(out_file),
        str(bam_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"aldy genotype failed for {gene}: {proc.stderr.strip()}")

    diplotype, ambiguous = _parse_aldy_output(out_file)
    return DiplotypeCall(
        gene=gene,
        diplotype=diplotype,
        method="aldy",
        ambiguous_candidates=ambiguous,
        raw_aldy_output=proc.stdout,
    )


def _parse_aldy_output(out_file: Path) -> tuple[str, list[str]]:
    if not out_file.exists():
        return "unknown", []
    lines = [line.strip() for line in out_file.read_text().splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return "unknown", []
    # Aldy's .aldy output lists one solution per line; first column is the
    # diplotype string. Multiple lines => multiple equally-likely solutions.
    candidates = [line.split("\t")[0] for line in lines]
    return candidates[0], candidates[1:] if len(candidates) > 1 else []


def call_diplotype_from_vcf(
    gene: str,
    variants: list[VariantCall],
    allele_definitions: dict[str, dict] | None = None,
) -> DiplotypeCall:
    """Best-effort rule-based diplotype call from VCF variants alone.

    Matches called variants against defining SNVs/indels for the gene and
    reports a diplotype like "*4/*10" when exactly two defining alleles are
    observed (one homozygous defining variant is reported as e.g. "*4/*4").
    Falls back to "*1/*1" (wild type / normal function) when no defining
    variants are observed -- this reflects "no evidence of a variant allele
    in this table", not a guarantee of true wild type, since structural
    variants (deletions, duplications, repeat expansions) are invisible to
    this path. See confidence.py, which flags this explicitly.
    """
    defs = allele_definitions or BOOTSTRAP_ALLELE_DEFINITIONS.get(gene, {})
    observed: list[tuple[str, str]] = []  # (allele_name, zygosity)

    by_rsid = {info["rsid"]: name for name, info in defs.items()}

    for v in variants:
        if v.rsid in by_rsid:
            allele_name = by_rsid[v.rsid]
            zygosity = "hom" if v.genotype in ("1/1", "1|1") else "het"
            observed.append((allele_name, zygosity))

    if not observed:
        return DiplotypeCall(gene=gene, diplotype="*1/*1", method="rule_based_vcf")

    if len(observed) == 1:
        name, zygosity = observed[0]
        diplotype = f"{name}/{name}" if zygosity == "hom" else f"{name}/*1"
        return DiplotypeCall(gene=gene, diplotype=diplotype, method="rule_based_vcf", matched_variants=variants)

    # More than one defining variant observed -- without phasing information
    # a VCF cannot tell us whether they're in cis (same allele, e.g. a
    # multi-SNV star allele) or trans (two different variant alleles). Report
    # the naive combination as the top candidate and list alternatives.
    names = [n for n, _ in observed]
    diplotype = "/".join(sorted(set(names))[:2]) if len(set(names)) > 1 else f"{names[0]}/{names[0]}"
    ambiguous = [f"phase unresolved among: {', '.join(sorted(set(names)))}"]
    return DiplotypeCall(
        gene=gene,
        diplotype=diplotype,
        method="rule_based_vcf",
        matched_variants=variants,
        ambiguous_candidates=ambiguous,
    )
