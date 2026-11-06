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

import contextlib
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import cyvcf2

TARGET_GENES = ["CYP2D6", "CYP2C19", "DPYD", "UGT1A1", "TPMT", "CYP2C9"]

# Genes the VCF-only path does not call. CYP2D6's common alleles include
# whole-gene deletion (*5), duplications and CYP2D6/CYP2D7 hybrids (*36,
# *13), none of which a VCF's SNV/indel records show, and *10's defining
# variant (rs1065852) also lies on *4 -- the rule-based path matched 0 of 6
# GeT-RM samples. A phenotype from it would look precise and be wrong, so
# these genes are reported as not called; use run_aldy on aligned reads.
VCF_UNSUPPORTED_GENES = {
    "CYP2D6": (
        "CYP2D6 is not called from VCF input: its common structural alleles "
        "(*5 deletion, duplications, CYP2D6/2D7 hybrids) are not visible in a VCF. "
        "Call it with Aldy from aligned reads (BAM/CRAM)."
    ),
}

# Approximate gene regions (GRCh38), used only to subset the VCF to the
# relevant loci before matching. Widened slightly beyond the coding sequence
# to catch upstream/intronic defining variants (e.g. splice-site alleles).
# Checked against Ensembl GRCh38 positions of every defining variant below:
# an earlier DPYD entry used GRCh37 coordinates (missing *2A at
# chr1:97450058 and c.2846A>T at chr1:97082391), and the CYP2C19 window
# started just after the *17 promoter variant (chr10:94761900).
GENE_REGIONS_GRCH38 = {
    "CYP2D6": "chr22:42125000-42132000",
    "CYP2C19": "chr10:94760000-94856000",
    "DPYD": "chr1:97075000-97998000",
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
#
# Optional keys beyond rsid/ref/alt/function:
#   alias_rsids          other IDs for the same variant (merged rsIDs, or the
#                        ID a specific callset uses, e.g. 1000 Genomes phase 3
#                        records the UGT1A1 TA repeat as rs34983651).
#   proxy_rsids          markers in strong linkage disequilibrium with the
#                        allele; a match on one of these alone reports the
#                        allele, with proxy_note added as an ambiguity.
#   indel_length_change  required len(ALT) - len(REF), for repeat alleles that
#                        share one rsID (UGT1A1 *36/*28/*37 are TA5/TA7/TA8;
#                        the reference genome carries TA6, i.e. *1).
#   cis                  a haplotype defined by other alleles in this table
#                        lying on the same chromosome (TPMT *3A = *3B + *3C).
#
# ref/alt are GRCh38 forward-strand alleles, for reference only; matching
# is by rsID (plus indel_length_change where given).
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
        "*13": {"rsid": "rs55886062", "ref": "A", "alt": "C", "function": "no function"},
        "c.2846A>T": {"rsid": "rs67376798", "ref": "T", "alt": "A", "function": "decreased function"},
        # HapB3 is tagged by c.1236G>A and its causal intronic variant
        # c.1129-5923C>G (rs75017182); either one identifies the haplotype.
        "c.1236G>A(HapB3)": {
            "rsid": "rs56038477", "ref": "C", "alt": "T", "function": "decreased function",
            "proxy_rsids": ["rs75017182"],
        },
    },
    "UGT1A1": {
        "*6": {"rsid": "rs4148323", "ref": "G", "alt": "A", "function": "decreased function"},
        # Promoter TA repeat. Matched on the insertion/deletion length, since
        # one rsID covers every repeat length.
        "*28": {
            "rsid": "rs3064744", "ref": "A(TA)6", "alt": "A(TA)7", "function": "decreased function",
            "alias_rsids": ["rs8175347", "rs34983651"],
            "indel_length_change": 2,
            "proxy_rsids": ["rs887829"],
            "proxy_note": (
                "*28 inferred from the *80 proxy (rs887829); the TA repeat was "
                "not observed directly -- confirm *28 before dosing irinotecan"
            ),
        },
        "*36": {
            "rsid": "rs3064744", "ref": "A(TA)6", "alt": "A(TA)5", "function": "increased function",
            "alias_rsids": ["rs8175347"],
            "indel_length_change": -2,
        },
        "*37": {
            "rsid": "rs3064744", "ref": "A(TA)6", "alt": "A(TA)8", "function": "decreased function",
            "alias_rsids": ["rs8175347"],
            "indel_length_change": 4,
        },
    },
    "TPMT": {
        "*2": {"rsid": "rs1800462", "ref": "C", "alt": "G", "function": "no function"},
        "*3B": {"rsid": "rs1800460", "ref": "C", "alt": "T", "function": "no function"},
        "*3C": {"rsid": "rs1142345", "ref": "T", "alt": "C", "function": "no function"},
        "*3A": {"rsid": None, "cis": ["*3B", "*3C"], "function": "no function"},
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
    alt_index: int = 1  # which ALT allele (1-based) this call represents


@dataclass
class DiplotypeCall:
    gene: str
    diplotype: str
    method: str  # "aldy" or "rule_based_vcf"
    matched_variants: list[VariantCall] = field(default_factory=list)
    ambiguous_candidates: list[str] = field(default_factory=list)
    raw_aldy_output: str | None = None
    # Alleles inferred only from an LD proxy (UGT1A1 *28 from *80). CPIC
    # reports such a result as indeterminate (see pipeline.py).
    proxy_only_alleles: list[str] = field(default_factory=list)


def load_variants(
    vcf_path: str | Path,
    genes: list[str] | None = None,
    min_depth: int = 10,
    require_pass: bool = True,
    dropped: dict[str, list[VariantCall]] | None = None,
) -> dict[str, list[VariantCall]]:
    """Parse a VCF with cyvcf2, apply quality filters, and return variant
    calls grouped by target gene.

    Records are dropped (not silently kept) if they fail the PASS filter
    (when ``require_pass``) or fall below ``min_depth``; both are common
    causes of an unreliable star-allele call and are surfaced later via
    confidence.py rather than hidden here.

    If ``dropped`` is given, alt-carrying records excluded by those filters
    are collected there per gene, so a defining variant lost to low depth
    can be reported instead of silently becoming a *1/*1 call.
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
        except (AssertionError, OSError, ValueError):
            # No tabix/CSI index (cyvcf2 raises AssertionError), or region
            # string didn't match this VCF's contig naming (e.g. "22" vs
            # "chr22") -- fall back to a full scan.
            records = _fallback_scan(path, region)

        for rec in records:
            depth = _extract_depth(rec)
            gt_str = _genotype_string(rec)  # per-sample; assume single-sample VCF
            if not _has_alt_allele(gt_str):
                continue
            excluded = (require_pass and rec.FILTER not in (None, "PASS", ".")) or (
                depth is not None and depth < min_depth
            )
            if excluded and dropped is None:
                continue
            target = dropped.setdefault(gene, []) if excluded else results[gene]
            for alt_index, alt in enumerate(rec.ALT, start=1):
                # At a multi-allelic site, only keep the ALT(s) this sample
                # actually carries (e.g. GT 0/2 carries ALT 2, not ALT 1).
                if not _alt_copies(gt_str, alt_index):
                    continue
                target.append(
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
                        alt_index=alt_index,
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
        except OSError:
            continue
        out.extend(rec for rec in vcf if rec.CHROM == candidate and start <= rec.POS <= end)
        if out:
            break
    return out


def _extract_depth(rec) -> int | None:
    # INFO/DP can be absent (None), multi-valued (a tuple) or non-numeric;
    # FORMAT/DP raises KeyError when the header doesn't declare it.
    with contextlib.suppress(TypeError, ValueError):
        dp = rec.INFO.get("DP")
        if dp is not None:
            return int(dp)
    with contextlib.suppress(KeyError, IndexError, TypeError, ValueError):
        fmt_dp = rec.format("DP")
        if fmt_dp is not None and len(fmt_dp):
            return int(fmt_dp[0][0])
    return None


def _genotype_string(rec) -> str:
    try:
        gt = rec.genotypes[0]
        alleles = "/".join(str(a) for a in gt[:-1])
        phased = "|" if len(gt) > 2 and gt[-1] else "/"
        return alleles.replace("/", phased) if phased == "|" else alleles
    except (IndexError, TypeError):  # no sample columns
        return "./."


def _has_alt_allele(genotype: str) -> bool:
    """True if a genotype string carries at least one non-reference,
    non-missing allele (e.g. "0/1", "1/1" -> True; "0/0", "./." -> False).

    A gene region VCF slice lists every polymorphic site seen across an
    entire population, not just sites where this one sample carries the
    alternate allele -- without this check, a homozygous-reference (0/0)
    call at a defining position would be misread as evidence the sample
    carries that star allele. Caught by running the real GeT-RM
    validation harness against actual 1000 Genomes population data,
    where it silently produced the same wrong diplotype for every
    sample; invisible against the hand-authored synthetic test fixture,
    which only ever included rows where the variant was truly present.
    """
    alleles = genotype.replace("|", "/").split("/")
    return any(a not in ("0", ".") for a in alleles)


def _alt_copies(genotype: str, alt_index: int = 1) -> int:
    """Number of copies of ALT allele ``alt_index`` in a genotype string
    ("0/1" -> 1, "1|1" -> 2, "1/2" -> 1 for either ALT, "0/0" -> 0)."""
    return genotype.replace("|", "/").split("/").count(str(alt_index))


def _alt_haplotypes(genotype: str, alt_index: int = 1) -> tuple[bool, ...] | None:
    """For a phased genotype, which haplotypes carry the ALT allele
    ("1|0" -> (True, False)); None when the genotype is unphased."""
    if "|" not in genotype:
        return None
    return tuple(a == str(alt_index) for a in genotype.split("|"))


def _split_rsids(rsid: str | None) -> list[str]:
    # Some callsets join several IDs with ";" (e.g. "rs123;rs456").
    return [r for r in (rsid or "").split(";") if r and r != "."]


def _allele_sort_key(allele: str) -> tuple:
    if allele.startswith("*"):
        digits = ""
        for ch in allele[1:]:
            if not ch.isdigit():
                break
            digits += ch
        if digits:
            return (0, int(digits), allele)
    return (1, 0, allele)


def normalize_diplotype(diplotype: str) -> str:
    """Put a diplotype's two alleles in a canonical order ("*2A/*1" ->
    "*1/*2A"), so calls can be compared with consensus genotypes regardless
    of which allele was written first."""
    alleles = diplotype.split("/")
    if len(alleles) != 2:
        return diplotype
    return "/".join(sorted(alleles, key=_allele_sort_key))


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
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
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

    # rsID -> [(allele_name, definition, matched_via_proxy)]. Several alleles
    # can share an rsID (UGT1A1 repeat lengths), and one allele can have
    # several IDs (aliases, LD proxies).
    lookup: dict[str, list[tuple[str, dict, bool]]] = {}
    for name, info in defs.items():
        for rsid in [info.get("rsid"), *info.get("alias_rsids", [])]:
            if rsid:
                lookup.setdefault(rsid, []).append((name, info, False))
        for rsid in info.get("proxy_rsids", []):
            lookup.setdefault(rsid, []).append((name, info, True))

    copies: dict[str, int] = {}
    haplotypes: dict[str, tuple[bool, ...] | None] = {}
    seen_directly: set[str] = set()

    for v in variants:
        n = _alt_copies(v.genotype, v.alt_index)
        if not n:
            continue
        for rsid in _split_rsids(v.rsid):
            for name, info, via_proxy in lookup.get(rsid, []):
                # The length check applies to the repeat itself, not to an
                # SNV proxy for it (UGT1A1 *80 for *28).
                change = info.get("indel_length_change")
                if change is not None and not via_proxy and len(v.alt) - len(v.ref) != change:
                    continue
                # Several markers for one allele count once, not as two
                # separate allele copies.
                copies[name] = max(copies.get(name, 0), n)
                haplotypes.setdefault(name, _alt_haplotypes(v.genotype, v.alt_index))
                if not via_proxy:
                    seen_directly.add(name)

    ambiguous: list[str] = []
    proxy_only: list[str] = []

    for name, info in defs.items():
        if name in copies and name not in seen_directly and info.get("proxy_note"):
            ambiguous.append(info["proxy_note"])
            proxy_only.append(name)

    # Haplotypes made of other alleles in cis (TPMT *3A = *3B + *3C).
    for name, info in defs.items():
        parts = info.get("cis")
        if not parts or not all(p in copies for p in parts):
            continue
        n_cis = min(copies[p] for p in parts)
        if all(copies[p] == 1 for p in parts):
            phases = [haplotypes.get(p) for p in parts]
            if all(phases):
                if len(set(phases)) > 1:
                    n_cis = 0  # phased onto different chromosomes: in trans
            else:
                ambiguous.append(
                    f"{name} assumed ({' + '.join(parts)} in cis, the common "
                    f"arrangement); {'/'.join(parts)} in trans is also possible "
                    "and would be a poor metabolizer -- phasing needed to confirm"
                )
        if n_cis:
            copies[name] = n_cis
            for p in parts:
                copies[p] -= n_cis
                if not copies[p]:
                    del copies[p]

    alleles = sorted(a for name, n in copies.items() for a in [name] * n)

    if not alleles:
        return DiplotypeCall(gene=gene, diplotype="*1/*1", method="rule_based_vcf")

    if len(alleles) == 1:
        diplotype = f"{alleles[0]}/*1"
    else:
        diplotype = "/".join(alleles[:2])
        distinct = sorted(set(alleles))
        in_trans = len(distinct) == 2 and all(haplotypes.get(a) for a in distinct) and (
            haplotypes[distinct[0]] != haplotypes[distinct[1]]
        )
        if len(alleles) > 2 or (len(distinct) > 1 and not in_trans and not any(
            defs[a].get("cis") for a in distinct
        )):
            # Without phasing a VCF cannot tell whether variants are in cis
            # (same allele) or trans (two different variant alleles). Report
            # the naive combination as the top candidate and list alternatives.
            ambiguous.append(f"phase unresolved among: {', '.join(distinct)}")

    return DiplotypeCall(
        gene=gene,
        diplotype=diplotype,
        method="rule_based_vcf",
        matched_variants=variants,
        ambiguous_candidates=ambiguous,
        proxy_only_alleles=proxy_only,
    )
