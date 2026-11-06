"""Diplotype call confidence flagging.

Wraps a ``caller.DiplotypeCall`` and flags ``confidence: low`` with
specific, actionable reasons whenever the call should not be treated as
clinical-grade, covering the three failure modes named in the project
spec: low coverage at defining positions, a novel variant not present in
the allele definition table, and phasing ambiguity. A fourth reason is
added that's specific to this pipeline's two calling paths: a VCF-only
rule-based call is always flagged for genes with known structural
variation, since it cannot see copy number.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .caller import DEFINING_SITES_GRCH38, DiplotypeCall, VariantCall, known_rsids

# Genes in this panel with clinically important structural variation
# (whole-gene deletions/duplications, gene-conversion/hybrid alleles) that
# a VCF's variant calls cannot represent.
GENES_REQUIRING_CN_CALLING = {"CYP2D6", "UGT1A1"}

DEFAULT_MIN_DEPTH = 20


@dataclass
class ConfidenceResult:
    level: str  # "high" or "low"
    reasons: list[str] = field(default_factory=list)
    # True when an allele could be present but unseen (a defining site with
    # no genotype call, or a defining variant dropped for depth/FILTER). The
    # pipeline then reports the phenotype as indeterminate instead of
    # letting a missing allele read as normal function.
    phenotype_blocked: bool = False

    @property
    def is_low(self) -> bool:
        return self.level == "low"


def assess_confidence(
    call: DiplotypeCall,
    all_gene_variants: list[VariantCall] | None = None,
    min_depth: int = DEFAULT_MIN_DEPTH,
    dropped_variants: list[VariantCall] | None = None,
    site_status: dict[str, str] | None = None,
) -> ConfidenceResult:
    """Assess confidence in a single gene's diplotype call.

    ``all_gene_variants`` should be the full (pre-defining-variant-match)
    variant list for the gene, if available, so novel/unmatched variants
    can be detected -- not just the ones that happened to match a known
    allele. ``dropped_variants`` are records ``load_variants`` excluded for
    low depth or a failed FILTER; any at a defining site is reported, since
    the call cannot see that allele. ``site_status`` (from load_variants)
    maps each defining site's rsID to "called" or "no-call"; sites missing
    from it had no VCF record.
    """
    reasons: list[str] = []

    if call.method == "rule_based_vcf":
        reasons.append(
            "VCF-only rule-based call: no copy-number or phasing information available"
        )
        if call.gene in GENES_REQUIRING_CN_CALLING:
            reasons.append(
                f"{call.gene} has clinically important structural variation "
                "(deletions/duplications/hybrid alleles) that cannot be detected "
                "from variant calls alone; run Aldy on aligned reads (BAM/CRAM) "
                "for a clinical-grade call"
            )

    if call.ambiguous_candidates:
        reasons.append(
            "phasing ambiguity: " + "; ".join(call.ambiguous_candidates)
        )

    low_depth = [v for v in call.matched_variants if v.depth is not None and v.depth < min_depth]
    if low_depth:
        positions = ", ".join(f"{v.chrom}:{v.pos}" for v in low_depth)
        reasons.append(f"low coverage (<{min_depth}x) at defining position(s): {positions}")

    known = known_rsids(call.gene)
    blocked = False

    lost = [
        v for v in dropped_variants or []
        if v.rsid and any(r in known for r in v.rsid.split(";"))
    ]
    if lost:
        positions = ", ".join(
            f"{v.chrom}:{v.pos} ({v.rsid}, DP={v.depth}, FILTER={v.filter or 'PASS'})" for v in lost
        )
        reasons.append(
            f"defining variant(s) present but excluded for low depth (<{min_depth}x) or a failed "
            f"FILTER: {positions} -- the diplotype cannot include this allele; review before use"
        )
        blocked = True

    if site_status is not None:
        sites = [rsid for *_, rsid in DEFINING_SITES_GRCH38.get(call.gene, [])]
        nocall = [r for r in sites if site_status.get(r) == "no-call"]
        absent = [r for r in sites if r not in site_status]
        if nocall:
            reasons.append(
                f"no genotype call at defining site(s) {', '.join(nocall)} -- the allele "
                "cannot be ruled out; phenotype reported as indeterminate"
            )
            blocked = True
        if absent:
            reasons.append(
                f"no VCF record at {len(absent)} of {len(sites)} defining site(s) ({', '.join(absent)}); "
                "in a variant-only VCF this usually means reference, but coverage there is not shown"
            )

    if all_gene_variants is not None:
        novel = [
            v for v in all_gene_variants
            if v.rsid not in (None, ".")
            and not any(r in known for r in v.rsid.split(";"))
        ]
        unnamed = [v for v in all_gene_variants if v.rsid in (None, ".")]
        if novel:
            positions = ", ".join(f"{v.chrom}:{v.pos} ({v.rsid})" for v in novel)
            reasons.append(f"variant(s) not in the allele definition table: {positions}")
        if unnamed:
            positions = ", ".join(f"{v.chrom}:{v.pos}" for v in unnamed)
            reasons.append(
                f"variant(s) with no rsID (possible novel/private variant not in "
                f"PharmVar): {positions}"
            )

    level = "low" if reasons else "high"
    return ConfidenceResult(level=level, reasons=reasons, phenotype_blocked=blocked)
