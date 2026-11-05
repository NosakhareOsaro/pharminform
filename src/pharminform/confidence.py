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

from .caller import BOOTSTRAP_ALLELE_DEFINITIONS, DiplotypeCall, VariantCall

# Genes in this panel with clinically important structural variation
# (whole-gene deletions/duplications, gene-conversion/hybrid alleles) that
# a VCF's variant calls cannot represent.
GENES_REQUIRING_CN_CALLING = {"CYP2D6", "UGT1A1"}

DEFAULT_MIN_DEPTH = 20


@dataclass
class ConfidenceResult:
    level: str  # "high" or "low"
    reasons: list[str] = field(default_factory=list)

    @property
    def is_low(self) -> bool:
        return self.level == "low"


def assess_confidence(
    call: DiplotypeCall,
    all_gene_variants: list[VariantCall] | None = None,
    min_depth: int = DEFAULT_MIN_DEPTH,
    dropped_variants: list[VariantCall] | None = None,
) -> ConfidenceResult:
    """Assess confidence in a single gene's diplotype call.

    ``all_gene_variants`` should be the full (pre-defining-variant-match)
    variant list for the gene, if available, so novel/unmatched variants
    can be detected -- not just the ones that happened to match a known
    allele. ``dropped_variants`` are records ``load_variants`` excluded for
    low depth or a failed FILTER; any at a defining site is reported, since
    the call cannot see that allele.
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

    known_rsids = {
        rsid
        for info in BOOTSTRAP_ALLELE_DEFINITIONS.get(call.gene, {}).values()
        for rsid in [info.get("rsid"), *info.get("alias_rsids", []), *info.get("proxy_rsids", [])]
        if rsid
    }

    lost = [
        v for v in dropped_variants or []
        if v.rsid and any(r in known_rsids for r in v.rsid.split(";"))
    ]
    if lost:
        positions = ", ".join(
            f"{v.chrom}:{v.pos} ({v.rsid}, DP={v.depth}, FILTER={v.filter or 'PASS'})" for v in lost
        )
        reasons.append(
            f"defining variant(s) present but excluded for low depth (<{min_depth}x) or a failed "
            f"FILTER: {positions} -- the diplotype cannot include this allele; review before use"
        )

    if all_gene_variants is not None:
        novel = [
            v for v in all_gene_variants
            if v.rsid not in (None, ".")
            and not any(r in known_rsids for r in v.rsid.split(";"))
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
    return ConfidenceResult(level=level, reasons=reasons)
