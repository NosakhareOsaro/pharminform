"""End-to-end pipeline: VCF in, ReportContext out.

Shared by the CLI, FastAPI endpoint, and Streamlit demo so the same
calling -> phenotyping -> guideline-lookup -> composite-risk logic isn't
duplicated three times across the interface layers.
"""

from __future__ import annotations

from pathlib import Path

from . import caller
from .composite_risk import NoPathwayDataError, load_drug_pathway_weights, score_drug
from .confidence import assess_confidence
from .cpic_dpwg import GuidelineLookup
from .pharmgkb_client import diplotype_to_phenotype
from .report import GeneReportRow, ReportContext


def run_pipeline(
    vcf_path: str | Path,
    min_depth: int = 20,
    patient_id: str | None = None,
    genes: list[str] | None = None,
) -> ReportContext:
    genes = genes or caller.TARGET_GENES
    dropped_by_gene: dict[str, list[caller.VariantCall]] = {}
    variants_by_gene = caller.load_variants(vcf_path, genes=genes, min_depth=min_depth, dropped=dropped_by_gene)

    gene_rows: list[GeneReportRow] = []
    gene_phenotypes: dict[str, str] = {}

    for gene in genes:
        if gene in caller.VCF_UNSUPPORTED_GENES:
            gene_rows.append(
                GeneReportRow(
                    gene=gene, diplotype="not called", phenotype=None, calling_method="none",
                    confidence_level="low", confidence_reasons=[caller.VCF_UNSUPPORTED_GENES[gene]],
                )
            )
            continue
        gene_variants = variants_by_gene.get(gene, [])
        call = caller.call_diplotype_from_vcf(gene, gene_variants)
        conf = assess_confidence(
            call, all_gene_variants=gene_variants, min_depth=min_depth,
            dropped_variants=dropped_by_gene.get(gene),
        )

        allele_functions = {
            name: info["function"]
            for name, info in caller.BOOTSTRAP_ALLELE_DEFINITIONS.get(gene, {}).items()
        }
        try:
            phenotype, _ = diplotype_to_phenotype(gene, call.diplotype, allele_functions)
        except ValueError:
            phenotype = None
        if call.proxy_only_alleles:
            # CPIC lists UGT1A1 *80 without the TA repeat itself as
            # indeterminate; no guideline is applied to it.
            phenotype = "indeterminate"

        if phenotype and phenotype != "indeterminate":
            gene_phenotypes[gene] = phenotype

        gene_rows.append(
            GeneReportRow(
                gene=gene,
                diplotype=call.diplotype,
                phenotype=phenotype,
                calling_method=call.method,
                confidence_level=conf.level,
                confidence_reasons=conf.reasons,
            )
        )

    guideline_lookup = GuidelineLookup()
    cpic_recommendations = []
    for drug in guideline_lookup.available_drugs():
        rec = guideline_lookup.lookup(drug, gene_phenotypes)
        if rec is not None:
            cpic_recommendations.append(
                {
                    "drug": rec.drug,
                    "phenotypes": rec.phenotypes,
                    "recommendation": rec.recommendation,
                    "classification": rec.classification,
                    "citations": rec.citations,
                    "assumed_phenotypes": rec.assumed_phenotypes,
                }
            )

    weights_data = load_drug_pathway_weights()
    composite_results = []
    for entry in weights_data["drugs"]:
        drug = entry["drug"]
        relevant_genes = [g["gene"] for g in entry["contributing_genes"] if g["gene"] in gene_phenotypes]
        if not relevant_genes:
            continue
        try:
            result = score_drug(drug, {g: gene_phenotypes[g] for g in relevant_genes}, weights_data)
        except NoPathwayDataError:
            continue
        composite_results.append(result.to_dict()["composite_risk"] | {"drug": drug})

    return ReportContext(
        patient_id=patient_id,
        gene_results=gene_rows,
        cpic_recommendations=cpic_recommendations,
        composite_results=composite_results,
    )
