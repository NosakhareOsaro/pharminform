---
title: "PharmInform: A Python Platform for Star Allele Calling, CPIC Dosing Guidance, and Multi-Gene Composite Pharmacogenomic Risk Scoring"
tags:
  - Python
  - pharmacogenomics
  - bioinformatics
  - clinical genomics
  - CPIC
authors:
  - name: Nosakhare Osaro
    orcid: 0009-0004-9778-1161
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 05 July 2026
bibliography: paper.bib
---

# Summary

PharmInform is a Python platform for clinical pharmacogenomics (PGx) that
takes a patient VCF (or, for copy-number-aware calling, aligned reads) and
produces star allele calls, CPIC-standard metabolizer phenotypes, and
CPIC dosing recommendations for six pharmacogenes commonly implicated in
adverse drug reactions and treatment failure: _CYP2D6_, _CYP2C19_,
_DPYD_, _UGT1A1_, _TPMT_, and _CYP2C9_. It is delivered as an installable
Python package, a command-line tool, a FastAPI REST endpoint, and a
Streamlit web demo, all sharing one pipeline
(`pharminform.pipeline.run_pipeline`). Beyond reproducing standard
per-gene PGx interpretation, PharmInform introduces an experimental
multi-gene composite risk score for drugs whose metabolism depends on
more than one of these genes -- a scenario existing open-source PGx
tools do not address (see Statement of Need).

# Statement of Need

Star allele calling from sequencing data is a solved-enough problem to
have several mature open-source tools: Aldy performs copy-number-aware
allelic decomposition from aligned reads and explicitly models structural
variation such as _CYP2D6_ deletions, duplications, and gene-conversion
hybrids [@numanagic2018aldy]; PharmCAT annotates VCF-derived diplotypes
against CPIC guidelines and generates a prescribing report
[@sangkuhl2020pharmcat]; StellarPGx uses genome-graph-based variant
detection specifically tuned for cytochrome P450 genes and reports the
highest published CYP2D6/GeT-RM concordance among compared tools (99%)
[@twesigomwe2021stellarpgx]. PharmInform does not attempt to out-perform
these tools at their core task -- it wraps Aldy directly for BAM/CRAM
input rather than reinventing copy-number calling, and is explicit in its
own documentation that its VCF-only fallback path (used when only a VCF,
not aligned reads, is available) cannot detect structural variation and
is always reported at reduced confidence.

What none of Aldy, PharmCAT, or StellarPGx do is combine per-gene
phenotypes into a single risk signal for drugs metabolized by _more than
one_ pharmacogene. Amitriptyline is a representative case: CYP2D6
hydroxylates it to a less-active metabolite while CYP2C19 demethylates it
to the active metabolite nortriptyline, so a patient's combined
CYP2D6/CYP2C19 status -- not either gene alone -- determines their actual
exposure risk. CPIC's own guideline for tricyclic antidepressants
discusses this qualitatively [@hicks2017tca] but, like the tools above,
does not publish a combined quantitative score. PharmInform's
`composite_risk` module fills this specific, narrow gap: an experimental,
transparently-documented combination of per-gene activity into one score,
never presented as a replacement for per-gene CPIC guidance.

A second, related gap PharmInform addresses is confidence flagging:
whether a given star allele call should be trusted at all. `confidence.py`
flags low-coverage defining positions, variants absent from the allele
definition table, phasing ambiguity, and -- specific to this project's
VCF-only fallback path -- genes with known structural variation that a
VCF cannot represent.

# Implementation

The composite score is deliberately simple and auditable rather than a
black-box model. Each gene's CPIC phenotype is mapped to a normalized
activity value (poor = 0, intermediate = 0.5, normal = 1.0, capped at 1.0
-- modeling reduced-function toxicity accumulation only). Each
contributing gene is weighted by its published fraction of the drug's
metabolism, curated in `data/drug_pathway_weights.json` with a PMID or
FDA-label citation for every non-null value (`DATA_SOURCES.md`); when no
published fraction exists for at least one gene, the model falls back to
equal weighting across the genes present, and this fallback is always
reported via a `weighting_method` field rather than silently assumed. The
weighted per-gene scores combine into a single `composite_function_score`
(0-1, lower = higher toxicity risk) and a three-tier `risk_tier`
(Low/Moderate/High), with thresholds chosen so that a single
poor-metabolizer gene at ≥50% effective weight can never resolve to
"Low." A separate `directional_flag` fires whenever any contributing gene
is an ultrarapid metabolizer (under-exposure or dangerous prodrug
over-activation risk, e.g. codeine to morphine via CYP2D6) --
deliberately not blended into the composite score, since it reflects a
distinct physiological mechanism. Three drugs in the curated dataset
(codeine, tramadol, clopidogrel) are prodrug-activation pathways where
this composite framing's directionality is inverted; this is flagged
per-drug rather than silently mis-scored, and is a known limitation (see
below). Every result carries a mandatory disclaimer identifying it as
exploratory and not a CPIC/DPWG guideline output.

# Validation

Monotonicity -- the property that a dominant poor-metabolizer gene must
never resolve to "Low" risk -- is validated by a parametrized test across
every drug in the reference dataset. Face validity is checked against
three real, published case reports: compound-heterozygous DPYD deficiency
causing grade 4 fluorouracil toxicity [@mukherji2019dpydcase], and
TPMT\*3A and TPMT\*3C homozygous deficiency causing azathioprine-induced
myelosuppression [@budhiraja2011tpmt3acase; @lee2013tpmt3ccase]; the
model correctly flags Moderate/High risk in all three.

Star allele calling accuracy was validated against GeT-RM (Genetic
Testing Reference Materials Coordination Program) consensus genotypes
[@pratt2016getrm; @gaedigk2019cyp2d6getrm] using real 1000 Genomes
Project sequencing data. Of 11 published CYP2D6 consensus samples
[@gaedigk2019cyp2d6getrm], 6 are present in the 1000 Genomes phase 3
panel (confirmed against both the official sample panel file and the
sequencing data's own sample list); the remaining 5 -- including
NA17244, the most widely-cited CYP2D6 duplication/hybrid reference
sample -- are not, and were excluded rather than guessed at. For the 6
panel-covered samples, real genotypes were fetched via `bcftools` and
passed through the same `caller.call_diplotype_from_vcf` function used
throughout the pipeline: 0/6 matched their published consensus
diplotype (0.0% concordance). Every mismatch is individually
attributable to one of two already-documented limitations, not an
unexplained failure: the illustrative bootstrap allele table covers
only 5 of CYP2D6's 150+ named haplotypes, and VCF-only matching cannot
represent the copy-number and gene-conversion alleles (e.g. `*36`, `x2`
duplication) that dominate this gene's structural variation --
precisely the scenario for which this project recommends Aldy against
aligned reads instead. Running this validation for real also surfaced
and fixed a genuine defect in the caller itself, where homozygous-
reference genotypes were misread as evidence of the variant allele --
invisible against hand-authored test fixtures, but immediately apparent
against real population data. Because the CYP2D6 causes are structural,
CYP2D6 is no longer called from VCF input. For the
other five genes, consensus genotypes from the CDC GeT-RM consolidated
table were translated to the calls the allele table can make and run
against the 1000 Genomes 30x high-coverage panel
[@byrskabishop2022highcoverage]: DPYD matched 94/94 diplotypes, TPMT
87/87, UGT1A1 59/60 (60/60 phenotypes), CYP2C9 81/81 and CYP2C19 88/88.
Most DPYD consensus calls come from the 2016 study [@pratt2016getrm],
which did not test c.2846A>T or HapB3, so all four DPYD variants are
compared only for the 9 samples from the DPYD study
[@gaedigk2024dpydgetrm], including all 7 carriers. The UGT1A1 mismatch is
a \*28/\*37 sample called \*28/\*28 because the \*80 proxy overrides the
directly observed \*37. Carriers remain few (7 DPYD, 9 TPMT, 35 UGT1A1,
19 CYP2C9, 50 CYP2C19). Separately, the phenotype model was checked
against CPIC's diplotype-to-phenotype tables for all 82 diplotypes the
allele table can call, which corrected the CYP2D6 \*10 and \*41 activity
values and the UGT1A1 \*36 function. Full per-sample results are in
`validation/getrm_concordance_report.md`.

# Availability

PharmInform is available under the MIT license at
`https://github.com/NosakhareOsaro/pharminform`, with a test suite
of over 100 tests, continuous integration, a Dockerfile, and a Bioconda
recipe skeleton (not yet submitted -- see README). The v0.1.0 release is
archived on Zenodo at `https://doi.org/10.5281/zenodo.21210588`.

# References
