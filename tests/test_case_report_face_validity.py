"""Face-validity tests: does the composite risk model flag Moderate/High
risk for patients whose phenotypes match real, published clinical case
reports of severe toxicity?

Every PMID below was verified against a real PubMed/PMC record via live
web search before being hardcoded here. A literature search specifically
for a combined CYP2D6+CYP2C19 poor-metabolizer amitriptyline case report
(the multi-gene scenario this project's composite model was built for)
did not turn up a citable case with both genotypes documented -- the
Medical Genetics Summary covering amitriptyline/CYP2D6/CYP2C19
(PMID:28520380) explicitly notes that "scientific literature has scarce
data" on combined-genotype patients. That absence is itself documented
here rather than papered over with an invented case; the amitriptyline
check below instead validates the model against the single-gene
(CYP2D6 PM) scenario that IS well-documented -- per the Medical Genetics
Summary at PMID:28520380, which itself summarizes FDA label content --
run through the same multi-gene composite pathway.
"""

from pharminform.composite_risk import load_drug_pathway_weights, score_drug

WEIGHTS_DATA = load_drug_pathway_weights()


def test_dpyd_compound_heterozygous_grade4_toxicity_case():
    """PMID:30898145 -- a patient heterozygous for three DPYD
    deficiency-associated variants (DPYD*2A / c.1601G>A / c.2194G>A)
    experienced grade 4 mucositis, febrile neutropenia, and skin
    desquamation after a single FOLFIRINOX (5-FU-containing) cycle.
    Modeled here as a poor metabolizer given the severity of observed
    toxicity from multiple deficiency alleles.
    """
    result = score_drug("fluorouracil", {"DPYD": "poor metabolizer"}, WEIGHTS_DATA)
    assert result.risk_tier in ("Moderate", "High")


def test_tpmt_homozygous_pancytopenia_case():
    """PMID:21712851 -- a patient homozygous for TPMT*3A (no-function/no-
    function, i.e. poor metabolizer) developed pancytopenia on standard-dose
    azathioprine. A second independent case with the same genotype-outcome
    pattern (TPMT*3C/*3C leukopenia) is reported at PMID:24142665.
    """
    result = score_drug("azathioprine", {"TPMT": "poor metabolizer"}, WEIGHTS_DATA)
    assert result.risk_tier in ("Moderate", "High")

    result_3c = score_drug("mercaptopurine", {"TPMT": "poor metabolizer"}, WEIGHTS_DATA)
    assert result_3c.risk_tier in ("Moderate", "High")


def test_amitriptyline_cyp2d6_poor_metabolizer_case():
    """PMID:28520380 (Dean L, Kane M -- NCBI Medical Genetics Summary,
    summarizing FDA label content) -- CYP2D6 poor metabolizers have
    higher-than-expected tricyclic plasma concentrations at usual doses
    and an increased risk of adverse effects; CPIC-summarized guidance
    recommends a 50% starting-dose reduction with therapeutic drug
    monitoring in this scenario. No combined CYP2D6+CYP2C19 case report
    with a citable PMID was found (see module docstring), so CYP2C19 is
    held at normal metabolizer to isolate the well-documented CYP2D6 effect.
    """
    result = score_drug(
        "amitriptyline",
        {"CYP2D6": "poor metabolizer", "CYP2C19": "normal metabolizer"},
        WEIGHTS_DATA,
    )
    assert result.risk_tier in ("Moderate", "High")
