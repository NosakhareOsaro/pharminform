"""Every diplotype the bootstrap allele table can call, checked against
CPIC's own diplotype-to-phenotype table (phenotype and, where CPIC gives
one, activity score).

The GeT-RM harness derives both sides' phenotypes with PharmInform's
model, so it cannot catch a wrong rule; this test can. The CPIC rows are
bundled in src/pharminform/data/cpic_diplotype_phenotype_panel.tsv (source
and retrieval date in its header) so the check runs offline.
"""

import csv

import pytest
from tests.getrm_oncology_consensus import DATA_DIR

from pharminform import caller
from pharminform.pharmgkb_client import diplotype_to_phenotype

CPIC_PANEL_TSV = DATA_DIR / "cpic_diplotype_phenotype_panel.tsv"


def _rows():
    with CPIC_PANEL_TSV.open() as fh:
        return list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))


def test_table_covers_every_callable_diplotype():
    expected = sum(
        (n + 1) * (n + 2) // 2 for n in (len(d) for d in caller.BOOTSTRAP_ALLELE_DEFINITIONS.values())
    )
    assert len(_rows()) == expected == 82


@pytest.mark.parametrize("row", _rows(), ids=lambda r: f"{r['gene']}:{r['pharminform_diplotype']}")
def test_phenotype_and_score_match_cpic(row):
    gene = row["gene"]
    functions = {n: i["function"] for n, i in caller.BOOTSTRAP_ALLELE_DEFINITIONS[gene].items()}
    phenotype, score = diplotype_to_phenotype(gene, row["pharminform_diplotype"], functions)
    assert phenotype == row["cpic_phenotype"].lower()
    if row["cpic_activity_score"] not in ("", "n/a"):
        assert score == pytest.approx(float(row["cpic_activity_score"]))
