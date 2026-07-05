"""PharmGKB API client and diplotype -> phenotype conversion.

The PharmGKB REST API (https://api.pharmgkb.org) is publicly reachable
without an API key for the read endpoints used here (confirmed live, e.g.
``GET https://api.pharmgkb.org/v1/data/gene?symbol=CYP2D6`` and
``GET https://api.pharmgkb.org/v1/data/chemical?name=amitriptyline``).

PharmGKB does not expose a single "diplotype in, phenotype out" endpoint --
gene-specific allele function assignments are published as CPIC allele
definition/function tables (spreadsheets), not a queryable REST resource.
``diplotype_to_phenotype`` therefore implements the CPIC activity-score
model directly: each allele is assigned a numeric activity value from its
CPIC-published function ("no function" = 0, "decreased function" = 0.5,
"normal function" = 1.0, "increased function" = 2.0), the two allele
values are summed into a diplotype activity score, and the score is mapped
to a phenotype using threshold bands.

The exact activity-score cutoffs used clinically differ slightly per gene
per CPIC's gene-specific standardization publications (e.g. CYP2D6's
bands are not identical to CYP2C19's). The DEFAULT_ACTIVITY_BANDS below is
a single generic band set used as a reasonable approximation across all
target genes -- it is documented here as a simplification, not asserted as
the exact CPIC value for every gene. Gene-specific overrides can be added
to ACTIVITY_BAND_OVERRIDES as they are verified against CPIC's published
tables.
"""

from __future__ import annotations

import requests

BASE_URL = "https://api.pharmgkb.org/v1/data"

ALLELE_FUNCTION_TO_ACTIVITY_VALUE = {
    "no function": 0.0,
    "decreased function": 0.5,
    "normal function": 1.0,
    "increased function": 2.0,
    "unknown function": 1.0,  # conservative default: treat as normal
}

# (lower_bound_inclusive, upper_bound_inclusive, phenotype)
# Bounded by the max reachable score (two "increased function" alleles at
# 2.0 each = 4.0) so every band is actually reachable under this model.
DEFAULT_ACTIVITY_BANDS: list[tuple[float, float, str]] = [
    (0.0, 0.0, "poor metabolizer"),
    (0.25, 1.0, "intermediate metabolizer"),
    (1.25, 2.0, "normal metabolizer"),
    (2.25, 3.0, "rapid metabolizer"),
    (3.25, float("inf"), "ultrarapid metabolizer"),
]

# Per-gene overrides where the generic bands above are known to diverge
# from CPIC's published gene-specific standardization. TPMT and DPYD are
# conventionally reported as a three-tier (PM/IM/NM) system in CPIC
# guidance rather than the five-tier CYP2D6-style system.
ACTIVITY_BAND_OVERRIDES: dict[str, list[tuple[float, float, str]]] = {
    "TPMT": [
        (0.0, 0.0, "poor metabolizer"),
        (0.25, 1.0, "intermediate metabolizer"),
        (1.25, 4.0, "normal metabolizer"),
    ],
    "DPYD": [
        (0.0, 0.0, "poor metabolizer"),
        (0.25, 1.0, "intermediate metabolizer"),
        (1.25, 4.0, "normal metabolizer"),
    ],
}


class PharmGKBClient:
    def __init__(self, base_url: str = BASE_URL, timeout: int = 15):
        self.base_url = base_url
        self.timeout = timeout

    def get_gene(self, symbol: str) -> dict | None:
        """GET /gene?symbol=... -- gene metadata (coordinates, whether CPIC
        publishes guidance, PharmVar coverage, citation)."""
        resp = requests.get(
            f"{self.base_url}/gene", params={"symbol": symbol}, timeout=self.timeout
        )
        resp.raise_for_status()
        payload = resp.json()
        data = payload.get("data") or []
        return data[0] if data else None

    def get_chemical(self, name: str) -> dict | None:
        """GET /chemical?name=... -- drug metadata (PharmGKB accession id,
        cross-references)."""
        resp = requests.get(
            f"{self.base_url}/chemical", params={"name": name}, timeout=self.timeout
        )
        resp.raise_for_status()
        payload = resp.json()
        data = payload.get("data") or []
        return data[0] if data else None


def diplotype_activity_score(
    diplotype: str, allele_functions: dict[str, str]
) -> float:
    """Sum activity values for both haplotypes in a diplotype string like
    "*4/*10". Alleles not present in ``allele_functions`` (including "*1",
    which is reference/normal function by convention) default to 1.0.
    """
    alleles = diplotype.split("/")
    if len(alleles) != 2:
        raise ValueError(f"Expected a diplotype like '*1/*4', got: {diplotype!r}")
    total = 0.0
    for allele in alleles:
        function = allele_functions.get(allele, "normal function" if allele == "*1" else "unknown function")
        total += ALLELE_FUNCTION_TO_ACTIVITY_VALUE.get(function, 1.0)
    return total


def diplotype_to_phenotype(
    gene: str, diplotype: str, allele_functions: dict[str, str]
) -> tuple[str, float]:
    """Convert a diplotype to a CPIC-style phenotype using the activity
    score model. Returns (phenotype, activity_score)."""
    score = diplotype_activity_score(diplotype, allele_functions)
    bands = ACTIVITY_BAND_OVERRIDES.get(gene, DEFAULT_ACTIVITY_BANDS)
    for low, high, phenotype in bands:
        if low <= score <= high:
            return phenotype, score
    return "indeterminate", score
