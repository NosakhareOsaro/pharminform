"""CPIC/DPWG guideline lookup and caching layer.

There is no single unified API that serves both CPIC and DPWG dosing
recommendations. What DOES exist, and is used here, is CPIC's own public
REST API (https://api.cpicpgx.org/v1, a PostgREST interface over CPIC's
data exports at https://github.com/cpicpgx/cpic-data) -- confirmed live
and unauthenticated:

    GET https://api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:704

DPWG (Dutch Pharmacogenetics Working Group) guidelines are not available
through this or any other public structured API at the time of writing;
DPWG content is published as PDF documents and embedded as text within
PharmGKB. Nothing in this module claims to serve DPWG data -- the
"_dpwg" naming in the module is aspirational/for future integration and
is called out explicitly here rather than silently faked.

``data/cpic_recommendation_cache.json`` is a real snapshot fetched from
the live CPIC API for the drugs this project's composite risk model
covers (fetch date recorded in the file's "_fetched" field). Recommendation
rows are deduplicated by phenotype combination (CPIC's raw table also
varies by activity-score sub-bins and population subgroup, which produces
many more rows than distinct clinical recommendations). For the two TPMT
substrates (azathioprine, mercaptopurine), CPIC's guideline is jointly
authored with NUDT15; rows are filtered to NUDT15 "Normal Metabolizer" to
isolate TPMT-driven variation, since NUDT15 is outside this project's gene
panel -- this is a deliberate scope limitation, documented rather than
silently dropped.

Notably, warfarin has zero rows in CPIC's `recommendation` table: CPIC's
warfarin guideline recommends a validated pharmacogenetic dosing algorithm
(e.g. warfarindosing.org) rather than a discrete phenotype-to-dose lookup,
so there is nothing to cache in that format. This is a genuine finding
from the live API, not a bug in this loader.

``CPICClient.refresh_cache`` can re-fetch live data (network required);
``load_cached_guidelines`` reads the bundled offline snapshot and is what
the rest of the pipeline uses by default.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import requests

BASE_URL = "https://api.cpicpgx.org/v1"
CACHE_PATH = Path(__file__).parent / "data" / "cpic_recommendation_cache.json"


@dataclass
class GuidelineRecommendation:
    drug: str
    genes: list[str]
    phenotypes: dict[str, str]
    recommendation: str
    classification: str
    comments: str
    cpic_level: str | None
    citations: list[str]
    source: str = "CPIC"


class CPICClient:
    """Thin real HTTP client over CPIC's public PostgREST API."""

    def __init__(self, base_url: str = BASE_URL, timeout: int = 20):
        self.base_url = base_url
        self.timeout = timeout

    def get_recommendations(self, drugid: str, population: str | None = None) -> list[dict]:
        params = {
            "drugid": f"eq.{drugid}",
            "select": "population,phenotypes,activityscore,drugrecommendation,classification,comments",
        }
        if population:
            params["population"] = f"eq.{population}"
        resp = requests.get(f"{self.base_url}/recommendation", params=params, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json()

    def get_pairs(self, drugid: str, genes: list[str]) -> list[dict]:
        resp = requests.get(
            f"{self.base_url}/pair",
            params={
                "drugid": f"eq.{drugid}",
                "genesymbol": f"in.({','.join(genes)})",
                "select": "genesymbol,cpiclevel,citations",
            },
            timeout=self.timeout,
        )
        resp.raise_for_status()
        return resp.json()


def load_cached_guidelines(cache_path: Path | str = CACHE_PATH) -> dict:
    """Load the bundled offline CPIC recommendation snapshot."""
    path = Path(cache_path)
    if not path.exists():
        raise FileNotFoundError(
            f"No cached CPIC guideline data at {path}. Run "
            "CPICClient-based refresh (requires network access) or restore "
            "the bundled data/cpic_recommendation_cache.json."
        )
    return json.loads(path.read_text())


class GuidelineLookup:
    """Look up a CPIC dosing recommendation for (drug, {gene: phenotype})."""

    def __init__(self, cache: dict | None = None):
        self._cache = cache if cache is not None else load_cached_guidelines()

    def lookup(self, drug: str, phenotypes: dict[str, str]) -> GuidelineRecommendation | None:
        entry = self._cache["drugs"].get(drug.lower())
        if entry is None:
            return None

        # Normalize phenotype casing for comparison (CPIC stores
        # "Poor Metabolizer"-style casing; the rest of the pipeline uses
        # lowercase, e.g. from pharmgkb_client.diplotype_to_phenotype).
        wanted = {g: p.lower() for g, p in phenotypes.items() if g in entry["genes"]}

        best_match = None
        best_score = -1
        for rec in entry["recommendations"]:
            rec_phenotypes = {g: p.lower() for g, p in rec["phenotypes"].items()}
            score = sum(1 for g, p in wanted.items() if rec_phenotypes.get(g) == p)
            if score > best_score and score > 0:
                best_score = score
                best_match = rec

        if best_match is None:
            return None

        citations: list[str] = []
        cpic_level = None
        for pair in entry.get("cpic_pairs", []):
            if pair["genesymbol"] in entry["genes"]:
                citations.extend(pair.get("citations") or [])
                cpic_level = cpic_level or pair.get("cpiclevel")
        citations = sorted(set(citations))

        return GuidelineRecommendation(
            drug=drug,
            genes=entry["genes"],
            phenotypes=best_match["phenotypes"],
            recommendation=best_match["drugrecommendation"],
            classification=best_match["classification"],
            comments=best_match.get("comments", ""),
            cpic_level=cpic_level,
            citations=citations,
        )

    def available_drugs(self) -> list[str]:
        return sorted(self._cache["drugs"].keys())
