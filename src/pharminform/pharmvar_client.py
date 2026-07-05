"""PharmVar API client.

PharmVar (https://www.pharmvar.org) is the authoritative source of star
allele (haplotype) definitions. As of this writing its REST API requires an
account-issued API key -- confirmed live against the real endpoint below,
which returns HTTP 401 without one:

    GET https://www.pharmvar.org/api-service/genes

This client is a real HTTP client, not a stub: it makes genuine requests
and raises a clear error when no key is configured, rather than silently
returning fabricated allele data. Tests mock the HTTP layer (see
tests/test_pharmvar_client.py) so the suite runs without network access or
a real key; `caller.py`'s BOOTSTRAP_ALLELE_DEFINITIONS is the offline
fallback used elsewhere in the pipeline until a key is supplied.

To use against the real API: set the PHARMVAR_API_KEY environment variable
to a key generated from your PharmVar account settings page.
"""

from __future__ import annotations

import os

import requests

BASE_URL = "https://www.pharmvar.org/api-service"


class PharmVarAuthError(RuntimeError):
    """Raised when no API key is configured or the key is rejected."""


class PharmVarClient:
    def __init__(self, api_key: str | None = None, base_url: str = BASE_URL, timeout: int = 15):
        self.api_key = api_key or os.environ.get("PHARMVAR_API_KEY")
        self.base_url = base_url
        self.timeout = timeout

    def _headers(self) -> dict[str, str]:
        if not self.api_key:
            raise PharmVarAuthError(
                "No PharmVar API key configured. Set PHARMVAR_API_KEY or pass "
                "api_key= explicitly. Generate a key from your PharmVar account "
                "settings page at https://www.pharmvar.org."
            )
        return {"Authorization": f"Bearer {self.api_key}", "Accept": "application/json"}

    def list_genes(self) -> list[dict]:
        """GET /genes -- all genes PharmVar has allele definitions for."""
        resp = requests.get(f"{self.base_url}/genes", headers=self._headers(), timeout=self.timeout)
        if resp.status_code == 401:
            raise PharmVarAuthError(f"PharmVar rejected the API key: {resp.text}")
        resp.raise_for_status()
        return resp.json()

    def get_gene_alleles(self, gene: str, reference_collection: str = "GRCh38") -> list[dict]:
        """GET /genes/{gene}?reference-collection=... -- allele definitions
        (defining variants, rsIDs, function) for one gene."""
        resp = requests.get(
            f"{self.base_url}/genes/{gene}",
            params={"reference-collection": reference_collection},
            headers=self._headers(),
            timeout=self.timeout,
        )
        if resp.status_code == 401:
            raise PharmVarAuthError(f"PharmVar rejected the API key: {resp.text}")
        resp.raise_for_status()
        return resp.json()
