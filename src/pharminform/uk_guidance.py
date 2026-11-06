"""UK-specific dosing guidance, versioned alongside the CPIC layer.

CPIC is an international consortium; UK labs and oncology pharmacy teams
work to UK documents, and the two do not always agree on the starting
dose. CPIC's DPYD guideline, for example, expresses fluoropyrimidine
dosing through a DPYD activity score, while the UK Chemotherapy Board
states dose percentages and an escalation schedule directly. Rather than
overwrite one with the other, this module serves the UK rule *beside* the
CPIC rule, each labelled with its source, issuer, version and date, so a
report can show both and a reader can see which body said what.

Two sources are encoded in ``data/uk_guidance.json``:

* **UK Chemotherapy Board (2020)** -- DPYD and fluoropyrimidines. Only the
  four variants the UK document tables (*2A, *13, c.2846A>T, HapB3) are
  in scope; the document's own caveat that a negative result does not
  exclude DPD deficiency from other causes is carried through.
* **CERSI-PGx (2026)** -- UGT1A1 and irinotecan. Note its scope limits:
  epithelial malignancy, non-pegylated irinotecan, no recommendation for
  rare malignancies such as sarcoma, and no dose reduction for
  intermediate metabolisers.

Nothing here is extrapolated past those texts. A genotype or phenotype
the source does not table returns ``None`` rather than a guessed dose.
A later edition of either guidance should be added as a new source id
with its own rules, not edited into these.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

DATA_PATH = Path(__file__).parent / "data" / "uk_guidance.json"

# Allele names as pharminform.caller reports them for DPYD.
DPYD_NO_FUNCTION = ("*2A", "*13")
DPYD_DECREASED_FUNCTION = ("c.2846A>T", "c.1236G>A(HapB3)")


@dataclass
class GuidanceSource:
    id: str
    name: str
    issuer: str
    version: str
    date: str
    scope: str
    notes: str = ""
    url: str | None = None
    doi: str | None = None
    pmid: str | None = None

    @property
    def label(self) -> str:
        return f"{self.issuer} {self.version}"


@dataclass
class UKRecommendation:
    """One UK dosing decision for one drug, with its provenance."""

    drug: str
    gene: str
    genotype: str
    group_label: str
    advice: str
    source: GuidanceSource
    start_dose_percent: int | None = None
    max_dose_percent: int | None = None
    max_dose_percent_with_caution: int | None = None
    dpd_activity: str | None = None
    titration: str | None = None
    alternatives: list[str] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)
    snomed_code: str | None = None

    def to_dict(self) -> dict:
        return {
            "drug": self.drug,
            "gene": self.gene,
            "genotype": self.genotype,
            "group": self.group_label,
            "advice": self.advice,
            "start_dose_percent": self.start_dose_percent,
            "max_dose_percent": self.max_dose_percent,
            "max_dose_percent_with_caution": self.max_dose_percent_with_caution,
            "dpd_activity": self.dpd_activity,
            "titration": self.titration,
            "alternatives": list(self.alternatives),
            "caveats": list(self.caveats),
            "snomed_code": self.snomed_code,
            "source": {
                "id": self.source.id,
                "name": self.source.name,
                "issuer": self.source.issuer,
                "version": self.source.version,
                "date": self.source.date,
                "url": self.source.url,
                "doi": self.source.doi,
                "pmid": self.source.pmid,
            },
        }


def load_uk_guidance(path: Path | str = DATA_PATH) -> dict:
    return json.loads(Path(path).read_text())


def _source(data: dict, source_id: str) -> GuidanceSource:
    raw = data["sources"][source_id]
    return GuidanceSource(
        id=source_id,
        name=raw["name"],
        issuer=raw["issuer"],
        version=raw["version"],
        date=raw["date"],
        scope=raw["scope"],
        notes=raw.get("notes", ""),
        url=raw.get("url"),
        doi=raw.get("doi"),
        pmid=raw.get("pmid"),
    )


def available_sources(data: dict | None = None) -> list[GuidanceSource]:
    data = data if data is not None else load_uk_guidance()
    return [_source(data, sid) for sid in sorted(data["sources"])]


def _split_diplotype(diplotype: str) -> list[str]:
    return [a.strip() for a in diplotype.split("/") if a.strip()]


def classify_dpyd(diplotype: str) -> tuple[int, int] | None:
    """(no-function alleles, decreased-function alleles) in a DPYD diplotype.

    ``None`` if the diplotype contains an allele outside the four the UK
    guidance tables, since the document says nothing about it.
    """
    alleles = _split_diplotype(diplotype)
    if len(alleles) != 2:
        return None
    no_function = decreased = 0
    for allele in alleles:
        if allele in DPYD_NO_FUNCTION:
            no_function += 1
        elif allele in DPYD_DECREASED_FUNCTION:
            decreased += 1
        elif allele != "*1":
            return None
    return no_function, decreased


def _matches(group: dict, no_function: int, decreased: int) -> bool:
    match = group["match"]
    if "total_variant_alleles" in match:
        return no_function + decreased >= match["total_variant_alleles"]
    return match["no_function"] == no_function and match["decreased_function"] == decreased


def dpyd_recommendations(
    diplotype: str, data: dict | None = None, drugs: list[str] | None = None
) -> list[UKRecommendation]:
    """UK Chemotherapy Board fluoropyrimidine advice for a DPYD diplotype."""
    data = data if data is not None else load_uk_guidance()
    block = data["dpyd"]
    counts = classify_dpyd(diplotype)
    if counts is None:
        return []
    no_function, decreased = counts

    # Two or more variant alleles first: the guidance treats loss of
    # activity as additive, so a compound genotype is its own group and
    # must not fall through to a heterozygous rule.
    groups = sorted(
        block["groups"],
        key=lambda g: -g["match"].get("total_variant_alleles", 0),
    )
    group = next((g for g in groups if _matches(g, no_function, decreased)), None)
    if group is None:
        return []

    source = _source(data, block["source"])
    titration = block["titration"]["text"] if group.get("titration") else None
    caveats = [source.notes] if source.notes else []
    if group["id"] == "normal":
        caveats.append(
            "Only *2A, *13, c.2846A>T and c.1236G>A (HapB3) are assessed; "
            "a negative result does not exclude DPD deficiency from other causes."
        )
    if no_function + decreased >= 2:
        caveats.append(block["combination_rule"])

    wanted = [d.lower() for d in drugs] if drugs else block["drugs"]
    return [
        UKRecommendation(
            drug=drug,
            gene="DPYD",
            genotype=diplotype,
            group_label=group["label"],
            advice=group["advice"],
            source=source,
            start_dose_percent=group.get("start_dose_percent"),
            max_dose_percent=group.get("max_dose_percent"),
            max_dose_percent_with_caution=group.get("max_dose_percent_with_caution"),
            dpd_activity=group.get("dpd_activity"),
            titration=titration,
            alternatives=list(group.get("alternatives", [])),
            caveats=caveats,
        )
        for drug in block["drugs"]
        if drug in wanted
    ]


def ugt1a1_recommendation(
    phenotype: str, diplotype: str | None = None, data: dict | None = None
) -> UKRecommendation | None:
    """CERSI-PGx irinotecan advice for a UGT1A1 phenotype."""
    data = data if data is not None else load_uk_guidance()
    block = data["ugt1a1"]
    entry = block["phenotypes"].get((phenotype or "").lower())
    if entry is None:
        return None
    source = _source(data, block["source"])
    caveats = [source.notes, block["testing"]]
    if diplotype and any(a == "*36" for a in _split_diplotype(diplotype)):
        caveats.append(f"UGT1A1 *36: {block['allele_notes']['*36']}")
    return UKRecommendation(
        drug="irinotecan",
        gene="UGT1A1",
        genotype=diplotype or "",
        group_label=entry["label"],
        advice=entry["advice"],
        source=source,
        start_dose_percent=entry.get("start_dose_percent"),
        caveats=[c for c in caveats if c],
        snomed_code=entry.get("snomed"),
    )


def uk_recommendations(
    gene_phenotypes: dict[str, str],
    gene_diplotypes: dict[str, str] | None = None,
    data: dict | None = None,
) -> list[UKRecommendation]:
    """Every UK recommendation the called genes support, in drug order."""
    data = data if data is not None else load_uk_guidance()
    gene_diplotypes = gene_diplotypes or {}
    out: list[UKRecommendation] = []
    if "DPYD" in gene_diplotypes:
        out.extend(dpyd_recommendations(gene_diplotypes["DPYD"], data))
    if "UGT1A1" in gene_phenotypes:
        rec = ugt1a1_recommendation(
            gene_phenotypes["UGT1A1"], gene_diplotypes.get("UGT1A1"), data
        )
        if rec is not None:
            out.append(rec)
    return sorted(out, key=lambda r: (r.drug, r.gene))
