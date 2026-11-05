"""GeT-RM consensus genotypes for DPYD, TPMT and UGT1A1, translated into
the call PharmInform's allele table is expected to make.

Source: CDC GeT-RM "Consolidated Table of all GeT-RM Pharmacogenetic and
HLA Reference Material Genotypes", V4 (1-22-25),
https://www.cdc.gov/lab-quality/media/files/2025/08/Consolidated_PGx-HLA_table_1-22-25-V4.xlsx
(listed on https://www.cdc.gov/lab-quality/php/get-rm/reference-materials.html).
``getrm_consolidated_v4_oncology.tsv`` is an unmodified extract of its
Coriell ID, CYP2D6, DPYD, NUDT15, TPMT and UGT1A1 columns (rows with any
DPYD/TPMT/UGT1A1 value; cells trimmed, tabs/newlines collapsed to spaces),
SHA-256 3a329a683efde66af1d1ea58eeb84b58508879aace4b9a697ecae78f55e0cb5f,
checked against the same extract taken from the downloaded workbook.

The consensus uses the full star-allele / HGVS nomenclature; the bootstrap
table only knows a few alleles per gene. Each sample is therefore either
translated to the expected *panel* call, or marked out of scope with a
reason -- never guessed:

* DPYD: the four CPIC actionable variants map to *2A, *13, c.2846A>T and
  HapB3; c.85T>C (*9A), c.496A>G, c.1601G>A (*4), c.1627A>G (*5),
  c.1896T>C and c.2194G>A (*6) are normal function in CPIC and are
  invisible to the panel (read as *1). Any other variant, deletion, or
  ambiguous star name puts the sample out of scope.
* TPMT: in scope only when every allele is *1, *2, *3A, *3B or *3C.
* UGT1A1: *60 (-3279T>G, not in the panel) reads as *1, and the *28 + *60
  haplotype as *28. Tentative calls in parentheses other than "(*60)", and
  alleles outside *1/*6/*28/*36/*37, put the sample out of scope.

Comparing the panel call is a test of what the caller claims to detect.
It says nothing about the alleles it cannot see; the out-of-scope count is
reported alongside so that limit stays visible.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

TSV_PATH = Path(__file__).parent.parent / "src" / "pharminform" / "data" / "getrm_consolidated_v4_oncology.tsv"

DPYD_ACTIONABLE = {
    "c.1905+1G>A": "*2A",
    "c.1679T>G": "*13",
    "c.2846A>T": "c.2846A>T",
    "c.1129-5923C>G": "c.1236G>A(HapB3)",
    "c.1236G>A": "c.1236G>A(HapB3)",
}
DPYD_PANEL_NORMAL = {"c.85T>C", "c.496A>G", "c.1601G>A", "c.1627A>G", "c.1896T>C", "c.2194G>A"}
DPYD_NORMAL_STARS = {"*1", "*9", "*4", "(*4)"}

TPMT_PANEL = {"*1", "*2", "*3A", "*3B", "*3C"}
UGT1A1_PANEL = {"*1", "*6", "*28", "*36", "*37"}


@dataclass
class ConsensusEntry:
    coriell_id: str
    gene: str
    consensus_raw: str
    consensus_diplotype: str | None  # expected panel call, None if out of scope
    scope_note: str = ""


def _diplotype(alleles: list[str]) -> str:
    # Same shape as caller.call_diplotype_from_vcf output; the harness
    # normalises allele order before comparing.
    alleles = sorted(a for a in alleles if a != "*1")
    if not alleles:
        return "*1/*1"
    if len(alleles) == 1:
        return f"{alleles[0]}/*1"
    return "/".join(alleles[:2])


def _dpyd_haplotype(text: str) -> tuple[list[str], str]:
    """Panel alleles on one DPYD haplotype (or an unphased variant list)."""
    text = text.strip()
    if text in ("no variants detected", "no variant detected"):
        return [], ""
    found: list[str] = []
    for item in [t.strip() for t in text.split(",") if t.strip()]:
        variant = item.replace(" het", "").strip()
        if variant in DPYD_ACTIONABLE:
            if DPYD_ACTIONABLE[variant] not in found:
                found.append(DPYD_ACTIONABLE[variant])
        elif variant not in DPYD_PANEL_NORMAL:
            return [], f"variant outside the panel: {variant}"
    return found, ""


def translate_dpyd(raw: str) -> tuple[str | None, str]:
    if re.fullmatch(r"[*()\d]+/[*()\d]+", raw.replace(" ", "")):
        stars = raw.replace(" ", "").split("/")
        if all(s in DPYD_NORMAL_STARS for s in stars):
            return "*1/*1", ""
        return None, f"star allele outside the panel: {raw}"
    if "/" in raw:
        alleles: list[str] = []
        for side in raw.split("/"):
            found, problem = _dpyd_haplotype(side)
            if problem:
                return None, problem
            # HapB3 markers on one haplotype are one allele; the same
            # allele on both haplotypes is two.
            alleles.extend(found)
        return _diplotype(alleles), ""
    found, problem = _dpyd_haplotype(raw)
    if problem:
        return None, problem
    return _diplotype(found), ""


def translate_tpmt(raw: str) -> tuple[str | None, str]:
    alleles = raw.replace(" ", "").split("/")
    if len(alleles) != 2 or not all(a in TPMT_PANEL for a in alleles):
        return None, f"allele outside the panel: {raw}"
    return _diplotype(alleles), ""


def _ugt1a1_allele(text: str) -> str | None:
    text = text.replace(" ", "")
    if text in ("*60", "(*60)"):
        return "*1"
    m = re.fullmatch(r"\((\*\d+)\+\*60\)", text)
    if m:
        text = m.group(1)
    return text if text in UGT1A1_PANEL else None


def translate_ugt1a1(raw: str) -> tuple[str | None, str]:
    sides = raw.split("/")
    if len(sides) != 2:
        return None, f"unexpected format: {raw}"
    alleles = [_ugt1a1_allele(s) for s in sides]
    if any(a is None for a in alleles):
        return None, f"tentative or out-of-panel allele: {raw}"
    return _diplotype(alleles), ""


TRANSLATORS = {"DPYD": translate_dpyd, "TPMT": translate_tpmt, "UGT1A1": translate_ugt1a1}


def load_oncology_consensus(path: Path | str = TSV_PATH) -> list[ConsensusEntry]:
    entries: list[ConsensusEntry] = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for gene, translate in TRANSLATORS.items():
                raw = row[gene].strip()
                if not raw:
                    continue
                expected, note = translate(raw)
                entries.append(ConsensusEntry(row["coriell_id"], gene, raw, expected, note))
    return entries


def in_scope_samples(entries: list[ConsensusEntry] | None = None) -> list[dict]:
    """In-scope entries in the shape getrm_harness.run_concordance_check takes."""
    entries = entries if entries is not None else load_oncology_consensus()
    return [
        {"coriell_id": e.coriell_id, "gene": e.gene, "consensus_diplotype": e.consensus_diplotype}
        for e in entries
        if e.consensus_diplotype is not None
    ]
