# PharmInform

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21210588.svg)](https://doi.org/10.5281/zenodo.21210588)

Clinical pharmacogenomics (PGx) platform: star allele calling, CPIC/DPWG
dosing guidance, and an experimental multi-gene composite risk score for
drugs affected by more than one pharmacogene -- a gap in existing
open-source PGx tools (see [Composite Risk Scoring](#composite-risk-scoring-experimental)
below).

**Status:** research-stage, not a clinical decision support tool. See
[Disclaimer](#disclaimer).

## Progress checklist

- [x] Repo scaffold, packaging, license
- [x] VCF ingestion + Aldy star allele calling wrapper
- [x] PharmVar + PharmGKB API clients
- [x] CPIC/DPWG guideline lookup + cache
- [x] OMOP CDM drug_exposure mapping
- [x] Diplotype call confidence flagging
- [x] Curated drug-pathway weight reference data
- [x] Composite multi-gene risk scoring engine
- [x] Directional ultrarapid-metabolizer flag
- [x] Composite risk test suite (monotonicity, fallback)
- [x] Face-validity case-report tests
- [x] GeT-RM validation: live concordance runs completed -- DPYD 63/64, TPMT 60/60, UGT1A1 37/41 diplotypes (41/41 phenotypes), CYP2C9 57/57, CYP2C19 62/62; CYP2D6 0/6 (see below)
- [x] Jinja2 HTML clinical report + PDF export
- [x] CPIC truth-set pytest suite
- [x] CLI
- [x] FastAPI endpoint
- [x] Streamlit web demo
- [x] NetworkX drug-gene visualization
- [x] Dockerfile
- [x] GitHub Actions CI
- [x] Bioconda recipe skeleton
- [x] Sphinx docs skeleton
- [x] Full README (this file)
- [x] JOSS paper draft

## Architecture

```
                 ┌─────────────────────┐
   VCF (or BAM) →│ caller.py            │→ diplotype per gene
                 │  - cyvcf2 (VCF path) │   (rule-based VCF matcher, or
                 │  - Aldy (BAM path)   │    real Aldy CLI for BAM/CRAM)
                 └──────────┬───────────┘
                            ↓
                 ┌─────────────────────┐
                 │ confidence.py        │→ confidence: high/low + reasons
                 └──────────┬───────────┘
                            ↓
                 ┌─────────────────────┐
                 │ pharmgkb_client.py   │→ phenotype (CPIC activity-score model)
                 └──────────┬───────────┘
                            ↓
        ┌───────────────────┴────────────────────┐
        ↓                                         ↓
┌───────────────────┐                  ┌─────────────────────────┐
│ cpic_dpwg.py       │                  │ composite_risk.py        │
│ per-gene CPIC       │                  │ multi-gene composite     │
│ dosing recs         │                  │ risk score (experimental)│
└──────────┬─────────┘                  └────────────┬─────────────┘
           └───────────────────┬───────────────────────┘
                                ↓
                     ┌────────────────────┐
                     │ report.py           │→ HTML / PDF clinical report
                     │ pipeline.py          │
                     └──────────┬───────────┘
                                ↓
              ┌─────────────────┼─────────────────┐
              ↓                 ↓                 ↓
          cli.py             api.py         streamlit_app/app.py
     (Typer CLI)        (FastAPI JSON)      (interactive web demo)
```

`omop_mapping.py` and `network.py` sit alongside this pipeline: the former
maps report output onto OMOP CDM `drug_exposure`-shaped records, the
latter builds a NetworkX gene↔drug graph from the same reference data used
by `composite_risk.py`.

## Install

```bash
git clone https://github.com/NosakhareOsaro/pharminform.git
cd pharminform
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

PDF export (`weasyprint`) needs native Pango/cairo/GObject libraries not
installed by `pip`:

- **macOS (Homebrew):** `brew install pango`, then run with
  `DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib` set (Apple Silicon) --
  confirmed necessary on this project's own development machine; a plain
  `brew install pango` was not sufficient on its own.
- **Linux (Debian/Ubuntu):** `apt-get install libpango-1.0-0
libpangoft2-1.0-0 libgdk-pixbuf-2.0-0 libcairo2` (see `Dockerfile` and
  `.github/workflows/ci.yml` for exact package names, including the
  `libgdk-pixbuf2.0-0` → `libgdk-pixbuf-2.0-0` rename on Debian trixie that
  broke this project's own first Docker build attempt).

Aldy (real BAM/CRAM-based, copy-number-aware star allele calling) is an
optional extra: `pip install -e ".[aldy]"`. It also needs an ILP solver
(Gurobi or Google OR-Tools/CBC) -- see [Aldy's docs](https://github.com/0xTCG/aldy).

## Usage

### CLI

```bash
pharminform report --vcf sample.vcf --out report.html --pdf report.pdf --min-depth 20
```

### FastAPI

```bash
uvicorn pharminform.api:app --reload
curl -X POST http://localhost:8000/report -F "vcf=@sample.vcf" -F "min_depth=20"
```

### Streamlit demo

```bash
streamlit run streamlit_app/app.py
```

### Docker

```bash
docker build -t pharminform .
docker run -p 8000:8000 pharminform          # FastAPI
docker run -p 8501:8501 pharminform streamlit run streamlit_app/app.py --server.address 0.0.0.0
```

## CPIC gene coverage

| Gene    | Star allele source (this project)                             | Structural variation?                                                                   |
| ------- | ------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| CYP2D6  | Aldy (BAM/CRAM) preferred; simplified rule-based VCF fallback | Yes -- deletions/duplications/hybrids; VCF-only calls are always flagged low-confidence |
| CYP2C19 | Aldy or rule-based VCF                                        | No                                                                                      |
| CYP2C9  | Aldy or rule-based VCF                                        | No                                                                                      |
| DPYD    | Aldy or rule-based VCF                                        | No                                                                                      |
| TPMT    | Aldy or rule-based VCF                                        | No                                                                                      |
| UGT1A1  | Aldy or rule-based VCF                                        | Yes (\*28 TA-repeat expansion; not representable as a simple SNV match)                 |

The rule-based VCF matcher (`caller.call_diplotype_from_vcf`) uses a
small, illustrative, non-exhaustive allele table
(`caller.BOOTSTRAP_ALLELE_DEFINITIONS`) -- not a substitute for PharmVar's
authoritative definitions (`pharmvar_client.py`, which requires an
account-issued API key) or for running Aldy directly against aligned reads.

## GeT-RM validation

**DPYD, TPMT, UGT1A1, CYP2C9 and CYP2C19.** Consensus genotypes from the CDC GeT-RM
consolidated table (V4, 1-22-25) were translated to the calls the
bootstrap allele table can make, and the samples that are also in 1000
Genomes phase 3 were run through the real pipeline:

| Gene | Attempted | Diplotype match | Phenotype match |
|---|---|---|---|
| DPYD | 64 | 63 | 63 |
| TPMT | 60 | 60 | 60 |
| UGT1A1 | 41 | 37 | 41 |
| CYP2C9 | 57 | 57 | 57 |
| CYP2C19 | 62 | 62 | 62 |

The one DPYD mismatch (NA07048) is a sample both 1000 Genomes call sets
show as heterozygous for c.2846A>T while the GeT-RM consensus says
\*1/\*1 -- the caller reports what the data contains, and the discrepancy
is unresolved. The four UGT1A1 mismatches are \*36/\*37 repeat lengths
that phase 3 does not represent; their phenotypes still match. Carriers
are few (5 DPYD, 7 TPMT, 28 UGT1A1, 10 CYP2C9, 33 CYP2C19), samples with alleles outside the
panel are excluded and counted, and this is analytical testing on public
reference data for research use, not clinical validation. The phenotype
column uses PharmInform's own model on both sides, so it does not test the
model itself; that is covered by the CPIC cases in
`tests/test_cpic_bands_and_alleles.py`. Details in
[`validation/getrm_concordance_report.md`](validation/getrm_concordance_report.md).

**CYP2D6.** A live concordance run against published GeT-RM consensus genotypes has
been completed for CYP2D6, using real 1000 Genomes sequencing data (11
consensus diplotypes transcribed from Gaedigk et al. 2019, PMID:31401124).
Of those 11 samples, 6 are present in the 1000 Genomes phase 3 panel
(verified against both the official IGSR panel file and the actual VCF
header) and were run through the real pipeline
(`bcftools` fetch -> `cyvcf2` parse -> `caller.call_diplotype_from_vcf`);
the other 5 were honestly skipped as not present in that panel.

**Result: 0/6 matched (0.0% concordance).** This is an expected,
attributable result, not a bug: every mismatch traces either to an
allele absent from the small illustrative `BOOTSTRAP_ALLELE_DEFINITIONS`
table (5 of CYP2D6's 150+ named haplotypes) or to copy-number/gene-conversion
alleles (`*36`, `x2` duplication, `+hybrid`) that VCF-only matching
cannot represent -- exactly the limitation this project has documented
throughout as requiring Aldy against aligned reads instead. Running this
for real also surfaced and fixed a genuine bug in `caller.py` (homozygous-
reference genotypes were being misread as evidence of the variant allele).
Full per-sample breakdown, panel-coverage table, and next steps (the
other five genes' consensus data, and a same-style Aldy/BAM comparison)
in [`validation/getrm_concordance_report.md`](validation/getrm_concordance_report.md).

## Composite Risk Scoring (Experimental)

Most PGx tools (Aldy, PharmCAT, StellarPGx) call star alleles and
phenotypes **per gene**, and CPIC/DPWG guidelines give **per-gene** dosing
recommendations. For drugs metabolized by more than one pharmacogene --
amitriptyline via both CYP2D6 and CYP2C19 is this project's running
example -- nothing in the existing open-source ecosystem combines those
per-gene signals into one risk indicator. `composite_risk.py` does that:

1. Each gene's phenotype is converted to a normalized activity score
   (0 = poor metabolizer/no function, 1 = normal, capped at 1 -- this
   models _reduced-function toxicity accumulation_ only).
2. Each contributing gene is weighted by its published fraction of the
   drug's metabolism (`fraction_metabolized` in
   `data/drug_pathway_weights.json`, cited in `DATA_SOURCES.md`) -- or, if
   no published fraction exists for at least one gene, by equal weighting
   across the genes present, **and this fallback is always reported**
   in the `weighting_method` field, never silently assumed.
3. The weighted scores combine into one `composite_function_score`
   (0-1, lower = more reduced capacity = higher toxicity risk) and a
   `risk_tier` (Low/Moderate/High), calibrated so a single poor-metabolizer
   gene at ≥50% effective weight can never resolve to "Low" -- validated by
   a parametrized monotonicity test across every drug in the reference
   dataset (`tests/test_composite_risk.py`).
4. A separate `directional_flag` fires independently whenever any
   contributing gene is an ultrarapid metabolizer (under-exposure or
   dangerous prodrug over-activation risk, e.g. codeine → morphine via
   CYP2D6 UM) -- deliberately **not blended** into the composite score,
   since it's a different physiological mechanism.

**This is exploratory, not a validated clinical scoring system.** Every
result carries a mandatory disclaimer field. Three drugs in the reference
dataset (codeine, tramadol, clopidogrel) are prodrug-activation pathways
where reduced gene function means _under-treatment_, not toxicity
accumulation -- the opposite of what the composite score's default
framing assumes; this is flagged per-drug (`directionality` field) and
discussed in full in `DATA_SOURCES.md` and `paper.md`'s limitations
section.

> **Disclaimer** (repeated from every generated report): Exploratory
> composite score. Not a CPIC/DPWG guideline output. For
> research/educational use only -- clinical decisions should rely on
> individual gene-level CPIC/DPWG recommendations.

## Testing

```bash
pip install -e ".[dev]"
pytest                       # full suite
ruff check src tests         # lint
DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib pytest tests/test_report.py  # macOS, to exercise real PDF export
```

101 tests as of this writing, covering VCF ingestion/filtering, the
CPIC-truth-set star allele matcher, confidence flagging, the CPIC/PharmGKB
API clients (mocked), composite risk (monotonicity, fallback, face
validity against 3 real published case reports), the CLI, FastAPI
endpoint, Streamlit app (via `AppTest`), and the NetworkX graph.

## Disclaimer

PharmInform is a research
project, not a certified clinical
decision support tool. Star allele calls, phenotypes, CPIC
recommendations, and especially the composite risk score should not be
used to make real prescribing decisions. Verify all outputs against
primary CPIC/DPWG/PharmVar sources.

---

Built PharmInform: clinical pharmacogenomics platform performing
CYP2D6/CYP2C19/DPYD star allele calling (Aldy + PharmVar), CPIC/DPWG
dosing recommendations, an original multi-gene composite risk scoring
model, and Streamlit web demo -- published to PyPI with Bioconda recipe
and submitted to JOSS. (GeT-RM concordance runs completed for DPYD,
TPMT, UGT1A1 and CYP2D6 -- see GeT-RM validation section above.)
