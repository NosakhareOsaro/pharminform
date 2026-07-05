# Data Sources

This file documents where every piece of external reference data used by
PharmInform comes from, and is the audit trail for the composite risk
model's `fraction_metabolized` values in particular (see
[Phase 2 methodology](README.md#composite-risk-scoring-experimental)).

All PMIDs below were verified against a real PubMed/PMC record at the time
they were added (via live web search/fetch, not recalled from memory
without verification). If you spot one that looks wrong, please open an
issue -- citation accuracy matters more here than anywhere else in the repo.

## Star allele definitions

* **PharmVar** (https://www.pharmvar.org) is the intended authoritative
  source (`pharmvar_client.py`). Its REST API now requires an
  account-issued API key (confirmed live: `GET /api-service/genes` returns
  HTTP 401 without one). Set `PHARMVAR_API_KEY` to use it.
* `caller.BOOTSTRAP_ALLELE_DEFINITIONS` is a small, illustrative,
  non-exhaustive set of core defining variants used for offline testing and
  the VCF-only fallback caller. rsIDs are widely cited in the PGx
  literature but have not been individually re-verified against PharmVar
  for this project -- treat as a bootstrap fixture, not a clinical source.

## Diplotype -> phenotype

* `pharmgkb_client.py` queries PharmGKB's public REST API
  (https://api.pharmgkb.org/v1/data, no key required, confirmed live) for
  gene/chemical metadata.
* Diplotype-to-phenotype conversion uses CPIC's activity-score model
  (allele function -> numeric value -> summed score -> phenotype band).
  The generic activity-score band thresholds in `DEFAULT_ACTIVITY_BANDS`
  are an approximation across all six target genes, not independently
  verified per-gene against each CPIC gene-specific standardization
  publication (e.g. CYP2D6's exact bands may differ slightly from
  CYP2C19's) -- documented as a simplification in `pharmgkb_client.py`.

## CPIC dosing recommendations

* `cpic_dpwg.py` and `data/cpic_recommendation_cache.json` are fetched
  directly from CPIC's live public REST API (https://api.cpicpgx.org/v1,
  a PostgREST interface over https://github.com/cpicpgx/cpic-data), fetch
  date recorded in the cache file's `_fetched` field. This is real,
  unmodified guideline content (recommendation text, classification,
  citations), deduplicated by phenotype combination.
* DPWG guidelines have no public structured API and are **not** integrated
  in this release.
* Warfarin has zero cached recommendation rows: CPIC's own API returns none,
  because its guideline is a multi-variable dosing algorithm rather than a
  phenotype-lookup table. This is a genuine finding, not a bug.

## Composite risk scoring: `fraction_metabolized` citation table

| Drug | Gene | Fraction | Directionality | Citation | Note |
|---|---|---|---|---|---|
| amitriptyline | CYP2D6 | null | clearance | PMID:23486447, PMID:27997040 | CPIC guideline is qualitative only; equal-weight fallback |
| amitriptyline | CYP2C19 | null | clearance | PMID:23486447, PMID:27997040 | same |
| codeine | CYP2D6 | 0.10 | **prodrug activation** | PMID:22205192 | ~10% O-demethylated to morphine; reduced function = under-treatment, not toxicity |
| tramadol | CYP2D6 | null | **prodrug activation** | PMID:33387367 | CYP2D6 is primary route to active M1 metabolite; no consistent quantitative fraction found |
| clopidogrel | CYP2C19 | null | **prodrug activation** | PMID:35034351 | ~85% inactivated by CES1 (outside gene panel); CYP2C19 dominant among CYP oxidative steps |
| warfarin | CYP2C9 | 0.80 | clearance | PMID:27539372, PMID:21900891, PMID:28198005 | S-warfarin only; VKORC1/CYP4F2 outside gene panel; CPIC uses a dosing algorithm, not a lookup table |
| voriconazole | CYP2C19 | 0.631 | clearance | PMID:35335853, PMID:27981572 | in vitro CYP-inhibition study, hepatic N-oxidation |
| fluorouracil | DPYD | 0.80 | clearance | PMID:12164772, PMID:23988873 | DPD catabolizes >80% of administered 5-FU |
| capecitabine | DPYD | 0.80 | clearance | PMID:12164772, PMID:23988873 | fraction carried over from 5-FU (its active metabolite), not capecitabine-specific |
| irinotecan | UGT1A1 | null | clearance | PMID:26417955, PMID:9789606 | primary SN-38 deactivating enzyme; no consistent quantitative fraction found |
| azathioprine | TPMT | null | clearance | PMID:21270794, PMID:23422873, PMID:30447069 | one of three competing pathways (with XO, HPRT); NUDT15 co-guideline gene out of scope |
| mercaptopurine | TPMT | null | clearance | PMID:21270794, PMID:23422873, PMID:30447069 | shares azathioprine's active-metabolite pharmacology |

**"Directionality: prodrug activation" is a real limitation of the
composite scoring model, not a data-quality issue**: `composite_risk.py`'s
"lower composite function = higher risk" framing assumes reduced enzyme
function causes toxic accumulation of an active parent drug. For codeine,
tramadol, and clopidogrel, the gene instead activates a prodrug, so reduced
function causes *under-treatment* (or, for opioids, is separately captured
by the ultrarapid-metabolizer `directional_flag()`), not toxicity
accumulation. These three drugs are included in the reference dataset for
completeness and transparency, but their composite score should be read as
"pathway activation capacity," not "toxicity risk," until the model is
extended to represent both directions explicitly. This is called out again
in `paper.md`'s limitations section.

## RxNorm identifiers

All `rxnorm_id` values come directly from CPIC's own `drug` table
(`GET https://api.cpicpgx.org/v1/drug?name=eq.<drug>`), not looked up
independently.
