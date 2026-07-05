# PharmInform

Clinical pharmacogenomics (PGx) platform: star allele calling, CPIC/DPWG dosing
guidance, and an experimental multi-gene composite risk score for drugs
affected by more than one pharmacogene.

**Status: under active development.** See the checklist below for progress.

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
- [ ] GeT-RM validation harness
- [ ] Jinja2 HTML clinical report + PDF export
- [ ] CPIC truth-set pytest suite
- [ ] CLI
- [ ] FastAPI endpoint
- [ ] Streamlit web demo
- [ ] NetworkX drug-gene visualization
- [ ] Dockerfile
- [ ] GitHub Actions CI
- [ ] Bioconda recipe skeleton
- [ ] Sphinx docs skeleton
- [ ] Full README (architecture, usage, GeT-RM results, composite risk methodology)
- [ ] JOSS paper draft

Full documentation will land in this README as each piece is built.
