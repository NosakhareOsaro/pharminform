# PharmInform

Clinical pharmacogenomics (PGx) platform: star allele calling, CPIC/DPWG dosing
guidance, and an experimental multi-gene composite risk score for drugs
affected by more than one pharmacogene.

**Status: under active development.** See the checklist below for progress.

## Progress checklist

- [x] Repo scaffold, packaging, license
- [ ] VCF ingestion + Aldy star allele calling wrapper
- [ ] PharmVar + PharmGKB API clients
- [ ] CPIC/DPWG guideline lookup + cache
- [ ] OMOP CDM drug_exposure mapping
- [ ] Diplotype call confidence flagging
- [ ] Curated drug-pathway weight reference data
- [ ] Composite multi-gene risk scoring engine
- [ ] Directional ultrarapid-metabolizer flag
- [ ] Composite risk test suite (monotonicity, fallback)
- [ ] Face-validity case-report tests
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
