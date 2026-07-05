# GeT-RM Concordance Report

**Status: harness implemented, live concordance not yet run.** This report
documents what exists, what was genuinely attempted, and exactly what's
needed to produce a real concordance table. No concordance percentage is
reported below because none has actually been measured -- see
`tests/getrm_harness.py` and `tests/test_getrm_validation.py`.

## What GeT-RM is

The Genetic Testing Reference Materials Coordination Program (GeT-RM), run
by the CDC in collaboration with the Coriell Cell Repositories, publishes
consensus pharmacogenetic genotypes for publicly available Coriell cell
line DNA samples, established by cross-laboratory concordance testing.
Two GeT-RM publications are relevant to this project's target genes:

* Pratt VM, et al. *Characterization of 137 Genomic DNA Reference Materials
  for 28 Pharmacogenetic Genes: A GeT-RM Collaborative Project.* Genet Med.
  2016. **PMID:26621017** -- covers CYP2C9, CYP2C19, DPYD, TPMT, UGT1A1
  (among others).
* Gaedigk A, et al. *Characterization of Reference Materials for Genetic
  Testing of CYP2D6 Alleles: A GeT-RM Collaborative Project.* J Mol Diagn.
  2019. **PMID:31401124** -- CYP2D6-specific, since CYP2D6's structural
  variation (deletions/duplications/hybrids) requires copy-number-aware
  characterization beyond the 2016 study's scope.

## What's in this repo right now

`src/pharminform/data/getrm_consensus_genotypes.json` has 11 real,
published CYP2D6 consensus diplotypes transcribed from the Gaedigk 2019
paper (PMID:31401124), including the well-known NA17244 sample
(consensus: `*2x2/*4x2+hybrid`). **The other five target genes'
consensus tables (published in Pratt 2016, PMID:26621017) have not yet
been transcribed into this repo** -- that's a to-do, not an omission of
data that doesn't exist.

## What was actually attempted (and what happened)

During development, this harness attempted a live connectivity check
against 1000 Genomes-hosted public sequencing data (the standard public
source for these Coriell samples' raw reads/variants, since many GeT-RM
Coriell IDs are also 1000 Genomes Project participants):

```
cyvcf2.VCF("https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr10...vcf.gz")
```

This failed with:

```
[E::hts_open_format] Failed to open file "https://...vcf.gz" : No such file or directory
OSError: Error opening https://...vcf.gz
```

This traces to the `cyvcf2` wheel installed in this environment bundling
an htslib build without remote-URL (libcurl) support -- a real,
reproducible environment limitation, not a stand-in for "we didn't try."
`tests/getrm_harness.py`'s `fetch_1000genomes_region` instead shells out to
`bcftools`, which supports remote indexed VCF access over HTTPS; `bcftools`
is also not installed in this environment (confirmed via
`shutil.which("bcftools")`), so every sample in the current run is
reported as `skipped_no_data` with that specific reason -- the harness
does not report a match/mismatch it hasn't actually verified.

The underlying 1000 Genomes URL itself **is** live and publicly reachable
(confirmed via `curl`, HTTP 200) -- the blocker is tooling in this
environment, not data access.

## What's needed to run this for real

1. Install `bcftools` (`brew install bcftools` or
   `conda install -c bioconda htslib bcftools`).
2. Confirm each Coriell sample ID in `getrm_consensus_genotypes.json` is
   actually present in the 1000 Genomes phase 3 sample panel (check against
   `igsr_samples.tsv` from https://www.internationalgenome.org) -- not
   every GeT-RM Coriell sample was sequenced by 1000 Genomes; some may
   need Coriell's own sequencing data or a targeted PGx panel instead.
3. Wire up the currently-stubbed step in `run_concordance_check`: fetch
   each gene's full region (not just the CYP2D6 CN-neutral demo window
   used for the connectivity smoke test), run `caller.call_diplotype_from_vcf`
   (or `caller.run_aldy` if a BAM/CRAM is fetched instead of a VCF -- strongly
   preferred for CYP2D6 given its structural variation), and compare the
   result to `consensus_diplotype`.
4. Transcribe the Pratt 2016 (PMID:26621017) consensus tables for
   CYP2C9, CYP2C19, DPYD, TPMT, and UGT1A1 into
   `getrm_consensus_genotypes.json` to extend coverage beyond CYP2D6.
5. Re-run `pytest tests/test_getrm_validation.py` and replace this report's
   "not yet run" status with a real per-gene concordance table generated
   from `getrm_harness.summarize()`.

## Concordance table (placeholder -- to be filled in once step 3 above is complete)

| Gene | Samples attempted | Matches | Mismatches | Concordance |
|---|---|---|---|---|
| CYP2D6 | -- | -- | -- | not yet run |
| CYP2C9 | -- | -- | -- | consensus data not yet transcribed |
| CYP2C19 | -- | -- | -- | consensus data not yet transcribed |
| DPYD | -- | -- | -- | consensus data not yet transcribed |
| TPMT | -- | -- | -- | consensus data not yet transcribed |
| UGT1A1 | -- | -- | -- | consensus data not yet transcribed |
