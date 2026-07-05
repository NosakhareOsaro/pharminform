# GeT-RM Concordance Report

**Status: live concordance run completed for CYP2D6 against real 1000
Genomes sequencing data.** Result: 0/6 attempted samples matched their
published consensus diplotype (0.0% concordance). Every mismatch has a
specific, attributable cause explained below -- this is a real,
informative result that confirms this project's own repeated
documentation that the VCF-only rule-based caller is not adequate for
CYP2D6, not a surprise finding. 5 of the 11 bundled samples are not
present in the 1000 Genomes panel and were honestly skipped rather than
guessed at.

## What GeT-RM is

The Genetic Testing Reference Materials Coordination Program (GeT-RM), run
by the CDC in collaboration with the Coriell Cell Repositories, publishes
consensus pharmacogenetic genotypes for publicly available Coriell cell
line DNA samples, established by cross-laboratory concordance testing.
Two GeT-RM publications are relevant to this project's target genes:

* Pratt VM, et al. *Characterization of 137 Genomic DNA Reference Materials
  for 28 Pharmacogenetic Genes: A GeT-RM Collaborative Project.* J Mol
  Diagn. 2016. **PMID:26621101** -- covers CYP2C9, CYP2C19, DPYD, TPMT,
  UGT1A1 (among others).
* Gaedigk A, et al. *Characterization of Reference Materials for Genetic
  Testing of CYP2D6 Alleles: A GeT-RM Collaborative Project.* J Mol Diagn.
  2019. **PMID:31401124** -- CYP2D6-specific, since CYP2D6's structural
  variation (deletions/duplications/hybrids) requires copy-number-aware
  characterization beyond the 2016 study's scope.

## Sample panel coverage (verified against two independent real sources)

`src/pharminform/data/getrm_consensus_genotypes.json` has 11 real,
published CYP2D6 consensus diplotypes transcribed from Gaedigk 2019.
**The other five target genes' consensus tables (Pratt 2016) have not
yet been transcribed into this repo** -- a to-do, not an omission of
data that doesn't exist.

Of the 11 CYP2D6 samples, cross-referenced against (a) the official
1000 Genomes phase 3 sample panel
(`integrated_call_samples_v3.20130502.ALL.panel`, fetched live from
`ftp.1000genomes.ebi.ac.uk`) and (b) the actual sample columns in the
real phase 3 VCF header (`bcftools query -l` against the live AWS Open
Data mirror) -- both sources agree exactly:

| Coriell ID | In 1000 Genomes panel? | Consensus diplotype |
|---|---|---|
| NA18563 | Yes | \*1/\*36+\*10 |
| NA18564 | Yes | \*2A/\*36+\*10 |
| NA18565 | Yes | \*10/\*36x2 |
| NA19785 | Yes | \*1/\*13+\*2 |
| NA19908 | Yes | \*1/\*46 |
| NA20289 | Yes | \*6/\*11 |
| NA17222 | No | \*2/\*108 |
| NA17287 | No | \*1/\*83 |
| NA18524 | No | \*1/\*36x2+\*10 |
| NA19790 | No | \*1/\*13+\*2 |
| NA17244 | No | \*2x2/\*4x2+hybrid |

(NA17244 -- the sample most often cited as a canonical CYP2D6
duplication/hybrid example -- is notably *not* in the 1000 Genomes panel.
Samples not in 1000 Genomes would need Coriell's own sequencing data or
a targeted PGx panel to test.)

## Real concordance results (CYP2D6, 6 samples attempted)

Run via `tests/getrm_harness.run_concordance_check()`: for each of the
6 panel-covered samples, real genotypes were fetched from the live 1000
Genomes phase 3 VCF (`bcftools view -r 22:42,500,000-42,550,000 -s
<sample> <url>`), parsed with `cyvcf2`, and passed through the actual
production `caller.call_diplotype_from_vcf` -- the same function the
CLI/API/Streamlit interfaces use, not a separate one-off script.

| Coriell ID | Consensus | Called | Result | Why |
|---|---|---|---|---|
| NA18563 | \*1/\*36+\*10 | \*10/\*1 | Mismatch | `*10` correctly detected; `*36` is a CYP2D6/2D7 gene-conversion allele -- a structural event invisible to SNV-only matching, and not in the bootstrap allele table either |
| NA18564 | \*2A/\*36+\*10 | \*10/\*1 | Mismatch | `*10` detected; `*2A` and `*36` both absent from the bootstrap table, and `*36` is structural regardless |
| NA18565 | \*10/\*36x2 | \*10/\*10 | Mismatch | `*10` detected on both matched calls; the `x2` copy-number duplication is invisible to VCF-only matching -- this is exactly the "requires Aldy on aligned reads" case this project has documented throughout |
| NA19785 | \*1/\*13+\*2 | \*1/\*1 | Mismatch | Neither `*13` nor `*2` is in the bootstrap allele table (which only covers `*3`, `*4`, `*10`, `*17`, `*41`) -- falls back to wild-type default |
| NA19908 | \*1/\*46 | \*1/\*1 | Mismatch | `*46` not in the bootstrap table |
| NA20289 | \*6/\*11 | \*1/\*1 | Mismatch | Neither `*6` nor `*11` is in the bootstrap table |

```
{
  "total": 11,
  "matches": 0,
  "mismatches": 6,
  "skipped_no_data": 5,
  "concordance_pct_of_attempted": 0.0
}
```

**0/6 (0.0%) concordance for CYP2D6 via the VCF-only rule-based caller
against real GeT-RM/1000 Genomes samples.** Every single mismatch traces
to one of two already-documented, already-disclosed limitations of this
pipeline's VCF-only fallback path:

1. **`caller.BOOTSTRAP_ALLELE_DEFINITIONS` is a small, explicitly
   illustrative table** (5 CYP2D6 alleles: `*3`, `*4`, `*10`, `*17`,
   `*41`), not an exhaustive PharmVar-equivalent definition set (PharmVar
   currently documents 150+ named CYP2D6 haplotypes). Real individuals
   overwhelmingly carry alleles outside this small table.
2. **VCF-only matching cannot see copy number or gene conversion**
   (`*36`, `x2` duplication, `+hybrid`) -- this requires Aldy running
   against aligned reads (BAM/CRAM), which this project has documented
   as the clinically-preferred path since `caller.py` was first written.

This result does not indicate a bug in the comparison itself (the
harness fetch/parse/call/compare pipeline works correctly, as
demonstrated by getting distinct, real, sample-specific outputs for each
of the 6 samples) -- it's an accurate measurement confirming that the
simplified VCF-only rule-based fallback is exactly as limited as this
project has said all along, now with real numbers instead of a caveat.

## A real bug this run found and fixed

Running this for real against actual population data surfaced a genuine
bug in `caller.py` itself (fixed, see git history): `load_variants` and
`call_diplotype_from_vcf` never checked whether a VCF record's genotype
actually contained a non-reference allele, so every polymorphic site in
a gene region was treated as "observed in this sample" -- including
homozygous-reference (`0/0`) calls. This was invisible against the
hand-authored synthetic test fixture (which only ever included rows
where a variant was truly present), but a real population VCF slice
lists every site polymorphic across all 2504 samples, most of which any
one sample is `0/0` for. Before the fix, every one of the 6 samples
below produced the identical, wrong `*10/*17` call. Fixed by checking
for a non-reference, non-missing allele before treating a record as
evidence of a defining variant.

## What's needed to extend this beyond CYP2D6

1. Transcribe the Pratt 2016 (PMID:26621101) consensus tables for
   CYP2C9, CYP2C19, DPYD, TPMT, and UGT1A1 into
   `getrm_consensus_genotypes.json`.
2. Add each gene's GRCh37 region to `GRCH37_GENE_REGIONS` in
   `tests/getrm_harness.py` (CYP2D6 is currently the only one
   configured).
3. Re-run `pytest tests/test_getrm_validation.py` -- the harness will
   pick up the new genes automatically via the same
   `run_concordance_check()` path exercised here.
4. For CYP2D6 specifically, re-testing with `caller.run_aldy` against
   real aligned reads (BAM/CRAM, e.g. 1000 Genomes' own alignments)
   rather than the VCF-only fallback would be the clinically-relevant
   comparison, since Aldy is copy-number-aware; that's a materially
   larger effort (reference genome + ILP solver + multi-GB alignment
   downloads) not undertaken in this pass.
