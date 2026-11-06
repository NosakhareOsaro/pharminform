# GeT-RM Concordance Report

**Status, 6 November 2026: the five genes PharmInform calls from a VCF
(DPYD, TPMT, UGT1A1, CYP2C9, CYP2C19) have been re-run on the 1000 Genomes
30x high-coverage panel, which covers more GeT-RM samples and carriers
than phase 3; the phenotype rules have been checked against CPIC's own
tables; the NA07048 DPYD discrepancy is explained; and CYP2D6 is no longer
called from VCF input. See "30x high-coverage run and CPIC table check"
at the end. The sections between here and there are the earlier phase 3
runs, kept as the record.**

Earlier phase 3 runs: all six genes, CYP2D6 (below), DPYD, TPMT and
UGT1A1, and CYP2C9 and CYP2C19 (sections further down). CYP2D6 result: 0/6 attempted samples matched their
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

## DPYD, TPMT and UGT1A1 (live run, 5 November 2026)

Code: branch `fix/cpic-phenotype-bands` (after commit 0378d4f). Data: 1000
Genomes phase 3 GRCh37 (`v5a`, AWS mirror), fetched with `bcftools` in
small windows around each defining variant. Consensus: CDC GeT-RM
consolidated table V4 (1-22-25), extracted to
`src/pharminform/data/getrm_consolidated_v4_oncology.tsv` and translated
to the expected panel call by `tests/getrm_oncology_consensus.py`. Raw
per-sample outcomes: `validation/getrm_oncology_live_results.json`.

| Gene | Consensus calls | Out of panel scope | Not in 1000 Genomes | Attempted | Diplotype match | Phenotype match |
|---|---|---|---|---|---|---|
| DPYD | 167 | 25 | 78 | 64 | 63/64 (98.4%) | 63/64 |
| TPMT | 147 | 13 | 74 | 60 | 60/60 (100.0%) | 60/60 |
| UGT1A1 | 197 | 5 | 151 | 41 | 37/41 (90.2%) | 41/41 |

"Out of panel scope" means the consensus contains an allele the bootstrap
table cannot detect (for example DPYD c.557A>G or exon deletions, TPMT
\*16 or \*40, tentative UGT1A1 calls); those samples are listed with the
reason by `load_oncology_consensus()` and are not counted either way.

### Every carrier and every mismatch

DPYD:

| Sample | Consensus (panel) | Called | Diplotype | Phenotype (consensus / called) |
|---|---|---|---|---|
| HG00118 | c.1236G>A(HapB3)/c.2846A>T | c.1236G>A(HapB3)/c.2846A>T | match | intermediate metabolizer / intermediate metabolizer |
| HG00129 | c.1236G>A(HapB3)/\*1 | c.1236G>A(HapB3)/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| HG00185 | \*2A/\*1 | \*2A/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| HG00332 | \*13/\*1 | \*13/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA07048 | \*1/\*1 | c.2846A>T/\*1 | mismatch | normal metabolizer / intermediate metabolizer |
| NA20901 | \*2A/\*1 | \*2A/\*1 | match | intermediate metabolizer / intermediate metabolizer |

TPMT:

| Sample | Consensus (panel) | Called | Diplotype | Phenotype (consensus / called) |
|---|---|---|---|---|
| HG00133 | \*2/\*1 | \*2/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| HG00589 | \*3C/\*1 | \*3C/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| HG01083 | \*2/\*1 | \*2/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA18966 | \*3C/\*1 | \*3C/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA19035 | \*3C/\*3C | \*3C/\*3C | match | poor metabolizer / poor metabolizer |
| NA19920 | \*3C/\*1 | \*3C/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA20296 | \*3C/\*1 | \*3C/\*1 | match | intermediate metabolizer / intermediate metabolizer |

UGT1A1:

| Sample | Consensus (panel) | Called | Diplotype | Phenotype (consensus / called) |
|---|---|---|---|---|
| HG00276 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA07056 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA07357 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA10847 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA11832 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA12003 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA12006 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA12156 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA12273 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA12717 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA12813 | \*28/\*28 | \*28/\*28 | match | poor metabolizer / poor metabolizer |
| NA12878 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA18519 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA18544 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA18868 | \*28/\*28 | \*28/\*28 | match | poor metabolizer / poor metabolizer |
| NA18942 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA18966 | \*6/\*1 | \*6/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA18973 | \*6/\*1 | \*6/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA18980 | \*6/\*6 | \*6/\*6 | match | poor metabolizer / poor metabolizer |
| NA19007 | \*6/\*1 | \*6/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA19095 | \*28/\*36 | \*28/\*1 | mismatch | intermediate metabolizer / intermediate metabolizer |
| NA19143 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA19147 | \*28/\*28 | \*28/\*28 | match | poor metabolizer / poor metabolizer |
| NA19207 | \*28/\*36 | \*28/\*1 | mismatch | intermediate metabolizer / intermediate metabolizer |
| NA19213 | \*36/\*1 | \*1/\*1 | mismatch | normal metabolizer / normal metabolizer |
| NA19239 | \*28/\*37 | \*28/\*28 | mismatch | poor metabolizer / poor metabolizer |
| NA19819 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |
| NA19908 | \*28/\*1 | \*28/\*1 | match | intermediate metabolizer / intermediate metabolizer |

### The mismatches

- **NA07048, DPYD.** GeT-RM consensus \*1/\*1; called c.2846A>T/\*1
  (intermediate). The call is what the data says: 1000 Genomes phase 3
  has NA07048 heterozygous (`0|1`) at rs67376798, and so does the 30x
  high-coverage release (`chr1:97082391 T>A`, `0|1`). The discrepancy is
  between the sequence data and the consensus, not in the caller.
  **Since explained** (6 November): the consensus for this sample comes
  from the 2016 study, whose assays did not test c.2846A>T -- see the 30x
  section below.
- **UGT1A1 \*36 and \*37 (NA19095, NA19207, NA19213, NA19239).** 1000
  Genomes phase 3 records the promoter TA repeat only as the TA7
  insertion (rs34983651, C>CAT), so TA5 (\*36) and TA8 (\*37) never
  reach the caller: \*36 reads as \*1 and \*37 as \*28. The phenotype is
  unaffected in all four (\*36 and \*1 are both normal function; \*37 and
  \*28 both decreased). A data limitation, not a caller fault; it would
  need a callset that represents every repeat length.

### Limits

- Real carriers are few: 5 DPYD, 7 TPMT and 28 UGT1A1
  consensus carriers among the attempted samples. Most matches are normal
  calls, which is why the synthetic tests in
  `tests/test_cpic_bands_and_alleles.py` still matter.
- Tested on real data: DPYD \*2A, \*13, c.2846A>T and HapB3 (including
  HG00118 with HapB3 and c.2846A>T together); TPMT \*2, \*3C
  heterozygous and homozygous; UGT1A1 \*6, \*28 heterozygous and
  homozygous. Not tested on real data: TPMT \*3A (no \*3A carrier in
  scope is in 1000 Genomes), so the cis/trans logic rests on the
  synthetic tests; UGT1A1 \*36 and \*37 (not represented in phase 3).
- Analytical concordance on public reference data, research use only --
  not clinical validation.

## CYP2C9 and CYP2C19 (live run, 5 November 2026)

Same method as the DPYD/TPMT/UGT1A1 run: CDC GeT-RM consolidated table V4,
CYP2C9 and CYP2C19 columns extracted unmodified to
`src/pharminform/data/getrm_consolidated_v4_cyp2c.tsv`, translated by
`tests/getrm_oncology_consensus.load_cyp2c_consensus`, run against 1000
Genomes phase 3. Raw outcomes: `validation/getrm_cyp2c_live_results.json`.

| Gene | Consensus calls | Out of panel scope | Not in 1000 Genomes | Attempted | Diplotype match | Phenotype match |
|---|---|---|---|---|---|---|
| CYP2C9 | 248 | 30 | 161 | 57 | 57/57 | 57/57 |
| CYP2C19 | 248 | 59 | 127 | 62 | 62/62 | 62/62 |

Carriers among the attempted samples -- CYP2C9: \*2/\*1 ×6, \*3/\*1 ×4;
CYP2C19: \*2/\*1 ×13, \*17/\*1 ×11, \*2/\*2 ×3, \*17/\*17 ×2, \*17/\*2 ×2, \*2/\*3 ×1, \*3/\*1 ×1. Out of scope: alleles outside the panel (for
example CYP2C9 \*5, \*8, \*9, \*11, \*13; CYP2C19 \*4, \*8, \*35) and every
call with a parenthetical alternative such as "\*1/\*1 (\*1/\*17)".

**A limit of the phenotype column.** The harness derives both phenotypes
with PharmInform's own model, so it shows that matching diplotypes give
matching phenotypes -- not that the model is right. Checking the model
directly against CPIC while loading these genes found three wrong rules,
fixed in the same change and pinned in
`tests/test_cpic_bands_and_alleles.py`: CYP2C9 \*1/\*2 was normal (CPIC:
intermediate) and \*2/\*3 intermediate (CPIC: poor); CYP2C19 \*2/\*17 and
\*3/\*17 were normal (CPIC: intermediate, since CYP2C19 is assigned by
allele function, not a summed score); CYP2D6 scores above 2.25 were
"rapid", a tier CYP2D6 does not have (CPIC: ultrarapid).

## 30x high-coverage run and CPIC table check (6 November 2026)

### Why a second data set

Phase 3 is low-coverage sequencing plus imputation, and has 2504 samples.
The 30x high-coverage panel (NYGC; Byrska-Bishop 2022, PMID 36055201) is
GRCh38, sequenced to 30x, and has 3202 samples -- phase 3 plus 698
relatives, several of them GeT-RM CEPH samples. It also calls the UGT1A1
TA repeat at every length (TA5, TA7, TA8), which phase 3 did not.

The 30x VCF has no rsIDs, so `tests/getrm_harness.GRCH38_DEFINING_SITES`
labels the defining sites by GRCh38 position and alleles, checked
against Ensembl GRCh38 and against the 30x records themselves. Run with
`run_concordance_check(samples, source="30x")`; one remote fetch per gene,
then per-sample calls locally. Raw outcomes:
`validation/getrm_30x_live_results.json`.

### Results

| Gene | Attempted (phase 3) | Attempted (30x) | Consensus carriers (30x) | Diplotype match | Phenotype match |
|---|---|---|---|---|---|
| DPYD | 64 | 94 | 7 | 94/94* | 94/94 |
| TPMT | 60 | 87 | 9 | 87/87 | 87/87 |
| UGT1A1 | 41 | 60 | 35 | 59/60 | 60/60 |
| CYP2C9 | 57 | 81 | 19 | 81/81 | 81/81 |
| CYP2C19 | 62 | 88 | 50 | 88/88 | 88/88 |

\*DPYD needs a caveat. 85 of the 94 DPYD samples have a star-name
consensus from the 2016 28-gene study (Pratt 2016, PMID 26621101), whose
assays did not include c.2846A>T or HapB3 (neither has a star name). For
those 85, only \*2A and \*13 are compared; c.2846A>T and HapB3 are left
out of the comparison (`ConsensusEntry.untested_alleles`). All four
variants are compared only for the 9 samples from the DPYD study (Gaedigk
2024, PMID 39032822), which include all 7 DPYD carriers: HG00118 (HapB3 +
c.2846A>T), HG00129 (HapB3), HG00185 and NA20901 (\*2A), HG00332 and
NA12248 (\*13), NA06991 (c.2846A>T). All 7 were called correctly.

**The one mismatch, NA19239, UGT1A1.** Consensus \*28/\*37, called
\*28/\*28 (the 30x data shows the TA7 and TA8 repeats directly). The
\*80 proxy marker is homozygous in this sample -- \*80 is in linkage
with both \*28 and \*37 -- and the proxy copy count overrides the
directly observed \*37. Phenotype is unaffected (both poor), but it is a
real caller fault, recorded as an open item below.

### NA07048 explained

Previously reported as unresolved. NA07048's DPYD consensus (\*1/\*1) is
a star-name call from the 2016 study, which did not test c.2846A>T, so
the consensus says nothing about that variant. The sequence evidence that
NA07048 carries it is now strong:

- Both 1000 Genomes call sets show it heterozygous at rs67376798.
- In the 30x panel's CEPH pedigree 1341, NA07048's mother NA07055 is also
  heterozygous and its father NA07034 is not -- a clean Mendelian
  transmission (NA07048 `0|1`, maternal allele).
- NA07055 has the same star-name \*1/\*1 consensus and the same call.
- NA06991, from the same pedigree but typed in the 2024 DPYD study, which
  did test c.2846A>T, has a consensus of c.2846A>T/\*1.

So there is no conflict to raise: the 2016 consensus did not cover the
variant. A note to GeT-RM could still be useful, suggesting NA07048 and
NA07055 as likely c.2846A>T carriers, but nothing here depends on it.

### Phenotype rules against CPIC's full tables

`tests/test_cpic_full_tables.py` checks every diplotype the panel can
call (82) against CPIC's diplotype-to-phenotype table (CPIC API, retrieved
6 November 2026, bundled as
`src/pharminform/data/cpic_diplotype_phenotype_panel.tsv`): phenotype for
all 82, and activity score where CPIC gives one. Phenotypes already
matched; three rules were corrected:

- CYP2D6 \*10 and \*41 have activity value 0.25 in CPIC, not the generic
  0.5 for decreased function (11 activity scores were wrong, e.g. \*1/\*10
  1.5 instead of 1.25; the phenotype happened to be the same in all 11).
- UGT1A1 \*36 is increased function in CPIC (was normal). UGT1A1 is now
  assigned by counting decreased alleles, as CPIC does, so \*36 cannot
  offset \*28.
- UGT1A1 \*80 seen without the TA repeat itself is "Indeterminate" in
  CPIC. A \*28 inferred only from the \*80 proxy is now reported as
  indeterminate, with no dosing guidance applied.

### CYP2D6 no longer called from VCF

0 of 6 matched on phase 3, and the causes are structural (\*5 deletion,
duplications, \*36 and other CYP2D6/2D7 hybrids) or in the allele table
(\*10's rs1065852 also lies on \*4), not fixable by adding alleles. The
pipeline now reports CYP2D6 as "not called" from VCF input, with the
reason, and applies no CYP2D6 phenotype. `caller.run_aldy` remains for
aligned reads.

Removing CYP2D6 exposed a guideline-lookup fault: for a two-gene drug the
lookup took the first recommendation matching either gene, so
amitriptyline (CYP2D6 + CYP2C19) silently assumed a CYP2D6 intermediate
metabolizer, and azathioprine and mercaptopurine assume a NUDT15
phenotype though NUDT15 is never called. The lookup now prefers the
recommendation that assumes normal metabolism for an untested gene, and
the report states it ("Assumes, not tested: CYP2D6 Normal Metabolizer").

## What's still needed

1. UGT1A1: a directly observed TA repeat should take precedence over the
   \*80 proxy (NA19239 above).
2. More carriers, especially DPYD (7) and TPMT (9). Public data is close
   to exhausted for these; the next source is a lab pilot with real
   clinical samples.
3. CYP2D6 with Aldy on real aligned reads (BAM/CRAM), if CYP2D6 is to be
   offered again; a materially larger effort (reference genome, ILP
   solver, multi-GB alignments).
4. NUDT15 is not called, so thiopurine guidance assumes a normal NUDT15;
   adding it would remove that assumption.
