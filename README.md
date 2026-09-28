# rnd-longmer-library-qc-nf

Nanopore QC for **direct-print long-oligo (750mer) paired VH+VL libraries** — the pathway
where each variant is printed as one full-length oligo rather than assembled from a VH half
and a VL half by overlap-extension PCR.

Two inputs: the **design reference CSV** and the run's **fastq_pass prefix(es)**. Everything
else is derived, so Ops runs the same command for any library on this construct architecture.

```bash
nextflow run Twistbioscience/rnd-longmer-library-qc-nf \
  --reference s3://.../AILK014-007-008-009_750mer_reference.csv \
  --input     samplesheet.csv \
  --outdir    s3://.../AILK014_007-009_QC \
  --library_name AILK014-007-008-009 \
  --account Absci --order_number Q-724749 --manufacture_date 2026-09-25
```

`samplesheet.csv`:

```csv
sublibrary,fastq_dir
AILK014-007,s3://twist-instrument-data/ngs/PromethION/<run-007>/fastq_pass
AILK014-008,s3://twist-instrument-data/ngs/PromethION/<run-008>/fastq_pass
AILK014-009,s3://twist-instrument-data/ngs/PromethION/<run-009>/fastq_pass
```

`sublibrary` must match the `order` column of the reference. See [docs/USAGE.md](docs/USAGE.md)
for the Ops runbook and [docs/METHODS.md](docs/METHODS.md) for what each number means.

---

## What it measures

| | |
|---|---|
| **VH/VL pairing** | Each read is aligned once against a combined index of designed VH-side and VL-side half-references; a half is assigned only if its best hit beats the runner-up by ≥15 matching bases over ≥120 nt. Reads are then classed `clean_designed` / `within_well` / `cross_well` / `cross_sublib` / `one_half` / `unmapped`. |
| **Internal block deletions** | Reads that still carry **both** primer sites are aligned with minimap2 splice preset against their **own** assigned construct; a single internal gap ≥50 nt is a call, cross-checked with a splice-model-free edlib alignment. |
| **Uniformity and dropout** | Per-variant and per-well read counts from correctly paired reads: dropout, CV, Gini, 95/5, within-2×-median. |
| **Full-length** | Four definitions against three denominators, because they disagree by >10 points on the same reads. |
| **Screening depth** | Poisson recovery per sub-library from each variant's own abundance and its own deletion rate. |
| **Calibration** | Classifier and deletion-caller controls are simulated from **this library's own constructs** on every run — a reference with more well-mates is a harder problem, so the operating point is re-measured rather than inherited. |

### Why content deconvolution, not soft-clips

In a conserved-framework antibody library a chimera aligns end-to-end to one designed gene
with no soft clip and no supplementary alignment, so an aligner-rescue caller scores it
clean. The undercount is library-specific and has been measured between ~9× and ~80× on
other Twist builds. This pipeline never uses a soft-clip heuristic for pairing.

---

## Outputs

```
<outdir>/
  reports/    <name>_qc_summary.md      headline table + denominators + calibration — read this first
              <name>_nanopore_QC.xlsx   summary, controls, read length, block deletion,
                                        full-length ladder, screening depth, per-variant, per-well
              <name>_QC_Summary.docx    one-page customer collateral (MGF format)
  plots/      fig01..fig17 (300 dpi PNG + SVG), figures_skipped.json
  metrics/    metrics_core.json, metrics_extra.json, blockdel_summary.json,
              blockdel_scope.json, controls.json, derived.json
  reference/  halves.fasta.gz, constructs.fasta.gz, gene_meta.json, reference_qc.json
  per_gene.csv, per_well.csv, blockdel_per_gene_all.csv, blockdel_events.csv
  pipeline_info/  timeline, report, trace, dag
```

**The narrative customer report is deliberately not automated.** The pipeline produces the
numbers, figures, workbook and the numeric one-page collateral; the interpretive report —
what the numbers mean for this build, what is and is not explained — is written by a
scientist from these outputs. A pipeline that also wrote the conclusions would be asserting
causes it cannot test.

### Reading the numbers

- Pairing percentages are shares of reads with **both halves assigned**, not of all reads.
- Deletion rates are shares of reads carrying **both primer sites** and assigned to a
  variant, on a 1-in-N subsample (default 8).
- Intact fractions are **upper bounds** (the caller cannot see deletions <50 nt or
  substitutions), so screening-depth cell counts are **minima**.
- `cross_well` reads are largely ligation-fused molecules formed during sequencing library
  prep, not chimeric constructs — check the read-length-by-class sheet before quoting a
  combined "total mispaired" figure.
- A metric that cannot be computed is reported as `n/a` with the reason, never as a
  placeholder number. A 95/5 ratio is undefined when the 5th percentile is zero.

---

## Reference CSV

One row per designed construct. Required columns:

| column | meaning |
|---|---|
| `gene` | unique construct id |
| `order` | sub-library id; must match `sublibrary` in the samplesheet |
| `well` | high-density synthesis well — the unit `within_well` chimera is defined against |
| `vh_len_nt` | VH length, used to cut the half-references |
| `construct_dna` | the full printed construct, flanks included |

Optional: `framework`, `pool`, `well_rep2`, `well_60perwell`, `length_nt`,
`intermediate_len_nt`, `gc_pct`.

`BUILD_REFERENCE` re-derives both FASTAs from `construct_dna` and `vh_len_nt` on every run,
so the half-references cannot drift out of step with the constructs. It checks flanks,
declared length, reading frame, and **half-reference uniqueness** — if two constructs share
a VH-side or VL-side half the classifier cannot tell them apart, and the run says so.
`--strict_reference` turns those warnings into a failure.

`bin/make_reference.py` in the design repo builds this CSV from the order files and verifies
it against the writer files that went to the printer.

---

## Validation

The Nextflow port was checked against the Databricks analysis of the same run (Q-724749,
AILK014-007/-008/-009, 93.5 M reads). Running the pipeline on a 720-construct slice of that
reference, with the reads that map to it:

| | pipeline (slice) | Databricks (full depth) |
|---|---|---|
| Correctly paired, -007 / -008 / -009 | 99.86 / 100.00 / 99.72% | 99.46 / 99.47 / 99.70% |
| ≥50 nt block deletion | 6.79 / 7.48 / 5.34% | 6.14 / 6.98 / 5.88% |
| Intact full-length (both-primer) | 93.2 / 92.5 / 94.7% | 93.8 / 93.0 / 94.1% |

Agreement is within what a few-hundred-read slice supports. Uniformity, dropout and
screening depth are **not** comparable at slice depth and are not listed. See
[docs/VALIDATION.md](docs/VALIDATION.md).

---

## Running the tests

```bash
nextflow run . -profile local        # end-to-end on the bundled fixtures, tools from PATH
nextflow run . -profile stub -stub   # wiring only, no tools, no data, no Wave
```

The fixtures are ~1,800 real reads against a 720-construct reference slice. They exist to
prove the pipeline runs end to end and to catch regressions — **the rates they produce are
not the library's rates**.

---

## Compute notes

- Wave builds each module's `conda "${moduleDir}/environment.yml"` into a cached ECR image.
  Do **not** set `conda.enabled = true` globally: that builds a native conda env in the work
  dir, which fails on an S3 work dir.
- `CLASSIFY_READS` and `BLOCKDEL_READS` are one task per FASTQ file. Both hold the
  half-reference index in memory (~1 GB for 136k halves); `BLOCKDEL_READS` additionally
  caches up to 512 per-construct splice indexes. Do not raise that cap — a large cache of
  live minimap2 indexes is what exhausts a node.
- `mappy` is pinned to 2.28. A different minimap2 recovers marginally different alignments,
  which moves the chimera rate in the fourth decimal place.
