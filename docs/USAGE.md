# Ops runbook

## What you need

1. **The design reference CSV** for the library, from the design team. One row per designed
   construct. It is the only design input — no FASTA to assemble, no barcode table to keep
   in step.
2. **The fastq_pass prefix for each sub-library.** One PromethION flowcell per sub-library
   is the usual layout, and the pipeline relies on that: each sub-library's reads are
   classified against the whole library, so reads that pair across sub-libraries are a
   built-in false-positive control. If several sub-libraries shared a flowcell, say so —
   that control is then not available and the numbers need a caveat.

## Running it

Write a samplesheet with one row per sub-library:

```csv
sublibrary,fastq_dir
AILK014-007,s3://twist-instrument-data/ngs/PromethION/20260924_1825_P2I-00342-B_260924_Q-724749_AILK0114-007/fastq_pass
AILK014-008,s3://twist-instrument-data/ngs/PromethION/20260924_1824_P2I-00341-A_260924_Q-724749_AILK014-008/fastq_pass
AILK014-009,s3://twist-instrument-data/ngs/PromethION/20260924_1824_P2I-00342-A_260924_Q-724749_AILK014-009/fastq_pass
```

`sublibrary` must match the `order` column in the reference exactly. Then launch from
Seqera (RND workspace) or the CLI:

```bash
nextflow run Twistbioscience/rnd-longmer-library-qc-nf -r main \
  --reference s3://.../AILK014-007-008-009_750mer_reference.csv \
  --input     s3://.../samplesheet.csv \
  --outdir    s3://.../AILK014_007-009_QC \
  --library_name AILK014-007-008-009 \
  --account Absci --order_number Q-724749 \
  --order_items 'AILK014-007, 008, 009' --manufacture_date 2026-09-25
```

The collateral header fields (`--account`, `--order_number`, `--order_items`,
`--manufacture_date`) only affect the one-page customer summary. Leave them out and those
cells are highlighted placeholders rather than guesses — but then the collateral is not
ready to send.

## Before you send anything

1. Open `reports/<name>_qc_summary.md`. It has the headline table, the denominators, and
   the calibration results.
2. Check `reference/reference_qc.json`: `n_problems` should be 0. Anything else means the
   reference and the printed library may not agree, and every number downstream inherits
   that.
3. Check `metrics/controls.json`. A run whose classifier produced false chimera calls, or
   whose deletion caller lost sensitivity, is not reportable — the rates are not defensible
   without these.
4. Check `plots/figures_skipped.json`. Each entry says which figure was skipped and why.
   Missing optional inputs (`--compare`, `--microhomology`) are expected; anything else is
   worth understanding before the numbers go out.
5. `reports/<name>_QC_Summary.docx` is the customer one-pager. `reports/<name>_nanopore_QC.xlsx`
   is the full data.

The **interpretive report** — what the numbers mean for this build — is written by the
project scientist from these outputs. The pipeline does not write conclusions.

## Timing and cost

One task per FASTQ file for the pairing pass, one per FASTQ file for the deletion pass.
A three-flowcell PromethION run (~400 files, ~93 M reads) is a few hundred tasks. Budget
from real reads, not simulated ones: the classifier runs ~400–800 reads/s per core on real
data, which includes junk, adapters, truncations and concatamers.

`--blockdel_subsample` defaults to 8. The deletion rate converges far below full depth;
raising it to 1 costs 8× the compute for no change in the reported rate.

## When something fails

| Symptom | Cause and fix |
|---|---|
| `No files matching <pattern> for sub-library X` | `fastq_dir` is wrong, or the run has no `.fastq.gz` (check for `fastq_fail` vs `fastq_pass`) |
| `reference is missing required columns` | Wrong CSV — check for the design-side order file rather than the built reference |
| `reference_qc.json` reports non-unique halves | Two constructs share a VH-side or VL-side half; the classifier cannot tell them apart. Design issue — raise it before reporting |
| Task killed, exit 137 | Out of memory. `CLASSIFY_READS`/`BLOCKDEL_READS` hold the half index; raise `process_low` memory |
| `no space left on device` | Use a `_16TB` compute environment |
| A sub-library shows `n/a` for 95/5 or screening depth | Its 5th-percentile variant has zero reads — that is a real result (dropout), not a bug. The summary names the reason |
| Wave `Missing credentials` locally | Wave needs Seqera credentials. Use `-profile local` (tools from PATH) or `-profile stub -stub` (wiring only) |
