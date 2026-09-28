# Deploying to Seqera (Twist BFX)

The pipeline runs on the self-hosted Seqera Platform at `seqera.tools.twistbio.co`, org
**BFX**, workspace **RND** (`73662236814646`). Kefan is a workspace member, not an owner —
registering a pipeline and launching works; creating compute environments or credentials
does not.

## One-time: register the pipeline

Seqera → RND workspace → **Pipelines** → *Add pipeline*:

| field | value |
|---|---|
| Name | `rnd-longmer-library-qc-nf` |
| Repository | `https://github.com/keyang-twist/rnd-longmer-library-qc-nf` |
| Revision | `main` (or a release tag once the repo is moved to the org) |
| Compute environment | `rnd_aws_batch_r6i-c6i-m6i` (`6zyQ3qfyeeQiwCqpOOMpqH`) |
| Work directory | `s3://twist-nonprod-rnd-nextflow-work` |

**Pre-run script** — pins the Nextflow version the pipeline is validated on:

```bash
export NXF_VER=25.10.4
```

**Nextflow config** — lets the compute env's account read the outputs it writes:

```groovy
aws { client { s3Acl = 'BucketOwnerFullControl'; uploadChunkSize = 10485760 } }
```

The repo is currently under the personal account `keyang-twist` and is **private**, so
Seqera needs a GitHub credential that can read it, or the repo needs moving to
`Twistbioscience` (which is where a production pipeline belongs — name it `bfx-` rather
than `rnd-` if it becomes customer-facing, and run `.github/setup-repo.sh` there for branch
protection and squash-merge).

## Launch parameters

```yaml
reference: s3://<bucket>/<path>/AILK014-007-008-009_750mer_reference.csv
input:     s3://<bucket>/<path>/samplesheet.csv
outdir:    s3://twist-nonprod-rnd-nextflow-work/kefan/AILK014_007-009_QC
library_name: AILK014-007-008-009
account: Absci
order_number: Q-724749
order_items: 'AILK014-007, 008, 009'
manufacture_date: '2026-09-25'
```

Give every run an explicit name (`kefan-ailk014-007-009-qc-YYYYMMDD`). Seqera otherwise
assigns a random pair like `suspicious_meitner`, which is unfindable a week later among
~2,600 runs.

## Before the first real launch

1. **Upload the reference and samplesheet to S3.** The reference is ~53 MB uncompressed;
   the pipeline reads it once in `BUILD_REFERENCE`.
2. **Dry-run the wiring**: launch with `stubRun: true`. It executes no real work and proves
   the params resolve and the compute env accepts the job.
3. **Then launch for real.** ~400 FASTQ files means ~800 tasks across the two per-file
   passes; the first run also pays for Wave building each module's image (cached in ECR
   afterwards).

## Checking the run

- `reference/reference_qc.json` → `n_problems` must be 0.
- `metrics/controls.json` → the run is not reportable without clean calibration.
- `reports/<name>_qc_summary.md` → headline table, denominators, calibration.
- `plots/figures_skipped.json` → each skipped figure and why.

`tower.yml` surfaces those files in the Seqera run page's Reports tab.

## If a task fails

Work from the failing task outward, not from the head-job log: `GET /workflow/<id>/tasks`,
find the `FAILED` task, note `exitStatus` and `workdir`, then read `.command.err` there.

| exit | cause |
|---|---|
| 137 | OOM — raise the label's memory; `CLASSIFY_READS`/`BLOCKDEL_READS` hold the half index |
| "no space left on device" | use a `_16TB` compute environment |
| `AccessDenied` on S3 | the compute env's credentials do not cover that bucket — a DevOps request, not something to work around |

Relaunch with `resume: true` and the failed run's `sessionId` after a transient failure so
completed tasks are not recomputed.
