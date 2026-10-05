#!/usr/bin/env nextflow

/*
 * rnd-longmer-library-qc-nf
 *
 * Nanopore QC for direct-print long-oligo (750mer) paired VH+VL libraries.
 *
 * Two inputs: the design reference CSV and the run's fastq_pass prefix(es).
 * Everything else is derived, so the same command works for any library built on this
 * construct architecture.
 *
 *   nextflow run . \
 *     --reference s3://.../AILK014-007-008-009_750mer_reference.csv \
 *     --input     samplesheet.csv \
 *     --outdir    s3://.../results
 */

nextflow.enable.dsl = 2

include { BUILD_REFERENCE } from './modules/local/build_reference/main.nf'
include { CLASSIFY_READS  } from './modules/local/classify_reads/main.nf'
include { BLOCKDEL_READS  } from './modules/local/blockdel_reads/main.nf'
include { RUN_CONTROLS    } from './modules/local/run_controls/main.nf'
include { MERGE_METRICS   } from './modules/local/merge_metrics/main.nf'
include { ANALYZE         } from './modules/local/analyze/main.nf'
include { MAKE_PLOTS      } from './modules/local/make_plots/main.nf'
include { MAKE_OUTPUTS    } from './modules/local/make_outputs/main.nf'

def helpMessage() {
    log.info """
    rnd-longmer-library-qc-nf — nanopore QC for direct-print long-oligo VH+VL libraries

    Required
      --reference   Design reference CSV. Columns: gene, order, well, vh_len_nt,
                    construct_dna (+ optional framework, pool, well_rep2, length_nt,
                    intermediate_len_nt, gc_pct). One row per designed construct.
      --input       Samplesheet CSV with header: sublibrary,fastq_dir
                    One row per sub-library; fastq_dir is a folder of *.fastq.gz
                    (local or s3://). `sublibrary` must match the `order` column of
                    the reference.
                    Launchpad alternative: --fastq_path plus --sublibrary (one
                    sub-library per run, no samplesheet).
      --outdir      Output directory (local or s3://).

    Common options
      --library_name          Prefix for the workbook / collateral / summary files.
      --blockdel_subsample    Keep every Nth read for the deletion pass (default 8).
      --skip_blockdel         Pairing pass only.
      --skip_controls         Skip calibration (NOT for a reportable run).
      --compare               Prior-build metrics JSON for the head-to-head figure.
      --account / --order_number / --order_items / --manufacture_date / --product
                              Customer collateral header fields.

    See docs/USAGE.md for the Ops runbook and README.md for what each output means.
    """.stripIndent()
}

workflow {

    if (params.help) { helpMessage(); return }
    if (!params.reference) { error "Missing --reference (design reference CSV). Run with --help." }
    if (!params.input && !(params.fastq_path && params.sublibrary)) {
        error "Missing --input (samplesheet CSV: sublibrary,fastq_dir), or --fastq_path plus --sublibrary."
    }
    if (!params.outdir)    { error "Missing --outdir." }

    ch_reference = Channel.fromPath(params.reference, checkIfExists: true).first()

    // Fan out one task per FASTQ file. Per-sub-library folders are listed rather than
    // globbed into one channel so a file can always be traced back to its sub-library —
    // the flowcells are separate, and the cross-sub-library class is the false-positive
    // control, so a mislabelled shard would quietly corrupt that control.
    def glob_fastqs = { sub, dir ->
        def pattern = dir.endsWith('/') ? "${dir}${params.fastq_pattern}" : "${dir}/${params.fastq_pattern}"
        def files = file(pattern)
        if (!files) { error "No files matching ${pattern} for sub-library ${sub}" }
        (files instanceof List ? files : [files]).collect { f -> tuple(sub, f) }
    }

    if (params.input) {
        ch_fastq = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                if (!row.sublibrary || !row.fastq_dir) {
                    error "Samplesheet needs columns 'sublibrary' and 'fastq_dir'; got: ${row}"
                }
                // __PROJDIR__ lets the bundled test samplesheet work from any launch directory
                tuple(row.sublibrary.trim(), row.fastq_dir.trim().replace('__PROJDIR__', "${projectDir}"))
            }
            .flatMap { sub, dir -> glob_fastqs(sub, dir) }
    } else {
        ch_fastq = Channel
            .of(tuple(params.sublibrary.trim(), params.fastq_path.trim()))
            .flatMap { sub, dir -> glob_fastqs(sub, dir) }
    }

    BUILD_REFERENCE(ch_reference)

    CLASSIFY_READS(ch_fastq, BUILD_REFERENCE.out.halves, BUILD_REFERENCE.out.meta)

    if (!params.skip_blockdel) {
        BLOCKDEL_READS(ch_fastq, BUILD_REFERENCE.out.halves, BUILD_REFERENCE.out.constructs,
                       BUILD_REFERENCE.out.meta)
        ch_bd     = BLOCKDEL_READS.out.json.collect()
        ch_events = BLOCKDEL_READS.out.events.collect()
    } else {
        ch_bd     = Channel.value([])
        ch_events = Channel.value([])
    }

    if (!params.skip_controls) {
        RUN_CONTROLS(BUILD_REFERENCE.out.halves, BUILD_REFERENCE.out.constructs, BUILD_REFERENCE.out.meta)
        ch_controls = RUN_CONTROLS.out.json
    } else {
        ch_controls = Channel.value([])
    }

    MERGE_METRICS(CLASSIFY_READS.out.json.collect(), ch_bd, ch_events, BUILD_REFERENCE.out.meta)

    // The controls live beside the metrics they calibrate, so no downstream step can read
    // a rate without its calibration being in the same directory.
    ch_metrics = MERGE_METRICS.out.core
        .mix(MERGE_METRICS.out.extra, MERGE_METRICS.out.per_gene, MERGE_METRICS.out.per_well,
             MERGE_METRICS.out.bd_summary.ifEmpty([]), MERGE_METRICS.out.bd_scope.ifEmpty([]),
             MERGE_METRICS.out.bd_per_gene.ifEmpty([]), MERGE_METRICS.out.bd_events.ifEmpty([]),
             ch_controls)
        .flatten()
        .collect()

    ANALYZE(ch_metrics, ch_reference, BUILD_REFERENCE.out.meta)
    MAKE_PLOTS(ch_metrics, ANALYZE.out.derived, BUILD_REFERENCE.out.meta)

    ch_logo = params.logo ? Channel.fromPath(params.logo, checkIfExists: true).first()
                          : Channel.fromPath("${projectDir}/assets/NO_LOGO").first()

    MAKE_OUTPUTS(ch_metrics, ANALYZE.out.derived, BUILD_REFERENCE.out.meta,
                 MAKE_PLOTS.out.plots.collect(), ch_logo)
}

workflow.onComplete = {
    if (workflow.success) {
        log.info "Done. Outputs in ${params.outdir}"
        log.info "  Read reports/<name>_qc_summary.md first — it carries the denominators and bounds."
    }
    else {
        log.info "Failed: ${workflow.errorMessage}"
    }
}
