process BLOCKDEL_READS {
    tag "${sublibrary}/${fastq.name}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sublibrary), path(fastq)
    path halves
    path constructs
    path meta

    output:
    path "*.blockdel.json",  emit: json
    path "*.events.tsv",     emit: events

    script:
    """
    blockdel_reads.py \\
        --fastq ${fastq} \\
        --halves ${halves} \\
        --constructs ${constructs} \\
        --meta ${meta} \\
        --sublibrary ${sublibrary} \\
        --subsample ${params.blockdel_subsample} \\
        --min-del ${params.min_del} \\
        --call-del ${params.call_del} \\
        --flank5 ${params.flank5} \\
        --flank3 ${params.flank3} \\
        --cpus ${task.cpus} \\
        --out ${sublibrary}.${fastq.simpleName}.blockdel.json \\
        --events-out ${sublibrary}.${fastq.simpleName}.events.tsv
    """

    stub:
    """
    echo '{"sublibrary":"${sublibrary}","fastq":"${fastq.name}","subsample":8,"n_subsampled":0,"primer_ladder":{},"strand":{},"eligible":{},"calls":{},"per_gene":{},"gap_size_hist_10nt":{},"gap_start_hist_10nt":{},"gap_end_hist_10nt":{},"n_gap_per_read":{},"full_length":{},"identity_median":{"gap50":null,"no_gap":null},"identity_n":{"gap50":0,"no_gap":0},"edlib_confirm":{"n":0,"confirmed":0},"params":{}}' > ${sublibrary}.${fastq.simpleName}.blockdel.json
    printf 'sublib\\tgene\\tgap_start\\tgap_len\\tqlen\\treflen\\n' > ${sublibrary}.${fastq.simpleName}.events.tsv
    """
}
