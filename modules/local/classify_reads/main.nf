process CLASSIFY_READS {
    tag "${sublibrary}/${fastq.name}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sublibrary), path(fastq)
    path halves
    path meta

    output:
    path "*.classify.json", emit: json

    script:
    def max_reads = params.max_reads_per_shard ? "--max-reads ${params.max_reads_per_shard}" : ''
    """
    classify_reads.py \\
        --fastq ${fastq} \\
        --halves ${halves} \\
        --meta ${meta} \\
        --sublibrary ${sublibrary} \\
        --margin ${params.margin} \\
        --min-mlen ${params.min_mlen} ${max_reads} \\
        --out ${sublibrary}.${fastq.simpleName}.classify.json
    """

    stub:
    """
    echo '{"sublibrary":"${sublibrary}","fastq":"${fastq.name}","n_reads":0,"class_counts":{},"per_gene":{},"per_well":{},"rlen_hist_25bp":{},"rlen_hist_25bp_by_cls":{},"length_classes":{},"params":{}}' > ${sublibrary}.${fastq.simpleName}.classify.json
    """
}
