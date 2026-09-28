process MERGE_METRICS {
    label 'process_merge'
    conda "${moduleDir}/environment.yml"

    input:
    path classify_json, stageAs: 'classify/*'
    path blockdel_json, stageAs: 'blockdel/*'
    path events,        stageAs: 'events/*'
    path meta,          stageAs: 'ref_gene_meta.json'

    // Outputs are named explicitly rather than globbed. A `*.json` glob would also match
    // the staged reference meta, which is a symlink into another task's work directory —
    // publishing that silently drops the whole publish set.
    output:
    path "metrics_core.json",          emit: core
    path "metrics_extra.json",         emit: extra
    path "per_gene.csv",               emit: per_gene
    path "per_well.csv",               emit: per_well
    path "blockdel_summary.json",      emit: bd_summary,  optional: true
    path "blockdel_scope.json",        emit: bd_scope,    optional: true
    path "blockdel_per_gene_all.csv",  emit: bd_per_gene, optional: true
    path "blockdel_events.csv",        emit: bd_events,   optional: true

    script:
    def bd = blockdel_json ? "--blockdel blockdel/*.json" : ""
    """
    merge_metrics.py --classify classify/*.json ${bd} --meta ${meta} --outdir .
    if ls events/*.tsv >/dev/null 2>&1; then
        merge_events.py --events events/*.tsv --max ${params.max_events} --out blockdel_events.csv
    fi
    """

    stub:
    """
    echo '{"total_reads":0,"n_files":0,"overall":{},"per_sublib":{}}' > metrics_core.json
    echo '{}' > metrics_extra.json
    echo 'gene,n' > per_gene.csv
    echo 'order,well,n_assigned,n_clean,n_within,n_cross' > per_well.csv
    """
}
