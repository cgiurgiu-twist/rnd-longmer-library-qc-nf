process BUILD_REFERENCE {
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    path reference

    output:
    path "halves.fasta.gz",     emit: halves
    path "constructs.fasta.gz", emit: constructs
    path "gene_meta.json",      emit: meta
    path "reference_qc.json",   emit: qc

    script:
    def strict = params.strict_reference ? '--strict' : ''
    """
    build_reference.py \\
        --reference ${reference} \\
        --flank5 ${params.flank5} \\
        --flank3 ${params.flank3} \\
        --intermediate-len ${params.intermediate_len} ${strict}
    """

    stub:
    """
    printf '>H_g1\\nACGT\\n>L_g1\\nACGT\\n' | gzip > halves.fasta.gz
    printf '>g1\\nACGTACGT\\n' | gzip > constructs.fasta.gz
    echo '{"g1":{"order":"O1","well":"A1","framework":"A","vh_len":4,"full_len":8,"h_len":4,"l_len":4}}' > gene_meta.json
    echo '{"n_constructs":1,"n_problems":0}' > reference_qc.json
    """
}
