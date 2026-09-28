process RUN_CONTROLS {
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    path halves
    path constructs
    path meta

    output:
    path "controls.json", emit: json

    script:
    """
    run_controls.py \\
        --halves ${halves} \\
        --constructs ${constructs} \\
        --meta ${meta} \\
        --n-negative ${params.n_control_negative} \\
        --n-positive ${params.n_control_positive} \\
        --n-del ${params.n_control_deletion} \\
        --flank5 ${params.flank5} \\
        --flank3 ${params.flank3} \\
        --out controls.json
    """

    stub:
    """
    echo '{"classifier":{"negative":{},"positive":{}},"deletion_caller":{}}' > controls.json
    """
}
