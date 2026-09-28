process MAKE_OUTPUTS {
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    path metrics, stageAs: 'metrics/*'
    path derived
    path meta
    path plots, stageAs: 'plots/*'
    path logo

    output:
    path "*.xlsx",        emit: workbook
    path "*_qc_summary.md", emit: summary
    path "*_QC_Summary.docx", emit: collateral, optional: true

    script:
    def mh   = params.microhomology ? "--microhomology ${params.microhomology}" : ""
    def logo_arg = logo.name != 'NO_LOGO' ? "--logo ${logo}" : ""
    """
    cp ${derived} metrics/derived.json
    make_workbook.py --metrics-dir metrics --meta ${meta} --outdir . --name ${params.library_name} ${mh}
    make_collateral.py \\
        --metrics-dir metrics --plots-dir plots \\
        --product '${params.product}' --account '${params.account}' \\
        --order-number '${params.order_number}' --order-items '${params.order_items}' \\
        --manufacture-date '${params.manufacture_date}' ${logo_arg} \\
        --out '${params.library_name}_QC_Summary.docx'
    """

    stub:
    """
    touch ${params.library_name}_nanopore_QC.xlsx ${params.library_name}_qc_summary.md
    """
}
