process MAKE_PLOTS {
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    path metrics, stageAs: 'metrics/*'
    path derived
    path meta

    output:
    path "plots/*", emit: plots

    script:
    def cmp = params.compare ? "--compare ${params.compare}" : ""
    """
    cp ${derived} metrics/derived.json
    make_plots.py --metrics-dir metrics --meta ${meta} --outdir plots ${cmp}
    """

    stub:
    """
    mkdir -p plots && touch plots/stub_fig.png
    """
}
