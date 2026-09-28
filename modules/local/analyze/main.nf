process ANALYZE {
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    path metrics, stageAs: 'metrics/*'
    path reference
    path meta

    output:
    path "derived.json", emit: derived

    script:
    def mh = params.microhomology ? "--microhomology ${params.microhomology}" : ""
    def ds = params.design_screen  ? "--design-screen ${params.design_screen}" : ""
    """
    analyze.py --metrics-dir metrics --reference ${reference} --meta ${meta} ${mh} ${ds} --out derived.json
    """

    stub:
    """
    echo '{}' > derived.json
    """
}
