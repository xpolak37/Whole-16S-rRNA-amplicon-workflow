process METASTANDARD_PLOTS {
    tag "${denoiser}/${classifier}"
    container 'quay.io/biocontainers/seaborn:0.13.2'
    cpus 2
    publishDir path: { "${params.outdir}/metastandard/${denoiser}/${classifier}" }, mode: 'copy'

    input:
    tuple val(denoiser), val(classifier), path(metastandard_tsv)

    output:
    path "*_barplot.png", emit: barplot
    path "*_heatmap.png", emit: heatmap

    script:
    def prefix = metastandard_tsv.baseName
    """
    python3 ${projectDir}/bin/plot_metastandard.py \\
        --input  ${metastandard_tsv} \\
        --top_n  ${params.metastandard_top_n} \\
        --prefix ${prefix}
    """
}
