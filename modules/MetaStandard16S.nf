process METASTANDARD {
    tag "${denoiser}/${classifier}"
    container 'quay.io/biocontainers/pandas:2.2.1'
    cpus 2
    publishDir path: { "${params.outdir}/metastandard/${denoiser}/${classifier}" }, mode: 'copy'

    input:
    tuple val(denoiser), val(classifier), path(taxa_table), path(asv_table)

    output:
    tuple val(denoiser), val(classifier), path("*.tsv"), emit: tsv

    script:
    """
    python3 ${projectDir}/bin/metastandard16S.py \\
        --asv_table  ${asv_table} \\
        --taxa_table ${taxa_table} \\
        --taxa_tool  ${classifier} \\
        --denoiser   ${denoiser} \\
        --level      ${params.tax_level} \\
        --run_id     ${params.run_id}
    """
}
