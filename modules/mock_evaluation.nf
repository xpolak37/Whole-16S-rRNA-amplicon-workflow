process MOCK_EVALUATION {
    tag "${denoiser}/${classifier}"
    publishDir path: { "${params.outdir}/mock_evaluation/${denoiser}/${classifier}" }, mode: 'copy'

    input:
    tuple val(denoiser), val(classifier), path(metastandard_tsv)
    path(mock_abundance)
    path(mock_taxa)
    path(mock_synonyms)

    output:
    path "*_barplot.png",     emit: barplot,     optional: true
    path "*_metrics.tsv",     emit: metrics,     optional: true
    path "*_composition.tsv", emit: composition, optional: true

    script:
    def workflow_id = "${denoiser}_${classifier}"
    def prefix = "${workflow_id}_mock"
    """
    python3 ${projectDir}/bin/mock_evaluation.py \\
        --metastandard_table ${metastandard_tsv} \\
        --mock_abundance ${mock_abundance} \\
        --mock_taxa ${mock_taxa} \\
        --synonyms ${mock_synonyms} \\
        --mock_pattern "${params.mock_pattern}" \\
        --output_prefix ${prefix} \\
        --top_n ${params.mock_top_n} \\
        --run_id ${params.run_id}
    """
}
