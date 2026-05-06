process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    cpus 2

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    """
    multiqc . \\
        --title "${params.multiqc_title}" \\
        --filename multiqc_report \\
        --force \\
        --interactive \\
        ${config} \\
        ${params.multiqc_extra_args}
    """
}
