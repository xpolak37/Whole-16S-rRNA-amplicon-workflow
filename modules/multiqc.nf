process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    # [stub] multiqc
    touch multiqc_report.html
    """
}
