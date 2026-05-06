process FASTQC {
    tag "${meta.id}"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    cpus 2

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    # [stub] FASTQC on ${meta.id}
    touch ${meta.id}_fastqc.html ${meta.id}_fastqc.zip
    """
}
