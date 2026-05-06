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
    // FastQC names outputs after the input file's basename, so the raw
    // (sample.fastq.gz) and trimmed (sample_trimmed.fastq.gz) passes can
    // share the same publishDir without their reports colliding.
    """
    fastqc \\
        --threads ${task.cpus} \\
        --quiet \\
        ${reads}
    """
}
