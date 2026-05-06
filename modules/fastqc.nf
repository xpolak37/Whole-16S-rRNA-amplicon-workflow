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
    // Mirror real FastQC: name outputs after the input file's basename so the
    // raw and trimmed passes don't collide downstream (e.g. in MULTIQC).
    def base = reads.getName().replace('.fastq.gz', '').replace('.fq.gz', '').replace('.fastq', '').replace('.fq', '')
    """
    # [stub] FASTQC on ${meta.id} (output base: ${base})
    touch ${base}_fastqc.html ${base}_fastqc.zip
    """
}
