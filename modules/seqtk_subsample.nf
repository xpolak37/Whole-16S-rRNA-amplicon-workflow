process SEQTK_SUBSAMPLE {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/subsampled" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.subsampled.fastq.gz"), emit: reads

    script:
    """
    seqtk sample -s100 ${reads} ${params.quick_depth} | gzip > ${meta.id}.subsampled.fastq.gz
    """
}
