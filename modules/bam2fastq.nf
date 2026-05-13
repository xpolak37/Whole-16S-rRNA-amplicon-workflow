// Convert a per-sample PacBio HiFi BAM to fastq.gz.
// Used when the samplesheet 'fastq' column points at a .bam file (already
// demultiplexed PacBio data). For the non-demuxed --bam mode, lima emits
// fastq.gz directly so this process is not used there.
//
// `cp -L` dereferences Nextflow's symlink staging so pbindex can write the
// .pbi sibling into the writeable work directory rather than the original
// (potentially read-only) data location. pbindex is idempotent — we run it
// unconditionally rather than threading optional .pbi inputs through.
process BAM2FASTQ {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/bam2fastq" }, mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads

    script:
    """
    cp -L ${bam} reads.bam
    pbindex reads.bam
    bam2fastq -o ${meta.id} -j ${task.cpus} reads.bam
    """
}
