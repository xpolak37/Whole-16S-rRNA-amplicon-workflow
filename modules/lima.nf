process LIMA_DEMUX {
    container 'quay.io/biocontainers/lima:2.9.0--h9ee0642_0'
    cpus 4

    input:
    tuple path(bam), path(barcodes)

    output:
    path "*.fastq.gz", emit: fastqs

    script:
    """
    # [stub] lima --isoseq --peek-guess ${bam} ${barcodes} demuxed.fastq.gz ${params.lima_extra_args ?: ''}
    # Stub emits two placeholder per-barcode fastqs to drive downstream channel cardinality.
    echo "@stub_read" | gzip > bc01.fastq.gz
    echo "@stub_read" | gzip > bc02.fastq.gz
    """
}
