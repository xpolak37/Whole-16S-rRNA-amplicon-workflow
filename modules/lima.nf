// Demultiplex a non-demuxed PacBio HiFi BAM into per-barcode fastq.gz.
// Lima emits .fastq.gz directly when the output filename ends in .fastq.gz,
// so no bam2fastq step is needed. --split-named --split produces one file
// per barcode pair named demux.<bcF>--<bcR>.fastq.gz; we strip the shared
// "demux." prefix so each file's simpleName (consumed as the sample id
// downstream) is the barcode-pair name.
//
// ASYMMETRIC matches PacBio's standard 16S kit (forward x reverse pairing).
// Override via params.lima_extra_args, e.g. --hifi-preset SYMMETRIC for
// symmetric kits. --peek-guess limits output to barcode pairs actually
// observed in the run, so unused entries in the kit FASTA don't produce
// phantom output files.
process LIMA_DEMUX {
    publishDir "${params.outdir}/lima", mode: 'copy'

    input:
    tuple path(bam), path(barcodes)

    output:
    path "*.fastq.gz",         emit: fastqs
    path "demux.lima.summary", emit: summary
    path "demux.lima.counts",  emit: counts

    script:
    """
    lima \\
        --hifi-preset ASYMMETRIC \\
        --peek-guess \\
        --split-named --split \\
        --num-threads ${task.cpus} \\
        ${bam} ${barcodes} demux.fastq.gz \\
        ${params.lima_extra_args}

    shopt -s nullglob
    files=( demux.*--*.fastq.gz )
    if [ \${#files[@]} -eq 0 ]; then
        echo "lima produced no per-barcode fastq output. Check demux.lima.summary." >&2
        exit 1
    fi
    for f in "\${files[@]}"; do
        mv "\$f" "\${f#demux.}"
    done
    """
}
