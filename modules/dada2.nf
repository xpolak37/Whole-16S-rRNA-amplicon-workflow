process DADA2_PACBIO {
    tag "dada2"
    container 'quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0'
    cpus 8
    publishDir "${params.outdir}/dada2/dada2", mode: 'copy'

    input:
    path '*'

    output:
    tuple val('dada2'), path('ASV_table.tsv'), path('ASV_sequences.fasta'), emit: asv
    path  'track_control.tsv',                                              emit: track
    path  'dada2.counts.tsv',                                               emit: counts

    script:
    """
    Rscript ${projectDir}/bin/dada2_pacbio.R \\
        --input . \\
        --nproc ${task.cpus} \\
        --minQ ${params.minQ} \\
        --minLen ${params.minLen} \\
        --maxLen ${params.maxLen} \\
        --maxN ${params.maxN} \\
        --maxEE ${params.maxEE}
    """
}

process DADA2_PACBIO_NODENOISE {
    tag "dada2_nodenoise"
    container 'quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0'
    cpus 8
    publishDir "${params.outdir}/dada2/dada2_nodenoise", mode: 'copy'

    input:
    path '*'

    output:
    tuple val('dada2_nodenoise'), path('ASV_table.tsv'), path('ASV_sequences.fasta'), emit: asv
    path  'track_control.tsv',                                                        emit: track
    path  'dada2_nodenoise.counts.tsv',                                               emit: counts

    script:
    """
    Rscript ${projectDir}/bin/dada2_pacbio_nodenoise.R \\
        --input . \\
        --nproc ${task.cpus} \\
        --minQ ${params.minQ} \\
        --minLen ${params.minLen} \\
        --maxLen ${params.maxLen} \\
        --maxN ${params.maxN} \\
        --maxEE ${params.maxEE} \\
        --minAbundance ${params.minAbundance}
    """
}
