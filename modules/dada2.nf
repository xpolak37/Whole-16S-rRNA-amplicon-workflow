process DADA2_PACBIO {
    tag "dada2"
    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'
    cpus 8
    publishDir "${params.outdir}/dada2/dada2", mode: 'copy'

    input:
    path '*'   // all oriented fastqs collected

    output:
    tuple val('dada2'), path('ASV_table.tsv'), path('ASV_sequences.fasta'), emit: asv
    path  'track_control.tsv',                                              emit: track
    path  'dada2.counts.tsv',                                               emit: counts

    script:
    """
    # [stub] real R script lives in bin/dada2_pacbio.R; Nextflow auto-runs in this task's work dir
    dada2_pacbio.R . . ${params.minQ} ${params.minLen} ${params.maxLen} ${params.maxN} ${params.maxEE}
    printf 'sample\\tcount\\n_cohort_\\t0\\n' > dada2.counts.tsv
    """
}

process DADA2_PACBIO_NODENOISE {
    tag "dada2_nodenoise"
    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'
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
    dada2_pacbio_nodenoise.R . . ${params.minQ} ${params.minLen} ${params.maxLen} ${params.maxN} ${params.maxEE} ${params.minAbundance}
    printf 'sample\\tcount\\n_cohort_\\t0\\n' > dada2_nodenoise.counts.tsv
    """
}
