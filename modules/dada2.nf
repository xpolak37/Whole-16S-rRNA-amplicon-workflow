process DADA2_PACBIO {
    tag "dada2"
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

// ============================================================
// Split DADA2 path: filter -> learn_errors -> denoise -> merge.
//
// Same statistics as DADA2_PACBIO, different scheduling. dada() runs with
// pool = FALSE, so samples are already denoised independently; splitting them
// into separate tasks changes when work happens, not what it computes. Verified
// byte-identical against DADA2_PACBIO on synthetic reads (dada2 1.38.0).
//
// The cohort-level script held a derep object -- unique sequences plus their
// quality matrix -- for every sample simultaneously, which is what pushed peak
// RSS to 39.4 GB at 30 samples and made larger cohorts run out of memory. Here
// each task holds one sample, and the merge stage sees only denoised ASVs.
// ============================================================

process DADA2_FILTER {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_filt.fastq.gz"), emit: filt, optional: true
    path "${meta.id}.filter_stats.tsv",                emit: stats

    script:
    """
    Rscript ${projectDir}/bin/dada2_filter.R \\
        --input ${reads} \\
        --sample ${meta.id} \\
        --nproc ${task.cpus} \\
        --minQ ${params.minQ} \\
        --minLen ${params.minLen} \\
        --maxLen ${params.maxLen} \\
        --maxN ${params.maxN} \\
        --maxEE ${params.maxEE}
    """
}

process DADA2_LEARN_ERRORS {
    tag "dada2/learn_errors"
    publishDir "${params.outdir}/dada2/dada2", mode: 'copy', pattern: 'qualBins.txt'

    input:
    path '*'

    output:
    path 'err.rds',     emit: err
    path 'qualBins.txt', emit: qualbins

    script:
    """
    Rscript ${projectDir}/bin/dada2_learn_errors.R \\
        --input . \\
        --nproc ${task.cpus} \\
        --nbases ${params.dada2_learn_nbases} \\
        --randomize ${params.dada2_learn_randomize ? 'TRUE' : 'FALSE'} \\
        --band_size ${params.dada2_band_size}
    """
}

process DADA2_DENOISE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(filt), path(err)

    output:
    path "${meta.id}.uniques.tsv",        emit: uniques
    path "${meta.id}.denoise_stats.tsv",  emit: stats

    script:
    """
    Rscript ${projectDir}/bin/dada2_denoise.R \\
        --input ${filt} \\
        --sample ${meta.id} \\
        --err ${err} \\
        --nproc ${task.cpus} \\
        --band_size ${params.dada2_band_size}
    """
}

process DADA2_MERGE {
    tag "dada2/merge"
    publishDir "${params.outdir}/dada2/dada2", mode: 'copy'

    input:
    path 'uniq/*'
    path 'fstats/*'
    path 'dstats/*'

    output:
    tuple val('dada2'), path('ASV_table.tsv'), path('ASV_sequences.fasta'), emit: asv
    path  'track_control.tsv',                                              emit: track
    path  'dada2.counts.tsv',                                               emit: counts

    script:
    """
    Rscript ${projectDir}/bin/dada2_merge.R \\
        --uniques uniq \\
        --filter_stats fstats \\
        --denoise_stats dstats \\
        --nproc ${task.cpus} \\
        --min_fold ${params.dada2_min_fold} \\
        --counts_name dada2.counts.tsv
    """
}

// Split no-denoise path: reuses DADA2_FILTER, then derep per sample and merge.
//
// This path merges raw uniques rather than dada() output, so the merge stage --
// not the derep stage -- is the scaling problem: makeSequenceTable builds a
// dense samples x sequences matrix over every distinct sequence in the cohort.
// dada2_nodenoise_merge.R prefilters sequences that provably cannot clear
// minAbundance before building that matrix. See the header of that script for
// the argument; the result is identical to the cohort-level table.

process DADA2_DEREP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(filt)

    output:
    path "${meta.id}.uniques.tsv.gz",   emit: uniques
    path "${meta.id}.derep_stats.tsv",  emit: stats

    script:
    """
    Rscript ${projectDir}/bin/dada2_derep.R \\
        --input ${filt} \\
        --sample ${meta.id}
    """
}

process DADA2_NODENOISE_MERGE {
    tag "dada2_nodenoise/merge"
    publishDir "${params.outdir}/dada2/dada2_nodenoise", mode: 'copy'

    input:
    path 'uniq/*'
    path 'fstats/*'
    path 'dstats/*'

    output:
    tuple val('dada2_nodenoise'), path('ASV_table.tsv'), path('ASV_sequences.fasta'), emit: asv
    path  'track_control.tsv',                                                        emit: track
    path  'dada2_nodenoise.counts.tsv',                                               emit: counts

    script:
    """
    Rscript ${projectDir}/bin/dada2_nodenoise_merge.R \\
        --uniques uniq \\
        --filter_stats fstats \\
        --derep_stats dstats \\
        --nproc ${task.cpus} \\
        --min_fold ${params.dada2_min_fold} \\
        --minAbundance ${params.minAbundance} \\
        --counts_name dada2_nodenoise.counts.tsv
    """
}

process DADA2_PACBIO_NODENOISE {
    tag "dada2_nodenoise"
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
