// Host + PhiX depletion via direct minimap2 calls with the map-hifi preset.
//
// hostile 1.1.0 hardcodes minimap2's `-x map-ont` preset, which assumes
// ~5–15 % error rate and under-calls alignments for HiFi data (~0.1 %).
// We bypass `hostile clean` and call minimap2 directly, using the same
// container (which ships minimap2 + samtools as runtime deps). Unmapped
// reads (FLAG 0x4) are extracted with `samtools fastq -f 4` and gzipped.

process HOST_REMOVAL {
    tag "${meta.id}"
    publishDir "${params.outdir}/hostile", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.clean.fastq.gz"), emit: reads
    path  "${meta.id}.host_removal.counts.tsv",         emit: counts

    script:
    """
    minimap2 -ax map-hifi -t ${task.cpus} \\
        ${params.hostile_index_dir}/human-t2t-hla-argos985-mycob140.mmi \\
        ${reads} \\
      | samtools fastq -n -f 4 - \\
      | gzip > ${meta.id}.clean.fastq.gz

    count=\$(zcat ${meta.id}.clean.fastq.gz | awk 'END{print NR/4}')
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.host_removal.counts.tsv
    """
}

process PHIX_REMOVAL {
    tag "${meta.id}"
    publishDir "${params.outdir}/hostile", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed_cleaned.fastq.gz"), emit: reads
    path  "${meta.id}.phix_removal.counts.tsv",                   emit: counts

    script:
    """
    minimap2 -ax map-hifi -t ${task.cpus} \\
        ${params.phix_index} \\
        ${reads} \\
      | samtools fastq -n -f 4 - \\
      | gzip > ${meta.id}_trimmed_cleaned.fastq.gz

    count=\$(zcat ${meta.id}_trimmed_cleaned.fastq.gz | awk 'END{print NR/4}')
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.phix_removal.counts.tsv
    """
}
