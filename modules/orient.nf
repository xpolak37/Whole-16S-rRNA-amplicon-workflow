process VSEARCH_ORIENT {
    tag "${meta.id}"
    publishDir "${params.outdir}/orient", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*-oriented.fq"),  emit: reads
    path  "${meta.id}.orient.counts.tsv",    emit: counts

    script:
    """
    vsearch \\
        --orient ${reads} \\
        --db ${params.silva_orient_db}/silva-27F-1492R-orient.fasta \\
        --fastqout ${meta.id}-oriented.fq \\
        --threads ${task.cpus}

    count=\$(awk 'END{print NR/4}' ${meta.id}-oriented.fq)
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.orient.counts.tsv
    """
}
