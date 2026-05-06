process VSEARCH_ORIENT {
    tag "${meta.id}"
    container 'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_1'
    cpus 4
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
        --db ${params.silva_orient_db} \\
        --fastqout ${meta.id}-oriented.fq \\
        --threads ${task.cpus}

    count=\$(awk 'END{print NR/4}' ${meta.id}-oriented.fq)
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.orient.counts.tsv
    """
}
