process VSEARCH_ORIENT {
    tag "${meta.id}"
    container 'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_1'
    cpus 4

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*-oriented.fq"),  emit: reads
    path  "${meta.id}.orient.counts.tsv",    emit: counts

    script:
    """
    # [stub] vsearch --orient ${reads} --db ${params.silva_orient_db} --fastqout ${meta.id}-oriented.fq
    touch ${meta.id}-oriented.fq
    printf 'sample\\tcount\\n%s\\t0\\n' "${meta.id}" > ${meta.id}.orient.counts.tsv
    """
}
