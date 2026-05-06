process HOST_REMOVAL {
    tag "${meta.id}"
    container 'quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0'
    cpus 4

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.clean.fastq.gz"),  emit: reads
    path  "${meta.id}.host_removal.counts.tsv", emit: counts

    script:
    """
    # [stub] hostile clean --index ${params.hostile_index_dir} ...
    touch ${meta.id}.clean.fastq.gz
    printf 'sample\\tcount\\n%s\\t0\\n' "${meta.id}" > ${meta.id}.host_removal.counts.tsv
    """
}

process PHIX_REMOVAL {
    tag "${meta.id}"
    container 'quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0'
    cpus 4

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_cleaned.fastq.gz"), emit: reads
    path  "${meta.id}.phix_removal.counts.tsv",  emit: counts

    script:
    """
    # [stub] hostile clean --index ${params.phix_fasta} ...
    touch ${meta.id}_trimmed_cleaned.fastq.gz
    printf 'sample\\tcount\\n%s\\t0\\n' "${meta.id}" > ${meta.id}.phix_removal.counts.tsv
    """
}
