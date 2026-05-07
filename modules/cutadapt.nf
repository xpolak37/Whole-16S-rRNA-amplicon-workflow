process CUTADAPT {
    tag "${meta.id}"
    container 'quay.io/biocontainers/cutadapt:4.6--py39hf95cd2a_1'
    cpus 4

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    path  "${meta.id}.cutadapt.counts.tsv",      emit: counts

    script:
    """
    cutadapt \\
        --quiet \\
        --cores ${task.cpus} \\
        --revcomp \\
        -g ${params.f_primer} \\
        -a ${params.r_primer} \\
        -a 'A{10}' -a 'G{10}' \\
        -o ${meta.id}_trimmed.fastq.gz \\
        ${reads}

    count=\$(zcat ${meta.id}_trimmed.fastq.gz | awk 'END{print NR/4}')
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.cutadapt.counts.tsv
    """
}
