process CUTADAPT {
    tag "${meta.id}"

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
        -e 0.2 \\
        -g "^${params.f_primer}" \\
        -a "${params.r_primer}\$" \\
        -a 'A{10}' -a 'G{10}' \\
        -o ${meta.id}_trimmed.fastq.gz \\
        ${reads}

    count=\$(zcat ${meta.id}_trimmed.fastq.gz | awk 'END{print NR/4}')
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.cutadapt.counts.tsv
    """
}
