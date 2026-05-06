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
    # [stub] cutadapt -g ${params.f_primer} -a ${params.r_primer} --revcomp ...
    touch ${meta.id}_trimmed.fastq.gz
    printf 'sample\\tcount\\n%s\\t0\\n' "${meta.id}" > ${meta.id}.cutadapt.counts.tsv
    """
}
