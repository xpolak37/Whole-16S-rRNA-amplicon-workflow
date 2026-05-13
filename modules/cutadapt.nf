process CUTADAPT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    path  "${meta.id}.cutadapt.counts.tsv",      emit: counts

    script:
    // Two-pass cutadapt. Pass 1 requires both primers and drops any read
    // missing one (--discard-untrimmed) — kills off-target / partial CCS reads
    // before they reach DADA2. Pass 2 trims residual polyA/polyG tails; it's
    // a no-op for clean HiFi but prevents tail variants from inflating ASVs
    // when present. Counts are post-pass-1 (pass 2 trims, does not drop).
    """
    cutadapt \\
        --quiet \\
        --cores ${task.cpus} \\
        --revcomp \\
        --discard-untrimmed \\
        -e 0.2 \\
        -g "^${params.f_primer}" \\
        -a "${params.r_primer}\$" \\
        -o ${meta.id}.primer.fastq.gz \\
        ${reads}

    cutadapt \\
        --quiet \\
        --cores ${task.cpus} \\
        -a 'A{10}' -a 'G{10}' \\
        -o ${meta.id}_trimmed.fastq.gz \\
        ${meta.id}.primer.fastq.gz

    count=\$(zcat ${meta.id}_trimmed.fastq.gz | awk 'END{print NR/4}')
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.cutadapt.counts.tsv
    """
}
