process CUTADAPT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    path  "${meta.id}.cutadapt.counts.tsv",      emit: counts

    script:
    // Two-pass cutadapt. Pass 1 uses an anchored linked adapter
    // "^27F...1492R\$" so BOTH primers are required (5' anchored explicitly,
    // 3' anchored with \$). cutadapt makes the 3' part of a linked adapter
    // optional by default — without the \$ anchor, --discard-untrimmed
    // would still keep reads missing 1492R. --revcomp flips reverse reads;
    // -e 0.1 = ~2 mismatches per primer (IUPAC codes don't count as errors).
    // Pass 2 trims residual polyA/polyG tails; a no-op for clean HiFi but
    // prevents tail variants from inflating ASVs when present. Counts
    // are post-pass-2 (pass 2 trims, does not drop).
    """
    cutadapt \\
        --quiet \\
        --cores ${task.cpus} \\
        --revcomp \\
        --discard-untrimmed \\
        -e 0.1 \\
        -g "^${params.f_primer}...${params.r_primer}\$" \\
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
