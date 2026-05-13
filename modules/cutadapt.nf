process CUTADAPT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    path  "${meta.id}.cutadapt.counts.tsv",      emit: counts

    script:
    // Two-pass cutadapt. Pass 1 uses an anchored linked adapter
    // "^27F...RC(1492R)\$" so BOTH primers are required (5' anchored
    // explicitly, 3' anchored with \$). The 3' half is the reverse complement
    // of the reverse primer — that's what actually appears at the 3' end of
    // a forward-oriented full-length amplicon. cutadapt does NOT auto-RC the
    // second half of a linked adapter, and --revcomp only flips the *read*,
    // not the pattern. cutadapt makes the 3' part of a linked adapter
    // optional by default — without the \$ anchor, --discard-untrimmed
    // would still keep reads missing 1492R.
    // -e 0.1 = ~2 mismatches per primer (IUPAC codes don't count as errors).
    // Pass 2 trims residual polyA/polyG tails; a no-op for clean HiFi but
    // prevents tail variants from inflating ASVs when present. Counts
    // are post-pass-2 (pass 2 trims, does not drop).
    def iupac_comp = [
        'A':'T','T':'A','G':'C','C':'G',
        'R':'Y','Y':'R','M':'K','K':'M','S':'S','W':'W',
        'B':'V','V':'B','D':'H','H':'D','N':'N'
    ]
    def r_primer_rc = params.r_primer.toUpperCase().reverse().collect { iupac_comp[it] ?: it }.join('')
    """
    cutadapt \\
        --quiet \\
        --cores ${task.cpus} \\
        --revcomp \\
        --discard-untrimmed \\
        -e 0.1 \\
        -g "^${params.f_primer}...${r_primer_rc}\$" \\
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
