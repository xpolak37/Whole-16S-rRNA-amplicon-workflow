process HOST_REMOVAL {
    tag "${meta.id}"
    container 'quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0'
    cpus 4
    publishDir "${params.outdir}/hostile", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.clean.fastq.gz"),  emit: reads
    path  "${meta.id}.host_removal.counts.tsv", emit: counts

    script:
    """
    # hostile's aligner module mkdir's \$XDG_DATA_HOME/hostile at import time.
    # Singularity --no-home leaves \$HOME pointing at a read-only path, so
    # redirect HOME to the writable Nextflow work dir.
    export HOME=\$PWD

    hostile clean \\
        --fastq1 ${reads} \\
        --index ${params.hostile_index_dir}/human-t2t-hla-argos985-mycob140.mmi \\
        --output . \\
        --threads ${task.cpus}

    out=\$(ls *.clean.fastq.gz | head -n1)
    count=\$(zcat "\$out" | awk 'END{print NR/4}')
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.host_removal.counts.tsv
    """
}

process PHIX_REMOVAL {
    tag "${meta.id}"
    container 'quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0'
    cpus 4
    publishDir "${params.outdir}/hostile", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_cleaned.fastq.gz"), emit: reads
    path  "${meta.id}.phix_removal.counts.tsv",  emit: counts

    script:
    """
    export HOME=\$PWD

    hostile clean \\
        --fastq1 ${reads} \\
        --index ${params.phix_fasta} \\
        --output . \\
        --threads ${task.cpus}

    # hostile appends .clean → reads is *.clean.fastq.gz, output is *.clean.clean.fastq.gz.
    # Rename to the bash pipeline's *_trimmed_cleaned.fastq.gz convention.
    raw=\$(ls *.clean.clean.fastq.gz | head -n1)
    final="${meta.id}_trimmed_cleaned.fastq.gz"
    mv "\$raw" "\$final"

    count=\$(zcat "\$final" | awk 'END{print NR/4}')
    printf 'sample\\tcount\\n%s\\t%s\\n' "${meta.id}" "\$count" > ${meta.id}.phix_removal.counts.tsv
    """
}
