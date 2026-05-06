process CUSTOM_SUMMARY_PARSE {
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/custom_summary", mode: 'copy'

    input:
    path fastqc_zips

    output:
    path "parsed.json",   emit: parsed
    path "top_seqs.fasta", emit: top_seqs

    script:
    """
    # [stub] parse fastqc zips for overrepresented sequences
    touch parsed.json top_seqs.fasta
    """
}

process CUSTOM_SUMMARY_BLAST {
    container 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'
    publishDir "${params.outdir}/custom_summary", mode: 'copy'

    input:
    path top_seqs

    output:
    path "blast_hits.tsv", emit: hits

    script:
    """
    # [stub] blast top_seqs against custom_summary_blast_db
    touch blast_hits.tsv
    """
}

process CUSTOM_SUMMARY_RENDER {
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/custom_summary", mode: 'copy'

    input:
    tuple path(parsed), path(blast_hits)

    output:
    path "custom_summary.html", emit: report

    script:
    """
    # [stub] render combined summary html
    touch custom_summary.html
    """
}
