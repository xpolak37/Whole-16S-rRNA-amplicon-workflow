process COLLATE_COUNTS {
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    path counts_files

    output:
    path 'read_counts_summary.tsv', emit: summary

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    script:
    """
    collate_counts.py read_counts_summary.tsv ${counts_files.join(' ')}
    """
}
