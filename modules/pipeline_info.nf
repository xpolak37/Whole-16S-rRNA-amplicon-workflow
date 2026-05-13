process COLLATE_COUNTS {

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
