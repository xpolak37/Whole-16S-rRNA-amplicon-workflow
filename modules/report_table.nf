/*
========================================================================================
    REPORT TABLE MODULE
========================================================================================
    Converts ONE canonical classifier's long output
        SeqID | Taxonomy(prefixed lineage) | Confidence
    into the WIDE, report-ready taxonomy table the downstream R reports consume
        SeqID | domain | phylum | class | order | family | genus | species

    This is ADDITIVE: the per-classifier long tables under
    taxonomy/<denoiser>/<classifier>/ are left untouched (they remain the
    provenance record and keep Confidence). The wide table is published to
    taxonomy/<denoiser>/report/taxa_table.tsv and keyed by the same ASV_<i>
    IDs as ASV_table.tsv, so it joins straight onto the ASV table.

    Which classifier feeds this is set by --report_classifier (default qblast).
----------------------------------------------------------------------------------------
*/

process REPORT_TABLE {
    tag "${denoiser}/${classifier}"
    publishDir path: { "${params.outdir}/taxonomy/${denoiser}/report" }, mode: 'copy'

    input:
    tuple val(denoiser), val(classifier), path(taxa_long)

    output:
    tuple val(denoiser), path('taxa_table.tsv'), emit: report

    script:
    """
    python3 ${projectDir}/bin/build_report_table.py \\
        --in ${taxa_long} \\
        --out taxa_table.tsv
    """
}
