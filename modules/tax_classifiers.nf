process QIIME_NAIVE_BAYES {
    tag "${denoiser}/qnb"
    container 'quay.io/qiime2/amplicon:2026.1'
    cpus 4

    input:
    tuple val(denoiser), path(asv_table), path(asv_fasta)

    output:
    tuple val(denoiser), val('qnb'), path('taxa_table_qnb.tsv'), emit: taxa

    script:
    def classifier = 'qnb'
    """
    # [stub] qiime feature-classifier classify-sklearn ...
    touch taxa_table_qnb.tsv
    """
}

process QIIME_BLAST {
    tag "${denoiser}/qblast"
    container 'quay.io/qiime2/amplicon:2026.1'
    cpus 4

    input:
    tuple val(denoiser), path(asv_table), path(asv_fasta)

    output:
    tuple val(denoiser), val('qblast'), path('taxa_table_qblast.tsv'), emit: taxa

    script:
    def classifier = 'qblast'
    """
    touch taxa_table_qblast.tsv
    """
}

process IDTAXA {
    tag "${denoiser}/idtaxa"
    container 'quay.io/biocontainers/bioconductor-decipher:3.6.0--r45h01b2380_0'
    cpus 4

    input:
    tuple val(denoiser), path(asv_table), path(asv_fasta)

    output:
    tuple val(denoiser), val('idtaxa'), path('taxa_table_idtaxa.tsv'), path('taxa_table_idtaxa_conf.tsv'), emit: taxa

    script:
    def classifier = 'idtaxa'
    """
    Rscript ${projectDir}/bin/idtaxa.R ${asv_fasta} ${params.classifiers_dir}/idtaxa.RData taxa_table_idtaxa.tsv taxa_table_idtaxa_conf.tsv
    """
}

process ASSIGNTAXONOMY {
    tag "${denoiser}/assigntaxonomy"
    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'
    cpus 4

    input:
    tuple val(denoiser), path(asv_table), path(asv_fasta)

    output:
    tuple val(denoiser), val('assigntaxonomy'), path('taxa_table_assigntaxonomy.tsv'), emit: taxa

    script:
    def classifier = 'assigntaxonomy'
    // bin/assigntaxonomy.R is real code from the sibling repo; it requires a
    // SILVA reference at classifiers_dir/silva_assigntaxonomy.fa.gz. Touch the
    // output instead when no reference dir is provided so scaffold validation
    // (and any --classifiers_dir-less run) still passes channel topology.
    if (params.classifiers_dir) {
        """
        Rscript ${projectDir}/bin/assigntaxonomy.R ${asv_fasta} ${params.classifiers_dir}/silva_assigntaxonomy.fa.gz taxa_table_assigntaxonomy.tsv
        """
    } else {
        """
        touch taxa_table_assigntaxonomy.tsv
        """
    }
}
