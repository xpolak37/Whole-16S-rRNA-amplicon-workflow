process QIIME_NAIVE_BAYES {
    tag "${denoiser}/qnb"
    container 'quay.io/qiime2/amplicon:2026.1'
    cpus 4
    publishDir path: { "${params.outdir}/taxonomy/${denoiser}/qnb" }, mode: 'copy'

    input:
    tuple val(denoiser), path(asv_table), path(asv_fasta)

    output:
    tuple val(denoiser), val('qnb'), path('taxa_table_qnb.tsv'), emit: taxa

    script:
    """
    export NUMBA_CACHE_DIR=\${PWD}/numba_cache
    export TMPDIR=\${PWD}/tmp
    mkdir -p \${NUMBA_CACHE_DIR} \${TMPDIR}

    qiime tools import \\
        --input-path ${asv_fasta} \\
        --output-path rep-seqs.qza \\
        --type 'FeatureData[Sequence]'

    qiime feature-classifier classify-sklearn \\
        --i-reads rep-seqs.qza \\
        --i-classifier ${params.classifiers_dir}/qnb_classifier.qza \\
        --p-n-jobs ${task.cpus} \\
        --p-confidence ${params.qiime_naive_bayes_confidence} \\
        --o-classification taxonomy.qza

    qiime tools export --input-path taxonomy.qza --output-path .

    python3 - <<'EOF'
import csv
with open("taxonomy.tsv", newline='') as fin, \\
     open("taxa_table_qnb.tsv", "w", newline='') as fout:
    reader = csv.reader(fin, delimiter="\\t")
    writer = csv.writer(fout, delimiter="\\t")
    header = next(reader)
    header[0], header[1], header[2] = "SeqID", "Taxonomy", "Confidence"
    writer.writerow(header)
    for row in reader:
        writer.writerow(row)
EOF
    """
}

process QIIME_BLAST {
    tag "${denoiser}/qblast"
    container 'quay.io/qiime2/amplicon:2026.1'
    cpus 4
    publishDir path: { "${params.outdir}/taxonomy/${denoiser}/qblast" }, mode: 'copy'

    input:
    tuple val(denoiser), path(asv_table), path(asv_fasta)

    output:
    tuple val(denoiser), val('qblast'), path('taxa_table_qblast.tsv'), emit: taxa

    script:
    """
    export NUMBA_CACHE_DIR=\${PWD}/numba_cache
    export TMPDIR=\${PWD}/tmp
    mkdir -p \${NUMBA_CACHE_DIR} \${TMPDIR}

    qiime tools import \\
        --input-path ${asv_fasta} \\
        --output-path rep-seqs.qza \\
        --type 'FeatureData[Sequence]'

    qiime feature-classifier classify-consensus-blast \\
        --i-query rep-seqs.qza \\
        --i-reference-reads ${params.classifiers_dir}/qblast_seqs.qza \\
        --i-reference-taxonomy ${params.classifiers_dir}/qblast_tax.qza \\
        --p-num-threads ${task.cpus} \\
        --p-perc-identity ${params.blast_percidentity} \\
        --p-strand both \\
        --p-min-consensus ${params.blast_minconsensus} \\
        --o-classification taxonomy.qza \\
        --o-search-results blast_results.qza

    qiime tools export --input-path taxonomy.qza --output-path .

    python3 - <<'EOF'
import csv
with open("taxonomy.tsv", newline='') as fin, \\
     open("taxa_table_qblast.tsv", "w", newline='') as fout:
    reader = csv.reader(fin, delimiter="\\t")
    writer = csv.writer(fout, delimiter="\\t")
    header = next(reader)
    header[0], header[1], header[2] = "SeqID", "Taxonomy", "Confidence"
    writer.writerow(header)
    for row in reader:
        writer.writerow(row)
EOF
    """
}

process IDTAXA {
    tag "${denoiser}/idtaxa"
    container 'quay.io/biocontainers/bioconductor-decipher:3.6.0--r45h01b2380_0'
    cpus 4
    publishDir path: { "${params.outdir}/taxonomy/${denoiser}/idtaxa" }, mode: 'copy'

    input:
    tuple val(denoiser), path(asv_table), path(asv_fasta)

    output:
    tuple val(denoiser), val('idtaxa'),
          path('taxa_table_idtaxa.tsv'),
          path('taxa_table_idtaxa_conf.tsv'), emit: taxa

    script:
    """
    Rscript ${projectDir}/bin/idtaxa.R \\
        ${asv_fasta} \\
        ${params.classifiers_dir}/idtaxa.RData \\
        ${task.cpus} \\
        taxa_table_idtaxa.tsv \\
        taxa_table_idtaxa_conf.tsv
    """
}

process ASSIGNTAXONOMY {
    tag "${denoiser}/assigntaxonomy"
    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'
    cpus 4
    publishDir path: { "${params.outdir}/taxonomy/${denoiser}/assigntaxonomy" }, mode: 'copy'

    input:
    tuple val(denoiser), path(asv_table), path(asv_fasta)

    output:
    tuple val(denoiser), val('assigntaxonomy'), path('taxa_table_assigntaxonomy.tsv'), emit: taxa

    script:
    // bin/assigntaxonomy.R is the sibling script; it writes
    // ${denoising_tool}_taxa_table.tsv. Rename to the scaffold's contract.
    if (params.classifiers_dir) {
        """
        Rscript ${projectDir}/bin/assigntaxonomy.R \\
            ${asv_fasta} \\
            ${params.classifiers_dir}/silva_assigntaxonomy.fa.gz \\
            ${task.cpus} \\
            ${denoiser}

        mv ${denoiser}_taxa_table.tsv taxa_table_assigntaxonomy.tsv
        """
    } else {
        """
        touch taxa_table_assigntaxonomy.tsv
        """
    }
}
