#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    16S PacBio HiFi Profiling Pipeline (single-end)
========================================================================================
    [optional LIMA_DEMUX] -> FASTQC_RAW -> CUTADAPT -> FASTQC_TRIMMED ->
    HOST_REMOVAL -> PHIX_REMOVAL -> VSEARCH_ORIENT ->
    {DADA2_PACBIO | DADA2_PACBIO_NODENOISE} ->
    {QIIME_NAIVE_BAYES | QIIME_BLAST | IDTAXA | ASSIGNTAXONOMY}
========================================================================================
*/

include { FASTQC as FASTQC_RAW }      from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED }  from './modules/fastqc'
include { MULTIQC }                   from './modules/multiqc'
include { CUSTOM_SUMMARY_PARSE; CUSTOM_SUMMARY_BLAST; CUSTOM_SUMMARY_RENDER } from './modules/custom_summary'
include { LIMA_DEMUX }                from './modules/lima'
include { CUTADAPT }                  from './modules/cutadapt'
include { HOST_REMOVAL; PHIX_REMOVAL } from './modules/hostile'
include { VSEARCH_ORIENT }            from './modules/orient'
include { DADA2_PACBIO; DADA2_PACBIO_NODENOISE } from './modules/dada2'
include { QIIME_NAIVE_BAYES; QIIME_BLAST; IDTAXA; ASSIGNTAXONOMY } from './modules/tax_classifiers'
include { COLLATE_COUNTS }            from './modules/pipeline_info'

// ============================================================
// Tool axis constants & resolution
// ============================================================
def VALID_DENOISERS   = ['dada2', 'dada2_nodenoise']
def VALID_CLASSIFIERS = ['qnb', 'qblast', 'idtaxa', 'assigntaxonomy']

def resolveTools(String csv, List valid, String label) {
    if (!csv) error "Empty ${label} list. Valid: ${valid}"
    def chosen = csv.tokenize(',').collect { it.trim() } as LinkedHashSet
    def bad = chosen - valid
    if (bad) error "Unknown ${label} tool(s): ${bad}. Valid: ${valid}"
    return chosen.toList()
}

def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf --input <samplesheet.csv> --outdir <dir>
      nextflow run main.nf --bam <hifi.bam> --barcodes <bc.fa> --outdir <dir>

    Mode 1 (samplesheet):
      --input              CSV with columns: sample,fastq

    Mode 2 (non-demuxed BAM):
      --bam                PacBio HiFi BAM
      --barcodes           Barcodes FASTA for lima
      --lima_extra_args    Extra args passed to lima (default: '')

    Common:
      --outdir             Output directory                            (default: ./results)
      --denoiser           Comma list: dada2,dada2_nodenoise           (default: dada2)
      --classifiers        Comma list: qnb,qblast,idtaxa,assigntaxonomy (default: all four)
      --all                Run every denoiser × classifier combination

    References (required at runtime by the relevant stages):
      --hostile_index_dir  hostile human index directory
      --phix_fasta         PhiX174 FASTA
      --silva_orient_db    Primer-anchored SILVA DB (vsearch --orient)
      --classifiers_dir    Directory with classifier reference files

    Primers (HiFi 16S — do not change):
      --f_primer           ${params.f_primer}
      --r_primer           ${params.r_primer}

    DADA2 (PacBio-tuned — do not change):
      --minQ ${params.minQ}  --minLen ${params.minLen}  --maxLen ${params.maxLen}  --maxN ${params.maxN}  --maxEE ${params.maxEE}  --minAbundance ${params.minAbundance}

    Custom summary:
      --custom_summary               (default: ${params.custom_summary})
      --custom_summary_blast         (default: ${params.custom_summary_blast})
      --custom_summary_blast_db      (default: ${params.custom_summary_blast_db})
      --custom_summary_top_overreps  (default: ${params.custom_summary_top_overreps})
      --min_reads_threshold          (default: ${params.min_reads_threshold})

    Resource ceilings:
      --max_cpus ${params.max_cpus}   --max_memory ${params.max_memory}   --max_time ${params.max_time}
    """.stripIndent()
}

// ============================================================
// Workflow
// ============================================================
workflow {
    if (params.help) {
        helpMessage()
        return
    }

    // mode validation
    def has_input = params.input != null
    def has_bam   = params.bam   != null

    if (has_input && has_bam)   error "Specify either --input <samplesheet.csv> or --bam <hifi.bam>, not both."
    if (!has_input && !has_bam) error "Must supply --input <samplesheet.csv> or --bam <hifi.bam> + --barcodes <bc.fasta>."
    if (has_bam && params.barcodes == null) error "--bam requires --barcodes <bc.fasta>."

    def denoisers   = params.all ? VALID_DENOISERS   : resolveTools(params.denoiser,    VALID_DENOISERS,   'denoiser')
    def classifiers = params.all ? VALID_CLASSIFIERS : resolveTools(params.classifiers, VALID_CLASSIFIERS, 'classifier')

    log.info """
    ========================================================
    16S PacBio HiFi Profiling Pipeline
    --------------------------------------------------------
    Input mode    : ${has_input ? 'samplesheet (' + params.input + ')' : 'bam (' + params.bam + ')'}
    Output dir    : ${params.outdir}
    Denoisers     : ${denoisers.join(', ')}
    Classifiers   : ${classifiers.join(', ')}
    ========================================================
    """.stripIndent()

    // ============================================================
    // Reads channel (both modes converge to [meta, fq])
    // ============================================================
    if (has_input) {
        ch_reads = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                if (!row.sample || !row.fastq) error "Samplesheet row missing 'sample' or 'fastq': ${row}"
                def sample_id = row.sample.trim()
                if (sample_id ==~ /.*\s.*/) error "Sample id contains whitespace: '${sample_id}'"
                def fq = file(row.fastq.trim(), checkIfExists: true)
                tuple([id: sample_id], fq)
            }
            .toList()
            .map { rows ->
                def ids  = rows.collect { it[0].id }
                def dups = ids.findAll { id -> ids.count(id) > 1 }.unique()
                if (dups) error "Duplicate sample ids in samplesheet: ${dups}"
                rows
            }
            .flatMap { it }
    } else {
        ch_bam      = Channel.fromPath(params.bam,      checkIfExists: true)
        ch_barcodes = Channel.fromPath(params.barcodes, checkIfExists: true)
        ch_reads    = LIMA_DEMUX(ch_bam.combine(ch_barcodes)).fastqs
            .flatten()
            .map { fq -> tuple([id: fq.simpleName], fq) }
    }

    // ============================================================
    // Per-sample QC + preprocessing
    // ============================================================
    FASTQC_RAW(ch_reads)
    CUTADAPT(ch_reads)
    FASTQC_TRIMMED(CUTADAPT.out.reads)
    HOST_REMOVAL(CUTADAPT.out.reads)
    PHIX_REMOVAL(HOST_REMOVAL.out.reads)
    VSEARCH_ORIENT(PHIX_REMOVAL.out.reads)

    // ============================================================
    // Cohort-level: collect oriented reads, fan to denoiser axis
    // ============================================================
    ch_oriented_pool = VSEARCH_ORIENT.out.reads.map { meta, fq -> fq }.collect()

    ch_asv = Channel.empty()
    if ('dada2' in denoisers) {
        DADA2_PACBIO(ch_oriented_pool)
        ch_asv = ch_asv.mix(DADA2_PACBIO.out.asv)
    }
    if ('dada2_nodenoise' in denoisers) {
        DADA2_PACBIO_NODENOISE(ch_oriented_pool)
        ch_asv = ch_asv.mix(DADA2_PACBIO_NODENOISE.out.asv)
    }

    // ============================================================
    // Classifier axis (Cartesian: each enabled classifier consumes ch_asv)
    // ============================================================
    if ('qnb' in classifiers)            QIIME_NAIVE_BAYES(ch_asv)
    if ('qblast' in classifiers)         QIIME_BLAST(ch_asv)
    if ('idtaxa' in classifiers)         IDTAXA(ch_asv)
    if ('assigntaxonomy' in classifiers) ASSIGNTAXONOMY(ch_asv)

    // ============================================================
    // Counts collation, MultiQC, custom summary
    // ============================================================
    ch_all_counts = CUTADAPT.out.counts
        .mix(HOST_REMOVAL.out.counts)
        .mix(PHIX_REMOVAL.out.counts)
        .mix(VSEARCH_ORIENT.out.counts)
        .collect()
    COLLATE_COUNTS(ch_all_counts)

    ch_multiqc_input = FASTQC_RAW.out.zip
        .map { meta, zip -> zip }
        .mix(FASTQC_TRIMMED.out.zip.map { meta, zip -> zip })
        .collect()
    MULTIQC(ch_multiqc_input)

    if (params.custom_summary) {
        ch_summary = CUSTOM_SUMMARY_PARSE(FASTQC_RAW.out.zip.map { meta, zip -> zip }.collect())
        if (params.custom_summary_blast) {
            CUSTOM_SUMMARY_BLAST(ch_summary.top_seqs)
            CUSTOM_SUMMARY_RENDER(ch_summary.parsed.combine(CUSTOM_SUMMARY_BLAST.out.hits))
        } else {
            CUSTOM_SUMMARY_RENDER(ch_summary.parsed.map { tuple(it, file('NO_BLAST')) })
        }
    }
}
