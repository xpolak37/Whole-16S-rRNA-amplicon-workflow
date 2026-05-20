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
include { BAM2FASTQ }                 from './modules/bam2fastq'
include { CUTADAPT }                  from './modules/cutadapt'
include { HOST_REMOVAL; PHIX_REMOVAL } from './modules/hostile'
include { VSEARCH_ORIENT }            from './modules/orient'
include { DADA2_PACBIO; DADA2_PACBIO_NODENOISE } from './modules/dada2'
include { QIIME_NAIVE_BAYES; QIIME_BLAST; IDTAXA; ASSIGNTAXONOMY } from './modules/tax_classifiers'
include { METASTANDARD }              from './modules/MetaStandard16S'
include { METASTANDARD_PLOTS }        from './modules/metastandard_plots'
include { MOCK_EVALUATION }           from './modules/mock_evaluation'
include { SEQTK_SUBSAMPLE }           from './modules/seqtk_subsample'
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
                           The 'fastq' column accepts .fastq[.gz], .fq[.gz],
                           or per-sample .bam (PacBio HiFi). BAMs are
                           auto-converted via pbtk bam2fastq.

    Mode 2 (non-demuxed BAM):
      --bam                PacBio HiFi BAM
      --barcodes           Barcodes FASTA for lima
      --lima_extra_args    Extra args passed to lima (default: '')

    Common:
      --outdir             Output directory                            (default: ./results)
      --denoiser           Comma list: dada2,dada2_nodenoise           (default: dada2)
      --classifiers        Comma list: qnb,qblast,idtaxa,assigntaxonomy (default: qnb,qblast,assigntaxonomy — idtaxa is opt-in, genus-only)
      --all                Run every denoiser × classifier combination
      --quick              Subsample reads to --quick_depth before FastQC (smoke test)
      --quick_depth        Reads per sample under --quick               (default: ${params.quick_depth})
      --skip_phix          Skip the PhiX depletion step                  (default: ${params.skip_phix})

    MetaStandard (cross-run unification + plots):
      --run_id             Run label baked into output filenames       (default: ${params.run_id})
      --tax_level          domain|phylum|...|species|asv               (default: ${params.tax_level})
      --metastandard_top_n Top-N taxa for stacked barplot              (default: ${params.metastandard_top_n})

    Mock community evaluation (off unless mock samples are present):
      --mock_evaluation    Enable mock evaluation                       (default: ${params.mock_evaluation})
      --mock_pattern       Regex matching mock sample IDs               (default: ${params.mock_pattern})
      --mock_top_n         Top-N genera in mock barplot                 (default: ${params.mock_top_n})
      --mock_abundance     Reference abundance CSV                      (default: bundled)
      --mock_taxa          Reference taxonomy CSV                       (default: bundled)
      --mock_synonyms      Genus synonym CSV                            (default: bundled)

    References (required at runtime by the relevant stages):
      --hostile_index_dir  hostile human index directory
      --phix_index         PhiX174 minimap2 index (.mmi)
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
        // Samplesheet 'fastq' column accepts .fastq[.gz], .fq[.gz] or .bam
        // (per-sample PacBio HiFi). BAMs are routed through BAM2FASTQ; the
        // rest pass straight through and the two streams are mixed back into
        // a single ch_reads of (meta, fastq.gz) tuples.
        ch_input_rows = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                if (!row.sample || !row.fastq) error "Samplesheet row missing 'sample' or 'fastq': ${row}"
                def sample_id = row.sample.trim()
                if (sample_id ==~ /.*\s.*/) error "Sample id contains whitespace: '${sample_id}'"
                def reads = file(row.fastq.trim(), checkIfExists: true)
                tuple([id: sample_id], reads)
            }
            .toList()
            .map { rows ->
                def ids  = rows.collect { it[0].id }
                def dups = ids.findAll { id -> ids.count(id) > 1 }.unique()
                if (dups) error "Duplicate sample ids in samplesheet: ${dups}"
                rows
            }
            .flatMap { it }

        ch_input_rows.branch {
            bam:   it[1].name.endsWith('.bam')
            fastq: true
        }.set { ch_branched }

        BAM2FASTQ(ch_branched.bam)
        ch_reads = BAM2FASTQ.out.reads.mix(ch_branched.fastq)
    } else {
        ch_bam      = Channel.fromPath(params.bam,      checkIfExists: true)
        ch_barcodes = Channel.fromPath(params.barcodes, checkIfExists: true)
        ch_reads    = LIMA_DEMUX(ch_bam.combine(ch_barcodes)).fastqs
            .flatten()
            .map { fq -> tuple([id: fq.simpleName], fq) }
    }

    // ============================================================
    // Optional --quick subsample (smoke-test mode; subsamples before FastQC)
    // ============================================================
    if (params.quick) {
        SEQTK_SUBSAMPLE(ch_reads)
        ch_reads = SEQTK_SUBSAMPLE.out.reads
    }

    // ============================================================
    // Per-sample QC + preprocessing
    // ============================================================
    FASTQC_RAW(ch_reads)
    CUTADAPT(ch_reads)
    FASTQC_TRIMMED(CUTADAPT.out.reads)
    HOST_REMOVAL(CUTADAPT.out.reads)

    // PhiX pass is optional: default on for safety, --skip_phix to bypass.
    if (params.skip_phix) {
        ch_after_phix  = HOST_REMOVAL.out.reads
        ch_phix_counts = Channel.empty()
    } else {
        PHIX_REMOVAL(HOST_REMOVAL.out.reads)
        ch_after_phix  = PHIX_REMOVAL.out.reads
        ch_phix_counts = PHIX_REMOVAL.out.counts
    }
    VSEARCH_ORIENT(ch_after_phix)

    // ============================================================
    // Cohort-level: collect oriented reads, fan to denoiser axis
    // ============================================================
    ch_oriented_pool = VSEARCH_ORIENT.out.reads.map { meta, fq -> fq }.collect()

    ch_asv             = Channel.empty()
    ch_denoiser_counts = Channel.empty()
    if ('dada2' in denoisers) {
        DADA2_PACBIO(ch_oriented_pool)
        ch_asv             = ch_asv.mix(DADA2_PACBIO.out.asv)
        ch_denoiser_counts = ch_denoiser_counts.mix(DADA2_PACBIO.out.counts)
    }
    if ('dada2_nodenoise' in denoisers) {
        DADA2_PACBIO_NODENOISE(ch_oriented_pool)
        ch_asv             = ch_asv.mix(DADA2_PACBIO_NODENOISE.out.asv)
        ch_denoiser_counts = ch_denoiser_counts.mix(DADA2_PACBIO_NODENOISE.out.counts)
    }

    // ============================================================
    // Classifier axis (Cartesian: each enabled classifier consumes ch_asv)
    // ============================================================
    ch_taxa = Channel.empty()
    if ('qnb' in classifiers) {
        QIIME_NAIVE_BAYES(ch_asv)
        ch_taxa = ch_taxa.mix(QIIME_NAIVE_BAYES.out.taxa)
    }
    if ('qblast' in classifiers) {
        QIIME_BLAST(ch_asv)
        ch_taxa = ch_taxa.mix(QIIME_BLAST.out.taxa)
    }
    if ('idtaxa' in classifiers) {
        IDTAXA(ch_asv)
        // IDTAXA emits an extra confidence-table path; drop it for MetaStandard
        ch_taxa = ch_taxa.mix(IDTAXA.out.taxa.map { d, c, t, conf -> tuple(d, c, t) })
    }
    if ('assigntaxonomy' in classifiers) {
        ASSIGNTAXONOMY(ch_asv)
        ch_taxa = ch_taxa.mix(ASSIGNTAXONOMY.out.taxa)
    }

    // ============================================================
    // MetaStandard: per (denoiser × classifier) unify ASV + taxa,
    // then plot. Mock evaluation runs against each MetaStandard TSV
    // when --mock_evaluation is set.
    // ============================================================
    ch_asv_keyed = ch_asv.map { d, asv, fa -> tuple(d, asv) }
    ch_meta_in   = ch_taxa.combine(ch_asv_keyed, by: 0)
        .map { d, c, t, asv -> tuple(d, c, t, asv) }
    METASTANDARD(ch_meta_in)
    METASTANDARD_PLOTS(METASTANDARD.out.tsv)

    if (params.mock_evaluation) {
        ch_mock_abundance = Channel.value(file(params.mock_abundance, checkIfExists: true))
        ch_mock_taxa      = Channel.value(file(params.mock_taxa,      checkIfExists: true))
        ch_mock_synonyms  = Channel.value(file(params.mock_synonyms,  checkIfExists: true))
        MOCK_EVALUATION(METASTANDARD.out.tsv,
                        ch_mock_abundance, ch_mock_taxa, ch_mock_synonyms)
    }

    // ============================================================
    // Counts collation, MultiQC, custom summary
    // ============================================================
    ch_all_counts = CUTADAPT.out.counts
        .mix(HOST_REMOVAL.out.counts)
        .mix(ch_phix_counts)
        .mix(VSEARCH_ORIENT.out.counts)
        .mix(ch_denoiser_counts)
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
