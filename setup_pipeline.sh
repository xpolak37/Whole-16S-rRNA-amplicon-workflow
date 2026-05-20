#!/usr/bin/env bash
#===============================================================================
# 16S PacBio HiFi PROFILING PIPELINE — SETUP SCRIPT
#===============================================================================
# Downloads/builds every reference artifact the Nextflow pipeline needs:
#   - Singularity container images for every process
#   - PhiX174 fasta for hostile second-pass decontamination
#   - hostile minimap2 index of human-t2t-hla-argos985-mycob140
#   - DADA2-format SILVA 138.2 training set (assignTaxonomy)
#   - QIIME naive-bayes classifier + dereplicated seqs/tax (qblast)
#   - vsearch --orient SILVA fasta (extracted from qiime build)
#   - DECIPHER IDTAXA training set (.RData)
#   - NCBI 16S_ribosomal_RNA BLAST DB (optional, for custom_summary)
#
# Usage:
#   ./setup_pipeline.sh [--rewrite] [INSTALL_DIR]
#
# Flags:
#   --rewrite, -r   Re-download/rebuild every artifact even if it already exists.
#                   Default behaviour is to skip steps whose output is present;
#                   use --rewrite for testing or to refresh stale resources.
#
# The script is idempotent: each step skips if its output already exists, so
# you can Ctrl-C and re-run safely. Total first-run time is dominated by the
# QIIME classifier build (~1-2h, downloads ~5GB of SILVA), which the user can
# bypass by placing pre-built classifiers under ${INSTALL_DIR}/classifiers/
# before running.
#
# After completion, run the pipeline with:
#   nextflow run main.nf -c ${INSTALL_DIR}/pipeline_paths.config \
#       --input <samplesheet.csv> --outdir <results>
#===============================================================================

set -euo pipefail

# Color codes
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

START_TIME=$(date +%s)

#-------------------------------------------------------------------------------
# Logging
#-------------------------------------------------------------------------------
log_info()    { echo -e "${GREEN}[INFO]${NC}    $(date '+%H:%M:%S') - $1" | tee -a "${LOGFILE:-/dev/null}"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC}    $(date '+%H:%M:%S') - $1" | tee -a "${LOGFILE:-/dev/null}"; }
log_error()   { echo -e "${RED}[ERROR]${NC}   $(date '+%H:%M:%S') - $1" | tee -a "${LOGFILE:-/dev/null}"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $(date '+%H:%M:%S') - $1" | tee -a "${LOGFILE:-/dev/null}"; }

# Remove the listed paths when --rewrite is in effect, so the existence check
# in each step falls through to the build/fetch branch.
maybe_rm() {
    [ "${REWRITE:-0}" = "1" ] || return 0
    for p in "$@"; do
        if [ -e "$p" ] || [ -L "$p" ]; then
            log_info "  --rewrite: removing $p"
            rm -rf "$p"
        fi
    done
}

print_header() {
    echo -e "${BLUE}======================================================================${NC}"
    echo -e "${BLUE}  16S PacBio HiFi PIPELINE — SETUP${NC}"
    echo -e "${BLUE}======================================================================${NC}"
    echo ""
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "$1 is not installed or not in PATH"
        exit 1
    fi
}

#-------------------------------------------------------------------------------
# Install directory
#-------------------------------------------------------------------------------
print_header

REWRITE=0
POSITIONAL_ARGS=()
for arg in "$@"; do
    case "$arg" in
        -h|--help)
            sed -n '2,32p' "$0"
            exit 0
            ;;
        --rewrite|-r)
            REWRITE=1
            ;;
        -*)
            log_error "Unknown flag: $arg (try --help)"
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$arg")
            ;;
    esac
done

if [ "${#POSITIONAL_ARGS[@]}" -eq 0 ]; then
    echo -e "${YELLOW}No installation directory provided.${NC}"
    echo "Enter full path where pipeline resources will live"
    echo "(creates: classifiers/ singularity_cache/ hostile_index/ phix/ silva_orient/ blast_db/ logs/)"
    read -rp "Installation directory: " INSTALL_DIR
    INSTALL_DIR=$(echo "$INSTALL_DIR" | xargs)
    [ -z "$INSTALL_DIR" ] && { log_error "No directory provided."; exit 1; }
else
    INSTALL_DIR="${POSITIONAL_ARGS[0]}"
fi

INSTALL_DIR=$(realpath -m "$INSTALL_DIR")
log_info "Installation directory: ${INSTALL_DIR}"

read -rp "Continue with this directory? (yes/no): " CONFIRM
[[ ! "$CONFIRM" =~ ^[Yy]([Ee][Ss])?$ ]] && { log_warn "Cancelled."; exit 0; }

#-------------------------------------------------------------------------------
# Directory layout
#-------------------------------------------------------------------------------
SING_DIR="${INSTALL_DIR}/singularity_cache"
CLASSIFIERS_DIR="${INSTALL_DIR}/classifiers"
HOSTILE_DIR="${INSTALL_DIR}/hostile_index"
PHIX_DIR="${INSTALL_DIR}/phix"
ORIENT_DIR="${INSTALL_DIR}/silva_orient"
BLAST_DB_DIR="${INSTALL_DIR}/blast_db"
LOG_DIR="${INSTALL_DIR}/logs"
# SIF -> sandbox extraction scratch. /tmp is often a small tmpfs that fills up
# during HOST_REMOVAL (the hostile image is the chunkiest), so we point
# Singularity at the install partition instead.
TMP_DIR="${INSTALL_DIR}/tmp"

mkdir -p "$SING_DIR" "$CLASSIFIERS_DIR" "$HOSTILE_DIR" "$PHIX_DIR" "$ORIENT_DIR" "$BLAST_DB_DIR" "$LOG_DIR" "$TMP_DIR"

LOGFILE="${LOG_DIR}/setup_$(date +%Y%m%d_%H%M%S).log"
touch "$LOGFILE"
log_info "Log file: ${LOGFILE}"

if [ "$REWRITE" = "1" ]; then
    log_warn "--rewrite enabled: existing artifacts will be removed and re-fetched"
fi

#-------------------------------------------------------------------------------
# Required host tools
#-------------------------------------------------------------------------------
log_info "Checking required tools..."
check_command wget
check_command tar
check_command singularity
log_success "Host tools OK"

#-------------------------------------------------------------------------------
# STEP 1 — Pull Singularity containers
#-------------------------------------------------------------------------------
echo ""
log_info "=== STEP 1 — Pulling Singularity containers ==="

# Filenames must match Nextflow's auto-generated cache name
# (docker URI with `://` and `/` replaced by `-`, plus `.img`),
# otherwise Nextflow re-pulls each container at runtime.
declare -A CONTAINERS=(
    ["quay.io-biocontainers-cutadapt-4.6--py39hf95cd2a_1.img"]="docker://quay.io/biocontainers/cutadapt:4.6--py39hf95cd2a_1"
    ["quay.io-biocontainers-hostile-1.1.0--pyhdfd78af_0.img"]="docker://quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0"
    ["quay.io-biocontainers-vsearch-2.27.0--h6a68c12_1.img"]="docker://quay.io/biocontainers/vsearch:2.27.0--h6a68c12_1"
    ["quay.io-biocontainers-bioconductor-dada2-1.38.0--r45ha27e39d_0.img"]="docker://quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"
    ["quay.io-biocontainers-bioconductor-decipher-3.6.0--r45h01b2380_0.img"]="docker://quay.io/biocontainers/bioconductor-decipher:3.6.0--r45h01b2380_0"
    ["quay.io-qiime2-amplicon-2026.1.img"]="docker://quay.io/qiime2/amplicon:2026.1"
    ["quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img"]="docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    ["quay.io-biocontainers-multiqc-1.21--pyhdfd78af_0.img"]="docker://quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
    ["quay.io-biocontainers-lima-2.12.0--h9ee0642_1.img"]="docker://quay.io/biocontainers/lima:2.12.0--h9ee0642_1"
    ["quay.io-biocontainers-pandas-2.2.1.img"]="docker://quay.io/biocontainers/pandas:2.2.1"
    ["quay.io-biocontainers-seaborn-0.13.2.img"]="docker://quay.io/biocontainers/seaborn:0.13.2"
    ["quay.io-biocontainers-seqtk-1.4--he4a0461_2.img"]="docker://quay.io/biocontainers/seqtk:1.4--he4a0461_2"
    ["quay.io-biocontainers-pbtk-3.5.0--h9ee0642_0.img"]="docker://quay.io/biocontainers/pbtk:3.5.0--h9ee0642_0"
    ["quay.io-biocontainers-blast-2.15.0--pl5321h6f7f691_1.img"]="docker://quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1"
    ["quay.io-biocontainers-entrez-direct-24.0--he881be0_0.img"]="docker://quay.io/biocontainers/entrez-direct:24.0--he881be0_0"
)

cd "$SING_DIR"
n=0; total=${#CONTAINERS[@]}
for img in "${!CONTAINERS[@]}"; do
    n=$((n+1))
    maybe_rm "$img"
    if [ -f "$img" ]; then
        log_warn  "[$n/$total] $img already present — skipping"
    else
        log_info  "[$n/$total] pulling $img"
        singularity pull --name "$img" "${CONTAINERS[$img]}" >> "$LOGFILE" 2>&1 \
            || { log_error "pull failed: $img"; exit 1; }
        log_success "[$n/$total] $img"
    fi
done

#-------------------------------------------------------------------------------
# STEP 2 — PhiX174 fasta (hostile second pass)
#-------------------------------------------------------------------------------
echo ""
log_info "=== STEP 2 — PhiX174 fasta ==="

PHIX_FASTA="${PHIX_DIR}/phiX174.fasta"
PHIX_MMI="${PHIX_DIR}/phiX174.mmi"

# When --rewrite is set, drop the fasta AND the .mmi; the .mmi must be rebuilt
# from the fresh fasta.
maybe_rm "$PHIX_FASTA" "$PHIX_MMI"

if [ -s "$PHIX_FASTA" ]; then
    log_warn "$PHIX_FASTA already present — skipping fetch"
else
    log_info "Fetching NC_001422.1 via entrez-direct..."
    singularity exec --bind "${PHIX_DIR}:${PHIX_DIR}" --pwd "${PHIX_DIR}" \
        "${SING_DIR}/quay.io-biocontainers-entrez-direct-24.0--he881be0_0.img" \
        bash -c "efetch -db nucleotide -id NC_001422.1 -format fasta > ${PHIX_FASTA}" \
        >> "$LOGFILE" 2>&1
    [ -s "$PHIX_FASTA" ] || { log_error "PhiX fetch failed"; exit 1; }
    log_success "PhiX174 fasta downloaded"
fi

# Pre-build the minimap2 index. Hostile in long-read mode requires .mmi —
# passing the fasta directly fails with "neither a valid custom index path
# nor a valid standard index name".
if [ -s "$PHIX_MMI" ]; then
    log_warn "$PHIX_MMI already present — skipping index build"
else
    log_info "Building PhiX minimap2 index..."
    singularity exec --bind "${PHIX_DIR}:${PHIX_DIR}" --pwd "${PHIX_DIR}" \
        "${SING_DIR}/quay.io-biocontainers-hostile-1.1.0--pyhdfd78af_0.img" \
        minimap2 -x map-ont -d "$PHIX_MMI" "$PHIX_FASTA" \
        >> "$LOGFILE" 2>&1 \
        || { log_error "PhiX minimap2 index build failed"; exit 1; }
    [ -s "$PHIX_MMI" ] || { log_error "PhiX .mmi missing after build"; exit 1; }
    log_success "PhiX minimap2 index ready"
fi

#-------------------------------------------------------------------------------
# STEP 3 — hostile minimap2 index of human-t2t-hla-argos985-mycob140
#-------------------------------------------------------------------------------
echo ""
log_info "=== STEP 3 — hostile minimap2 index ==="

HOSTILE_FA="${HOSTILE_DIR}/human-t2t-hla-argos985-mycob140.fa.gz"
HOSTILE_MMI="${HOSTILE_DIR}/human-t2t-hla-argos985-mycob140.mmi"

# When --rewrite is set, drop the fasta AND the .mmi; the .mmi must be rebuilt
# from the fresh fasta.
maybe_rm "$HOSTILE_FA" "$HOSTILE_MMI"

# 3a — fetch the source fasta (hostile fetch only downloads .fa.gz; it does not
# pre-build the .mmi).
if [ -s "$HOSTILE_FA" ]; then
    log_warn "$HOSTILE_FA already present — skipping fetch"
else
    log_info "Fetching hostile source fasta (~940 MB)..."
    # hostile streams via Python's tempfile.gettempdir() which defaults to /tmp
    # inside the container; /tmp is typically a small tmpfs, so redirect TMPDIR
    # to the install partition.
    mkdir -p "${HOSTILE_DIR}/tmp"
    singularity exec --bind "${HOSTILE_DIR}:${HOSTILE_DIR}" --pwd "${HOSTILE_DIR}" \
        --env HOSTILE_CACHE_DIR="${HOSTILE_DIR}" \
        --env TMPDIR="${HOSTILE_DIR}/tmp" \
        "${SING_DIR}/quay.io-biocontainers-hostile-1.1.0--pyhdfd78af_0.img" \
        hostile fetch --aligner minimap2 --name human-t2t-hla-argos985-mycob140 \
        >> "$LOGFILE" 2>&1 \
        || { log_error "hostile fetch failed"; exit 1; }
    [ -s "$HOSTILE_FA" ] || { log_error "hostile fasta missing after fetch"; exit 1; }
    log_success "hostile fasta downloaded"
fi

# 3b — pre-build the minimap2 .mmi once. Without this, every HOST_REMOVAL task
# would rebuild the index from .fa.gz (~5 min each).
if [ -s "$HOSTILE_MMI" ]; then
    log_warn "$HOSTILE_MMI already present — skipping index build"
else
    log_info "Pre-building minimap2 index (~5-10 min, output ~14 GB)..."
    singularity exec --bind "${HOSTILE_DIR}:${HOSTILE_DIR}" --pwd "${HOSTILE_DIR}" \
        --env TMPDIR="${HOSTILE_DIR}/tmp" \
        "${SING_DIR}/quay.io-biocontainers-hostile-1.1.0--pyhdfd78af_0.img" \
        minimap2 -x map-ont -d "$HOSTILE_MMI" "$HOSTILE_FA" \
        >> "$LOGFILE" 2>&1 \
        || { log_error "minimap2 index build failed"; exit 1; }
    [ -s "$HOSTILE_MMI" ] || { log_error "hostile .mmi missing after build"; exit 1; }
    log_success "hostile minimap2 index ready"
fi

#-------------------------------------------------------------------------------
# STEP 4 — AssignTaxonomy SILVA 138.2 training set
#-------------------------------------------------------------------------------
echo ""
log_info "=== STEP 4 — AssignTaxonomy SILVA training set ==="

ASSIGNTAX_FA="${CLASSIFIERS_DIR}/silva_assigntaxonomy.fa.gz"
ASSIGNTAX_URL="${ASSIGNTAX_URL:-https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz}"

maybe_rm "$ASSIGNTAX_FA"

if [ -s "$ASSIGNTAX_FA" ]; then
    log_warn "$ASSIGNTAX_FA already present — skipping"
else
    log_info "Downloading $ASSIGNTAX_URL ..."
    wget -O "$ASSIGNTAX_FA" "$ASSIGNTAX_URL" >> "$LOGFILE" 2>&1 \
        || { log_error "AssignTaxonomy download failed (try setting ASSIGNTAX_URL)"; exit 1; }
    log_success "AssignTaxonomy fasta ready"
fi

#-------------------------------------------------------------------------------
# STEP 5 — QIIME classifiers + orient DB (build via qiime2 container)
#-------------------------------------------------------------------------------
echo ""
log_info "=== STEP 5 — QIIME classifiers + vsearch orient DB ==="

QNB_QZA="${CLASSIFIERS_DIR}/qnb_classifier.qza"
QBLAST_SEQS_QZA="${CLASSIFIERS_DIR}/qblast_seqs.qza"
QBLAST_TAX_QZA="${CLASSIFIERS_DIR}/qblast_tax.qza"
ORIENT_FA="${ORIENT_DIR}/silva-27F-1492R-orient.fasta"
QBUILD_DIR="${CLASSIFIERS_DIR}/_build"

# Also wipe _build so its inner .qza idempotency checks don't short-circuit
# the rebuild.
maybe_rm "$QNB_QZA" "$QBLAST_SEQS_QZA" "$QBLAST_TAX_QZA" "$ORIENT_FA" "$QBUILD_DIR"

if [ -s "$QNB_QZA" ] && [ -s "$QBLAST_SEQS_QZA" ] && [ -s "$QBLAST_TAX_QZA" ] && [ -s "$ORIENT_FA" ]; then
    log_warn "QIIME classifiers + orient fasta already present — skipping"
else
    log_info "Running QIIME RESCRIPt build (downloads SILVA ~5 GB, runs ~1-2 h)..."

    mkdir -p "$QBUILD_DIR"

    cat > "${QBUILD_DIR}/build.sh" <<'QBUILD'
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export NUMBA_CACHE_DIR="$PWD/numba_cache"
export TMPDIR="$PWD/tmp"
mkdir -p "$NUMBA_CACHE_DIR" "$TMPDIR"

# Mirrors bin/extract_silva.sh, with SILVA 138.2.
if [ ! -s silva-138.2-ssu-nr99-rna-seqs.qza ]; then
    qiime rescript get-silva-data \
        --p-version '138.2' \
        --p-target 'SSURef_NR99' \
        --p-include-species-labels \
        --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
        --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza
fi

if [ ! -s silva-138.2-ssu-nr99-seqs.qza ]; then
    qiime rescript reverse-transcribe \
        --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
        --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza
fi

if [ ! -s silva-138.2-ssu-nr99-seqs-cleaned.qza ]; then
    qiime rescript cull-seqs \
        --i-sequences silva-138.2-ssu-nr99-seqs.qza \
        --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza
fi

if [ ! -s silva-138.2-ssu-nr99-seqs-filt.qza ]; then
    qiime rescript filter-seqs-length-by-taxon \
        --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
        --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
        --p-labels Archaea Bacteria Eukaryota \
        --p-min-lens 900 1000 1400 \
        --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
        --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza
fi

if [ ! -s silva-138.2-ssu-nr99-seqs-derep-super.qza ]; then
    qiime rescript dereplicate \
        --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza \
        --i-taxa silva-138.2-ssu-nr99-tax.qza \
        --p-rank-handles 'domain' 'phylum' 'class' 'order' 'family' 'genus' 'species' \
        --p-mode 'uniq' \
        --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-super.qza \
        --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-super.qza
fi

if [ ! -s silva-138.2-ssu-nr99-seqs-27F-1492R.qza ]; then
    qiime feature-classifier extract-reads \
        --i-sequences silva-138.2-ssu-nr99-seqs-derep-super.qza \
        --p-f-primer AGRGTTYGATYMTGGCTCAG \
        --p-r-primer RGYTACCTTGTTACGACTT \
        --p-n-jobs 2 \
        --p-read-orientation 'forward' \
        --o-reads silva-138.2-ssu-nr99-seqs-27F-1492R.qza
fi

if [ ! -s qblast_seqs.qza ] || [ ! -s qblast_tax.qza ]; then
    qiime rescript dereplicate \
        --i-sequences silva-138.2-ssu-nr99-seqs-27F-1492R.qza \
        --i-taxa silva-138.2-ssu-nr99-tax-derep-super.qza \
        --p-mode 'uniq' \
        --o-dereplicated-sequences qblast_seqs.qza \
        --o-dereplicated-taxa qblast_tax.qza
fi

if [ ! -s qnb_classifier.qza ]; then
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads qblast_seqs.qza \
        --i-reference-taxonomy qblast_tax.qza \
        --o-classifier qnb_classifier.qza
fi

# Export the dereplicated 27F-1492R seqs to fasta for vsearch --orient
if [ ! -s silva-27F-1492R-orient.fasta ]; then
    qiime tools export --input-path qblast_seqs.qza --output-path orient_export
    mv orient_export/dna-sequences.fasta silva-27F-1492R-orient.fasta
    rm -rf orient_export
fi
QBUILD
    chmod +x "${QBUILD_DIR}/build.sh"

    singularity exec --bind "${INSTALL_DIR}:${INSTALL_DIR}" --pwd "${QBUILD_DIR}" \
        "${SING_DIR}/quay.io-qiime2-amplicon-2026.1.img" bash "${QBUILD_DIR}/build.sh" \
        2>&1 | tee -a "$LOGFILE"

    cp "${QBUILD_DIR}/qnb_classifier.qza"        "$QNB_QZA"
    cp "${QBUILD_DIR}/qblast_seqs.qza"           "$QBLAST_SEQS_QZA"
    cp "${QBUILD_DIR}/qblast_tax.qza"            "$QBLAST_TAX_QZA"
    cp "${QBUILD_DIR}/silva-27F-1492R-orient.fasta" "$ORIENT_FA"
    log_success "QIIME classifiers + orient fasta ready"
fi

#-------------------------------------------------------------------------------
# STEP 6 — IDTAXA training set (DECIPHER pre-built, renamed for our R script)
#-------------------------------------------------------------------------------
echo ""
log_info "=== STEP 6 — IDTAXA training set (.RData) ==="

IDTAXA_RDATA="${CLASSIFIERS_DIR}/idtaxa.RData"
# DECIPHER retired the direct-download links on www2.decipher.codes in 2025
# and now hosts SILVA_SSU_r138.2_v2.RData on Google Drive (file id below).
# `confirm=t` bypasses the >100 MB virus-scan interstitial.
IDTAXA_GDRIVE_ID="${IDTAXA_GDRIVE_ID:-1w3wdSCpSihntWkbP_zvXz7r3s-tNB8DV}"
IDTAXA_URL="${IDTAXA_URL:-https://drive.usercontent.google.com/download?id=${IDTAXA_GDRIVE_ID}&export=download&confirm=t}"

maybe_rm "$IDTAXA_RDATA"

if [ -s "$IDTAXA_RDATA" ]; then
    log_warn "$IDTAXA_RDATA already present — skipping"
else
    raw="${CLASSIFIERS_DIR}/_idtaxa_raw.RData"
    log_info "Downloading $IDTAXA_URL ..."
    wget -O "$raw" "$IDTAXA_URL" >> "$LOGFILE" 2>&1 \
        || { log_error "IDTAXA RData download failed (try setting IDTAXA_URL)"; exit 1; }

    log_info "Renaming variable trainingSet -> trainingSet_custom for bin/idtaxa.R compatibility..."
    singularity exec --bind "${CLASSIFIERS_DIR}:${CLASSIFIERS_DIR}" --pwd "${CLASSIFIERS_DIR}" \
        "${SING_DIR}/quay.io-biocontainers-bioconductor-decipher-3.6.0--r45h01b2380_0.img" \
        Rscript -e "load('${raw}'); \
                    if (!exists('trainingSet_custom')) trainingSet_custom <- trainingSet; \
                    save(trainingSet_custom, file='${IDTAXA_RDATA}')" \
        >> "$LOGFILE" 2>&1 \
        || { log_error "IDTAXA RData rename failed"; exit 1; }
    rm -f "$raw"
    log_success "IDTAXA RData ready"
fi

#-------------------------------------------------------------------------------
# STEP 7 — NCBI 16S BLAST DB (optional)
#-------------------------------------------------------------------------------
echo ""
if [ "${SKIP_BLAST_DB:-0}" = "1" ]; then
    log_warn "=== STEP 7 — BLAST 16S DB skipped (SKIP_BLAST_DB=1) ==="
else
    log_info "=== STEP 7 — BLAST 16S DB (set SKIP_BLAST_DB=1 to skip) ==="
    if [ "$REWRITE" = "1" ]; then
        log_info "  --rewrite: removing existing BLAST DB files"
        rm -f "$BLAST_DB_DIR"/16S_ribosomal_RNA.*
    fi
    if ls "$BLAST_DB_DIR"/16S_ribosomal_RNA.n* 1>/dev/null 2>&1; then
        log_warn "BLAST DB already present — skipping"
    else
        mkdir -p "${BLAST_DB_DIR}/tmp"
        singularity exec --bind "${BLAST_DB_DIR}:${BLAST_DB_DIR}" --pwd "${BLAST_DB_DIR}" \
            --env TMPDIR="${BLAST_DB_DIR}/tmp" \
            "${SING_DIR}/quay.io-biocontainers-blast-2.15.0--pl5321h6f7f691_1.img" \
            update_blastdb.pl --decompress 16S_ribosomal_RNA \
            >> "$LOGFILE" 2>&1 \
            || log_warn "BLAST DB download failed (non-fatal)"
        log_success "BLAST DB ready (or skipped on failure)"
    fi
fi

#-------------------------------------------------------------------------------
# STEP 8 — Emit pipeline_paths.config
#-------------------------------------------------------------------------------
echo ""
log_info "=== STEP 8 — Writing pipeline_paths.config ==="

CONFIG_FILE="${INSTALL_DIR}/pipeline_paths.config"
cat > "$CONFIG_FILE" <<EOF
/*
 * Generated by setup_pipeline.sh on $(date)
 * Use:  nextflow run main.nf -c ${CONFIG_FILE} ...
 */
params {
    singularity_cache_dir   = '${SING_DIR}'
    hostile_index_dir       = '${HOSTILE_DIR}'
    phix_index              = '${PHIX_MMI}'
    silva_orient_db         = '${ORIENT_FA}'
    classifiers_dir         = '${CLASSIFIERS_DIR}'
    blast_db_dir            = '${BLAST_DB_DIR}'
    custom_summary_blast_db = '16S_ribosomal_RNA'
}
EOF
log_success "Config: $CONFIG_FILE"

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
END_TIME=$(date +%s)
DUR=$((END_TIME - START_TIME))
echo ""
echo -e "${GREEN}======================================================================${NC}"
echo -e "${GREEN}  SETUP COMPLETE — ${DUR}s${NC}"
echo -e "${GREEN}======================================================================${NC}"
echo ""
echo "Resources at: ${INSTALL_DIR}"
echo "Config file : ${CONFIG_FILE}"
echo ""
echo "Before each pipeline run, redirect Singularity's SIF-extraction scratch"
echo "to the install partition (host /tmp is usually a small tmpfs that fills"
echo "up during HOST_REMOVAL):"
echo "  export SINGULARITY_TMPDIR=${TMP_DIR}"
echo "  export APPTAINER_TMPDIR=${TMP_DIR}"
echo ""
echo "Smoke test (one denoiser × one classifier, no qiime/blast):"
echo "  nextflow run main.nf -c ${CONFIG_FILE} \\"
echo "      --input test/samplesheet.csv \\"
echo "      --denoiser dada2 \\"
echo "      --classifiers idtaxa,assigntaxonomy"
echo ""
