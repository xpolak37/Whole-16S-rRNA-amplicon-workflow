#!/usr/bin/env bash
# TODO: real implementation deferred to its own spec.
# This pipeline expects the following reference artifacts to exist and have
# their paths supplied as Nextflow params (see nextflow.config):
#
#   hostile_index_dir   — hostile human-t2t-hla-argos985-mycob140 index directory
#   phix_fasta          — PhiX174 reference fasta for hostile second pass
#   silva_orient_db     — primer-anchored SILVA fasta for vsearch --orient (27F-1492R)
#   classifiers_dir     — directory containing:
#                            * QIIME naive-bayes classifier .qza
#                            * QIIME BLAST classifier .qza pair (sequences + taxonomy)
#                            * IDTAXA training set .RData
#                            * AssignTaxonomy SILVA reference fasta (.fa.gz)
#
# Once implemented, this script will download or build each of the above into
# the user-specified cache directories. For now it is a help stub only.

set -euo pipefail

print_help() {
    sed -n '2,17p' "$0"
    echo ""
    echo "Usage: setup_pipeline.sh [--help]"
    echo ""
    echo "(stub — real implementation pending)"
}

case "${1:-}" in
    -h|--help|"") print_help; exit 0 ;;
    *) echo "Unknown argument: $1"; print_help; exit 2 ;;
esac
