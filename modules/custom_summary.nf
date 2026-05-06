/*
========================================================================================
    CUSTOM SUMMARY MODULES
========================================================================================
    Builds a small per-run QC report from raw FastQC outputs.

    Three processes because no single biocontainer ships both Python and
    BLAST+:
      1. CUSTOM_SUMMARY_PARSE  — pandas image, parses FastQC zips into
                                  JSON + a FASTA of overrepresented seqs
      2. CUSTOM_SUMMARY_BLAST  — BLAST image, queries the local
                                  16S_ribosomal_RNA db one seq at a time
                                  (best-effort: failures land in a
                                  sentinel file so the pipeline keeps going)
      3. CUSTOM_SUMMARY_RENDER — pandas image, emits the HTML summary
----------------------------------------------------------------------------------------
*/

process CUSTOM_SUMMARY_PARSE {
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/custom_summary", mode: 'copy'

    input:
    path fastqc_zips

    output:
    path "parsed.json",    emit: parsed
    path "top_seqs.fasta", emit: top_seqs

    script:
    """
    mkdir -p fastqc_input
    for f in ${fastqc_zips}; do
        cp "\$f" fastqc_input/
    done

    python3 ${projectDir}/bin/build_custom_summary.py parse \\
        --fastqc-dir fastqc_input \\
        --out-json parsed.json \\
        --out-fasta top_seqs.fasta \\
        --top-overreps ${params.custom_summary_top_overreps}
    """
}

process CUSTOM_SUMMARY_BLAST {
    container 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'
    containerOptions "--bind ${params.blast_db_dir}"
    publishDir "${params.outdir}/custom_summary", mode: 'copy'
    // Best-effort: never sink the run if remote BLAST flakes.
    errorStrategy 'ignore'

    input:
    path top_seqs

    output:
    path "blast_hits.tsv", emit: hits

    script:
    """
    export BLASTDB="${params.blast_db_dir}"

    if [ ! -s ${top_seqs} ]; then
        echo "BLAST_SKIPPED: no overrepresented sequences" > blast_hits.tsv
        exit 0
    fi

    # BLAST one sequence at a time with a per-query timeout — a slow or
    # failing query can't sink the whole step, and any hits already
    # written are preserved.
    touch blast_hits.tsv

    while IFS= read -r header || [ -n "\$header" ]; do
        IFS= read -r seq || true
        qid=\${header#>}

        tmpfasta=\${qid}.fa
        printf '>%s\\n%s\\n' "\$qid" "\$seq" > "\$tmpfasta"

        set +e
        timeout ${params.custom_summary_blast_timeout} blastn \\
            -query "\$tmpfasta" \\
            -db ${params.custom_summary_blast_db} \\
            -outfmt '6 qseqid ssaccver stitle bitscore' \\
            -out "\${qid}_raw.tsv" 2>> blast.err
        rc=\$?
        set -e

        rm -f "\$tmpfasta"

        if [ \$rc -ne 0 ]; then
            echo "[WARN] blastn for \$qid exited with status \$rc — skipping" >&2
            continue
        fi

        # Keep best hit (highest bitscore) per query.
        awk -v qid="\$qid" '
        {
            bs=\$NF; acc=\$2
            stitle=""
            for(i=3;i<NF;i++) stitle = stitle (i==3?"":OFS) \$i
            if (!(qid in best) || bs > best[qid]) {
                best[qid] = bs
                line[qid] = qid "\\t" acc " " stitle
            }
        }
        END { if (qid in line) print line[qid] }
        ' "\${qid}_raw.tsv" >> blast_hits.tsv

        rm -f "\${qid}_raw.tsv"
    done < ${top_seqs}

    if [ ! -s blast_hits.tsv ]; then
        echo "BLAST_SKIPPED: all sequences failed or produced no hits" > blast_hits.tsv
    fi
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
    def attempted_flag = params.custom_summary_blast ? "--blast-attempted" : "--no-blast-attempted"
    """
    python3 ${projectDir}/bin/build_custom_summary.py render \\
        --in-json ${parsed} \\
        --blast-tsv ${blast_hits} \\
        --output-html custom_summary.html \\
        --output-txt raw_custom_summary.txt \\
        --threshold ${params.min_reads_threshold} \\
        --blast-db ${params.custom_summary_blast_db} \\
        ${attempted_flag}
    """
}
