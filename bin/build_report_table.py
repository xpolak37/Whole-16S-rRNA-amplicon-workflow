#!/usr/bin/env python3
"""
Build a report-ready WIDE taxonomy table from a classifier's long output.

The classifiers emit one row per ASV as:

    SeqID <tab> Taxonomy <tab> Confidence

where Taxonomy is a single prefixed, semicolon-joined lineage, e.g.

    d__Bacteria;p__Bacillota;...;g__Enterococcus;s__Enterococcus_faecalis

A no-hit row is the literal token "Unassigned"; a row may also stop early
(e.g. at genus, with no s__ field).

Downstream R reports (g354 / bile-acid family) expect a wide table with one
column per rank, plain values (no d__/s__ prefixes), and missing ranks filled
with the literal string "unassigned":

    SeqID <tab> domain <tab> phylum <tab> class <tab> order <tab> family <tab> genus <tab> species

Parsing is prefix-based (not positional), so a lineage that skips a rank still
lands each value in the right column. This mirrors analysis/build_phyloseq.R.

Confidence is intentionally dropped — the wide table is the analysis input; the
per-classifier long file remains the provenance record (and keeps Confidence).
"""

import argparse
import csv
import sys

# Ordered (prefix -> output column). Lowercase headers match what the existing
# g354 / bile-acid reports already read; do not reorder.
RANKS = [
    ("d__", "domain"),
    ("p__", "phylum"),
    ("c__", "class"),
    ("o__", "order"),
    ("f__", "family"),
    ("g__", "genus"),
    ("s__", "species"),
]
UNASSIGNED = "unassigned"


def split_lineage(taxonomy):
    """Map a prefixed semicolon lineage to one value per rank.

    Missing/blank ranks become UNASSIGNED. Unknown tokens (e.g. a bare
    "Unassigned") match no prefix and are ignored, leaving every rank
    UNASSIGNED.
    """
    parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
    out = []
    for prefix, _ in RANKS:
        value = UNASSIGNED
        for part in parts:
            if part.startswith(prefix):
                stripped = part[len(prefix):].strip()
                if stripped:
                    value = stripped
                break
        out.append(value)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile", required=True,
                    help="Long classifier table (SeqID, Taxonomy, Confidence)")
    ap.add_argument("--out", dest="outfile", required=True,
                    help="Wide report table to write")
    args = ap.parse_args()

    header_out = ["SeqID"] + [col for _, col in RANKS]

    n = 0
    with open(args.infile, newline="") as fin, \
         open(args.outfile, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")

        header_in = next(reader, None)
        if header_in is None:
            sys.exit(f"error: {args.infile} is empty")
        if header_in[:2] != ["SeqID", "Taxonomy"]:
            sys.exit(f"error: unexpected header {header_in!r}; "
                     "expected SeqID, Taxonomy, [Confidence]")

        writer.writerow(header_out)
        for row in reader:
            if not row:
                continue
            seqid, taxonomy = row[0], (row[1] if len(row) > 1 else "")
            writer.writerow([seqid] + split_lineage(taxonomy))
            n += 1

    sys.stderr.write(f"build_report_table: wrote {n} ASVs to {args.outfile}\n")


if __name__ == "__main__":
    main()
