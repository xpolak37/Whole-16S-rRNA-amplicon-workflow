#!/usr/bin/env python3
"""Merge per-chunk QIIME classify-consensus-blast results into one taxa table.

classify-consensus-blast assigns each query sequence independently: every ASV is
searched against the reference and its consensus taxonomy is derived from its own
hits. Nothing is normalised across the query set, and BLAST E-values depend on
the reference size rather than on how many queries share a batch. Splitting the
rep-seqs and merging here is therefore exactly equivalent to one large run -- it
only buys parallelism, since classify-consensus-blast was observed pinned to a
single core regardless of --p-num-threads.

Rows are re-sorted by the numeric part of the ASV id so the merged table matches
the order a single unchunked run would have produced, independent of the order
chunks happen to finish in.

Usage: merge_blast_chunks.py <output.tsv> <chunk.tsv> [<chunk.tsv> ...]
"""
import csv
import re
import sys

ASV_RE = re.compile(r"^ASV_(\d+)$")


def main(argv):
    if len(argv) < 3:
        sys.exit("usage: merge_blast_chunks.py <output.tsv> <chunk.tsv> [...]")
    out_path, chunk_paths = argv[1], argv[2:]

    rows = []
    for path in chunk_paths:
        with open(path, newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            if header is None:
                continue
            for row in reader:
                if row:
                    rows.append(row)

    if not rows:
        sys.exit("merge_blast_chunks.py: no classified sequences in any chunk")

    # Restore ASV_<i> order. Anything not matching the pattern keeps its relative
    # position at the end rather than being silently dropped or reordered.
    def sort_key(item):
        idx, row = item
        m = ASV_RE.match(row[0])
        return (0, int(m.group(1)), 0) if m else (1, 0, idx)

    rows = [row for _, row in sorted(enumerate(rows), key=sort_key)]

    unclassifiable = sum(1 for r in rows if not ASV_RE.match(r[0]))
    if unclassifiable:
        print(f"warning: {unclassifiable} row(s) had non-ASV_<i> ids and were "
              f"appended in chunk order", file=sys.stderr)

    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["SeqID", "Taxonomy", "Confidence"])
        writer.writerows(rows)

    print(f"merged {len(chunk_paths)} chunk(s) -> {len(rows)} classified ASVs")


if __name__ == "__main__":
    main(sys.argv)
