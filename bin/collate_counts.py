#!/usr/bin/env python3
"""Pivot per-stage *.counts.tsv files into one wide ledger.

Each input file has the schema `sample\\tcount` (header + one or more rows).
The stage is deduced from the filename suffix:

    <sample>.cutadapt.counts.tsv          -> Cutadapt        (per-sample, 1 row)
    <sample>.host_removal.counts.tsv      -> Host            (per-sample, 1 row)
    <sample>.phix_removal.counts.tsv      -> Phix            (per-sample, 1 row)
    <sample>.orient.counts.tsv            -> Oriented        (per-sample, 1 row)
    dada2.counts.tsv                      -> Denoised_dada2  (cohort, N rows)
    dada2_nodenoise.counts.tsv            -> Denoised_dada2_nodenoise (cohort, N rows)

Output column order: Sample, then stages in pipeline order. Missing values
appear as empty cells.

Usage: collate_counts.py <output_tsv> <counts_tsv> [<counts_tsv>...]
"""
import csv
import sys
from pathlib import Path

# Filename suffix -> wide-ledger column name. Order here defines column order.
STAGE_MAP = [
    ("cutadapt",       "Cutadapt"),
    ("host_removal",   "Host"),
    ("phix_removal",   "Phix"),
    ("orient",         "Oriented"),
    ("dada2",                "Denoised_dada2"),
    ("dada2_nodenoise",      "Denoised_dada2_nodenoise"),
]


def stage_for(path: Path) -> str:
    """Map filename to a stage column. Returns None if no stage matches."""
    # Strip trailing .counts.tsv, then look at the last dot-separated token
    # (per-sample files) or the whole basename (cohort files).
    name = path.name
    if not name.endswith(".counts.tsv"):
        return None
    stem = name[: -len(".counts.tsv")]
    # Cohort files: the whole stem matches a known stage key
    for key, col in STAGE_MAP:
        if stem == key:
            return col
    # Per-sample files: <sample>.<stage>
    if "." in stem:
        suffix = stem.rsplit(".", 1)[1]
        for key, col in STAGE_MAP:
            if suffix == key:
                return col
    return None


def main():
    if len(sys.argv) < 3:
        sys.exit("collate_counts.py: expected <output> <counts_tsv>...")
    out_path = Path(sys.argv[1])
    inputs = [Path(p) for p in sys.argv[2:]]

    # rows[sample][stage_col] = count
    rows: dict[str, dict[str, str]] = {}
    seen_stages: list[str] = []

    for path in inputs:
        col = stage_for(path)
        if col is None:
            print(f"[collate_counts] skipping unrecognized file: {path.name}",
                  file=sys.stderr)
            continue
        if col not in seen_stages:
            seen_stages.append(col)
        with path.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            if header != ["sample", "count"]:
                sys.exit(f"collate_counts.py: {path} has unexpected header "
                         f"{header!r} (expected ['sample', 'count'])")
            for row in reader:
                if not row:
                    continue
                sample, count = row[0], row[1]
                rows.setdefault(sample, {})[col] = count

    # Stable column order: pipeline order from STAGE_MAP, restricted to seen.
    ordered_cols = [col for _, col in STAGE_MAP if col in seen_stages]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Sample", *ordered_cols])
        for sample in sorted(rows):
            writer.writerow([sample, *(rows[sample].get(c, "") for c in ordered_cols)])

    print(f"[collate_counts] wrote {out_path} "
          f"({len(rows)} samples x {len(ordered_cols)} stages)")


if __name__ == "__main__":
    main()
