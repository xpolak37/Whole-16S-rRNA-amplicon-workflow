#!/usr/bin/env python3
"""Pivot per-stage *.counts.tsv files into one wide ledger plus a retention table.

Each input file has the schema `sample\\tcount` (header + one or more rows).
The stage is deduced from the filename suffix:

    <sample>.cutadapt.counts.tsv          -> Cutadapt        (per-sample, 1 row)
    <sample>.host_removal.counts.tsv      -> Host            (per-sample, 1 row)
    <sample>.phix_removal.counts.tsv      -> Phix            (per-sample, 1 row)
    <sample>.orient.counts.tsv            -> Oriented        (per-sample, 1 row)
    dada2.counts.tsv                      -> Denoised_dada2  (cohort, N rows)
    dada2_nodenoise.counts.tsv            -> Denoised_dada2_nodenoise (cohort, N rows)

Two output files are written:
  <out>.tsv                Sample x Stage absolute counts.
  <out>_retention.tsv      Sample x Stage fraction retained from the previous stage
                           (first stage column blank — no predecessor).

Samples that fall to zero reads at any stage are logged to stderr.

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


def stage_for(path: Path):
    """Map filename to a stage column. Returns None if no stage matches."""
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


def fmt_pct(num: int, den: int) -> str:
    if den == 0:
        return "n/a"
    return f"{100.0 * num / den:.1f}%"


def main():
    if len(sys.argv) < 3:
        sys.exit("collate_counts.py: expected <output> <counts_tsv>...")
    out_path = Path(sys.argv[1])
    inputs = [Path(p) for p in sys.argv[2:]]

    # rows[sample][stage_col] = count (int)
    rows: dict[str, dict[str, int]] = {}
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
                try:
                    rows.setdefault(sample, {})[col] = int(count)
                except ValueError:
                    sys.exit(f"collate_counts.py: non-integer count {count!r} "
                             f"for sample {sample!r} in {path}")

    # Stable column order: pipeline order from STAGE_MAP, restricted to seen.
    ordered_cols = [col for _, col in STAGE_MAP if col in seen_stages]

    # --- absolute counts table -------------------------------------------------
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Sample", *ordered_cols])
        for sample in sorted(rows):
            writer.writerow([
                sample,
                *(rows[sample].get(c, "") for c in ordered_cols),
            ])

    print(f"[collate_counts] wrote {out_path} "
          f"({len(rows)} samples x {len(ordered_cols)} stages)")

    # --- retention table -------------------------------------------------------
    # For each sample, compute the fraction retained relative to the previous
    # stage's count. Stages with no predecessor count (e.g. first stage, or
    # missing earlier stages because a sample was dropped) get blank.
    retention_path = out_path.with_name(out_path.stem + "_retention.tsv")
    with retention_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Sample", *ordered_cols])
        for sample in sorted(rows):
            sample_counts = rows[sample]
            cells = []
            prev_count = None
            for col in ordered_cols:
                cur = sample_counts.get(col)
                if cur is None or prev_count is None:
                    cells.append("")
                else:
                    cells.append(fmt_pct(cur, prev_count))
                if cur is not None:
                    prev_count = cur
            writer.writerow([sample, *cells])

    print(f"[collate_counts] wrote {retention_path}")

    # --- warnings --------------------------------------------------------------
    # Flag samples that fell to zero at any stage. PacBio HiFi 16S typically has
    # 1k–50k reads/sample; a drop to zero usually means primer mismatch or a
    # library that didn't yield 16S.
    for sample in sorted(rows):
        for col in ordered_cols:
            if rows[sample].get(col) == 0:
                print(f"[collate_counts] WARNING: sample {sample!r} has 0 reads "
                      f"after stage {col}", file=sys.stderr)
                break


if __name__ == "__main__":
    main()
