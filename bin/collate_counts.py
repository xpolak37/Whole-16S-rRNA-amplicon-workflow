#!/usr/bin/env python3
"""Stub — joins per-stage *.counts.tsv into one wide TSV.

Real implementation in a later spec. Stub writes a header-only TSV.

Usage: collate_counts.py <output_tsv> <counts_tsv> [<counts_tsv>...]
"""
import sys
from pathlib import Path

def main():
    if len(sys.argv) < 3:
        sys.exit("collate_counts.py: expected <output> <counts_tsv>...")
    out = Path(sys.argv[1])
    counts_files = sys.argv[2:]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("sample\t" + "\t".join(Path(p).stem.split('.')[1] for p in counts_files) + "\n")
    print(f"[stub] collate_counts.py wrote header-only {out}")

if __name__ == "__main__":
    main()
