#!/usr/bin/env python3
"""
plot_metastandard.py — Visualise MetaStandard16S output tables.

Produces two plots per input TSV:
  1. Stacked bar chart  — top N taxa by mean abundance, rest collapsed to "Other"
  2. Heatmap            — all taxa × samples (no top-N cap)

Label extraction rules
----------------------
- For taxonomic-level output (TaxID has no "|"):
    Walk ranks from deepest to shallowest; use the first non-"Unclassified" rank
    value with its prefix, e.g. "g__Limosilactobacillus", "f__Lactobacillaceae".
- For ASV-level output (TaxID contains "|"):
    Strip the sequence after "|", apply the same walk, then append "_N" to
    distinguish multiple ASVs resolving to the same base label.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


# Rank order shallowest → deepest (used for walking back up)
RANKS = ["d", "p", "c", "o", "f", "g", "s"]


# ---------------------------------------------------------------------------
# Label helpers
# ---------------------------------------------------------------------------

def extract_base_label(taxid: str) -> str:
    """Return the deepest non-Unclassified rank label, e.g. 'g__Bacillus'."""
    # Strip ASV sequence if present
    tax_part = taxid.split("|")[0]
    parts = [p.strip() for p in tax_part.split(";")]

    # Build prefix → value map
    rank_map = {}
    for part in parts:
        if "__" in part:
            prefix, value = part.split("__", 1)
            rank_map[prefix.strip()] = value.strip()

    # Walk deepest → shallowest, return first classified
    for rank in reversed(RANKS):
        value = rank_map.get(rank, "Unclassified")
        if value and value != "Unclassified":
            return f"{rank}__{value}"

    return "d__Unclassified"


def is_asv_level(df: pd.DataFrame) -> bool:
    return df["TaxID"].str.contains("|", regex=False).any()


def make_labels(df: pd.DataFrame) -> list:
    """
    Build display labels for every row.
    For ASV-level data, append _1/_2/... to duplicate base labels.
    """
    base = [extract_base_label(t) for t in df["TaxID"]]

    if not is_asv_level(df):
        return base

    # Count occurrences and assign incrementing suffixes
    seen: dict = defaultdict(int)
    labels = []
    for b in base:
        seen[b] += 1
        labels.append(f"{b}_{seen[b]}")
    return labels


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def prepare_plot_df(df: pd.DataFrame, labels: list) -> pd.DataFrame:
    """Return a copy of df with TaxID replaced by display labels."""
    plot_df = df.copy()
    plot_df["TaxID"] = labels
    plot_df = plot_df.set_index("TaxID")
    return plot_df


def collapse_to_top_n(plot_df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Keep top_n taxa by mean abundance; sum the rest as 'Other'."""
    mean_abund = plot_df.mean(axis=1).sort_values(ascending=False)
    top_taxa = mean_abund.index[:top_n].tolist()
    other_taxa = mean_abund.index[top_n:].tolist()

    top_df = plot_df.loc[top_taxa].copy()
    if other_taxa:
        other_row = plot_df.loc[other_taxa].sum(axis=0)
        other_row.name = "Other"
        top_df = pd.concat([top_df, other_row.to_frame().T])
    return top_df


def plot_stacked_bar(plot_df: pd.DataFrame, top_n: int, output_path: Path, title: str):
    """Stacked bar chart — samples on x-axis, taxa as stacked segments."""
    top_df = collapse_to_top_n(plot_df, top_n)

    n_taxa = len(top_df)
    palette = sns.color_palette("tab20", n_taxa)
    # Make "Other" grey if present
    if "Other" in top_df.index:
        other_idx = list(top_df.index).index("Other")
        palette[other_idx] = (0.75, 0.75, 0.75)

    fig, ax = plt.subplots(figsize=(max(6, len(plot_df.columns) * 1.2), 7))

    bottom = [0.0] * len(top_df.columns)
    for i, (taxon, row) in enumerate(top_df.iterrows()):
        values = row.values.astype(float)
        ax.bar(top_df.columns, values, bottom=bottom,
               color=palette[i], label=taxon, width=0.6)
        bottom = [b + v for b, v in zip(bottom, values)]

    ax.set_ylabel("Relative abundance")
    ax.set_xlabel("Sample")
    ax.set_title(title)
    ax.set_ylim(0, 1)
    plt.xticks(rotation=45, ha="right")

    # Legend outside the plot
    handles, labels_leg = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels_leg[::-1],
              bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, frameon=False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(plot_df: pd.DataFrame, output_path: Path, title: str):
    """Clustered heatmap — taxa × samples."""
    data = plot_df.astype(float).fillna(0.0).replace([float('inf'), float('-inf')], 0.0)

    # Avoid clustering when only one sample or one taxon
    row_cluster = len(data) > 1
    col_cluster = len(data.columns) > 1

    # Disable clustering if all values in a dimension are identical (zero variance)
    if row_cluster and (data.nunique(axis=1) <= 1).all():
        row_cluster = False
    if col_cluster and (data.nunique(axis=0) <= 1).all():
        col_cluster = False

    height = max(6, len(data) * 0.3)
    width  = max(5, len(data.columns) * 0.8)

    # Cap figure size to avoid memory errors with large sample counts
    max_dim = 60
    if width > max_dim or height > max_dim:
        scale = max_dim / max(width, height)
        width  = width * scale
        height = height * scale
    dpi = 150 if max(width, height) <= 40 else 100

    show_xlabels = len(data.columns) <= 100

    sys.setrecursionlimit(10000)
    g = sns.clustermap(
        data,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        cmap="YlOrRd",
        figsize=(width, height),
        xticklabels=show_xlabels,
        yticklabels=True,
        linewidths=0.3 if len(data.columns) <= 100 else 0,
        cbar_kws={"label": "Relative abundance"},
    )
    g.ax_heatmap.set_xlabel("Sample")
    g.ax_heatmap.set_ylabel("")
    g.fig.suptitle(title, y=1.01, fontsize=11)

    if show_xlabels:
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=7)

    g.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(g.fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Plot MetaStandard16S output")
    parser.add_argument("--input",    required=True, type=Path,
                        help="MetaStandard TSV file")
    parser.add_argument("--top_n",   type=int, default=10,
                        help="Top N taxa for stacked bar chart (default: 10)")
    parser.add_argument("--prefix",  default=None,
                        help="Output file prefix (default: input stem)")
    return parser.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(args.input, sep="\t")
    if "TaxID" not in df.columns:
        sys.exit(f"ERROR: no TaxID column found in {args.input}")

    sample_cols = [c for c in df.columns if c != "TaxID"]
    if not sample_cols:
        sys.exit(f"ERROR: no sample columns found in {args.input}")

    prefix = args.prefix or args.input.stem
    labels = make_labels(df)
    plot_df = prepare_plot_df(df, labels)

    stem = Path(args.input).stem
    title = stem.replace("_", " ")

    bar_path  = Path(f"{prefix}_barplot.png")
    heat_path = Path(f"{prefix}_heatmap.png")

    plot_stacked_bar(plot_df, args.top_n, bar_path,
                     title=f"{title} — top {args.top_n} taxa")
    plot_heatmap(plot_df, heat_path, title=title)

    print(f"Saved: {bar_path}, {heat_path}")


if __name__ == "__main__":
    main()
