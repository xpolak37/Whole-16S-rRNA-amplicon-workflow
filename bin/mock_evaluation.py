#!/usr/bin/env python3

r"""
Mock Community Evaluation Script
================================
Automatically evaluates mock community samples from the 16S pipeline output.

WHAT IT DOES
------------
1. Reads the MetaStandard output TSV file (genus-level relative abundances)
2. Identifies mock samples using regex pattern (default: Mock\d+)
3. Loads the reference mock community composition
4. Aggregates reference to genus level and compares with observed
5. Generates a stacked barplot showing reference vs observed compositions
6. Calculates and reports deviation metrics (Bray-Curtis, correlation)

INPUTS
------
--metastandard_table  Path to MetaStandard output TSV (TaxID + sample columns)
--mock_abundance      Path to mock reference abundance CSV (ASV, Mock Reference)
--mock_taxa           Path to mock reference taxonomy CSV (ASV, Domain...Species)
--mock_pattern        Regex pattern to identify mock samples (default: Mock\d+)
--output_prefix       Prefix for output files (default: mock_evaluation)
--top_n               Number of top genera to show (rest grouped as "Other")

OUTPUT
------
- <prefix>_barplot.png      : Stacked barplot comparing reference vs mock samples
- <prefix>_metrics.tsv      : Deviation metrics (Bray-Curtis, Pearson correlation)
- <prefix>_composition.tsv  : Full composition table (reference + samples)

USAGE
-----
python mock_evaluation.py \\
    --metastandard_table dada2PE_NaiveBayes_run01_genus.tsv \\
    --mock_abundance composition/mock_asv_abundance.csv \\
    --mock_taxa composition/mock_taxa.csv \\
    --output_prefix mock_evaluation

DEPENDENCIES
------------
pandas, matplotlib, numpy, scipy
"""

import argparse
import re
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import braycurtis
from scipy.stats import pearsonr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate mock community samples against reference composition"
    )
    parser.add_argument(
        "--metastandard_table",
        required=True,
        help="MetaStandard output TSV file with TaxID and sample columns"
    )
    parser.add_argument(
        "--mock_abundance",
        required=True,
        help="Mock reference abundance CSV (ASV, Mock Reference columns)"
    )
    parser.add_argument(
        "--mock_taxa",
        required=True,
        help="Mock reference taxonomy CSV (ASV, Domain...Genus columns)"
    )
    parser.add_argument(
        "--mock_pattern",
        default=r"Mock\d+",
        help="Regex pattern to identify mock samples (default: Mock\\d+)"
    )
    parser.add_argument(
        "--output_prefix",
        default="mock_evaluation",
        help="Prefix for output files"
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=15,
        help="Number of top genera to display (rest grouped as 'Other')"
    )
    parser.add_argument(
        "--run_id",
        default="run01",
        help="Run ID for labeling outputs"
    )
    parser.add_argument(
        "--synonyms",
        default=None,
        help="CSV file with old_name,canonical_name columns for genus name normalisation"
    )
    return parser.parse_args()


def load_synonyms(path):
    """
    Load a two-column CSV (old_name, canonical_name) and return a dict
    mapping every old name to its canonical equivalent.
    Returns an empty dict if path is None or missing.
    """
    if path is None:
        return {}
    p = Path(path)
    if not p.exists():
        print(f"WARNING: synonyms file not found: {path}")
        return {}
    synonyms = {}
    df = pd.read_csv(p)
    for _, row in df.iterrows():
        synonyms[row["name"].strip()] = row["canonical_name"].strip()
    return synonyms


def normalize_genus(genus, synonyms):
    """Return the canonical name for a genus, or the genus itself if unknown."""
    return synonyms.get(genus, genus)


def load_reference_composition(abundance_path, taxa_path, synonyms=None):
    """
    Load mock reference files and aggregate to genus level.

    Returns DataFrame with Genus as index and 'Reference' as column.
    """
    synonyms = synonyms or {}

    # Load abundance and taxonomy
    abundance = pd.read_csv(abundance_path)
    taxa = pd.read_csv(taxa_path)

    # Merge on ASV
    merged = pd.merge(abundance, taxa, on="ASV", how="left")

    # Normalize genus names using synonym dictionary
    merged["Genus"] = merged["Genus"].apply(lambda g: normalize_genus(g, synonyms))

    # Aggregate to genus level
    genus_abundance = merged.groupby("Genus")["Mock Reference"].sum()

    # Normalize to relative abundance (should already sum to 1, but ensure)
    genus_abundance = genus_abundance / genus_abundance.sum()

    return genus_abundance.rename("Reference")


def extract_genus_from_taxid(taxid):
    """
    Extract genus from MetaStandard TaxID string.
    Format: d__Bacteria;p__Firmicutes;c__Bacilli;o__...;f__...;g__Lactobacillus
    """
    if pd.isna(taxid):
        return "Unclassified"
    
    parts = taxid.split(";")
    for part in parts:
        if part.startswith("g__"):
            genus = part.replace("g__", "").strip()
            if genus and genus != "Unclassified":
                return genus
    return "Unclassified"


def identify_mock_samples(columns, pattern):
    """
    Identify mock sample columns using regex pattern.
    """
    regex = re.compile(pattern, re.IGNORECASE)
    mock_cols = [col for col in columns if regex.search(col)]
    return mock_cols


def load_metastandard_table(path, mock_pattern, synonyms=None):
    """
    Load MetaStandard output and filter to mock samples.
    Aggregate to genus level if needed.
    """
    synonyms = synonyms or {}
    df = pd.read_csv(path, sep="\t")

    # keep rows that do NOT start with d__Unclassified
    df = df[~df["TaxID"].str.startswith("d__Unclassified", na=False)]
    
    # Identify mock samples
    sample_cols = [col for col in df.columns if col != "TaxID"]
    mock_cols = identify_mock_samples(sample_cols, mock_pattern)

    if not mock_cols:
        print(f"WARNING: No mock samples found matching pattern '{mock_pattern}'")
        print(f"Available samples: {sample_cols}")
        return None, []

    print(f"Found {len(mock_cols)} mock samples: {mock_cols}")

    # Extract genus from TaxID and normalize via synonym dictionary
    df["Genus"] = df["TaxID"].apply(extract_genus_from_taxid).apply(
        lambda g: normalize_genus(g, synonyms)
    )

    # Aggregate to genus level (sum abundances)
    genus_df = df.groupby("Genus")[mock_cols].sum()

    # Normalize each sample to relative abundance
    genus_df = genus_df.div(genus_df.sum(axis=0), axis=1)

    return genus_df, mock_cols


def calculate_metrics(reference, observed):
    """
    Calculate deviation metrics between reference and observed compositions.
    
    Returns dict with:
    - bray_curtis: Bray-Curtis dissimilarity (0 = identical, 1 = completely different)
    - pearson_r: Pearson correlation coefficient
    - pearson_p: p-value for Pearson correlation
    - rmse: Root mean square error
    """
    # Align indices
    all_genera = reference.index.union(observed.index)
    ref_aligned = reference.reindex(all_genera, fill_value=0)
    obs_aligned = observed.reindex(all_genera, fill_value=0)
    
    # Bray-Curtis
    bc = braycurtis(ref_aligned.values, obs_aligned.values)
    
    # Pearson correlation
    r, p = pearsonr(ref_aligned.values, obs_aligned.values)
    
    # RMSE
    rmse = np.sqrt(np.mean((ref_aligned.values - obs_aligned.values) ** 2))
    
    return {
        "bray_curtis": bc,
        "pearson_r": r,
        "pearson_p": p,
        "rmse": rmse
    }


def prepare_plot_data(reference, observed_df, top_n):
    """
    Prepare data for stacked barplot.
    Combines reference and observed samples, keeps top_n genera.
    """
    # Combine reference and observed
    combined = pd.merge(observed_df, reference, how="outer", on=["Genus"])
    combined = combined.fillna(0)
    # Calculate mean abundance across all samples for ranking
    combined["mean_abundance"] = combined.mean(axis=1)
    combined = combined.sort_values("mean_abundance", ascending=False)
    
    # Keep top N genera, group rest as "Other"
    top_genera = combined.head(top_n).index.tolist()
    
    # Remove Unclassified from top if present, handle separately
    if "Unclassified" in top_genera:
        top_genera.remove("Unclassified")
        top_genera = top_genera[:top_n-1]  # Keep space for Unclassified
    
    # Create final dataframe
    sample_cols = [col for col in combined.columns if col != "mean_abundance"]
    
    plot_df = pd.DataFrame(index=top_genera + ["Other", "Unclassified"], columns=sample_cols)
    
    for genus in top_genera:
        if genus in combined.index:
            plot_df.loc[genus] = combined.loc[genus, sample_cols]
    
    # Sum "Other" genera
    other_genera = [g for g in combined.index if g not in top_genera and g != "Unclassified"]
    if other_genera:
        plot_df.loc["Other"] = combined.loc[other_genera, sample_cols].sum()
    else:
        plot_df.loc["Other"] = 0
    
    # Add Unclassified
    if "Unclassified" in combined.index:
        plot_df.loc["Unclassified"] = combined.loc["Unclassified", sample_cols]
    else:
        plot_df.loc["Unclassified"] = 0
    
    plot_df = plot_df.fillna(0).astype(float)
    
    # Remove rows that are all zeros
    plot_df = plot_df.loc[(plot_df != 0).any(axis=1)]
    
    return plot_df


def create_barplot(plot_df, output_path, run_id):
    """
    Create stacked barplot comparing reference vs mock samples.
    """
    # Color palette - distinct colors for genera
    n_genera = len(plot_df)
    cmap = plt.colormaps.get_cmap('tab20')
    colors = [cmap(i / max(n_genera, 1)) for i in range(n_genera)]
    
    # Special colors for Other and Unclassified
    genus_colors = {}
    color_idx = 0
    for genus in plot_df.index:
        if genus == "Other":
            genus_colors[genus] = "#808080"  # Gray
        elif genus == "Unclassified":
            genus_colors[genus] = "#D3D3D3"  # Light gray
        else:
            genus_colors[genus] = colors[color_idx]
            color_idx += 1
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(10, len(plot_df.columns) * 1.5), 8))
    
    # Transpose for plotting (samples as x-axis)
    plot_data = plot_df.T
    
    # Reorder columns to put Reference first
    cols = plot_data.columns.tolist()
    if "Reference" in plot_data.index.tolist():
        sample_order = ["Reference"] + [s for s in plot_data.index if s != "Reference"]
        plot_data = plot_data.reindex(sample_order)
    
    # Create stacked bar
    bottom = np.zeros(len(plot_data))
    
    for genus in plot_df.index:
        values = plot_data[genus].values
        ax.bar(
            plot_data.index,
            values,
            bottom=bottom,
            label=genus,
            color=genus_colors[genus],
            edgecolor="white",
            linewidth=0.5
        )
        bottom += values
    
    # Formatting
    ax.set_ylabel("Relative Abundance", fontsize=12)
    ax.set_xlabel("Sample", fontsize=12)
    ax.set_title(f"Mock Community Composition: Reference vs Observed\n(Run: {run_id})", fontsize=14)
    ax.set_ylim(0, 1.05)
    
    # Rotate x-axis labels
    plt.xticks(rotation=45, ha="right")
    
    # Legend outside plot
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=9,
        title="Genus"
    )
    
    # Add reference line
    ax.axhline(y=1.0, color="black", linestyle="--", alpha=0.3, linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"Barplot saved to: {output_path}")


def main():
    args = parse_args()
    
    print("=" * 60)
    print("Mock Community Evaluation")
    print("=" * 60)
    
    # Load synonym dictionary
    synonyms = load_synonyms(args.synonyms)
    if synonyms:
        print(f"Loaded {len(synonyms)} genus synonyms from: {args.synonyms}")

    # Load reference composition
    print(f"\nLoading reference composition from:")
    print(f"  Abundance: {args.mock_abundance}")
    print(f"  Taxonomy:  {args.mock_taxa}")

    reference = load_reference_composition(args.mock_abundance, args.mock_taxa, synonyms)
    print(f"\nReference composition (genus level):")
    print(reference.sort_values(ascending=False).to_string())

    # Load MetaStandard table and extract mock samples
    print(f"\nLoading MetaStandard table: {args.metastandard_table}")
    print(f"Mock sample pattern: {args.mock_pattern}")

    observed_df, mock_cols = load_metastandard_table(
        args.metastandard_table,
        args.mock_pattern,
        synonyms
    )
    
    if observed_df is None or len(mock_cols) == 0:
        print(
            "\nNOTE: No mock samples found in this dataset — skipping mock "
            "evaluation. This is not an error; the pipeline will continue."
        )
        sys.exit(0)
    
    # Calculate metrics for each mock sample
    print("\n" + "=" * 60)
    print("Deviation Metrics (Reference vs Observed)")
    print("=" * 60)
    
    metrics_list = []
    for sample in mock_cols:
        metrics = calculate_metrics(reference, observed_df[sample])
        metrics["sample"] = sample
        metrics_list.append(metrics)
        
        print(f"\n{sample}:")
        print(f"  Bray-Curtis dissimilarity: {metrics['bray_curtis']:.4f}")
        print(f"  Pearson correlation:       {metrics['pearson_r']:.4f} (p={metrics['pearson_p']:.2e})")
        print(f"  RMSE:                      {metrics['rmse']:.4f}")
    
    # Save metrics
    metrics_df = pd.DataFrame(metrics_list)
    metrics_df = metrics_df[["sample", "bray_curtis", "pearson_r", "pearson_p", "rmse"]]
    metrics_path = f"{args.output_prefix}_metrics.tsv"
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    print(f"\nMetrics saved to: {metrics_path}")
    
    # Prepare and save full composition table
    composition_df = prepare_plot_data(reference, observed_df, args.top_n)
    composition_path = f"{args.output_prefix}_composition.tsv"
    composition_df.to_csv(composition_path, sep="\t")
    print(f"Composition table saved to: {composition_path}")
    
    # Create barplot
    barplot_path = f"{args.output_prefix}_barplot.png"
    create_barplot(composition_df, barplot_path, args.run_id)
    
    print("\n" + "=" * 60)
    print("Mock community evaluation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
