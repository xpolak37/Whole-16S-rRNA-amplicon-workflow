#!/usr/bin/env python3

"""
MetaStandard 16S - ASV Taxonomic Profile Unifier
=============================================
Merges and standardises ASV-based taxonomic abundance profiles from a
denoising tool (DADA2 here; UNOISE/Deblur in the sibling Illumina pipeline)
into a single, unified TSV table.

WHAT IT DOES
------------
1. Auto-detects the denoising tool from input filenames, OR uses --denoiser
   when provided (preferred for callers that already know the tool).
2. Parses ASV and taxonomy tables, linking sequences to taxonomic annotations.
3. Aggregates ASV counts to the requested taxonomic level (default: genus).
   Use --level asv to skip aggregation and keep individual ASV rows.
4. Standardises taxonomy strings using rank prefixes (d__, p__, c__, o__, f__, g__).
5. Merges all samples into one wide-format table (taxa x samples).
6. Writes the result to a TSV file.

INPUTS
------
--asv_table   ASV count table (TSV). Rows are SeqIDs, columns are samples.

--taxa_table  Taxonomy table (TSV). Must contain SeqID, Taxonomy, Confidence.
              Taxonomy uses rank-prefixed semicolon format:
                d__Bacteria;p__Firmicutes;c__Bacilli;...;g__Lactobacillus

--taxa_tool   Name of the classifier (e.g. qnb, qblast, idtaxa, assigntaxonomy).

--denoiser    Optional. Denoising tool used to produce --asv_table
              (e.g. dada2, dada2_nodenoise). When set, overrides filename
              detection and is used as the asv_tool prefix in the output
              filename. PacBio Nextflow callers should always pass this.

--level       Taxonomic level: domain, phylum, class, order, family,
              genus, species, asv. Default: genus.

--run_id      Label appended to the output filename. Default: run01

OUTPUT
------
<asv_tool>_<taxa_tool>_<run_id>_<level>.tsv in the working directory.
"""

import argparse
import pandas as pd

RANK_PREFIXES = ["d", "p", "c", "o", "f", "g", "s"]

LEVEL_TO_PREFIX = {
    "domain":  "d",
    "phylum":  "p",
    "class":   "c",
    "order":   "o",
    "family":  "f",
    "genus":   "g",
    "species": "s",
}

def parse_args():
    parser = argparse.ArgumentParser(description="MetaStandard: unify taxonomic profiles")

    parser.add_argument(
        "--asv_table",
        required=True,
        help="Input ASV table"
    )

    parser.add_argument(
        "--taxa_table",
        required=True,
        help="Input taxa table"
    )

    parser.add_argument(
        "--taxa_tool",
        required=True,
        help="Tool used to generate the taxonomy"
    )

    parser.add_argument(
        "--denoiser",
        default=None,
        help="Denoiser tool (overrides filename detection used for output prefix)"
    )

    parser.add_argument(
        "--level",
        default="genus",
        help="Taxonomic level (domain, phylum, class, order, family, genus, species, asv)"
    )

    parser.add_argument(
        "--run_id",
        default="run01",
        help="ID of the run in order to recognize the parameters used"
    )

    return parser.parse_args()


def detect_tool(f):

    if "dada2_paired" in f.lower():
        return "dada2PE"

    if "dada2_single" in f.lower():
        return "dada2SE"

    if "unoise" in f.lower():
        return "unoise"

    if "deblur" in f.lower():
        return "deblur"

    return "unknown"


def parse_taxonomy(taxonomy_string):
    """Parse d__Bacteria;p__...;g__Genus style string into a dictionary"""
    ranks = {}
    if pd.isna(taxonomy_string):
        return ranks
    for part in taxonomy_string.split(";"):
        part = part.strip()
        if "__" in part:
            prefix, value = part.split("__", 1)
            ranks[prefix.strip()] = value.strip()
    return ranks


def build_tax_df(taxa_table):
    """Parse taxonomy column into a DataFrame of rank columns (d..s)."""
    tax_parsed = taxa_table["Taxonomy"].apply(parse_taxonomy)
    tax_df = pd.DataFrame(tax_parsed.tolist(), index=taxa_table["SeqID"])
    tax_df = tax_df.reindex(columns=RANK_PREFIXES, fill_value="Unclassified")
    tax_df = tax_df.replace(r'^\s*$', "Unclassified", regex=True)
    tax_df = tax_df.fillna("Unclassified")
    return tax_df


def aggregate_to_level(asv_table, taxa_table, level):
    """Aggregate ASV counts to the requested taxonomic level."""
    depth_prefix = LEVEL_TO_PREFIX[level]
    rank_cols = RANK_PREFIXES[: RANK_PREFIXES.index(depth_prefix) + 1]

    tax_df = build_tax_df(taxa_table)
    tax_df = tax_df.reindex(columns=rank_cols, fill_value="Unclassified")

    asv_indexed = asv_table.set_index("SeqID")
    merged = tax_df.join(asv_indexed, how="outer")

    sample_cols = [col for col in merged.columns if col not in rank_cols]

    merged[rank_cols] = merged[rank_cols].fillna("Unclassified")
    merged[sample_cols] = merged[sample_cols].fillna(0)

    grouped = merged.groupby(rank_cols)[sample_cols].sum().reset_index()

    grouped["TaxID"] = grouped[rank_cols].apply(
        lambda row: ";".join([f"{p}__{row[p]}" for p in rank_cols]),
        axis=1
    )

    grouped = grouped.drop(columns=rank_cols)
    grouped = grouped[["TaxID"] + sample_cols]

    return grouped


def aggregate_to_asv(asv_table, taxa_table):
    """Keep individual ASV rows; TaxID = full taxonomy string | ASV sequence."""
    tax_df = build_tax_df(taxa_table)

    asv_indexed = asv_table.set_index("SeqID")
    merged = tax_df.join(asv_indexed, how="outer")

    sample_cols = [col for col in merged.columns if col not in RANK_PREFIXES]

    merged[RANK_PREFIXES] = merged[RANK_PREFIXES].fillna("Unclassified")
    merged[sample_cols] = merged[sample_cols].fillna(0)

    tax_string = merged[RANK_PREFIXES].apply(
        lambda row: ";".join([f"{p}__{row[p]}" for p in RANK_PREFIXES]),
        axis=1
    )
    merged["TaxID"] = tax_string + "|" + merged.index

    result = merged[["TaxID"] + sample_cols].reset_index(drop=True)

    return result


def main():
    args = parse_args()

    level = args.level.lower()
    if level not in LEVEL_TO_PREFIX and level != "asv":
        raise ValueError(
            f"Unknown level '{args.level}'. "
            f"Choose from: {', '.join(list(LEVEL_TO_PREFIX.keys()) + ['asv'])}"
        )

    asv_table  = pd.read_csv(args.asv_table,  sep="\t")
    taxa_table = pd.read_csv(args.taxa_table, sep="\t")
    asv_tool   = args.denoiser if args.denoiser else detect_tool(args.taxa_table)

    if level == "asv":
        result_df = aggregate_to_asv(asv_table, taxa_table)
    else:
        result_df = aggregate_to_level(asv_table, taxa_table, level)

    sample_cols = [col for col in result_df.columns if col != "TaxID"]
    result_df[sample_cols] = result_df[sample_cols].div(
        result_df[sample_cols].sum(axis=0), axis=1
    )

    outfile = f"{asv_tool}_{args.taxa_tool}_{args.run_id}_{level}.tsv"
    result_df.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
