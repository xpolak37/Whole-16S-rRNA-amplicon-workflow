#!/usr/bin/env python3
import pandas as pd
import sys

# --- check arguments ---
if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} <fastqc_file> <outfile> <col_index>")

fastqc_file = sys.argv[1]   # table with "Sample" and "fastqc-total_sequences"
outfile     = sys.argv[2]   # table with columns: SampleID  Raw  Trimmed  Decontaminated
col_index   = int(sys.argv[3])  # 0-based index of the column to update


# --- load both tables ---
df_fastqc = pd.read_csv(fastqc_file, sep="\t")
df_out    = pd.read_csv(outfile, sep="\t")

# --- normalize SampleIDs (strip trailing _001 if present) ---
df_fastqc["Sample_norm"] = df_fastqc["Sample"].str.replace("_trimmed_cleaned.fastq.gz", "", regex=True)
df_out["Sample_norm"]    = df_out["Sample"].str.replace("_trimmed_cleaned.fastq.gz", "", regex=True)

df_fastqc["Sample_norm"] = df_fastqc["Sample_norm"].str.replace("_mergedpairs.fastq.gz", "_L001_R1", regex=True)
df_fastqc["Sample_norm"] = df_fastqc["Sample_norm"].str.replace(r"(S\d+)$", r"\1_L001_R1", regex=True)

df_fastqc["Sample_norm"] = df_fastqc["Sample_norm"].str.replace("_001$", "", regex=True)
df_out["Sample_norm"]    = df_out["Sample_norm"].str.replace("_001$", "", regex=True)

df_fastqc["Sample_norm"] = df_fastqc["Sample_norm"].str.replace("_", "-", regex=True)
df_out["Sample_norm"]    = df_out["Sample_norm"].str.replace("_", "-", regex=True)

# --- keep only rows with fastqc-total_sequences ---
df_fastqc = df_fastqc[df_fastqc["fastqc-total_sequences"].notnull()].copy()

# --- scale values ---
if col_index==2:
    df_fastqc["fastqc_scaled"] = (df_fastqc["fastqc-total_sequences"] * 1000000).round().astype(int)
else:
    df_fastqc["fastqc_scaled"] = (df_fastqc["fastqc-total_sequences"]).astype(int)
# --- merge (left join keeps all outfile rows) ---
df_merged = df_out.merge(
    df_fastqc[["Sample_norm", "fastqc_scaled"]],
    on="Sample_norm",
    how="left"
)

# --- update the given column index ---
colname = df_out.columns[col_index]
df_out[colname] = df_merged["fastqc_scaled"].fillna(df_out[colname])
df_out = df_out.drop(columns=["Sample_norm"])

# --- save back to the same file ---
df_out.to_csv(outfile, sep="\t", index=False)
