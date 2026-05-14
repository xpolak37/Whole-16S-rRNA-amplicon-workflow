#!/usr/bin/env python3
"""
Build a custom QC summary report (raw txt + HTML) for a 16S pipeline run.

Two modes:

  parse   Read raw FastQC zip files and emit:
            - a JSON blob with adapters, per-file read counts and the
              deduped overrepresented sequences with their FastQC sources
            - a FASTA file with the unique overrepresented sequences,
              ready to feed into a remote BLAST step

  render  Read the parse-mode JSON plus an (optionally empty) BLAST
          hits TSV and emit both `raw_custom_summary.txt` (legacy plain
          text format) and a styled `custom_summary.html` report.

Splitting into two modes lets the parse/render steps run in a Python
biocontainer while the BLAST step runs in a BLAST biocontainer — neither
upstream image happens to ship both interpreters at once.
"""

import argparse
import html
import json
import statistics
import sys
import zipfile
from collections import OrderedDict
from pathlib import Path


# ---------------------------------------------------------------------------
# FastQC parsing
# ---------------------------------------------------------------------------

ADAPTER_PRESENCE_THRESHOLD_PCT = 0.1  # %% along any position to count as "present"


def _slugify_adapter(name: str) -> str:
    return (
        name.strip()
        .lower()
        .replace("'", "")
        .replace("-", "_")
        .replace(" ", "_")
    )


def _read_fastqc_data(zip_path: Path) -> str:
    """Return the contents of fastqc_data.txt inside a FastQC .zip."""
    with zipfile.ZipFile(zip_path) as zf:
        target = next(
            (n for n in zf.namelist() if n.endswith("/fastqc_data.txt")),
            None,
        )
        if target is None:
            raise ValueError(f"No fastqc_data.txt in {zip_path}")
        with zf.open(target) as fh:
            return fh.read().decode("utf-8", errors="replace")


def _split_modules(data: str):
    """Yield (module_name, [lines]) tuples from a fastqc_data.txt blob."""
    current_name = None
    current_lines = []
    for line in data.splitlines():
        if line.startswith(">>") and not line.startswith(">>END_MODULE"):
            current_name = line[2:].rsplit("\t", 1)[0]
            current_lines = []
        elif line.startswith(">>END_MODULE"):
            if current_name is not None:
                yield current_name, current_lines
            current_name = None
            current_lines = []
        else:
            if current_name is not None:
                current_lines.append(line)


def parse_fastqc_zip(zip_path: Path):
    """Return a dict with filename, total_sequences, adapters, overreps."""
    data = _read_fastqc_data(zip_path)
    filename = zip_path.stem
    if filename.endswith("_fastqc"):
        filename = filename[: -len("_fastqc")]
    total_sequences = None
    adapters_present = set()
    overreps = []  # list of (sequence, count, pct, source)

    for module_name, lines in _split_modules(data):
        if module_name == "Basic Statistics":
            for ln in lines:
                if ln.startswith("#") or not ln.strip():
                    continue
                parts = ln.split("\t")
                if len(parts) >= 2 and parts[0] == "Filename":
                    filename = parts[1]
                if len(parts) >= 2 and parts[0] == "Total Sequences":
                    try:
                        total_sequences = int(parts[1])
                    except ValueError:
                        pass

        elif module_name == "Adapter Content":
            header = None
            for ln in lines:
                if ln.startswith("#"):
                    header = ln.lstrip("#").rstrip().split("\t")
                    continue
                if not ln.strip() or header is None:
                    continue
                parts = ln.split("\t")
                # parts[0] is Position, rest are per-adapter percentages
                for col_name, val in zip(header[1:], parts[1:]):
                    try:
                        if float(val) >= ADAPTER_PRESENCE_THRESHOLD_PCT:
                            adapters_present.add(_slugify_adapter(col_name))
                    except ValueError:
                        continue

        elif module_name == "Overrepresented sequences":
            for ln in lines:
                if ln.startswith("#") or not ln.strip():
                    continue
                parts = ln.split("\t")
                if len(parts) < 4:
                    continue
                seq, count, pct, source = parts[0], parts[1], parts[2], parts[3]
                try:
                    count_i = int(count)
                except ValueError:
                    count_i = 0
                try:
                    pct_f = float(pct)
                except ValueError:
                    pct_f = 0.0
                overreps.append((seq, count_i, pct_f, source))

    return {
        "filename": filename,
        "total_sequences": total_sequences,
        "adapters": adapters_present,
        "overreps": overreps,
    }


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------


def depth_summary(counts):
    """Return min, q1, median, mean, q3, max for a list of ints."""
    if not counts:
        return None
    sorted_counts = sorted(counts)
    n = len(sorted_counts)
    mn = sorted_counts[0]
    mx = sorted_counts[-1]
    mean = statistics.fmean(sorted_counts)
    median = statistics.median(sorted_counts)
    if n >= 2:
        # Use inclusive method to mimic R's default quantile() behaviour
        # closely enough for a summary report.
        qs = statistics.quantiles(sorted_counts, n=4, method="inclusive")
        q1, _, q3 = qs
    else:
        q1 = q3 = median
    return {
        "min": mn,
        "q1": q1,
        "median": median,
        "mean": mean,
        "q3": q3,
        "max": mx,
    }


# ---------------------------------------------------------------------------
# Mode: parse
# ---------------------------------------------------------------------------


def cmd_parse(args):
    zips = sorted(args.fastqc_dir.glob("*_fastqc.zip"))
    if not zips:
        sys.stderr.write(
            f"[custom_summary] No FastQC zips found in {args.fastqc_dir}\n"
        )

    all_adapters = set()
    per_file_counts = []          # list of [filename, total_sequences]
    overreps_source = {}          # seq -> first-seen FastQC source label
    overreps_counts = {}          # seq -> cumulative count across all files

    for z in zips:
        try:
            parsed = parse_fastqc_zip(z)
        except Exception as exc:
            sys.stderr.write(f"[custom_summary] Failed to parse {z}: {exc}\n")
            continue
        all_adapters.update(parsed["adapters"])
        if parsed["total_sequences"] is not None:
            per_file_counts.append(
                [parsed["filename"], parsed["total_sequences"]]
            )
        for seq, count, _pct, source in parsed["overreps"]:
            if seq not in overreps_source:
                overreps_source[seq] = source
            overreps_counts[seq] = overreps_counts.get(seq, 0) + count

    # Sort by cumulative count descending and keep only the top N.
    overreps_ranked = sorted(
        overreps_source.keys(),
        key=lambda s: overreps_counts[s],
        reverse=True,
    )[: args.top_overreps]

    if len(overreps_source) > args.top_overreps:
        sys.stderr.write(
            f"[custom_summary] {len(overreps_source)} overrepresented sequences found; "
            f"keeping top {args.top_overreps} by abundance.\n"
        )

    payload = {
        "adapters": sorted(all_adapters),
        "per_file_counts": per_file_counts,
        "overreps": [[seq, overreps_source[seq]] for seq in overreps_ranked],
    }
    args.out_json.write_text(json.dumps(payload, indent=2))

    with open(args.out_fasta, "w") as fh:
        for i, seq in enumerate(overreps_ranked):
            fh.write(f">orep_{i}\n{seq}\n")


# ---------------------------------------------------------------------------
# Mode: render
# ---------------------------------------------------------------------------


def render_raw_txt(adapters, depth, overreps_unique, low_samples, threshold, blast_hits=None):
    lines = []
    lines.append("ADAPTER CONTENT:")
    for a in sorted(adapters):
        lines.append(a)
    lines.append("SEQUENCING DEPTH SUMMARY")
    if depth is None:
        lines.append("(no samples)")
    else:
        lines.append("   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. ")
        lines.append(
            f"{depth['min']:>7d} {depth['q1']:>7.0f} {depth['median']:>7.0f} "
            f"{depth['mean']:>7.0f} {depth['q3']:>7.0f} {depth['max']:>7d} "
        )
    lines.append("OVERREPRESENTED SEQUENCES:")
    for seq, fastqc_source in overreps_unique:
        # Prefer the BLAST-derived header (accession + description) when
        # available, falling back to the FastQC source label.
        header = (blast_hits or {}).get(seq) or fastqc_source
        lines.append(f">{header}")
        lines.append(seq)
    lines.append(f"SAMPLES BELOW {threshold // 1000}k:")
    for name, count in low_samples:
        lines.append(f"{name}\t{count}")
    return "\n".join(lines) + "\n"


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>16S Custom Run Summary{title_suffix}</title>
<style>
  :root {{
    --fg: #1f2933;
    --muted: #52606d;
    --accent: #0b7285;
    --bg: #f7f9fb;
    --card: #ffffff;
    --border: #d9e2ec;
    --warn: #b45309;
    --ok: #047857;
  }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    color: var(--fg);
    background: var(--bg);
    margin: 0;
    padding: 2rem;
    line-height: 1.5;
  }}
  h1 {{ margin-top: 0; color: var(--accent); }}
  h2 {{ margin-top: 0; color: var(--accent); font-size: 1.15rem; }}
  .meta {{ color: var(--muted); font-size: 0.9rem; margin-bottom: 1.5rem; }}
  .card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    margin-bottom: 1.25rem;
    box-shadow: 0 1px 2px rgba(0,0,0,0.03);
  }}
  table {{
    border-collapse: collapse;
    width: 100%;
    font-size: 0.9rem;
  }}
  th, td {{
    text-align: left;
    padding: 0.4rem 0.6rem;
    border-bottom: 1px solid var(--border);
  }}
  th {{ background: #f0f4f8; }}
  td.num {{ text-align: left; font-variant-numeric: tabular-nums; }}
  .seq {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 0.8rem; word-break: break-all; }}
  .pill {{
    display: inline-block;
    padding: 0.15rem 0.55rem;
    margin: 0.15rem 0.25rem 0.15rem 0;
    border-radius: 999px;
    background: #e3f2f7;
    color: var(--accent);
    font-size: 0.85rem;
  }}
  .empty {{ color: var(--muted); font-style: italic; }}
  .warn  {{ color: var(--warn); }}
  .ok    {{ color: var(--ok); }}
</style>
</head>
<body>
<h1>16S Custom Run Summary</h1>
<div class="meta">{meta_line}</div>

<div class="card">
  <h2>Adapter content</h2>
  {adapter_block}
</div>

<div class="card">
  <h2>Sequencing depth (raw reads)</h2>
  {depth_block}
</div>

<div class="card">
  <h2>Samples below {threshold_disp} reads</h2>
  {low_block}
</div>

<div class="card">
  <h2>Overrepresented sequences</h2>
  <p class="meta">{blast_status}</p>
  {overrep_block}
</div>

</body>
</html>
"""


def render_html(adapters, depth, overreps_unique, low_samples, threshold,
                blast_hits, blast_attempted, run_id, blast_db="16S_ribosomal_RNA"):
    title_suffix = f" — {run_id}" if run_id else ""
    meta_line = html.escape(f"Run ID: {run_id}" if run_id else "")

    if adapters:
        adapter_block = "".join(
            f'<span class="pill">{html.escape(a)}</span>' for a in sorted(adapters)
        )
    else:
        adapter_block = '<p class="empty ok">No adapters detected above threshold.</p>'

    if depth is None:
        depth_block = '<p class="empty">No samples available.</p>'
    else:
        depth_block = (
            "<table><thead><tr>"
            "<th>Min</th><th>1st Qu.</th><th>Median</th>"
            "<th>Mean</th><th>3rd Qu.</th><th>Max</th>"
            "</tr></thead><tbody><tr>"
            f'<td class="num">{depth["min"]:,}</td>'
            f'<td class="num">{depth["q1"]:,.0f}</td>'
            f'<td class="num">{depth["median"]:,.0f}</td>'
            f'<td class="num">{depth["mean"]:,.0f}</td>'
            f'<td class="num">{depth["q3"]:,.0f}</td>'
            f'<td class="num">{depth["max"]:,}</td>'
            "</tr></tbody></table>"
        )

    if low_samples:
        rows = "".join(
            f"<tr><td>{html.escape(name)}</td>"
            f'<td class="num warn">{count:,}</td></tr>'
            for name, count in low_samples
        )
        low_block = (
            "<table><thead><tr><th>Sample</th><th>Reads</th></tr></thead>"
            f"<tbody>{rows}</tbody></table>"
        )
    else:
        low_block = '<p class="empty ok">All samples are above the threshold.</p>'

    if not overreps_unique:
        overrep_block = '<p class="empty ok">No overrepresented sequences detected.</p>'
    else:
        rows = []
        show_blast_col = blast_hits is not None
        for seq, _fastqc_source in overreps_unique:
            blast_cell = ""
            if show_blast_col:
                hit = blast_hits.get(seq, "")
                blast_cell = f"<td>{html.escape(hit) if hit else '<span class=\"empty\">no hit</span>'}</td>"
            rows.append(
                f"<tr><td class='seq'>{html.escape(seq)}</td>"
                f"{blast_cell}</tr>"
            )
        header_cells = "<th>Sequence</th>"
        if show_blast_col:
            header_cells += f"<th>Top BLAST hit ({html.escape(blast_db)})</th>"
        overrep_block = (
            f"<table><thead><tr>{header_cells}</tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table>"
        )

    if not blast_attempted:
        blast_status = "BLAST disabled."
    elif blast_hits is None:
        blast_status = (
            f"BLAST against {blast_db} was attempted but failed. "
            "Step skipped — report still produced."
        )
    elif not overreps_unique:
        blast_status = "Nothing to BLAST — no overrepresented sequences."
    else:
        n_hit = sum(1 for s, _ in overreps_unique if blast_hits.get(s))
        blast_status = (
            f"BLAST against {blast_db}: {n_hit}/{len(overreps_unique)} "
            "sequences have a top hit shown below."
        )

    return HTML_TEMPLATE.format(
        title_suffix=html.escape(title_suffix),
        meta_line=meta_line,
        adapter_block=adapter_block,
        depth_block=depth_block,
        threshold_disp=f"{threshold:,}",
        low_block=low_block,
        overrep_block=overrep_block,
        blast_status=html.escape(blast_status),
    )


def _load_blast_tsv(path: Path, expected_seqs):
    """
    Load a blastn -outfmt '6 qseqid stitle' style TSV.

    Returns a dict {sequence: top_hit_title}, or None if the file is
    missing / contains the sentinel marker meaning the BLAST step was
    skipped (no internet, etc.).
    """
    if path is None or not path.exists() or path.stat().st_size == 0:
        return None
    text = path.read_text()
    # The blast process writes "BLAST_SKIPPED" if the remote call failed
    # so we can distinguish "BLAST ran and produced zero hits" from
    # "BLAST never ran successfully".
    if text.strip().startswith("BLAST_SKIPPED"):
        return None
    id_to_seq = {f"orep_{i}": seq for i, seq in enumerate(expected_seqs)}
    hits = {}
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t", 1)
        if len(parts) != 2:
            continue
        qid, title = parts
        seq = id_to_seq.get(qid)
        if seq and seq not in hits:
            hits[seq] = title
    return hits


def cmd_render(args):
    payload = json.loads(args.in_json.read_text())
    adapters = set(payload.get("adapters", []))
    per_file_counts = [tuple(x) for x in payload.get("per_file_counts", [])]
    overreps_unique = [tuple(x) for x in payload.get("overreps", [])]

    counts_only = [c for _, c in per_file_counts]
    depth = depth_summary(counts_only)
    low_samples = sorted(
        [(name, c) for name, c in per_file_counts if c < args.threshold],
        key=lambda x: (x[1], x[0]),
    )

    expected_seqs = [s for s, _ in overreps_unique]
    blast_hits = _load_blast_tsv(args.blast_tsv, expected_seqs) if args.blast_attempted else None

    raw_txt = render_raw_txt(
        adapters, depth, overreps_unique, low_samples, args.threshold, blast_hits
    )
    args.output_txt.write_text(raw_txt)

    html_doc = render_html(
        adapters,
        depth,
        overreps_unique,
        low_samples,
        args.threshold,
        blast_hits,
        blast_attempted=args.blast_attempted,
        run_id=args.run_id,
        blast_db=args.blast_db,
    )
    args.output_html.write_text(html_doc)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="mode", required=True)

    pp = sub.add_parser("parse", help="Parse FastQC zips into JSON + FASTA")
    pp.add_argument("--fastqc-dir", required=True, type=Path)
    pp.add_argument("--out-json", required=True, type=Path)
    pp.add_argument("--out-fasta", required=True, type=Path)
    pp.add_argument("--top-overreps", type=int, default=10,
                    help="Top N overrepresented sequences to include in report and BLAST (most abundant first).")
    pp.set_defaults(func=cmd_parse)

    pr = sub.add_parser("render", help="Render the HTML + raw txt summary")
    pr.add_argument("--in-json", required=True, type=Path)
    pr.add_argument("--blast-tsv", type=Path, default=None)
    pr.add_argument("--output-html", required=True, type=Path)
    pr.add_argument("--output-txt", required=True, type=Path)
    pr.add_argument("--threshold", type=int, default=10000)
    pr.add_argument("--run-id", default="")
    pr.add_argument("--blast-db", default="16S_ribosomal_RNA",
                    help="Name of the BLAST database used (shown in the report).")
    blast_group = pr.add_mutually_exclusive_group()
    blast_group.add_argument("--blast-attempted", dest="blast_attempted",
                             action="store_true")
    blast_group.add_argument("--no-blast-attempted", dest="blast_attempted",
                             action="store_false")
    pr.set_defaults(blast_attempted=True, func=cmd_render)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
