#!/usr/bin/env python3
"""
report_to_html.py
Converts pipeline_read_counts_report.tsv to styled HTML.
"""
import csv
import re
import sys
import pandas as pd
from pathlib import Path

BASE_CSS = """
body{font-family:Arial,sans-serif;margin:24px;background:#f8f9fa;font-size:12px}
h1{color:#2c3e50;font-size:18px}
h2{color:#34495e;margin-top:20px;font-size:14px}
p.note{font-size:11px;color:#888;margin-top:4px}
table{border-collapse:collapse;margin-bottom:20px;width:100%}
th{background:#2c3e50;color:#fff;padding:6px 10px;text-align:right;white-space:nowrap;font-size:11px}
th:first-child{text-align:left}
td{padding:4px 10px;border-bottom:1px solid #dde;text-align:right;white-space:nowrap}
td:first-child{text-align:left;font-weight:bold}
tr:nth-child(even){background:#ecf0f1}
tr:hover{background:#d5dbdb}
"""

# Cycled through in order for however many groups are present.
GROUP_COLOR_POOL = [
    "#0070b8", "#cc0000", "#1a8a3d", "#7d3c98",
    "#c9760d", "#0e8a8a", "#a83279", "#5d4037",
]


def load_samplesheet(path: Path) -> dict:
    """Return {sample: group} from a sample,group,fastq_1,fastq_2 sample sheet."""
    mapping = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mapping[row["sample"].strip()] = row["group"].strip()
    return mapping


def group_css_class(group: str) -> str:
    """A CSS-safe class name derived from an arbitrary group label."""
    return "grp-" + "".join(c if c.isalnum() else "-" for c in group.lower())

COL_ORDER = [
    "Sample",
    "Raw reads",
    "Trimmed",
    "STAR mapped %",
    "STAR unmapped",
    # "BBSplit: {ref} matched %" columns are dynamic (one per reference
    # genome BBSplit was run with) - spliced in here at render time.
    "MAGs matched %",
    "Candidate transRNAs total reads",
    "Kraken classified %",
    "BLAST classified %",
    "Classified total %",
    "Total transRNAs after all filters %",
    "Total transRNAs after all filters (RPM)",
    "Transmissible RNA representative sequences %",
    "Transmissible RNA representative sequences (RPM)",
    "% of transRNAs with >=5 duplicates",
]
BBSPLIT_RE = re.compile(r"^BBSplit: .+ matched %$")

def fmt(x):
    try:
        v = int(x)
        return f"{v:,}" if v >= 0 else ""
    except Exception:
        return "" if (str(x) in ("nan","")) else str(x)

def make_table(df, sample_to_group, group_css):
    static_cols = [c for c in COL_ORDER if c in df.columns]
    bbsplit_cols = sorted(c for c in df.columns if BBSPLIT_RE.match(c))
    insert_at = static_cols.index("STAR unmapped") + 1 if "STAR unmapped" in static_cols else len(static_cols)
    present = static_cols[:insert_at] + bbsplit_cols + static_cols[insert_at:]
    th = "".join(f"<th>{c}</th>" for c in present)
    rows_html = []
    for _, row in df.iterrows():
        sample = str(row.get("Sample",""))
        grp = sample_to_group.get(sample)
        cls = group_css.get(grp, "")
        cells = []
        for col in present:
            val = fmt(row.get(col,""))
            td_cls = f' class="{cls}"' if col=="Sample" and cls else ""
            cells.append(f"<td{td_cls}>{val}</td>")
        rows_html.append("<tr>" + "".join(cells) + "</tr>")
    return (
        f"<table><thead><tr>{th}</tr></thead>"
        f"<tbody>{''.join(rows_html)}</tbody></table>"
    )

def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} report.tsv report.html samplesheet.csv")
    df = pd.read_csv(sys.argv[1], sep="\t")
    sample_to_group = load_samplesheet(sys.argv[3])
    groups = sorted(set(sample_to_group.values()))
    group_css = {g: group_css_class(g) for g in groups}

    # Sort samples by the order their group first appears in the sample sheet,
    # then by sample name within a group.
    group_order = {g: i for i, g in enumerate(groups)}
    def sort_key(s):
        grp = sample_to_group.get(str(s))
        return (group_order.get(grp, len(groups)), str(s))
    df = df.iloc[df["Sample"].map(sort_key).argsort()]

    css_rules = "\n".join(
        f".{cls}{{color:{GROUP_COLOR_POOL[i % len(GROUP_COLOR_POOL)]}}}"
        for i, (g, cls) in enumerate(group_css.items())
    )

    table_html = make_table(df, sample_to_group, group_css)
    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8">
<title>transRNA Pipeline — Read Count Report</title>
<style>{BASE_CSS}
{css_rules}
</style>
</head><body>
<h1>transRNA Pipeline — Read Count Report</h1>
<p class="note">
  Decon-a/b/MAGs columns show reads REMAINING after each decontamination step.<br>
  Candidate totals use re-inflated counts (duplicate weights applied).<br>
  STAR mapped = Trimmed − STAR unmapped.
</p>
{table_html}
</body></html>"""
    Path(sys.argv[2]).write_text(html, encoding="utf-8")
    print(f"HTML written: {sys.argv[2]}", file=sys.stderr)

if __name__ == "__main__":
    main()