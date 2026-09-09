#!/usr/bin/env python3
"""
plot_report_summary.py
──────────────────────
Reads pipeline_read_counts_report.tsv and produces a single-panel read
filtering cascade: stacked horizontal bars showing what fraction of trimmed
reads is removed/passes at each filtering step, per sample.

Six segments, all as % of Trimmed reads: bee (STAR mapped), contaminants
(BBSplit, summed across whatever reference genomes were actually run),
metagenome (MAGs), invertebrates (Kraken+BLAST bad-hit discards), Not
classified (remainder), transRNAs (final survivors).

Group (for sample-label colour) comes from the real sample sheet, not a
naming-convention guess -- general-purpose across any dataset.

Usage:
    python3 plot_report_summary.py <pipeline_read_counts_report.tsv> --samplesheet <sheet.csv> [output.png]
"""

import argparse
import csv
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

TODAY = date.today().strftime("%d_%m_%Y")

CASCADE_COLORS = {
    "bee":            "#d62728",
    "contaminants":   "#ff7f0e",
    "metagenome":     "#2ca02c",
    "invertebrates":  "#9467bd",
    "Not classified": "#bcbd22",
    "transRNAs":      "#1f77b4",
}


def load_samplesheet(path: Path) -> dict:
    mapping = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mapping[row["sample"].strip()] = row["group"].strip()
    return mapping


def group_palette(groups: list) -> dict:
    cmap = plt.get_cmap("tab10")
    return {g: cmap(i % 10) for i, g in enumerate(groups)}


def load_report(path: Path, sample_to_group: dict) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["Group"] = df["Sample"].map(sample_to_group).fillna("unknown")
    return df


def bbsplit_ref_cols(df: pd.DataFrame):
    return [c for c in df.columns if c.startswith("BBSplit: ") and c.endswith(" matched %")]


def build_cascade(df: pd.DataFrame) -> pd.DataFrame:
    ref_cols = bbsplit_ref_cols(df)

    cas = pd.DataFrame(index=df.index)
    cas["Sample"] = df["Sample"]
    cas["Group"] = df["Group"]

    cas["bee"] = df["STAR mapped %"]
    star_unmapped_pct = 100 - cas["bee"]

    bbsplit_pct_of_unmapped = df[ref_cols].sum(axis=1) if ref_cols else 0.0
    cas["contaminants"] = bbsplit_pct_of_unmapped / 100 * star_unmapped_pct

    cas["metagenome"] = df["MAGs matched %"]
    cas["invertebrates"] = df.get("Invertebrates matched %", 0.0)

    accounted = cas["bee"] + cas["contaminants"] + cas["metagenome"] + cas["invertebrates"]
    cas["Not classified"] = (100 - accounted - df["Total transRNAs after all filters %"]).clip(lower=0)
    cas["transRNAs"] = df["Total transRNAs after all filters %"]

    return cas


def plot_cascade(ax, cas: pd.DataFrame, palette: dict):
    steps = ["bee", "contaminants", "metagenome", "invertebrates", "Not classified", "transRNAs"]
    samples = cas["Sample"].tolist()
    y_pos = np.arange(len(samples))
    bar_h = 0.65

    lefts = np.zeros(len(samples))
    for step in steps:
        vals = cas[step].clip(lower=0).to_numpy()
        bars = ax.barh(y_pos, vals, left=lefts, height=bar_h,
                       color=CASCADE_COLORS[step], label=step, edgecolor="white", linewidth=0.4)
        for bar, val, left in zip(bars, vals, lefts):
            if val > 1.5:
                ax.text(left + val / 2, bar.get_y() + bar.get_height() / 2,
                        f"{val:.1f}%", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
        lefts += vals

    ax.set_yticks(y_pos)
    ax.set_yticklabels(samples, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("% of trimmed reads", fontsize=10)
    ax.set_xlim(0, 102)
    ax.set_title("Read filtering cascade\n(% of trimmed reads)", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    for label in ax.get_yticklabels():
        sample = label.get_text()
        grp = cas.loc[cas["Sample"] == sample, "Group"].values
        label.set_color(palette.get(grp[0], "black") if len(grp) else "black")

    legend_patches = [mpatches.Patch(color=CASCADE_COLORS[s], label=s) for s in steps]
    ax.legend(handles=legend_patches, loc="upper center",
              bbox_to_anchor=(0.5, -0.12), ncol=6, fontsize=9,
              framealpha=0.9, edgecolor="grey")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("report", help="pipeline_read_counts_report.tsv")
    ap.add_argument("--samplesheet", required=True, help="sample,group,fastq_1,fastq_2 sample sheet")
    ap.add_argument("output", nargs="?", default=None)
    args = ap.parse_args()

    report_path = Path(args.report)
    out_path = Path(args.output) if args.output else \
               report_path.parent / f"{TODAY}_pipeline_report_summary.png"

    sample_to_group = load_samplesheet(Path(args.samplesheet))
    df = load_report(report_path, sample_to_group)
    cas = build_cascade(df)

    groups_present = sorted(set(df["Group"]) - {"unknown"})
    palette = group_palette(groups_present)
    palette["unknown"] = "grey"

    fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.55)))

    plot_cascade(ax, cas, palette)

    fig.suptitle("transRNA Pipeline — Read Count & Classification Summary",
                 fontsize=13, fontweight="bold", y=1.02)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
