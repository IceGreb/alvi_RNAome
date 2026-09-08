#!/usr/bin/env python3
"""
plot_report_summary.py
──────────────────────
Reads pipeline_read_counts_report.tsv and produces a two-panel summary figure:

  Panel 1 — Read cascade: stacked horizontal bars showing what fraction of
             trimmed reads is removed/passes at each filtering step, per sample.

  Panel 2 — Final metrics heatmap: Kraken%, BLAST%, Total transRNAs%,
             Representative seqs%, and % with >=5 dups, per sample.

Group (for bar/label colour and legend) comes from the real sample sheet,
not a naming-convention guess -- general-purpose across any dataset, not
just samples named RJ*/T*.

Usage:
    python3 plot_report_summary.py <pipeline_read_counts_report.tsv> --samplesheet <sheet.csv> [output.png]
"""

import argparse
import csv
import sys
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
    "STAR mapped":     "#d62728",
    "BBSplit removed": "#ff7f0e",
    "MAGs":            "#8c564b",
    "Not classified":  "#bcbd22",
    "transRNAs":       "#2ca02c",
}

HEATMAP_COLS = [
    "Kraken classified %",
    "BLAST classified %",
    "Classified total %",
    "Total transRNAs after all filters %",
    "Transmissible RNA representative sequences %",
    "% of transRNAs with >=5 duplicates",
]

HEATMAP_LABELS = [
    "Kraken\nclassified %",
    "BLAST\nclassified %",
    "Classified\ntotal %",
    "Total\ntransRNAs %",
    "Representative\nseqs %",
    "≥5 dups\n% of transRNAs",
]


def load_samplesheet(path: Path) -> dict[str, str]:
    """Return {sample: group} from a sample,group,fastq_1,fastq_2 sample sheet."""
    mapping = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mapping[row["sample"].strip()] = row["group"].strip()
    return mapping


def group_palette(groups: list[str]) -> dict[str, str]:
    """Cycle a qualitative colormap across however many distinct groups
    actually appear -- not a hardcoded two-colour RJ/ST scheme."""
    cmap = plt.get_cmap("tab10")
    return {g: cmap(i % 10) for i, g in enumerate(groups)}


def load_report(path: Path, sample_to_group: dict[str, str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["Group"] = df["Sample"].map(sample_to_group).fillna("unknown")
    return df


def bbsplit_ref_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith("BBSplit: ") and c.endswith(" matched %")]


def build_cascade(df: pd.DataFrame) -> pd.DataFrame:
    """Express each filtering step as % of Trimmed reads.

    Real report schema (aggregate_report.py): "STAR mapped %" is given
    directly; decontamination is a dynamic set of "BBSplit: {ref} matched %"
    columns (one per reference genome BBSplit was run against, not a fixed
    Host+Human/viral split), summed here into one "BBSplit removed" segment
    so this works for any number/choice of BBSplit references, not just the
    two this pipeline happened to use historically. "MAGs matched %" is the
    real metagenome-assignment fraction (of Trimmed). The remainder down to
    "Total transRNAs after all filters %" is "Not classified".
    """
    t = df["Trimmed"]
    ref_cols = bbsplit_ref_cols(df)

    cas = pd.DataFrame(index=df.index)
    cas["Sample"] = df["Sample"]
    cas["Group"] = df["Group"]

    cas["STAR mapped"] = df["STAR mapped %"]
    star_unmapped_pct = 100 - cas["STAR mapped"]

    # Each BBSplit ref's % is of STAR-unmapped reads (per aggregate_report.py);
    # sum them (of STAR-unmapped), then express that sum as % of Trimmed.
    bbsplit_pct_of_unmapped = df[ref_cols].sum(axis=1) if ref_cols else 0.0
    cas["BBSplit removed"] = bbsplit_pct_of_unmapped / 100 * star_unmapped_pct

    cas["MAGs"] = df["MAGs matched %"]

    accounted = cas["STAR mapped"] + cas["BBSplit removed"] + cas["MAGs"]
    cas["Not classified"] = (100 - accounted - df["Total transRNAs after all filters %"]).clip(lower=0)
    cas["transRNAs"] = df["Total transRNAs after all filters %"]

    return cas


def plot_cascade(ax, cas: pd.DataFrame, palette: dict[str, str]):
    steps = ["STAR mapped", "BBSplit removed", "MAGs", "Not classified", "transRNAs"]
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
    ax.set_title("Read filtering cascade\n(% of trimmed reads)", fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    for label in ax.get_yticklabels():
        sample = label.get_text()
        grp = cas.loc[cas["Sample"] == sample, "Group"].values
        label.set_color(palette.get(grp[0], "black") if len(grp) else "black")

    legend_patches = [mpatches.Patch(color=CASCADE_COLORS[s], label=s) for s in steps]
    ax.legend(handles=legend_patches, loc="upper center",
              bbox_to_anchor=(0.5, -0.08), ncol=3, fontsize=8,
              framealpha=0.9, edgecolor="grey")


def plot_heatmap(ax, df: pd.DataFrame, palette: dict[str, str]):
    present = [c for c in HEATMAP_COLS if c in df.columns]
    labels = [HEATMAP_LABELS[HEATMAP_COLS.index(c)] for c in present]

    matrix = df[present].to_numpy(dtype=float)
    samples = df["Sample"].tolist()
    groups = df["Group"].tolist()

    im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0,
                   vmax=max(matrix.max(), 1))

    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=9)
    ax.set_title("Classification & transRNA metrics (%)", fontsize=11, fontweight="bold")

    for label in ax.get_yticklabels():
        sample = label.get_text()
        grp = groups[samples.index(sample)]
        label.set_color(palette.get(grp, "black"))

    for i in range(len(samples)):
        for j in range(len(present)):
            val = matrix[i, j]
            text_col = "white" if val > matrix.max() * 0.6 else "black"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=8, color=text_col, fontweight="bold")

    plt.colorbar(im, ax=ax, shrink=0.8, label="%")

    legend_patches = [mpatches.Patch(color=c, label=g) for g, c in palette.items()]
    ax.legend(handles=legend_patches, loc="upper center",
              bbox_to_anchor=(0.5, -0.08), ncol=min(len(legend_patches), 4), fontsize=8,
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

    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(18, max(6, len(df) * 0.9)),
        gridspec_kw={"width_ratios": [1.4, 1]}
    )

    plot_cascade(ax1, cas, palette)
    plot_heatmap(ax2, df, palette)

    fig.suptitle("transRNA Pipeline — Read Count & Classification Summary",
                 fontsize=13, fontweight="bold", y=1.01)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
