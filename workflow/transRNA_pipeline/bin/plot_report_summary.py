#!/usr/bin/env python3
"""
plot_report_summary.py
──────────────────────
Reads pipeline_read_counts_report.tsv and produces two separate figures:

  1. Read filtering cascade (left panel of the original 2-panel design,
     now its own file) -- stacked horizontal bars showing what fraction of
     trimmed reads is removed/passes at each filtering step, per sample.
     Six segments, all as % of Trimmed reads: bee (STAR mapped),
     contaminants (BBSplit, summed across whatever reference genomes were
     actually run), metagenome (MAGs), taxonomically excluded (Kraken+BLAST
     bad-hit discards: virus/host/invertebrate), alignment excluded
     (remainder -- length/mismatch/gap/etc. failures), transRNAs (final
     survivors).

  2. Classification & transRNA metrics heatmap (the original design's right
     panel) -- Kraken/BLAST classified %, Total transRNAs %, representative
     seqs %, % with >=5 duplicates, per sample.

Group (for sample-label colour) comes from the real sample sheet, not a
naming-convention guess -- general-purpose across any dataset.

Usage:
    python3 plot_report_summary.py <pipeline_read_counts_report.tsv> --samplesheet <sheet.csv> [output_base]
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

# Real, editable <text> elements in SVG output (not outlined into paths --
# matplotlib's default) -- same convention used by every other SVG-producing
# script in this project (e.g. plot_size_distr_v2.py's own comment on this).
plt.rcParams["svg.fonttype"] = "none"

TODAY = date.today().strftime("%d_%m_%Y")

CASCADE_COLORS = {
    "bee":                    "#d62728",
    "contaminants":           "#ff7f0e",
    "metagenome":             "#2ca02c",
    "taxonomically excluded": "#9467bd",
    "alignment excluded":     "#bcbd22",
    "transRNAs":              "#1f77b4",
}

HEATMAP_COLS = [
    "Kraken classified %",
    "BLAST classified %",
    "Classified total %",
    "% of transRNAs with >=5 duplicates",
    "Transmissible RNA representative sequences %",
    "Total transRNAs after all filters %",
]

HEATMAP_LABELS = [
    "Kraken\nclassified %",
    "BLAST\nclassified %",
    "Classified\ntotal %",
    "≥5 dups\n% of transRNAs",
    "Representative\nseqs %",
    "Total\ntransRNAs %",
]


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


def build_display_labels(df: pd.DataFrame) -> dict:
    """{Sample: display label}, e.g. "TB1" -> "SLT4" -- samples within a
    group get a clean sequential label ("{group}{i}", i from a sort of
    their real sample IDs), so the underlying report table keeps real,
    file-traceable sample IDs while the plots show a friendlier name. A
    group with only one sample keeps its real ID (nothing to number)."""
    labels = {}
    for group, sub in df.groupby("Group"):
        samples_sorted = sorted(sub["Sample"].tolist())
        if len(samples_sorted) <= 1:
            labels.update({s: s for s in samples_sorted})
        else:
            labels.update({s: f"{group}{i}" for i, s in enumerate(samples_sorted, start=1)})
    return labels


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
    cas["taxonomically excluded"] = df.get("Invertebrates matched %", 0.0)

    accounted = cas["bee"] + cas["contaminants"] + cas["metagenome"] + cas["taxonomically excluded"]
    cas["alignment excluded"] = (100 - accounted - df["Total transRNAs after all filters %"]).clip(lower=0)
    cas["transRNAs"] = df["Total transRNAs after all filters %"]

    return cas


def plot_cascade(ax, cas: pd.DataFrame, palette: dict, display_labels: dict, warp_power: float = 0.5):
    """Segment DISPLAY width is a power-law warp of its true value
    (width_display ~ value**warp_power, renormalised per sample so the
    warped widths still sum to 100), not the raw value -- so a small
    segment (e.g. contaminants at 1.9%) gets disproportionately more
    visual width than a linear stacked bar would give it, regardless of
    where it falls in the stack (unlike a fixed position-based axis
    break, this works no matter how big "bee" happens to be for a given
    sample). True percentages are unchanged and still shown in the
    labels -- only the drawn width is warped."""
    steps = ["bee", "contaminants", "metagenome", "taxonomically excluded", "alignment excluded", "transRNAs"]
    samples = cas["Sample"].tolist()
    y_pos = np.arange(len(samples))
    bar_h = 0.65

    raw = cas[steps].clip(lower=0).to_numpy()               # true % values, one row per sample
    warped = np.power(raw, warp_power)
    warped = warped / warped.sum(axis=1, keepdims=True) * 100.0  # renormalise to still sum to 100 per row

    lefts = np.zeros(len(samples))
    for j, step in enumerate(steps):
        vals_true   = raw[:, j]
        vals_display = warped[:, j]
        bars = ax.barh(y_pos, vals_display, left=lefts, height=bar_h,
                       color=CASCADE_COLORS[step], label=step, edgecolor="white", linewidth=0.4)
        for bar, val_true, val_disp, left in zip(bars, vals_true, vals_display, lefts):
            if val_disp > 1.0:
                ax.text(left + val_disp / 2, bar.get_y() + bar.get_height() / 2,
                        f"{val_true:.1f}%", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
        lefts += vals_display

    ax.set_yticks(y_pos)
    ax.set_yticklabels([display_labels.get(s, s) for s in samples], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Segment width (√-scaled for visibility -- labels show true % of trimmed reads)", fontsize=9)
    ax.set_xticks([])
    ax.set_xlim(0, 102)
    ax.set_title("Read filtering cascade\n(segment widths \u221a-scaled; labels = true %)", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    # Zipped against the real sample list (not label.get_text(), which is
    # now the display name, not the real Sample ID the Group lookup needs).
    for label, sample in zip(ax.get_yticklabels(), samples):
        grp = cas.loc[cas["Sample"] == sample, "Group"].values
        label.set_color(palette.get(grp[0], "black") if len(grp) else "black")

    legend_patches = [mpatches.Patch(color=CASCADE_COLORS[s], label=s) for s in steps]
    ax.legend(handles=legend_patches, loc="upper center",
              bbox_to_anchor=(0.5, -0.12), ncol=6, fontsize=9,
              framealpha=0.9, edgecolor="grey")


def plot_heatmap(ax, df: pd.DataFrame, palette: dict, display_labels: dict):
    present = [c for c in HEATMAP_COLS if c in df.columns]
    labels  = [HEATMAP_LABELS[HEATMAP_COLS.index(c)] for c in present]

    matrix  = df[present].to_numpy(dtype=float)
    samples = df["Sample"].tolist()
    groups  = df["Group"].tolist()

    im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0,
                   vmax=max(matrix.max(), 1))

    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels([display_labels.get(s, s) for s in samples], fontsize=9)
    # imshow already places row 0 at the top -- no invert_yaxis() needed here
    ax.set_title("Classification & transRNA metrics (%)", fontsize=12, fontweight="bold")

    # Zipped against the real sample list (not label.get_text(), which is
    # now the display name, not the real Sample ID the Group lookup needs).
    for label, sample in zip(ax.get_yticklabels(), samples):
        grp = df.loc[df["Sample"] == sample, "Group"].values
        label.set_color(palette.get(grp[0], "black") if len(grp) else "black")

    for i in range(len(samples)):
        for j in range(len(present)):
            val = matrix[i, j]
            text_col = "white" if val > matrix.max() * 0.6 else "black"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=8, color=text_col, fontweight="bold")

    plt.colorbar(im, ax=ax, shrink=0.8, label="%")

    groups_present = sorted(set(groups))
    legend_patches = [mpatches.Patch(color=palette.get(g, "grey"), label=g) for g in groups_present]
    ax.legend(handles=legend_patches, loc="upper center",
              bbox_to_anchor=(0.5, -0.1), ncol=min(len(groups_present), 6), fontsize=8,
              framealpha=0.9, edgecolor="grey")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("report", help="pipeline_read_counts_report.tsv")
    ap.add_argument("--samplesheet", required=True, help="sample,group,fastq_1,fastq_2 sample sheet")
    ap.add_argument("output", nargs="?", default=None)
    args = ap.parse_args()

    report_path = Path(args.report)
    # Base path, no extension -- both real, editable vector formats get
    # written below (svg.fonttype="none" above keeps SVG text as live
    # <text>, not outlined paths; PDF is matplotlib-native, no extra step).
    out_base = Path(args.output) if args.output else \
               report_path.parent / f"{TODAY}_pipeline_report_summary"
    out_base = Path(str(out_base).removesuffix(".svg").removesuffix(".pdf"))
    heatmap_base = out_base.parent / out_base.name.replace("_summary", "_heatmap")

    sample_to_group = load_samplesheet(Path(args.samplesheet))
    df = load_report(report_path, sample_to_group)
    cas = build_cascade(df)

    groups_present = sorted(set(df["Group"]) - {"unknown"})
    palette = group_palette(groups_present)
    palette["unknown"] = "grey"
    display_labels = build_display_labels(df)

    def save(fig, base: Path):
        svg_path = base.with_suffix(".svg")
        pdf_path = base.with_suffix(".pdf")
        fig.savefig(svg_path, dpi=200, bbox_inches="tight")
        fig.savefig(pdf_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {svg_path}")
        print(f"Saved: {pdf_path}")

    # ── Figure 1: read filtering cascade (own file) ───────────────────────────
    fig1, ax1 = plt.subplots(figsize=(10, max(4, len(df) * 0.55)))
    plot_cascade(ax1, cas, palette, display_labels)
    fig1.suptitle("transRNA Pipeline — Read Count & Classification Summary",
                  fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    save(fig1, out_base)

    # ── Figure 2: classification & transRNA metrics heatmap (own file) ────────
    fig2, ax2 = plt.subplots(figsize=(8, max(4, len(df) * 0.5)))
    plot_heatmap(ax2, df, palette, display_labels)
    plt.tight_layout()
    save(fig2, heatmap_base)


if __name__ == "__main__":
    main()
