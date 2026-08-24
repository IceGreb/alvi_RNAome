#!/usr/bin/env python3
"""
plot_lengths.py
───────────────
Reads *_1_ge5_detected.hist files (weighted length histograms produced by the
ge5 script) and plots the average reinflated length distribution per sample
group over the 15–200 nt window.

Hist file format (no header):
    length<TAB>weighted_count

Only mate-1 hist files are used (*_1_ge5_detected.hist); both mates carry
identical content in the merged-pool pipeline, so mate 2 is skipped to avoid
double-counting.

Sample grouping is read from the sample sheet (--samplesheet), not inferred
from the sample name — any number of groups, with any labels, is supported.

Output (one plot):
    <date>_avg_reinflated_length_distribution_15to100nt.png
"""

import csv
import matplotlib
matplotlib.use("Agg")

import sys
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLOT_MIN = 15
PLOT_MAX = 100
HIGHLIGHT_TICKS = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]

# Cycled through in order for however many groups are present.
GROUP_COLOR_POOL = [
    "#00a4ff", "#ff0000", "#2ca02c", "#9467bd",
    "#ff7f0e", "#17becf", "#e377c2", "#8c564b",
]


def load_samplesheet(path: Path) -> dict:
    """Return {sample: group} from a sample,group,fastq_1,fastq_2 sample sheet."""
    mapping = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mapping[row["sample"].strip()] = row["group"].strip()
    return mapping


def load_hist(path: Path) -> pd.DataFrame:
    """Load a length-histogram file and return a percent-normalised DataFrame."""
    df = pd.read_csv(path, sep="\t", header=None, names=["length", "count"])
    df = df.dropna()
    df["length"] = df["length"].astype(int)
    df["count"]  = pd.to_numeric(df["count"], errors="coerce").fillna(0)
    total = df["count"].sum()
    if total <= 0:
        return pd.DataFrame(columns=["length", "percent"])
    df["percent"] = df["count"] / total * 100.0
    return df[["length", "percent"]]


def sample_name(path: Path) -> str:
    """Strip '_ge5_detected.hist' suffix to get the base sample name."""
    return path.name.replace("_ge5_detected.hist", "")


def align_and_mean(dfs: list, bins: np.ndarray) -> np.ndarray:
    """Align each sample to a shared bin axis and return the mean curve."""
    if not dfs:
        return np.zeros(len(bins))
    aligned = []
    for df in dfs:
        # Restrict to window and re-normalise within it
        sub = df[df["length"].between(PLOT_MIN, PLOT_MAX)].copy()
        sub_total = sub["percent"].sum()
        if sub_total > 0:
            sub["percent"] = sub["percent"] / sub_total * 100.0
        s = sub.set_index("length")["percent"].reindex(bins, fill_value=0.0)
        aligned.append(s.to_numpy())
    return np.vstack(aligned).mean(axis=0)


def make_plot(means: dict, bins: np.ndarray, grouped: dict, group_colors: dict, outfile: Path):
    fig, ax = plt.subplots(figsize=(12, 7))

    # Alternating column shading
    for i, b in enumerate(bins):
        if i % 2 == 1:
            ax.axvspan(b - 0.5, b + 0.5, color="gray", alpha=0.08, zorder=0)

    # Vertical dotted guide lines
    for x in HIGHLIGHT_TICKS:
        if PLOT_MIN <= x <= PLOT_MAX:
            ax.axvline(x=x, color="gray", linestyle=":", linewidth=1.8,
                       alpha=1.0, zorder=1)

    for group in sorted(grouped):
        if grouped[group]:
            ax.plot(bins, means[group], color=group_colors[group], lw=2,
                     label=group, zorder=3)

    ax.set_xlabel("Sequence length (nt)", fontsize=13)
    ax.set_ylabel("Weighted percentage of reads (%)", fontsize=13)
    ax.set_title(f"Average Length Distribution ({PLOT_MIN}–{PLOT_MAX} nt)",
                 fontsize=15)
    ax.legend(frameon=True, facecolor="white")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.set_ylim(bottom=0)
    ax.set_xlim(PLOT_MIN - 0.5, PLOT_MAX + 0.5)

    ticks = sorted({PLOT_MIN} | {t for t in HIGHLIGHT_TICKS if PLOT_MIN <= t <= PLOT_MAX})
    ax.set_xticks(ticks)
    ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved: {outfile}", file=sys.stderr)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("hist_dir", nargs="?", default=".")
    parser.add_argument("--samplesheet", required=True,
                         help="sample,group,fastq_1,fastq_2 sample sheet")
    args = parser.parse_args()

    hist_dir = Path(args.hist_dir)
    today    = date.today().strftime("%d_%m_%Y")

    sample_to_group = load_samplesheet(args.samplesheet)
    groups = sorted(set(sample_to_group.values()))
    group_colors = {g: GROUP_COLOR_POOL[i % len(GROUP_COLOR_POOL)]
                     for i, g in enumerate(groups)}

    hist_files = sorted(hist_dir.glob("*_ge5_detected.hist"))
    if not hist_files:
        print(f"No *_1_ge5_detected.hist files found in {hist_dir}", file=sys.stderr)
        sys.exit(1)

    grouped = {g: [] for g in groups}

    for f in hist_files:
        name  = sample_name(f)
        group = sample_to_group.get(name)
        if group is None:
            print(f"Skipping {f.name}: sample '{name}' not found in {args.samplesheet}",
                  file=sys.stderr)
            continue

        df = load_hist(f)
        if df.empty or df["percent"].sum() == 0:
            print(f"Skipping {f.name}: empty or zero-weight histogram",
                  file=sys.stderr)
            continue

        grouped[group].append(df)
        print(f"Loaded {f.name} → group={group}, "
              f"n_bins={len(df)}, sum%={df['percent'].sum():.4f}",
              file=sys.stderr)

    if not any(grouped.values()):
        print("No valid hist files found for any group.", file=sys.stderr)
        sys.exit(1)

    bins  = np.arange(PLOT_MIN, PLOT_MAX + 1)
    means = {label: align_and_mean(grouped[label], bins) for label in groups}

    for label in groups:
        if grouped[label]:
            print(f"{label}: n={len(grouped[label])}, "
                  f"mean curve sum={means[label].sum():.4f}%",
                  file=sys.stderr)

    outfile = hist_dir / f"{today}_avg_reinflated_length_distribution_15to100nt.png"
    make_plot(means, bins, grouped, group_colors, outfile)


if __name__ == "__main__":
    main()
