#!/usr/bin/env python3
import csv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import colorsys
from datetime import date
import matplotlib.ticker as mticker
from matplotlib import gridspec

# ---------- Fixed special colors ----------
OTHER_COLOR = "#c7c7c7"         # Other
UNCLASSIFIED_COLOR = "#7f7f7f"  # Unclassified (a.k.a. Unassigned)

# Base palette pool for non-special taxa (no greys)
BASE_COLOR_POOL = [
    "#ff7f0e",
    "#1f77b4", "#2b9c2b", "#dbdb8d", "#c49c94",
    "#aec7e8", "#c5b0d5", "#bcbd22", "#d62728", "#8c564b", "#f7b6d2",
    "#2ca02c", "#ffbb78",
    "#ff9896", "#9edae5",
    "#17becf", "#e377c2", "#98df8a"
]

# Suggestive base colors for some common high-level taxa (kept if present in that rank)
BASE_TAXA_COLORS = {
    "eukaryota": "#1f77b4",
    "bacteria": "#ff7f0e",
    "archaea": "#2ca02c",
    "viruses": "#d62728",
}

RANKS = ["Domain", "Kingdom", "Order", "Species"]

# ------------------ Utilities ------------------

def _norm_key(name: str) -> str:
    return str(name).strip().lower()

def _norm_series(series):
    return series.astype(str).str.strip().str.lower()

def _display_name(taxon: str) -> str:
    low = _norm_key(taxon)
    if low == "other":
        return "Other"
    if low in ("unassigned", "unclassified"):
        return "Unclassified"
    return str(taxon)

def hex_to_rgb01(hex_str):
    hex_str = hex_str.lstrip("#")
    r = int(hex_str[0:2], 16)/255.0
    g = int(hex_str[2:4], 16)/255.0
    b = int(hex_str[4:6], 16)/255.0
    return (r, g, b)

def rgb01_to_hex(rgb):
    r = max(0, min(255, int(round(rgb[0]*255))))
    g = max(0, min(255, int(round(rgb[1]*255))))
    b = max(0, min(255, int(round(rgb[2]*255))))
    return f"#{r:02x}{g:02x}{b:02x}"

def vary_color(hex_color, step_idx, total_needed):
    r, g, b = hex_to_rgb01(hex_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    if total_needed <= 1:
        return hex_color
    l_offsets = [0.15, -0.10, 0.25, -0.20, 0.05, -0.05, 0.30, -0.25]
    s_offsets = [0.10, -0.10, 0.05, -0.05, 0.15, -0.15, 0.08, -0.08]
    l_new = max(0.20, min(0.85, l + l_offsets[step_idx % len(l_offsets)]))
    s_new = max(0.35, min(1.00, s + s_offsets[step_idx % len(s_offsets)]))
    r2, g2, b2 = colorsys.hls_to_rgb(h, l_new, s_new)
    return rgb01_to_hex((r2, g2, b2))

def expand_palette(required_n, base_pool):
    out = []
    if required_n <= len(base_pool):
        return base_pool[:required_n]
    out.extend(base_pool)
    remaining = required_n - len(base_pool)
    idx = 0
    while remaining > 0:
        base = base_pool[idx % len(base_pool)]
        variant = vary_color(base, step_idx=idx, total_needed=required_n)
        if variant.lower() not in [c.lower() for c in out]:
            out.append(variant)
            remaining -= 1
        idx += 1
        if idx > 10000:
            break
    return out

# ------------------ IO & Merge ------------------

def read_top10_files(top10_dir):
    files = glob.glob(os.path.join(top10_dir, "*_top10.tsv"))
    data = {}
    for f in files:
        base = os.path.basename(f).replace("_top10.tsv", "")
        parts = base.split("_")
        dataset = "_".join(parts[:-1])   # e.g., RJ_long... or ST_short...
        rank = parts[-1]                 # e.g., Domain/Kingdom/Order/Species
        df = pd.read_csv(f, sep="\t")
        if dataset not in data:
            data[dataset] = {}
        data[dataset][rank] = df
    return data

def load_samplesheet_groups(path):
    """Return the distinct group labels from a sample,group,... sample sheet,
    longest-first so prefix matching below is unambiguous."""
    groups = set()
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            groups.add(row["group"].strip())
    return sorted(groups, key=len, reverse=True)

def simplify_dataset_name(ds_name, known_groups=None):
    """Recover the group label a top10 filename's leading tokens belong to.

    Filenames are "{group}_{tool}_{rank}_top10.tsv"; ds_name is everything
    before "_{rank}". known_groups (from the sample sheet) lets this match
    group labels that themselves contain underscores; without it, fall back
    to the first underscore-delimited token.
    """
    if known_groups:
        for g in known_groups:
            if ds_name == g or ds_name.startswith(g + "_"):
                return g
    return ds_name.split("_", 1)[0]

def merge_long_short(data, known_groups=None):
    merged = {}
    buckets = {}
    for ds_full, ranks_map in data.items():
        ds_simple = simplify_dataset_name(ds_full, known_groups)
        for rank, df in ranks_map.items():
            if {"Taxonomy", "Count"} - set(df.columns):
                continue
            buckets.setdefault((ds_simple, rank), []).append(df[["Taxonomy", "Count"]].copy())

    for (ds_simple, rank), df_list in buckets.items():
        if not df_list:
            continue
        cat = pd.concat(df_list, ignore_index=True)
        agg = cat.groupby("Taxonomy", as_index=False)["Count"].sum()
        total = agg["Count"].sum()
        agg["Percentage"] = (agg["Count"] / total * 100.0) if total > 0 else 0.0
        agg["Display"] = agg["Taxonomy"].apply(_display_name)
        merged.setdefault(ds_simple, {})[rank] = agg
    return merged

# ------------------ Threshold + Truncate ------------------

def apply_threshold_other_and_truncate(df, rank, pct_threshold=0.1, max_normals=10):
    if df is None or df.empty:
        return df
    df = df.copy()
    tiny_mask = df["Percentage"] < pct_threshold
    tiny_sum = df.loc[tiny_mask, "Count"].sum()
    major = df.loc[~tiny_mask, ["Taxonomy", "Count"]].copy()
    if tiny_sum > 0:
        if (_norm_series(major["Taxonomy"]) == "other").any():
            major.loc[_norm_series(major["Taxonomy"]) == "other", "Count"] += tiny_sum
        else:
            major = pd.concat([major, pd.DataFrame({"Taxonomy": ["other"], "Count": [tiny_sum]})],
                              ignore_index=True)
    if rank == "Domain":
        if (_norm_series(major["Taxonomy"]) == "other").any():
            other_count = major.loc[_norm_series(major["Taxonomy"]) == "other", "Count"].sum()
            if (_norm_series(major["Taxonomy"]) == "archaea").any():
                major.loc[_norm_series(major["Taxonomy"]) == "archaea", "Count"] += other_count
            else:
                major = pd.concat([major, pd.DataFrame({"Taxonomy": ["Archaea"], "Count": [other_count]})],
                                  ignore_index=True)
            major = major.loc[_norm_series(major["Taxonomy"]) != "other"].reset_index(drop=True)
    total = major["Count"].sum()
    major["Percentage"] = (major["Count"] / total * 100.0) if total > 0 else 0.0
    major["Display"] = major["Taxonomy"].apply(_display_name)

    is_other = _norm_series(major["Taxonomy"]) == "other"
    is_uncls = _norm_series(major["Taxonomy"]).isin(["unassigned", "unclassified"])
    specials = major.loc[is_other | is_uncls].copy()
    normals = major.loc[~(is_other | is_uncls)].copy()

    normals = normals.sort_values(by="Display", key=lambda s: s.str.lower())
    if len(normals) > max_normals:
        normals = normals.head(max_normals)

    def special_key(t):
        t = _norm_key(t)
        if t == "other": return 1
        if t in ("unassigned", "unclassified"): return 2
        return 0
    specials = specials.sort_values(by="Taxonomy", key=lambda s: s.map(special_key))

    out = pd.concat([normals, specials], ignore_index=True)
    total2 = out["Count"].sum()
    out["Percentage"] = (out["Count"] / total2 * 100.0) if total2 > 0 else 0.0

    return out.reset_index(drop=True)

# ------------------ Rank-specific palettes ------------------

def build_rank_palettes(data_by_ds_rank):
    rank_palettes = {}
    for rank in RANKS:
        taxa = []
        for ds in data_by_ds_rank.keys():
            df = data_by_ds_rank[ds].get(rank)
            if df is None or df.empty:
                continue
            taxa.extend(df["Taxonomy"].tolist())

        taxa_set = sorted({t for t in taxa}, key=lambda s: s.lower() if isinstance(s, str) else str(s))
        specials = [t for t in taxa_set if _norm_key(t) in ("other", "unassigned", "unclassified")]
        normals = [t for t in taxa_set if _norm_key(t) not in ("other", "unassigned", "unclassified")]

        assigned = {}
        used_colors = set()
        if any(_norm_key(t) == "other" for t in specials):
            assigned["other"] = OTHER_COLOR; used_colors.add(OTHER_COLOR)
        if any(_norm_key(t) in ("unassigned", "unclassified") for t in specials):
            assigned["unclassified"] = UNCLASSIFIED_COLOR; used_colors.add(UNCLASSIFIED_COLOR)

        for t in normals:
            key = _norm_key(t)
            if key in BASE_TAXA_COLORS:
                col = BASE_TAXA_COLORS[key]
                if col not in used_colors:
                    assigned[key] = col; used_colors.add(col)

        remaining_normals = [t for t in normals if _norm_key(t) not in assigned]
        need = len(remaining_normals)
        expanded_pool = expand_palette(need, [c for c in BASE_COLOR_POOL if c not in used_colors])
        for t, col in zip(remaining_normals, expanded_pool):
            key = _norm_key(t)
            assigned[key] = col
            used_colors.add(col)

        rank_palette = dict(assigned)
        if "other" in assigned:
            rank_palette["Other"] = assigned["other"]
        if "unclassified" in assigned:
            rank_palette["Unclassified"] = assigned["unclassified"]

        rank_palettes[rank] = rank_palette
    return rank_palettes

# ------------------ Plot (aligned & VERTICAL rank titles) ------------------

def plot_multi_panel(data, rank_palettes, datasets, ranks, outfile, x_break_point=2.5):
    """
    Generates a multi-panel plot with horizontal bar charts and a broken x-axis.
    Ensures within each rank, color uniqueness across datasets.
    Row titles (Domain/Kingdom/Order/Species) are aligned left and drawn VERTICALLY.
    """
    n_rows = len(ranks)
    n_cols = len(datasets)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 3.5 * n_rows), squeeze=False)

    for i, rank in enumerate(ranks):
        palette = rank_palettes.get(rank, {})
        for j, dataset in enumerate(datasets):
            placeholder_ax = axes[i, j]
            placeholder_ax.axis('off')

            df = data.get(dataset, {}).get(rank)
            if df is None or df.empty:
                continue

            gs = gridspec.GridSpecFromSubplotSpec(
                1, 2, subplot_spec=placeholder_ax.get_subplotspec(),
                width_ratios=[1, 4], wspace=0.06
            )
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1], sharey=ax1)

            d = 0.003
            kw = dict(transform=ax1.transAxes, color='k', clip_on=False, zorder=1)
            ax1.plot((1 - d, 1 + d), (-d, +d), **kw)
            ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kw)
            kw.update(transform=ax2.transAxes)
            ax2.plot((-d, +d), (-d, +d), **kw)
            ax2.plot((-d, +d), (1 - d, 1 + d), **kw)

            df = df.reset_index(drop=True)
            for idx, row in df.iterrows():
                taxon = row["Taxonomy"]
                pct = float(row["Percentage"])
                key = _norm_key(taxon)

                if key == "other":
                    color = OTHER_COLOR
                elif key in ("unassigned", "unclassified"):
                    color = UNCLASSIFIED_COLOR
                else:
                    color = palette.get(key, "#999999")

                ax1.barh(idx, min(pct, x_break_point), color=color, height=0.8, zorder=2)
                if pct > x_break_point:
                    ax2.barh(idx, pct - x_break_point, left=x_break_point, color=color, height=0.8, zorder=2)

                ax2.text(min(max(pct, x_break_point) + 0.5, 99.0), idx, f"{pct:.2f}%",
                         va="center", fontsize=7, zorder=10, clip_on=False)

            ax1.set_xlim(0, x_break_point)
            ax2.set_xlim(x_break_point, 100)
            ax1.invert_yaxis()

            ax1.spines['right'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            ax2.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

            ax1.grid(False); ax2.grid(False)
            ax2.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=5, prune='lower'))

            ax1.set_yticks(range(len(df)))
            ax1.set_yticklabels(df["Display"], fontsize=9)
            ax1.tick_params(axis='y', length=0)

            # (do not set per-axes ylabel)

            if i == 0:
                ax1.set_title(dataset, x=2.5, y=1.05, fontsize=14)
            if i == len(ranks) - 1:
                ax1.set_xlabel(" ")
                ax2.set_xlabel(" ")

    # Draw one VERTICAL rank label per row at a fixed left x, centered vertically
    fig.canvas.draw()
    left_x = 0.03  # figure coordinates
    for i, rank in enumerate(ranks):
        pos = axes[i, 0].get_position()
        y_center = pos.y0 + pos.height / 2.0
        # rotation=90 to keep the original vertical look
        fig.text(left_x, y_center, rank, va='center', ha='center', fontsize=12, rotation=90)

    fig.supxlabel("Percentage (%)", fontsize=12, y=0.03)
    # leave room on the left for vertical row titles
    plt.tight_layout(rect=[0.085, 0.05, 1, 0.97])
    plt.savefig(outfile, dpi=300)
    plt.close()

# ------------------ Main ------------------

def main(top10_dir, outdir, pct_threshold=0.1, max_normals=10, x_break_point=2.5, samplesheet=None):
    os.makedirs(outdir, exist_ok=True)

    known_groups = load_samplesheet_groups(samplesheet) if samplesheet else None

    data_raw = read_top10_files(top10_dir)
    data_merged = merge_long_short(data_raw, known_groups)

    for ds in list(data_merged.keys()):
        for rank in list(data_merged[ds].keys()):
            df = data_merged[ds][rank]
            data_merged[ds][rank] = apply_threshold_other_and_truncate(
                df, rank, pct_threshold=pct_threshold, max_normals=max_normals
            )

    rank_palettes = build_rank_palettes(data_merged)

    desired_datasets = sorted(data_merged.keys())

    today = date.today().strftime("%d_%m_%Y")
    outfile = os.path.join(outdir, f"{today}_ge5_taxonomy_top10.png")
    plot_multi_panel(data_merged, rank_palettes, desired_datasets, RANKS, outfile, x_break_point=x_break_point)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Plot merged (long+short) multi-panel taxa with broken X-axis; "
                    "<0.1% → Other; Domain 'Other' folded into Archaea; "
                    "max 10 taxa + Other + Unclassified; unique colors per rank."
    )
    parser.add_argument("--top10_dir", required=True, help="Directory with *_top10.tsv files")
    parser.add_argument("--outdir", required=True, help="Directory to save plots")
    parser.add_argument("--threshold", type=float, default=0.1, help="Percentage threshold (default: 0.1)")
    parser.add_argument("--max-normals", type=int, default=10, help="Max normal taxa per panel (default: 10)")
    parser.add_argument("--xbreak", type=float, default=2.5, help="X-axis break point (default: 2.5)")
    parser.add_argument("--samplesheet", required=False,
                         help="sample,group,fastq_1,fastq_2 sample sheet; used to recover group "
                              "labels from top10 filenames when a group name contains underscores")
    args = parser.parse_args()
    main(args.top10_dir, args.outdir, pct_threshold=args.threshold, max_normals=args.max_normals,
         x_break_point=args.xbreak, samplesheet=args.samplesheet)
