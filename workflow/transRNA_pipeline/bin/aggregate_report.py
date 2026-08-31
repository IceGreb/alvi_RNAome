#!/usr/bin/env python3
"""
aggregate_report.py
────────────────────
Produces pipeline_read_counts_report.tsv, dataset_summary_report.tsv,
and virus_exclusion_report.tsv from staged pipeline output files.

Column sources (per sample, in report order):
  Raw reads                                <- raw seqkit stats (both mates summed)
  Trimmed                                  <- trimmed seqkit stats (both mates summed)
  STAR mapped %                            <- Log.final.out (unique + multi-mapped %)
  STAR unmapped                            <- star_unmapped seqkit stats (both mates)
  BBSplit: {ref} matched %                 <- one column per reference genome BBSplit was run
                                              with, read straight from its own refstats.txt
                                              (%unambiguousReads + %ambiguousReads for that
                                              reference, exactly as BBSplit reports them)
  MAGs matched %                           <- (virus_survivors - mags_survivors) / virus_survivors * 100
  Candidate transRNAs total reads          <- mags_survivors (reads entering BLAST/Kraken)
  Kraken classified %                      <- kraken C-reads / mags_survivors * 100
  BLAST classified %                       <- blast unique query IDs / (mags_survivors - kraken_classified_reads) * 100
                                              denominator = Kraken-unclassified individual reads passed to BLAST
                                              (kraken counts doubled at ingestion: pairs → individual reads)
  Classified total %                       <- sum of above two
  Total transRNAs after all filters %      <- merged_reads / t_both * 100
                                              merged_reads = reads in blast+kraken pool before dedup
  Total transRNAs after all filters (RPM)  <- merged_reads / t_both * 1e6
  Transmissible RNA representative seqs %  <- reinflated_ge5 / t_both * 100
                                              reinflated_ge5 = sum of dup weights for >=5-dup sequences
  Transmissible RNA representative seqs (RPM) <- reinflated_ge5 / t_both * 1e6
  % of transRNAs with >=5 duplicates       <- reinflated_ge5 / merged_reads * 100
"""

import csv
import glob
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

SEQKIT_PATTERNS = {
    "raw":          re.compile(r"^(.+)_raw_stats\.tsv$"),
    "trimmed":      re.compile(r"^(.+)_trimmed_stats\.tsv$"),
    "star_unmapped":re.compile(r"^(.+)_star_unmapped_stats\.tsv$"),
    "mags":         re.compile(r"^(.+)_noMAGs_stats\.tsv$"),
}
COLLAPSE_RE  = re.compile(r"^(.+)_collapse_stats\.tsv$")
WIDS_RE      = re.compile(r"^(.+)_ge5_detected_weighted_ids\.tsv$")
STARLOG_RE   = re.compile(r"^(.+)_Log\.final\.out$")
REFSTATS_RE  = re.compile(r"^(.+)_bbsplit_refstats\.txt$")
MATE_RE      = re.compile(r"_[12]$")


def load_samplesheet(path: str) -> dict:
    """Return {sample: group} from a sample,group,fastq_1,fastq_2 sample sheet."""
    mapping = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mapping[row["sample"].strip()] = row["group"].strip()
    return mapping


def infer_group(sample: str, sample_to_group: dict) -> str:
    # Some callers pass a mate-suffixed field (e.g. "RJ1_1"); strip it before lookup.
    base = MATE_RE.sub("", sample)
    return sample_to_group.get(base, "unknown")


def pct(numerator, denominator, decimals=2):
    if denominator and denominator > 0:
        return round(numerator / denominator * 100, decimals)
    return 0.0


def rpm(numerator, denominator, decimals=2):
    if denominator and denominator > 0:
        return round(numerator / denominator * 1_000_000, decimals)
    return 0.0


def seqkit_total_reads(path: Path) -> int:
    """Return total read count across all rows (both mates) from a seqkit stats TSV."""
    try:
        df = pd.read_csv(path, sep="\t")
        return int(df["num_seqs"].sum())
    except Exception as e:
        print(f"Warning seqkit {path}: {e}", file=sys.stderr)
        return 0


def parse_star_log(path: Path) -> float:
    """Return total mapped % (unique + multi-mapped) from STAR Log.final.out."""
    unique_pct = 0.0
    multi_pct  = 0.0
    try:
        with path.open() as f:
            for line in f:
                if "Uniquely mapped reads %" in line:
                    unique_pct = float(line.split("|")[1].strip().rstrip("%"))
                elif "% of reads mapped to multiple loci" in line:
                    multi_pct  = float(line.split("|")[1].strip().rstrip("%"))
    except Exception as e:
        print(f"Warning STAR log {path}: {e}", file=sys.stderr)
    return round(unique_pct + multi_pct, 2)


def parse_refstats(path: Path) -> dict:
    """Return {ref_name: matched_pct} from a BBSplit refstats.txt.

    matched_pct = %unambiguousReads + %ambiguousReads, i.e. BBSplit's own
    reported % of input reads assigned to that reference — used as-is, no
    re-derivation. ref_name is BBSplit's own '#name' column, unmodified.
    """
    matched = {}
    try:
        df = pd.read_csv(path, sep="\t")
        df.columns = [c.lstrip("#") for c in df.columns]
        for _, row in df.iterrows():
            ref = str(row["name"])
            matched[ref] = float(row["%unambiguousReads"]) + float(row["%ambiguousReads"])
    except Exception as e:
        print(f"Warning refstats {path}: {e}", file=sys.stderr)
    return matched


def weighted_total(path: Path) -> int:
    total = 0
    try:
        with path.open() as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        total += int(parts[2])
                    except ValueError:
                        pass
    except Exception as e:
        print(f"Warning wids {path}: {e}", file=sys.stderr)
    return total


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--output",         required=True)
    ap.add_argument("--summary-output", default="dataset_summary_report.tsv")
    ap.add_argument("--virus-output",   default="virus_exclusion_report.tsv")
    ap.add_argument("--samplesheet",    required=True,
                     help="sample,group,fastq_1,fastq_2 sample sheet")
    args = ap.parse_args()

    sample_to_group = load_samplesheet(args.samplesheet)

    seqkit_data     = defaultdict(dict)
    trimmed_both    = {}   # total reads both mates (denominator for % columns)
    star_mapped_pct = {}
    collapse_total  = {}
    collapse_ge5    = {}
    merged_reads    = {}   # reads in merged.fq before rmdup
    kraken_counts   = defaultdict(int)
    blast_unique    = defaultdict(int)
    final_trans     = defaultdict(int)
    all_samples     = set()
    virus_rows      = []
    refstats_data   = {}   # {sample: {ref_name: matched_pct}}
    all_refs        = set()

    # ── seqkit stats (all counts use both-mate totals for consistency) ────────
    for fname in sorted(glob.glob("*.tsv")):
        for step, pattern in SEQKIT_PATTERNS.items():
            m = pattern.match(fname)
            if m:
                sample = m.group(1)
                all_samples.add(sample)
                seqkit_data[sample][step] = seqkit_total_reads(Path(fname))
                if step == "trimmed":
                    trimmed_both[sample] = seqkit_total_reads(Path(fname))
                break

    # ── STAR Log.final.out → mapped % ────────────────────────────────────────
    for fname in sorted(glob.glob("*_Log.final.out")):
        m = STARLOG_RE.match(fname)
        if m:
            sample = m.group(1)
            all_samples.add(sample)
            star_mapped_pct[sample] = parse_star_log(Path(fname))

    # ── BBSplit refstats → per-reference matched % ────────────────────────────
    for fname in sorted(glob.glob("*_bbsplit_refstats.txt")):
        m = REFSTATS_RE.match(fname)
        if m:
            sample = m.group(1)
            all_samples.add(sample)
            refs = parse_refstats(Path(fname))
            refstats_data[sample] = refs
            all_refs.update(refs.keys())

    # ── collapse stats ────────────────────────────────────────────────────────
    for fname in sorted(glob.glob("*_collapse_stats.tsv")):
        m = COLLAPSE_RE.match(fname)
        if m:
            sample = m.group(1)
            all_samples.add(sample)
            try:
                df = pd.read_csv(fname, sep="\t")
                row = df[df["mate"].astype(str) == "1"]
                if not row.empty:
                    merged_reads[sample]   = int(row["merged_reads"].iloc[0])
                    collapse_total[sample] = int(row["total_groups"].iloc[0])
                    collapse_ge5[sample]   = int(row["ge5_groups"].iloc[0])
            except Exception as e:
                print(f"Warning collapse {fname}: {e}", file=sys.stderr)

    # ── filter stats: kraken + blast unique queries + virus exclusions ────────
    for fname in sorted(glob.glob("*_filter_stats.tsv")):
        try:
            df_s = pd.read_csv(fname, sep="\t")
            for _, row in df_s.iterrows():
                raw_sample = str(row.get("sample", "")).strip()
                tool       = str(row.get("tool",   "")).strip()
                rm_v       = row.get("rm_virus", 0)

                if tool == "kraken":
                    sample = raw_sample
                    all_samples.add(sample)
                    # Each Kraken line = 1 read pair; ×2 converts to individual reads
                    # to match mags (seqkit, both mates) and blast (unique_queries per mate).
                    kraken_counts[sample] += int(row.get("classified", 0)) * 2
                elif tool == "blast":
                    # Blast filter_stats sample field is "RJ1_1" / "RJ1_2";
                    # strip mate suffix to match the base sample name.
                    sample = MATE_RE.sub("", raw_sample)
                    all_samples.add(sample)
                    blast_unique[sample] += int(row.get("unique_queries", 0))

                virus_rows.append({
                    "Dataset": infer_group(raw_sample, sample_to_group),
                    "Sample":  raw_sample,
                    "Tool":    tool,
                    "Virus reads excluded": rm_v,
                })
        except Exception as e:
            print(f"Warning filter stats {fname}: {e}", file=sys.stderr)

    # ── weighted_ids: reinflated transRNA total (informational) ───────────────
    for fname in sorted(glob.glob("*_ge5_detected_weighted_ids.tsv")):
        m = WIDS_RE.match(fname)
        if m:
            sample = m.group(1)
            all_samples.add(sample)
            final_trans[sample] += weighted_total(Path(fname))

    # ── Build main report ─────────────────────────────────────────────────────
    sorted_refs = sorted(all_refs)
    rows = []
    for sample in sorted(all_samples):
        sd       = seqkit_data.get(sample, {})
        raw      = sd.get("raw",           0)
        trimmed  = sd.get("trimmed",       0)  # both mates total
        t_both   = trimmed_both.get(sample, trimmed)
        unmapped = sd.get("star_unmapped", 0)
        mags     = sd.get("mags",  0)
        krak     = kraken_counts.get(sample, 0)
        blast    = blast_unique.get(sample,  0)
        c_total    = collapse_total.get(sample, 0)
        c_ge5      = collapse_ge5.get(sample,   0)
        mreads     = merged_reads.get(sample,   0)
        reinflated = final_trans.get(sample,    0)
        refs       = refstats_data.get(sample, {})

        row = {
            "Sample":                                        sample,
            "Group":                                         infer_group(sample, sample_to_group),
            "Raw reads":                                     raw,
            "Trimmed":                                       trimmed,
            "STAR mapped %":                                 star_mapped_pct.get(sample, 0.0),
            "STAR unmapped":                                 unmapped,
        }
        for ref in sorted_refs:
            row[f"BBSplit: {ref} matched %"] = refs.get(ref, 0.0)
        row.update({
            # Denominator is STAR-unmapped (not a "virus survivors" count, which
            # no longer exists now that host/virus stats come from BBSplit's own
            # refstats above) -- % of the STAR-unmapped pool removed as MAGs.
            "MAGs matched %":                                pct(unmapped - mags, unmapped),
            # All three variables are in individual-read units:
            # mags  : seqkit stats on *_noMAGs_stats.tsv (both mates summed)
            # krak  : *_kraken_filter_stats.tsv → 'classified' × 2  (pairs converted to reads at ingestion)
            # blast : *_blast_filter_stats.tsv → 'unique_queries' mate1 + mate2
            # mags - krak = Kraken-unclassified individual reads passed to BLAST
            "Candidate transRNAs total reads":               mags,
            "Kraken classified %":                           pct(krak,         mags),
            "BLAST classified %":                            pct(blast,        mags - krak),
            "Classified total %":                            pct(krak + blast, mags),
            "Total transRNAs after all filters %":               pct(mreads,     t_both),
            "Total transRNAs after all filters (RPM)":           rpm(mreads,     t_both),
            "Transmissible RNA representative sequences %":      pct(reinflated, t_both),
            "Transmissible RNA representative sequences (RPM)":  rpm(reinflated, t_both),
            "% of transRNAs with >=5 duplicates":               pct(reinflated, mreads),
        })
        rows.append(row)

    df = pd.DataFrame(rows)
    df.drop(columns=["Group"]).to_csv(args.output, sep="\t", index=False)
    print(f"Report: {args.output} ({len(df)} samples)", file=sys.stderr)

    # ── Dataset summary ───────────────────────────────────────────────────────
    summary_rows = []
    for grp in sorted(df["Group"].unique()):
        sub = df[df["Group"] == grp]
        if sub.empty:
            continue
        trimmed_tot = sub["Trimmed"].sum()
        summary_row = {
            "Dataset":                                       grp,
            "n_samples":                                     len(sub),
            "Trimmed (both mates total)":                     int(trimmed_tot),
            "Avg STAR mapped %":                             round(sub["STAR mapped %"].mean(), 2),
        }
        for ref in sorted_refs:
            col = f"BBSplit: {ref} matched %"
            summary_row[f"Avg {col}"] = round(sub[col].mean(), 2)
        summary_row.update({
            "Avg MAGs matched %":                            round(sub["MAGs matched %"].mean(), 2),
            "Avg Kraken classified %":                       round(sub["Kraken classified %"].mean(), 2),
            "Avg BLAST classified %":                        round(sub["BLAST classified %"].mean(), 2),
            "Avg Classified total %":                        round(sub["Classified total %"].mean(), 2),
            "Avg Total transRNAs after all filters %":           round(sub["Total transRNAs after all filters %"].mean(), 2),
            "Avg Total transRNAs after all filters (RPM)":       round(sub["Total transRNAs after all filters (RPM)"].mean(), 2),
            "Avg Transmissible RNA representative seqs %":       round(sub["Transmissible RNA representative sequences %"].mean(), 2),
            "Avg Transmissible RNA representative seqs (RPM)":   round(sub["Transmissible RNA representative sequences (RPM)"].mean(), 2),
            "Avg % of transRNAs with >=5 duplicates":            round(sub["% of transRNAs with >=5 duplicates"].mean(), 2),
        })
        summary_rows.append(summary_row)
    pd.DataFrame(summary_rows).to_csv(args.summary_output, sep="\t", index=False)
    print(f"Summary: {args.summary_output}", file=sys.stderr)

    # ── Virus exclusion report ────────────────────────────────────────────────
    if virus_rows:
        pd.DataFrame(virus_rows).to_csv(args.virus_output, sep="\t", index=False)
        print(f"Virus report: {args.virus_output}", file=sys.stderr)


if __name__ == "__main__":
    main()
