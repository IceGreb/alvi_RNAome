#!/usr/bin/env python3

import argparse
import gzip
import re
from pathlib import Path
from collections import Counter, defaultdict

import pandas as pd


RANKS = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

IDS_RE = re.compile(r"^(.*?)_ge5_detected_ids\.txt(?:\.gz)?$")

# Terms from taxonkit that mean "we don't know" at this rank
# These are replaced with "unclassified" in all ranks except Species
NOISE_TERMS = frozenset({
    "unassigned",
    "unclassified",
    "",
})

# At Domain level, these mean the read is from a known domain but
# the term itself should be normalised
DOMAIN_SYNONYMS = {
    "acellular root": "Viruses",
}

def normalise_taxon(value: str, rank: str) -> str:
    """
    Normalise a single taxon value for a given rank.

    Rules:
      - Empty string → "unclassified"
      - Known noise terms → "unclassified"
      - Domain synonyms (acellular root) → canonical name
      - Everything else → kept as-is (including "uncultured bacterium"
        at species rank, which is a valid NCBI species)
    """
    value = value.strip()

    # Empty → unclassified
    if not value:
        return "unclassified"

    # Known noise terms → unclassified (case-insensitive)
    if value.lower() in NOISE_TERMS:
        return "unclassified"

    # Domain-level synonyms
    if rank == "Domain":
        canonical = DOMAIN_SYNONYMS.get(value.lower())
        if canonical:
            return canonical

    return value


def open_maybe_gzip(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def clean_read_id(x: str) -> str:
    x = str(x).strip()
    if x.startswith("@"):
        x = x[1:]
    return x


def load_id_file(path: Path) -> set[str]:
    ids = set()
    with open_maybe_gzip(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rid = clean_read_id(line.split()[0])
            if rid:
                ids.add(rid)
    return ids


def discover_samples(ids_dir: Path) -> dict[str, Path]:
    samples = {}

    for p in list(ids_dir.glob("*_ge5_detected_ids.txt")) + list(ids_dir.glob("*_ge5_detected_ids.txt.gz")):
        m = IDS_RE.match(p.name)
        if not m:
            continue
        sample = m.group(1)
        samples[sample] = p

    return samples


def read_total_reads(reads_table: Path) -> dict[str, int]:
    df = pd.read_csv(reads_table, sep="\t")

    required = {"Sample", "TotalReads"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"{reads_table} is missing columns: {missing}")

    return {
        str(row["Sample"]).strip(): int(row["TotalReads"])
        for _, row in df.iterrows()
    }


def resolve_file(path: Path) -> Path | None:
    if path.exists():
        return path

    gz = Path(str(path) + ".gz")
    if gz.exists():
        return gz

    return None


def load_samplesheet(path: Path) -> dict[str, str]:
    """Return {sample: group} from a sample,group,fastq_1,fastq_2 sample sheet."""
    import csv
    mapping = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mapping[row["sample"].strip()] = row["group"].strip()
    return mapping


def sample_group(sample: str, sample_to_group: dict[str, str]) -> str | None:
    return sample_to_group.get(sample)


def init_counts():
    return defaultdict(Counter)


def add_counts(dest, src):
    for rank in RANKS:
        dest[rank].update(src[rank])


def rpm_normalise_counts(counts, total_reads: int | None):
    """
    Convert raw reinflated counts to RPM counts per rank.

    These RPM counts can then be safely summed across samples.
    """
    out = init_counts()

    if not total_reads or total_reads <= 0:
        for rank in RANKS:
            for taxon, count in counts[rank].items():
                out[rank][taxon] += float(count)
        return out

    for rank in RANKS:
        for taxon, count in counts[rank].items():
            out[rank][taxon] += float(count) / total_reads * 1_000_000.0

    return out


def parse_dup_detail_file(path: Path | None) -> dict[str, int]:
    """
    Parses duplicated.detail files of this format:

        occurrence_count<TAB>seqID1, seqID2, seqID3, ...

    Returns:

        seqID -> duplicate weight

    No filtering is done here, because selected IDs are already filtered upstream.
    """
    weights = defaultdict(int)

    if path is None:
        return {}

    with open_maybe_gzip(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue

            try:
                occ = int(parts[0].strip())
            except ValueError:
                continue

            seqids = [clean_read_id(x.strip()) for x in parts[1].split(",") if x.strip()]

            for rid in seqids:
                if rid:
                    weights[rid] += occ

    return dict(weights)


def get_weight(read_id: str, weights: dict[str, int] | None) -> int:
    if not weights:
        return 1
    return int(weights.get(read_id, 1))


def parse_semicolon_lineage(lineage: str) -> dict:
    """
    Parse a semicolon-delimited lineage string:
      Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species

    Classifies each read to the highest available rank.
    Ranks beyond the last known rank are "unclassified".

    Examples:
      "Eukaryota;Streptophyta;;;" → Domain=Eukaryota, Kingdom=Streptophyta,
                                    Phylum=unclassified, ...
      "Bacteria;;;;;;uncultured bacterium" → Domain=Bacteria,
                                    Kingdom=unclassified, ...,
                                    Species=uncultured bacterium
    """
    parts = [x.strip() for x in str(lineage).split(";")]

    # Pad to 8 ranks
    while len(parts) < 8:
        parts.append("")

    parsed = {}
    for rank, value in zip(RANKS, parts[:8]):
        parsed[rank] = normalise_taxon(value, rank)

    return parsed


def parse_blast_lineage(line: str) -> dict | None:
    """
    Parse lineage from an annotated BLAST TSV line.

    The annotated file (output of filter_blast_all_conditions.py, which
    passes through the lineage column from annotate_blast_lineage.py) has
    the lineage as its LAST column, semicolon-delimited:
      Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species

    This matches the Kraken lineage format so both use the same parser.
    """
    fields = line.rstrip("\n").split("\t")

    if len(fields) < 2:
        return None

    lineage = fields[-1].strip()

    # Must look like a lineage (has semicolons or is a single taxon)
    # Reject lines where the last field is a numeric taxid (not yet annotated)
    if lineage.isdigit():
        return None

    return parse_semicolon_lineage(lineage)


def parse_kraken_lineage(line: str) -> dict | None:
    """
    Parse lineage from an annotated Kraken TSV line.

    The annotated file has the lineage as its LAST column, semicolon-delimited:
      Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species
    """
    fields = line.rstrip("\n").split("\t")

    if not fields:
        return None

    lineage = fields[-1].strip()

    # Sanity check: should contain semicolons
    if not lineage:
        return None

    return parse_semicolon_lineage(lineage)


def extract_blast_read_id(line: str) -> str | None:
    if not line.strip():
        return None

    if "\t" in line:
        rid = line.split("\t", 1)[0]
    else:
        rid = line.split(maxsplit=1)[0]

    return clean_read_id(rid)


def extract_kraken_read_id(line: str) -> str | None:
    if not line.strip():
        return None

    fields = line.rstrip("\n").split("\t")

    if len(fields) >= 2:
        return clean_read_id(fields[1])

    return None


def process_blast_file(path: Path | None, wanted_ids: set[str], weights: dict[str, int] | None):
    counts = init_counts()
    matched_ids = set()

    stats = {
        "lines_seen": 0,
        "wanted_hits": 0,
        "unique_ids_assigned": 0,
        "weighted_hits": 0,
        "missing_lineage": 0,
    }

    if path is None or not path.exists():
        return counts, matched_ids, stats

    with open_maybe_gzip(path) as f:
        for line in f:
            if not line.strip():
                continue

            stats["lines_seen"] += 1

            rid = extract_blast_read_id(line)
            if not rid or rid not in wanted_ids:
                continue

            stats["wanted_hits"] += 1
            matched_ids.add(rid)

            parsed = parse_blast_lineage(line)
            if parsed is None:
                stats["missing_lineage"] += 1
                parsed = {rank: "unassigned" for rank in RANKS}

            weight = get_weight(rid, weights)
            stats["weighted_hits"] += weight

            for rank in RANKS:
                counts[rank][parsed[rank]] += weight

    stats["unique_ids_assigned"] = len(matched_ids)

    return counts, matched_ids, stats


def process_kraken_file(
    path: Path | None,
    wanted_1: set[str],
    wanted_2: set[str],
    weights_1: dict[str, int] | None,
    weights_2: dict[str, int] | None,
):
    counts = init_counts()
    matched_1 = set()
    matched_2 = set()

    stats = {
        "lines_seen": 0,
        "wanted_hits": 0,
        "unique_ids_assigned": 0,
        "weighted_hits": 0,
        "missing_lineage": 0,
        "ambiguous_mate_ids": 0,
    }

    if path is None or not path.exists():
        return counts, matched_1, matched_2, stats

    wanted_all = wanted_1 | wanted_2

    with open_maybe_gzip(path) as f:
        for line in f:
            if not line.strip():
                continue

            stats["lines_seen"] += 1

            rid = extract_kraken_read_id(line)
            if not rid or rid not in wanted_all:
                continue

            stats["wanted_hits"] += 1

            parsed = parse_kraken_lineage(line)
            if parsed is None:
                stats["missing_lineage"] += 1
                parsed = {rank: "unassigned" for rank in RANKS}

            in_1 = rid in wanted_1
            in_2 = rid in wanted_2

            if in_1 and in_2:
                stats["ambiguous_mate_ids"] += 1
                matched_1.add(rid)
                matched_2.add(rid)
                weight = get_weight(rid, weights_1) + get_weight(rid, weights_2)
            elif in_1:
                matched_1.add(rid)
                weight = get_weight(rid, weights_1)
            else:
                matched_2.add(rid)
                weight = get_weight(rid, weights_2)

            stats["weighted_hits"] += weight

            for rank in RANKS:
                counts[rank][parsed[rank]] += weight

    stats["unique_ids_assigned"] = len(matched_1 | matched_2)

    return counts, matched_1, matched_2, stats


def counts_to_df(counts_for_rank: Counter, total_reads: int | None = None, already_rpm: bool = False):
    rows = []

    for taxon, value in counts_for_rank.items():
        value = float(value)

        if already_rpm:
            raw_count = float("nan")
            rpm = value
        else:
            raw_count = value
            if total_reads and total_reads > 0:
                rpm = raw_count / total_reads * 1_000_000.0
            else:
                rpm = raw_count

        rows.append({
            "Taxonomy": taxon,
            "RawCount": raw_count,
            "RPM": rpm,
            "Count": rpm,
        })

    df = pd.DataFrame(rows)

    if df.empty:
        return pd.DataFrame(columns=["Taxonomy", "RawCount", "RPM", "Count", "Percentage"])

    total_rpm = df["RPM"].sum()
    df["Percentage"] = (df["RPM"] / total_rpm * 100.0) if total_rpm > 0 else 0.0

    return df.sort_values("RPM", ascending=False).reset_index(drop=True)


def save_summary(counts, outprefix: Path, total_reads: int | None = None, already_rpm: bool = False):
    for rank in RANKS:
        df = counts_to_df(counts[rank], total_reads=total_reads, already_rpm=already_rpm)
        df.to_csv(f"{outprefix}_{rank}_summary.tsv", sep="\t", index=False)


def make_top10(counts, outprefix: Path, total_reads: int | None = None, already_rpm: bool = False):
    for rank in RANKS:
        df = counts_to_df(counts[rank], total_reads=total_reads, already_rpm=already_rpm)

        if df.empty:
            out = pd.DataFrame(columns=["Taxonomy", "RawCount", "RPM", "Count", "Percentage"])
            out.to_csv(f"{outprefix}_{rank}_top10.tsv", sep="\t", index=False)
            continue

        unassigned_mask = df["Taxonomy"].astype(str).str.lower().isin(["unassigned", "unclassified"])

        unassigned_raw = df.loc[unassigned_mask, "RawCount"].sum(skipna=True)
        unassigned_rpm = df.loc[unassigned_mask, "RPM"].sum()

        df_nonun = df.loc[~unassigned_mask].copy()
        df_top = df_nonun.head(10).copy()

        other_raw = df_nonun.iloc[10:]["RawCount"].sum(skipna=True)
        other_rpm = df_nonun.iloc[10:]["RPM"].sum()

        rows = df_top[["Taxonomy", "RawCount", "RPM", "Count"]].to_dict("records")

        rows.append({
            "Taxonomy": "unassigned",
            "RawCount": unassigned_raw,
            "RPM": unassigned_rpm,
            "Count": unassigned_rpm,
        })

        rows.append({
            "Taxonomy": "other",
            "RawCount": other_raw,
            "RPM": other_rpm,
            "Count": other_rpm,
        })

        top_df = pd.DataFrame(rows)

        total_rpm = top_df["RPM"].sum()
        top_df["Percentage"] = (top_df["RPM"] / total_rpm * 100.0) if total_rpm > 0 else 0.0

        top_df.to_csv(f"{outprefix}_{rank}_top10.tsv", sep="\t", index=False)


def main(args):
    ids_dir = Path(args.ids_dir)
    blast_dir = Path(args.blast_dir) if args.blast_dir else None
    kraken_dir = Path(args.kraken_dir) if args.kraken_dir else None
    dup_dir = Path(args.dup_dir)
    outdir = Path(args.outdir)
    reads_table = Path(args.reads_table)

    outdir.mkdir(parents=True, exist_ok=True)

    total_reads_by_sample = read_total_reads(reads_table)
    samples = discover_samples(ids_dir)

    if not samples:
        raise ValueError(f"No *_ge5_detected_ids.txt files found in {ids_dir}")

    sample_to_group = load_samplesheet(args.samplesheet)
    groups = sorted(set(sample_to_group.values()))

    # A group with only one sample can share its label with that sample
    # (group defaults to the sample's own ID when unset) — its combined
    # output would otherwise be written to the exact same path as the
    # per-sample output, silently overwriting it. Count actual group
    # membership among samples being processed so those can be skipped.
    group_sizes = Counter(
        sample_group(s, sample_to_group) for s in samples
        if sample_group(s, sample_to_group) is not None
    )

    group_blast_counts    = {g: init_counts() for g in groups}
    group_kraken_counts   = {g: init_counts() for g in groups}
    group_combined_counts = {g: init_counts() for g in groups}

    run_report_rows = []

    for sample in sorted(samples):
        group = sample_group(sample, sample_to_group)

        if group is None:
            print(f"Skipping {sample}: sample not found in {args.samplesheet}")
            continue

        print(f"\nProcessing sample: {sample} -> group {group}")

        ids_file = samples[sample]
        wanted_ids_raw = load_id_file(ids_file)

        # COLLAPSE_READS_{PRE,POST}_TAXONOMY append /1 and /2 to every read ID
        # so mates stay distinguishable after merging (fixed 2026-09 -- prior
        # to that, mate 1 and mate 2 of the same fragment shared one bare ID
        # and could collide once collapsed). Split them here so each mate's
        # BLAST/Kraken lookup uses only its own pool. BLAST/Kraken TSVs still
        # carry the original unsuffixed read IDs, so strip the suffix before
        # matching.
        wanted_1_original = {rid[:-2] for rid in wanted_ids_raw if rid.endswith('/1')}
        wanted_2_original = {rid[:-2] for rid in wanted_ids_raw if rid.endswith('/2')}

        print(f"  wanted IDs: {len(wanted_ids_raw)} total "
              f"({len(wanted_1_original)} mate-1, {len(wanted_2_original)} mate-2)")

        total_reads = total_reads_by_sample.get(sample)

        if total_reads is None:
            print(f"  WARNING: {sample} not found in {reads_table}; RPM will use raw counts for this sample")

        # COLLAPSE_READS_{PRE,POST}_TAXONOMY produce a single merged detail
        # file, not separate per-mate files. The /1 and /2 suffixes on every
        # rep_id let us reconstruct per-mate weight dicts with unsuffixed
        # keys, matching the original IDs in the BLAST/Kraken TSVs.
        dup_file = resolve_file(dup_dir / f"{sample}_merged_duplicated.detail.txt")
        merged_weights = parse_dup_detail_file(dup_file)

        weights_1 = {k[:-2]: v for k, v in merged_weights.items() if k.endswith('/1')}
        weights_2 = {k[:-2]: v for k, v in merged_weights.items() if k.endswith('/2')}

        if dup_file:
            print(f"  duplicate detail (merged): {dup_file} "
                  f"({len(weights_1)} mate-1 clusters, {len(weights_2)} mate-2 clusters)")
        else:
            print("  duplicate detail: not found, using weight 1 for all reads")

        blast1_file = resolve_file(blast_dir / f"{sample}_1_blast_all_lengths_filtered.tsv") if blast_dir else None
        blast2_file = resolve_file(blast_dir / f"{sample}_2_blast_all_lengths_filtered.tsv") if blast_dir else None
        kraken_file = resolve_file(kraken_dir / f"{sample}_kraken_invertebrates_filtered.tsv") if kraken_dir else None

        sample_blast_counts = init_counts()
        sample_kraken_counts = init_counts()
        sample_combined_counts = init_counts()

        assigned_1 = set()
        assigned_2 = set()

        if args.priority == "blast":
            # ---------------- BLAST first ----------------
            if blast_dir:
                counts, matched_1, stats = process_blast_file(blast1_file, wanted_1_original, weights_1)
                add_counts(sample_blast_counts, counts)
                add_counts(sample_combined_counts, counts)
                assigned_1.update(matched_1)

                run_report_rows.append({
                    "Sample": sample,
                    "Group": group,
                    "Tool": "blast",
                    "Mate": "1",
                    "Priority": args.priority,
                    "File": str(blast1_file) if blast1_file else "missing",
                    "WantedIDs": len(wanted_1_original),
                    "IDsAlreadyAssignedBeforeThisTool": 0,
                    "IDsAvailableToThisTool": len(wanted_1_original),
                    **stats,
                })

                counts, matched_2, stats = process_blast_file(blast2_file, wanted_2_original, weights_2)
                add_counts(sample_blast_counts, counts)
                add_counts(sample_combined_counts, counts)
                assigned_2.update(matched_2)

                run_report_rows.append({
                    "Sample": sample,
                    "Group": group,
                    "Tool": "blast",
                    "Mate": "2",
                    "Priority": args.priority,
                    "File": str(blast2_file) if blast2_file else "missing",
                    "WantedIDs": len(wanted_2_original),
                    "IDsAlreadyAssignedBeforeThisTool": 0,
                    "IDsAvailableToThisTool": len(wanted_2_original),
                    **stats,
                })

            # ---------------- Kraken only remaining IDs ----------------
            if kraken_dir:
                remaining_1 = wanted_1_original - assigned_1
                remaining_2 = wanted_2_original - assigned_2

                counts, matched_k1, matched_k2, stats = process_kraken_file(
                    kraken_file,
                    remaining_1,
                    remaining_2,
                    weights_1,
                    weights_2,
                )

                add_counts(sample_kraken_counts, counts)
                add_counts(sample_combined_counts, counts)
                assigned_1.update(matched_k1)
                assigned_2.update(matched_k2)

                run_report_rows.append({
                    "Sample": sample,
                    "Group": group,
                    "Tool": "kraken",
                    "Mate": "both",
                    "Priority": args.priority,
                    "File": str(kraken_file) if kraken_file else "missing",
                    "WantedIDs": len(wanted_1_original | wanted_2_original),
                    "IDsAlreadyAssignedBeforeThisTool": len((wanted_1_original & assigned_1) | (wanted_2_original & assigned_2)),
                    "IDsAvailableToThisTool": len(remaining_1 | remaining_2),
                    **stats,
                })

        else:
            # ---------------- Kraken first ----------------
            if kraken_dir:
                counts, matched_k1, matched_k2, stats = process_kraken_file(
                    kraken_file,
                    wanted_1_original,
                    wanted_2_original,
                    weights_1,
                    weights_2,
                )

                add_counts(sample_kraken_counts, counts)
                add_counts(sample_combined_counts, counts)
                assigned_1.update(matched_k1)
                assigned_2.update(matched_k2)

                run_report_rows.append({
                    "Sample": sample,
                    "Group": group,
                    "Tool": "kraken",
                    "Mate": "both",
                    "Priority": args.priority,
                    "File": str(kraken_file) if kraken_file else "missing",
                    "WantedIDs": len(wanted_1_original | wanted_2_original),
                    "IDsAlreadyAssignedBeforeThisTool": 0,
                    "IDsAvailableToThisTool": len(wanted_1_original | wanted_2_original),
                    **stats,
                })

            # ---------------- BLAST only remaining IDs ----------------
            if blast_dir:
                remaining_1 = wanted_1_original - assigned_1
                remaining_2 = wanted_2_original - assigned_2

                counts, matched_1, stats = process_blast_file(blast1_file, remaining_1, weights_1)
                add_counts(sample_blast_counts, counts)
                add_counts(sample_combined_counts, counts)
                assigned_1.update(matched_1)

                run_report_rows.append({
                    "Sample": sample,
                    "Group": group,
                    "Tool": "blast",
                    "Mate": "1",
                    "Priority": args.priority,
                    "File": str(blast1_file) if blast1_file else "missing",
                    "WantedIDs": len(wanted_1_original),
                    "IDsAlreadyAssignedBeforeThisTool": len(wanted_1_original - remaining_1),
                    "IDsAvailableToThisTool": len(remaining_1),
                    **stats,
                })

                counts, matched_2, stats = process_blast_file(blast2_file, remaining_2, weights_2)
                add_counts(sample_blast_counts, counts)
                add_counts(sample_combined_counts, counts)
                assigned_2.update(matched_2)

                run_report_rows.append({
                    "Sample": sample,
                    "Group": group,
                    "Tool": "blast",
                    "Mate": "2",
                    "Priority": args.priority,
                    "File": str(blast2_file) if blast2_file else "missing",
                    "WantedIDs": len(wanted_2_original),
                    "IDsAlreadyAssignedBeforeThisTool": len(wanted_2_original - remaining_2),
                    "IDsAvailableToThisTool": len(remaining_2),
                    **stats,
                })

        # ---------------- Per-sample outputs ----------------
        if blast_dir:
            save_summary(sample_blast_counts, outdir / f"{sample}_blast", total_reads=total_reads, already_rpm=False)
            make_top10(sample_blast_counts, outdir / f"{sample}_blast", total_reads=total_reads, already_rpm=False)

        if kraken_dir:
            save_summary(sample_kraken_counts, outdir / f"{sample}_kraken", total_reads=total_reads, already_rpm=False)
            make_top10(sample_kraken_counts, outdir / f"{sample}_kraken", total_reads=total_reads, already_rpm=False)

        if blast_dir and kraken_dir:
            save_summary(sample_combined_counts, outdir / f"{sample}_combined", total_reads=total_reads, already_rpm=False)
            make_top10(sample_combined_counts, outdir / f"{sample}_combined", total_reads=total_reads, already_rpm=False)

        # ---------------- Add sample-normalised RPM to its group ----------------
        sample_blast_rpm = rpm_normalise_counts(sample_blast_counts, total_reads)
        sample_kraken_rpm = rpm_normalise_counts(sample_kraken_counts, total_reads)
        sample_combined_rpm = rpm_normalise_counts(sample_combined_counts, total_reads)

        add_counts(group_blast_counts[group], sample_blast_rpm)
        add_counts(group_kraken_counts[group], sample_kraken_rpm)
        add_counts(group_combined_counts[group], sample_combined_rpm)

        total_requested_unique = len(wanted_1_original | wanted_2_original)
        total_assigned_unique = len(assigned_1 | assigned_2)

        run_report_rows.append({
            "Sample": sample,
            "Group": group,
            "Tool": "assignment_summary",
            "Mate": "both",
            "Priority": args.priority,
            "File": "NA",
            "WantedIDs": total_requested_unique,
            "IDsAlreadyAssignedBeforeThisTool": "NA",
            "IDsAvailableToThisTool": "NA",
            "lines_seen": "NA",
            "wanted_hits": "NA",
            "unique_ids_assigned": total_assigned_unique,
            "weighted_hits": "NA",
            "missing_lineage": "NA",
            "ambiguous_mate_ids": "NA",
            "unassigned_selected_ids": total_requested_unique - total_assigned_unique,
        })

    # ---------------- Per-group combined outputs ----------------
    for group in groups:
        if group_sizes.get(group, 0) < 2:
            # Solo group: its label may equal its one sample's ID (default
            # group), which would make the combined-output path below
            # collide with, and silently overwrite, that sample's own
            # per-sample output. The per-sample output already covers this
            # group fully, so skip writing a redundant/colliding one.
            print(f"Skipping combined output for group '{group}': only one sample "
                  f"in this group; its per-sample output already covers it")
            continue

        if blast_dir:
            save_summary(group_blast_counts[group], outdir / f"{group}_blast", already_rpm=True)
            make_top10(group_blast_counts[group], outdir / f"{group}_blast", already_rpm=True)

        if kraken_dir:
            save_summary(group_kraken_counts[group], outdir / f"{group}_kraken", already_rpm=True)
            make_top10(group_kraken_counts[group], outdir / f"{group}_kraken", already_rpm=True)

        if blast_dir and kraken_dir:
            save_summary(group_combined_counts[group], outdir / f"{group}_combined", already_rpm=True)
            make_top10(group_combined_counts[group], outdir / f"{group}_combined", already_rpm=True)

    report_df = pd.DataFrame(run_report_rows)
    report_df.to_csv(outdir / "fetch_reinflate_report.tsv", sep="\t", index=False)

    print("\nDone.")
    print(f"Output directory: {outdir}")
    print(f"Report: {outdir / 'fetch_reinflate_report.tsv'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Fetch selected BLAST/Kraken taxonomy results from sample_ge5_detected_ids.txt (mate-suffixed "
            "/1 /2 IDs from a single merged candidate pool), assign each read ID to only one tool, "
            "reinflate using the merged duplicate detail file, normalise by reads-per-million, and write "
            "per-sample plus per-group combined summaries."
        )
    )

    parser.add_argument("--ids-dir", required=True, help="Directory with sample_ge5_detected_ids.txt files")
    parser.add_argument("--blast-dir", required=False, help="Directory with corrected BLAST filtered TSV files")
    parser.add_argument("--kraken-dir", required=False, help="Directory with corrected Kraken filtered TSV files")
    parser.add_argument("--dup-dir", required=True, help="Directory with SAMPLE_merged_duplicated.detail.txt files (one per sample, mate-suffixed IDs)")
    parser.add_argument("--reads-table", required=True, help="reads_posttrim_tab.tsv with Sample and TotalReads columns")
    parser.add_argument("--samplesheet", required=True, help="sample,group,fastq_1,fastq_2 sample sheet")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--priority",
        choices=["blast", "kraken"],
        default="blast",
        help="If a selected read ID appears in both tools, which tool gets it. Default: blast",
    )

    args = parser.parse_args()

    if not args.blast_dir and not args.kraken_dir:
        raise ValueError("Provide at least one of --blast-dir or --kraken-dir")

    main(args)