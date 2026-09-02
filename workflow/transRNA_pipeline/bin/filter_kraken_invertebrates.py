#!/usr/bin/env python3
"""
filter_kraken_invertebrates.py
Only CLASSIFIED (C) rows are kept. Unclassified (U) rows are silently
skipped — they are NOT lost, they are handled by the BLAST path.

Among classified rows, removes:
  • Viruses/acellular organisms (Domain rank)
  • All invertebrate phyla (Phylum rank, from invertebrate_phyla.txt)
  • Host genera: Homo, Mus, Canis, Felis, Apis (Genus rank)
  • Reads where either mate length < min_len nt

Now also writes a filter_stats TSV with exclusion counts per category,
used by aggregate_report.py for the virus exclusion report.

Lineage col 6: Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species
"""
from pathlib import Path
import argparse, re, sys

COMMON_NAMES = {
    "sponges","coral","jellyfish","anemones","comb jellies","flatworms",
    "jaw worms","mesozoa","proboscis worms","gastrotrichs","rotifers",
    "roundworms","horsehair worms","mud dragons","spiny-crown worms",
    "acanthocephalans","spiny-headed worms","brush heads","pandora",
    "cycliophorans","goblet worms","marine mats","moss animals","bryozoans",
    "horseshoe worms","brachipods","lampshells","molluscs","slugs","snails",
    "squid","peanut worms","segmented worms","earthworms","ragworms",
    "spoon worms","beard worms","water bears","velvet worms","insects",
    "spiders","crabs","etc","starfish","urchins","arrow worms",
    "acorn worms","vertebrates","invertebrates",
}
HOST_GENERA   = {"homo","mus","canis","felis"}  # apis dropped: always caught
                                                  # by the invertebrate check
                                                  # first (see main(), below)
VIRUS_DOMAINS = {"viruses","virus","acellular root","acellular organisms"}

def log(msg): print(msg, file=sys.stderr, flush=True)

def clean_name(s):
    s = s.strip().strip(" ,;:.()[]{}")
    return re.sub(r"\s+", " ", s)

def extract_names(line):
    parts = re.split(r"[()]", line.strip())
    names = set()
    for part in parts:
        part = part.replace(",", " or ")
        for token in re.split(r"\s+or\s+", part):
            name = clean_name(token)
            if name and name.lower() not in COMMON_NAMES:
                names.add(name)
    return names

def load_excluded_phyla(path):
    names = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line: names.update(extract_names(line))
    return names

def build_pattern(names):
    escaped = sorted((re.escape(x) for x in names), key=len, reverse=True)
    return re.compile(r"(?:^|;)(?:" + "|".join(escaped) + r")(?:;|$)")

def both_mates_ok(length_field, min_len):
    parts = length_field.split("|")
    try: return all(int(p) >= min_len for p in parts if p)
    except ValueError: return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exclude_file")
    ap.add_argument("input_tsv")
    ap.add_argument("output_tsv")
    ap.add_argument("--min-len", type=int, default=18)
    ap.add_argument("--stats-out", default=None,
                    help="Optional path to write per-category exclusion counts TSV")
    ap.add_argument("--unclassified-out", default=None,
                    help="Optional path to write U-row read IDs, one per line — "
                         "the real BLAST query list under collapse_mode=pre_taxonomy "
                         "(vs. the pre-computed params.blast_dir fallback otherwise). "
                         "Previously U rows were just counted and dropped; this makes "
                         "the 'U rows -> BLAST path' handoff a real, wired output "
                         "instead of an unenforced assumption about how an external "
                         "BLAST run was set up.")
    args = ap.parse_args()

    excluded = load_excluded_phyla(args.exclude_file)
    if not excluded: sys.exit("ERROR: No phyla names parsed")
    pattern = build_pattern(excluded)

    total = classified = kept = 0
    rm_u = rm_short = rm_virus = rm_host = rm_inv = bad = 0

    funclass = open(args.unclassified_out, "w") if args.unclassified_out else None
    try:
        with open(args.input_tsv) as fin, open(args.output_tsv, "w") as fout:
            for line in fin:
                total += 1
                raw = line.rstrip("\n")
                fields = raw.split("\t")
                if len(fields) < 6: bad += 1; continue

                status       = fields[0].strip()
                read_id      = fields[1].strip()
                length_field = fields[3].strip()
                lineage_col  = fields[5].strip()

                # U rows → BLAST path, not lost
                if status != "C":
                    rm_u += 1
                    if funclass: funclass.write(read_id + "\n")
                    continue
                classified += 1

                if not both_mates_ok(length_field, args.min_len): rm_short += 1; continue

                parts = [x.strip() for x in lineage_col.split(";")]
                while len(parts) < 8: parts.append("unassigned")
                domain = parts[0].lower()
                genus  = parts[6].lower()

                # Invertebrate check before host-genus: Apis is itself an
                # invertebrate (Arthropoda), so it's always caught here first —
                # HOST_GENERA below no longer needs "apis" listed at all, and
                # true off-target host contamination (Homo/Mus/Canis/Felis)
                # still gets its own accurate rm_host count.
                if domain in VIRUS_DOMAINS:     rm_virus += 1; continue
                if pattern.search(lineage_col): rm_inv   += 1; continue
                if genus in HOST_GENERA:        rm_host  += 1; continue

                fout.write(raw + "\n")
                kept += 1
    finally:
        if funclass: funclass.close()

    summary = (
        f"{Path(args.input_tsv).name}: total={total} "
        f"classified={classified} kept={kept} | "
        f"skipped_unclassified(→BLAST)={rm_u} rm_short={rm_short} "
        f"rm_virus={rm_virus} rm_host={rm_host} rm_invertebrate={rm_inv} "
        f"bad_format={bad} | excluded_phyla={len(excluded)}"
    )
    log(summary)

    # Write stats file for the aggregate report
    stats_path = args.stats_out
    if stats_path is None:
        # Default: same stem as output but _filter_stats.tsv
        stats_path = str(Path(args.output_tsv).with_suffix("")) + "_filter_stats.tsv"
    with open(stats_path, "w") as sf:
        sf.write("sample\ttool\ttotal\tclassified\tkept\t"
                 "rm_unclassified\trm_short\trm_virus\trm_host\trm_invertebrate\tbad\n")
        sample = Path(args.input_tsv).name.replace("_kraken_annotated.tsv", "")
        sf.write(f"{sample}\tkraken\t{total}\t{classified}\t{kept}\t"
                 f"{rm_u}\t{rm_short}\t{rm_virus}\t{rm_host}\t{rm_inv}\t{bad}\n")

if __name__ == "__main__": main()