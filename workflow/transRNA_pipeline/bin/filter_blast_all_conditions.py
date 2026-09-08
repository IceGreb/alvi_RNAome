#!/usr/bin/env python3
"""
filter_blast_all_conditions.py
───────────────────────────────
Filters an annotated BLAST TSV (output of annotate_blast_lineage.py).

TOP-HIT-ONLY filtering (changed 04_09_2026, was candidate-level/all-rows
before): a query is kept purely on the strength of its own top hit — the
other hit rows in its group are never even inspected. A tied or lower-
ranked bad hit no longer disqualifies a candidate whose best hit is clean.
(Old behavior: kept only if EVERY hit row passed. Real-data check on
02_09_2026's strict-filter run found 57.3% of the rows causing a
disqualification were merely TIED in bitscore with a passing top hit, not
an obviously worse alternative — motivating this change.)

BLAST's own outfmt-6 output is already sorted best-hit-first (bitscore
non-increasing) within a query's block — verified against the real data
(60.76M rows, zero violations) rather than assumed — so the first row of a
group IS its top hit, and every other row in the group can be skipped
without even parsing it.

Filtering criteria (top hit only):
  1.  Domain not in: Viruses, Virus, Acellular root, Acellular organisms
  2.  Phylum not in: excluded set from invertebrate_phyla.txt
  3.  Species contains no host fragment:
        homo sapiens, mus musculus, canis lupus, felis catus, apis mellifera
  4.  mismatch == 0
  5.  gapopen == 0
  6.  abs(qend - qstart) + 1 == alignment_length  (full-span)
  7.  alignment_length >= min_len  (default 18)
  8.  row must parse (>=9 tab-separated fields, numeric BLAST columns) —
      an unparseable top row can't be verified to pass, so it disqualifies
      its candidate (see row_outcome's "bad" case).

Relies on rows for the same query being contiguous in the input (verified:
zero interleaved query IDs across all 60.76M rows of RJ1 mate 1) — this is
what lets a single forward streaming pass do top-hit-only filtering
without loading the file or sorting it.

The lineage column is passed through to the output unchanged so the
taxonomy parser can re-use it directly.

Usage:
  python filter_blast_all_conditions.py [--min-len N] [--stats-out file]
      phyla.txt input_annotated.tsv output.tsv
"""

from pathlib import Path
import argparse
import re
import sys

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

VIRUS_DOMAINS = {"viruses","virus","acellular root","acellular organisms"}

HOST_FRAGS = [
    "homo sapiens","mus musculus","canis lupus","felis catus","apis mellifera"
]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


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
            if line:
                names.update(extract_names(line))
    return names


def is_host(species: str) -> bool:
    s = species.strip().lower()
    return any(f in s for f in HOST_FRAGS)


def parse_lineage(lineage_col: str) -> list:
    """
    Parse the semicolon-delimited lineage column (last column of annotated file).
    Returns list of 8 rank strings: [Domain, Kingdom, Phylum, ..., Species]
    """
    parts = [x.strip() for x in lineage_col.split(";")]
    while len(parts) < 8:
        parts.append("")
    return parts[:8]


def row_outcome(fields: list, excluded_phyla: set, min_len: int) -> str:
    """
    Classify a single row. Returns "kept" or the specific reason it fails
    (matching the stats file's rm_* column names, minus the "rm_" prefix):
    "bad", "virus", "phylum", "host", "mismatch", "gap", "not_fullspan", "short".

    Pulled out of the old main()'s inline loop unchanged in logic — same
    checks, same order, same thresholds — just callable per-row now instead
    of being embedded in a single streaming loop.
    """
    if len(fields) < 9:
        return "bad"
    try:
        aln_len  = int(fields[3])
        mismatch = int(fields[4])
        gapopen  = int(fields[5])
        qstart   = int(fields[6])
        qend     = int(fields[7])
    except ValueError:
        return "bad"

    qspan = abs(qend - qstart) + 1
    ranks   = parse_lineage(fields[-1])
    domain  = ranks[0].lower()
    phylum  = ranks[2]
    species = ranks[7]

    if domain in VIRUS_DOMAINS:      return "virus"
    if phylum in excluded_phyla:     return "phylum"
    if is_host(species):             return "host"
    if mismatch != 0:                return "mismatch"
    if gapopen  != 0:                return "gap"
    if qspan != aln_len:             return "not_fullspan"
    if aln_len < min_len:            return "short"
    return "kept"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exclude_file")
    ap.add_argument("input_tsv")
    ap.add_argument("output_tsv")
    ap.add_argument("--min-len",   type=int, default=18)
    ap.add_argument("--stats-out", default=None)
    args = ap.parse_args()

    excluded = load_excluded_phyla(args.exclude_file)
    if not excluded:
        sys.exit("ERROR: No phyla names parsed from exclude file")

    total = kept = unique_queries = 0
    rm_v = rm_p = rm_h = rm_mm = rm_gap = rm_span = rm_short = bad = 0

    def tally(reason):
        nonlocal rm_v, rm_p, rm_h, rm_mm, rm_gap, rm_span, rm_short, bad
        if reason == "virus":         rm_v    += 1
        elif reason == "phylum":      rm_p    += 1
        elif reason == "host":        rm_h    += 1
        elif reason == "mismatch":    rm_mm   += 1
        elif reason == "gap":         rm_gap  += 1
        elif reason == "not_fullspan":rm_span += 1
        elif reason == "short":       rm_short+= 1
        elif reason == "bad":         bad     += 1

    # ── State for the CURRENT group only. The decision is made the moment
    # a group's first (top) row is seen — every later row of the same
    # group is skipped without even being parsed, since only the top hit
    # matters now. "Flushing" (writing the kept top row, or tallying its
    # fail reason) happens whenever the query ID changes (the previous
    # group is complete) and once more after the loop ends (the last
    # group never triggers an ID-change).
    current_qid = None
    top_row     = None   # top row's raw text, to write out if kept
    decision    = None   # "kept" or a row_outcome() fail reason

    with open(args.input_tsv) as fin, open(args.output_tsv, "w") as fout:
        for line in fin:
            total += 1
            fields_probe = line.split("\t", 1)
            qid = fields_probe[0].strip() if fields_probe else ""

            if qid != current_qid:
                # ── close out the PREVIOUS group (if any) ──────────────────
                if top_row is not None:
                    if decision == "kept":
                        fout.write(top_row)
                        kept += 1
                    else:
                        tally(decision)
                # ── this row IS the new group's top hit — decide now ────────
                current_qid = qid
                raw = line.rstrip("\n")
                fields = raw.split("\t")
                top_row = raw + "\n"
                decision = row_outcome(fields, excluded, args.min_len)
                unique_queries += 1
            # else: not this group's top row — skip, it no longer matters

        # ── flush the LAST group in the file ──────────────────────────────
        if top_row is not None:
            if decision == "kept":
                fout.write(top_row)
                kept += 1
            else:
                tally(decision)

    summary = (
        f"{Path(args.input_tsv).name}: total_rows={total} "
        f"unique_queries={unique_queries} kept_candidates={kept} | "
        f"rm_virus={rm_v} rm_phylum={rm_p} rm_host={rm_h} "
        f"rm_mismatch={rm_mm} rm_gap={rm_gap} "
        f"rm_not_fullspan={rm_span} rm_short={rm_short} "
        f"bad_format={bad} | excluded_phyla={len(excluded)}"
    )
    log(summary)

    # Write stats file — same columns as before. "kept" now means kept
    # CANDIDATES (== rows written, since it's 1:1 now), not kept rows.
    stats_path = args.stats_out
    if stats_path is None:
        stats_path = str(Path(args.output_tsv).with_suffix("")) + "_filter_stats.tsv"
    stem = Path(args.input_tsv).name.replace("_blast_annotated.tsv", "")
    with open(stats_path, "w") as sf:
        sf.write("sample\ttool\ttotal\tunique_queries\tkept\t"
                 "rm_virus\trm_phylum\trm_host\trm_mismatch\t"
                 "rm_gap\trm_not_fullspan\trm_short\tbad\n")
        sf.write(
            f"{stem}\tblast\t{total}\t{unique_queries}\t{kept}\t"
            f"{rm_v}\t{rm_p}\t{rm_h}\t{rm_mm}\t"
            f"{rm_gap}\t{rm_span}\t{rm_short}\t{bad}\n"
        )


if __name__ == "__main__":
    main()
