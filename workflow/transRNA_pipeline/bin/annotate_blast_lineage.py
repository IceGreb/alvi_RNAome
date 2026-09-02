#!/usr/bin/env python3
"""
annotate_blast_lineage.py
─────────────────────────
Annotates BLAST tabular output with an 8-rank lineage column.

--input accepts a glob pattern (a single real file path also works, since
glob() matches it trivially) — e.g. "RJ1_1_*_top_hits.tsv" for a sample-mate
whose BLAST run was split into many chunk files. Order the chunks are read
in doesn't matter: each is a self-contained, non-overlapping slice of the
original run (no query's rows are split across two chunks), and this script
only aggregates a taxid vocabulary and re-annotates rows — it doesn't need
same-query rows to be contiguous the way the candidate-level filter does.

Two passes over the matched files, deliberately NOT buffering line content
between them (only Pass 1's small taxid set is kept in memory):
  Pass 1: read every matched file, collect the set of unique non-zero
          taxids only (not the lines themselves).
  Pass 2: re-read the same files a second time, resolving each row's
          lineage from the now-built taxid->lineage map and writing
          directly to output.
This trades one extra sequential read pass for a memory footprint bounded
by the taxid vocabulary size, not by row count -- for a 60M-row sample-mate
that's the difference between tens of GB and a rounding error.

Writes the temp taxid file to the current working directory (the Nextflow
work dir on RDS) instead of /tmp, to avoid filling the login node's /tmp.
"""

import argparse
import glob
import re
import subprocess
import sys
from pathlib import Path

TAXID_SEP = re.compile(r"[;,]")

LINEAGE_FMT = (
    "{domain|acellular root|superkingdom}"
    ";{kingdom};{phylum};{class};{order};{family};{genus};{species}"
)

FALLBACK = ";;;;;;;"   # 7 semicolons = 8 empty fields


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def first_taxid(field: str) -> int:
    for part in TAXID_SEP.split(field.strip()):
        part = part.strip()
        if part.isdigit():
            v = int(part)
            if v > 0:
                return v
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",     required=True,
                     help="A file path, or a glob pattern matching multiple "
                          "chunk files (e.g. 'RJ1_1_*_top_hits.tsv')")
    ap.add_argument("--output",    required=True)
    ap.add_argument("--taxid-col", type=int, default=13)
    args = ap.parse_args()

    taxid_idx = args.taxid_col - 1

    input_files = sorted(glob.glob(args.input))   # sorted: reproducible logs
                                                    # only -- order doesn't
                                                    # affect correctness here
    if not input_files:
        sys.exit(f"ERROR: no files matched {args.input!r}")
    log(f"Pass 1: scanning {len(input_files)} file(s) matching {args.input!r} "
        f"(taxid col={args.taxid_col})")

    unique_taxids = set()
    total_lines = 0
    for fpath in input_files:
        with open(fpath) as fh:
            for line in fh:
                total_lines += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) > taxid_idx:
                    tid = first_taxid(fields[taxid_idx])
                    if tid:
                        unique_taxids.add(tid)
    unique_taxids = sorted(unique_taxids)
    log(f"  {total_lines} lines across {len(input_files)} file(s), "
        f"{len(unique_taxids)} unique non-zero taxids")

    # ── Write temp taxid file to current dir (Nextflow work dir on RDS) ───────
    # Avoids filling /tmp on the login/compute node
    tmp_path = Path(".") / f"_taxids_tmp_{Path(args.output).stem}.txt"
    log(f"Running taxonkit reformat2 (temp file: {tmp_path})...")

    try:
        with tmp_path.open("w") as tmp:
            for tid in unique_taxids:
                tmp.write(f"{tid}\n")

        proc = subprocess.run(
            ["taxonkit", "reformat2", "-I", "1", "-f", LINEAGE_FMT, str(tmp_path)],
            capture_output=True,
            text=True,
        )
    finally:
        tmp_path.unlink(missing_ok=True)

    if proc.returncode != 0:
        log(f"[ERROR] taxonkit reformat2 failed (exit {proc.returncode}):")
        log(proc.stderr[:2000])
        sys.exit(1)

    if proc.stderr.strip():
        log(proc.stderr.strip())

    # ── Parse output: {taxid: lineage_string} ────────────────────────────────
    lineage_map = {}
    for out_line in proc.stdout.splitlines():
        parts = out_line.strip().split("\t")
        if len(parts) < 2:
            continue
        try:
            tid = int(parts[0])
        except ValueError:
            continue
        lineage_map[tid] = parts[1].strip()

    log(f"  {len(lineage_map)} lineages resolved")

    # ── Pass 2: re-read the same files, write annotated output ───────────────
    # Deliberately re-reads from disk rather than replaying a buffered line
    # list from Pass 1 -- see module docstring.
    log(f"Pass 2: writing {args.output}")
    no_lineage = 0

    with open(args.output, "w") as out:
        for fpath in input_files:
            with open(fpath) as fh:
                for line in fh:
                    raw = line.rstrip("\n")
                    fields = raw.split("\t")
                    tid = first_taxid(fields[taxid_idx]) if len(fields) > taxid_idx else 0
                    lineage = lineage_map.get(tid, FALLBACK)
                    if not lineage:
                        lineage = FALLBACK
                        no_lineage += 1
                    out.write(raw + "\t" + lineage + "\n")

    if no_lineage:
        log(f"  {no_lineage} lines used fallback (taxid 0 or not resolved)")
    log("Done.")


if __name__ == "__main__":
    main()
