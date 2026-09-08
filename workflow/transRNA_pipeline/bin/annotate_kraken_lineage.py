#!/usr/bin/env python3
"""
annotate_kraken_lineage.py
──────────────────────────
Annotates a Kraken2 output file with an 8-rank lineage column.

Kraken field 3 looks like: "Eukaryota (taxid 2759)"
taxonkit needs a plain integer taxid, so this script:

  1. Reads the file once, extracts taxid from the "(taxid N)" pattern in
     field 3, and keeps only the set of unique non-zero taxids -- not the
     lines themselves (a multi-GB Kraken output shouldn't need a matching
     multi-GB buffer just to re-emit the same lines with one extra column).
  2. Writes a temp file of unique taxids (one per line, plain integers)
  3. Runs: taxonkit reformat2 -I 1 -f '{...}' <tmpfile>
  4. Re-reads the file a second time, mapping each row's taxid to its
     resolved lineage and writing directly to output.

Output lineage column (semicolon-delimited):
  Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species

Empty ranks from taxonkit are left as empty strings — downstream normalisation
in the taxonomy parser handles "unclassified" labelling so each read is
classified to the highest available rank.

Usage:
  python3 annotate_kraken_lineage.py --input kraken.txt --output annotated.tsv
"""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

TAXID_RE    = re.compile(r"\(taxid\s+(\d+)\)")
LINEAGE_FMT = (
    "{domain|acellular root|superkingdom};{kingdom};{phylum};{class};{order};{family};{genus};{species}"
)
FALLBACK = ";;;;;;;"   # 8 empty semicolon-separated fields (7 separators)


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def extract_taxid(field3: str) -> int:
    m = TAXID_RE.search(field3)
    return int(m.group(1)) if m else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  required=True, help="Kraken2 output file")
    ap.add_argument("--output", required=True, help="Annotated output TSV")
    args = ap.parse_args()

    # ── Pass 1: collect only the unique taxid vocabulary ──────────────────────
    log(f"Pass 1: reading {args.input}")
    unique_taxids = set()
    total_lines = 0

    with open(args.input) as fh:
        for line in fh:
            total_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                tid = extract_taxid(fields[2])
                if tid:
                    unique_taxids.add(tid)

    unique_taxids = sorted(unique_taxids)
    log(f"  {total_lines} lines, {len(unique_taxids)} unique non-zero taxids")

    # ── Run taxonkit reformat2 on unique taxids ───────────────────────────────
    log("Running taxonkit reformat2...")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".taxids.txt",
                                     delete=False) as tmp:
        tmp_path = tmp.name
        for tid in unique_taxids:
            tmp.write(f"{tid}\n")

    proc = subprocess.run(
        [
            "taxonkit", "reformat2",
            "-I", "1",
            "-f", LINEAGE_FMT,
            tmp_path,
        ],
        capture_output=True,
        text=True,
    )
    Path(tmp_path).unlink(missing_ok=True)

    if proc.returncode != 0:
        log(f"[ERROR] taxonkit reformat2 failed (exit {proc.returncode}):")
        log(proc.stderr[:2000])
        sys.exit(1)

    if proc.stderr.strip():
        log(proc.stderr.strip())

    # ── Parse taxonkit output: {taxid: lineage_string} ───────────────────────
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

    # ── Pass 2: re-read the file, write annotated output ──────────────────────
    # Re-reads from disk rather than replaying a buffered line list from
    # Pass 1 -- see module docstring.
    log(f"Pass 2: writing {args.output}")
    no_lineage = 0

    with open(args.input) as fh, open(args.output, "w") as out:
        for line in fh:
            raw = line.rstrip("\n")
            fields = raw.split("\t")
            tid = extract_taxid(fields[2]) if len(fields) >= 3 else 0
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