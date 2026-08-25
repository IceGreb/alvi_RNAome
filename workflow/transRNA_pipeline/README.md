# transRNA Pipeline v1.2

Paired-end RNA-seq pipeline: raw reads → filtering → taxonomy → plots.  
Datasets: **RJ** (Royal Jelly) and **ST** (Systemic larval tissue), *Apis mellifera* by default or to be selected by user.

---

## How to run

```bash
# From the CSD3 login node — do NOT sbatch this:
bash run.sh           # fresh run from scratch
bash run.sh --resume  # resume after failure (Nextflow skips completed tasks)
```

To change partition or project, edit the top two lines of `run.sh`:
```bash
PARTITION="icelake"
PROJECT="<your-slurm-project>"
```

---

## Samples

Samples, groups, and raw read paths are all defined in `samples.csv`
(nf-core-style sample sheet). Any number of groups, with any labels, is supported,
and `group` is optional.

```csv
sample,group,fastq_1,fastq_2
RJ1,RJ,/path/to/raw_reads/RJ1_1.fq.gz,/path/to/raw_reads/RJ1_2.fq.gz
RJ2,RJ,/path/to/raw_reads/RJ2_1.fq.gz,/path/to/raw_reads/RJ2_2.fq.gz
T1GMN,ST,/path/to/raw_reads/T1GMN_1.fq.gz,/path/to/raw_reads/T1GMN_2.fq.gz
Extra1,,/path/to/raw_reads/Extra1_1.fq.gz,/path/to/raw_reads/Extra1_2.fq.gz
```

The example above (RJ = Royal Jelly, ST = Systemic larval tissue) reflects
the dataset used in the accompanying manuscript, but the pipeline itself
makes no assumption about group names or count. Edit `samples.csv` to
add, remove, or relabel samples/groups for a different dataset. Leaving
`group` blank (like `Extra1` above) is fine — a sample with no group is
treated as its own, single-sample group in every per-group summary/plot.

The sample sheet is parsed and validated against `assets/schema_input.json`
via the [nf-schema](https://nextflow-io.github.io/nf-schema/) plugin — the
same mechanism current nf-core pipelines (e.g. nf-core/mag) use for their
`--input` sample sheet. It checks required columns, that `fastq_1`/`fastq_2`
point to existing `.fastq.gz`/`.fq.gz` files, and that sample IDs are
unique, and fails fast with a clear error if not. Each sample's resolved
`id`/`group` then travels through the whole pipeline as a small Groovy map
(`meta`) attached to every channel item, rather than as separate
positional fields — again the same pattern nf-core pipelines use.

Trimming (TrimGalore) is run by the pipeline itself — see Step 2 below.
Every other pre-computed input (STAR unmapped reads, BBsplit outputs,
Kraken2 output, BLAST output) is still located via the directory params in
`params.yml` plus a fixed per-sample filename convention — see the comments
in `params.yml` for the exact expected filenames.

---

## Pipeline steps and outputs

```
Step 0   RAW-READ QC (FastQC + MultiQC, run via Singularity)
           Skippable via skip_fastqc. First thing done to the raw reads.
           → reports/00_fastqc_raw/{sample}_{1,2}_fastqc.{html,zip}
           → reports/00_fastqc_raw/multiqc_report.html

Step 1   Raw reads (fq.gz, paths from samples.csv fastq_1/fastq_2)
           → reports/01_raw/{sample}_raw_stats.tsv

Step 2   ADAPTER/QUALITY TRIMMING (TrimGalore, run by the pipeline)
           Runs with TrimGalore's own defaults (quality 20, stringency 1,
           length 20 — see params.yml to adjust, or trim_extra_args for
           anything else, e.g. explicit adapters). Skippable via
           skip_trimming (falls back to pre-trimmed reads in trimmed_dir).
           Also runs TrimGalore's own --fastqc, aggregated by a second
           MultiQC report (skippable via skip_fastqc, same as Step 0).
           → trimmed/{sample}_{1,2}_trimmed.fq.gz
           → reports/02_trimmed/{sample}_{1,2}.fq.gz_trimming_report.txt
           → reports/02_trimmed/{sample}_{1,2}_val_{1,2}_fastqc.{html,zip}
           → reports/02_trimmed/multiqc_report.html
           → reports/02_trimmed/{sample}_trimmed_stats.tsv
           → reports/reads_posttrim_tab.tsv            ← produced fresh here

Step 3   STAR unmapped reads (pre-computed)
           → reports/03_star/{sample}_star_unmapped_stats.tsv
           → reports/03_star/{sample}_Log.final.out

Step 4a  BBsplit host/human filtered (pre-computed)
           → reports/04a_bbsplit_host/{sample}_bbsplit_host_stats.tsv

Step 4b  BBsplit antiviral filtered (pre-computed)
           → reports/04b_bbsplit_virus/{sample}_bbsplit_virus_stats.tsv

Step 5   BBsplit MAG-clean reads (pre-computed)
           → reports/05_no_mags/{sample}_noMAGs_stats.tsv

Step 6   KRAKEN annotation + filtering
           annotate_kraken_lineage.py  →  taxonkit appends lineage col 6
           filter_kraken_invertebrates.py:
             - KEEPS C (classified) rows only
             - U rows silently skipped → handled by BLAST path (not lost)
             - Removes: viruses, invertebrates, host genera, length < 18 nt
           → kraken_annotated/{sample}_kraken_annotated.tsv
           → kraken_filtered/{sample}_kraken_invertebrates_filtered.tsv

Step 7   BLAST annotation + filtering (mates 1 and 2 independently)
           annotate_blast_lineage.py   →  taxonkit appends 8 taxonomy cols
           filter_blast_all_conditions.py:
             - Removes: viruses, invertebrates, host species
             - Removes: mismatch≠0, gapopen≠0, not-full-span, length<18 nt
           → blast_annotated/{sample}_{1,2}_blast_annotated.tsv
           → blast_filtered/{sample}_{1,2}_blast_all_lengths_filtered.tsv

Step 8   COLLAPSE (extract from trimmed reads → merge → seqkit rmdup)
           Passing IDs from BLAST and Kraken filtered TSVs are used to fetch
           sequences from trimmed FASTQ files. Mate 1 and mate 2 BLAST-passing
           reads are fetched separately; Kraken-passing IDs are fetched from
           mate 1 only. All three are cat-merged and collapsed with seqkit rmdup.
           → collapsed/{sample}_merged.fq
           → collapsed/{sample}_merged_collapsed_clean.fq
           → collapsed/{sample}_merged_duplicated.detail.txt
           → collapsed/{sample}_collapse_stats.tsv

Step 9   CANDIDATE SELECTION (ge5_18nt_filtering_pipeline_for_both_mates.py)
           Candidate = ≥5 duplicate occurrences AND length ≥ 18 nt.
           Lengths read from collapsed_clean.fq; weights from detail.txt.
           → candidates/{sample}_ge5_detected_ids.txt
           → candidates/{sample}_ge5_detected.hist          (weighted by dup count)
           → candidates/{sample}_ge5_detected_weighted_ids.tsv

Step 10  EXTRACT FASTA (seqkit grep + fq2fa on collapsed_clean.fq)
           → final_transRNAs_fasta/{sample}_final_transRNAs.fasta

Step 11  LENGTH HISTOGRAMS + PLOT
           → final_transRNAs_length_hists/{sample}_ge5_detected.hist
           → plots/{date}_avg_reinflated_length_distribution_15to100nt.png

Step 12  TAXONOMY PARSER  (called ONCE GLOBALLY — all samples together)
           04_05_2026_transRNA_taxonomy_parser.py:
             - Re-inflates using duplicate weights
             - RPM-normalises using reads_posttrim_tab.tsv
             - Priority: blast (configurable)
             - Produces per-sample AND per-group tables (groups from samples.csv)
           → final_transRNAs_taxonomies/{sample}_{blast,kraken,combined}_{Rank}_summary.tsv
           → final_transRNAs_taxonomies/{sample}_{blast,kraken,combined}_{Rank}_top10.tsv
           → final_transRNAs_taxonomies/{group}_{blast,kraken,combined}_{Rank}_summary.tsv
           → final_transRNAs_taxonomies/{group}_{blast,kraken,combined}_{Rank}_top10.tsv
             (skipped for a group with only one sample — its per-sample
             output above already covers it; avoids the two colliding
             since a solo group's label defaults to that sample's own ID)
           → final_transRNAs_taxonomies/fetch_reinflate_report.tsv

Step 13  TAXONOMY PLOT (plot_top10_taxa_global_colors.py)
           Reads each group's {group}_*_top10.tsv files.
           Multi-panel broken-axis horizontal bar plot; ranks: Domain/Kingdom/Order/Species
           → plots/*.png  (one column per group)

Step 14  AGGREGATE READ-COUNT REPORT
           → reports/pipeline_read_counts_report.tsv
           → reports/pipeline_read_counts_report.html
           → reports/dataset_summary_report.tsv
           → reports/virus_exclusion_report.tsv
```

---

## Filtering criteria — complete specification

### Kraken (`filter_kraken_invertebrates.py`)

| | Criterion | How checked |
|---|---|---|
| **SKIP** | Unclassified reads (U) | C/U flag col 1 — these go to BLAST path |
| **REMOVE** | Either mate < 18 nt | length field `"35\|35"`, both sides checked |
| **REMOVE** | Virus/acellular domain | Lineage[0] (Domain rank) |
| **REMOVE** | Host genera | Lineage[6] (Genus) ∈ {Homo, Mus, Canis, Felis, Apis} |
| **REMOVE** | Invertebrate phyla | Regex on full lineage vs invertebrate_phyla.txt names |

Lineage appended as col 6: `Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species`  
(taxonkit `{k};{K};{p};{c};{o};{f};{g};{s}` with `--fill-miss-rank`)

### BLAST (`filter_blast_all_conditions.py`)

| Criterion | Field |
|---|---|
| Virus/acellular domain | last-8 col 0 (Domain) |
| Invertebrate phylum | last-8 col 2 (Phylum), exact match |
| Host species | last-8 col 7 (Species), substring: homo/mus/canis lupus/felis/apis mellifera |
| mismatch ≠ 0 | col 5 |
| gapopen ≠ 0 | col 6 |
| Not full-span: abs(qend−qstart)+1 ≠ length | cols 4,7,8 |
| Alignment length < 18 nt | col 4 |

---

## Taxonomy lineage — both tools

```bash
echo "$taxid" \
  | taxonkit lineage --data-dir $DB \
  | taxonkit reformat --data-dir $DB \
      --format "{k};{K};{p};{c};{o};{f};{g};{s}" \
      --fill-miss-rank --miss-taxid-repl unassigned
```

8 ranks always present: **Domain ; Kingdom ; Phylum ; Class ; Order ; Family ; Genus ; Species**  
Kingdom is rank index 1. All 8 ranks appear in all taxonomy tables and plots.

---

## taxonkit setup (once only)

```bash
mkdir -p ~/.taxonkit && cd ~/.taxonkit
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
echo "9606" | taxonkit lineage   # test: should print Homo sapiens lineage
```

---

## Software requirements

```
nextflow/24.04.4
conda or mamba          (mamba strongly recommended — see below)
singularity/apptainer   (only for the FASTQC/MULTIQC processes)
```

That's it — every actual bioinformatics tool (seqkit, taxonkit, trim-galore,
cutadapt, pigz, pandas, numpy, matplotlib, seaborn...) is *not* a manual
install. Each process declares its own conda env via a `conda "..."`
directive pointing at a small, version-pinned file in `config/envs/`:

| Env file | Used by | Contains |
|---|---|---|
| `trim_galore.yaml` | `TRIMGALORE` | trim-galore, cutadapt, pigz |
| `seqkit.yaml` | `COUNT_*`, `EXTRACT_FASTA` | seqkit |
| `collapse.yaml` | `COLLAPSE_READS` | seqkit, python (calls both) |
| `taxonkit.yaml` | `ANNOTATE_KRAKEN`, `ANNOTATE_BLAST` | taxonkit, python (scripts shell out to taxonkit) |
| `python_core.yaml` | `FILTER_KRAKEN`, `FILTER_BLAST`, `SELECT_CANDIDATES` | python only (pure-stdlib scripts) |
| `python_analysis.yaml` | `MAKE_READS_POSTTRIM_TAB`, `TAXONOMY_PARSER`, `PLOT_LENGTH_DIST`, `PLOT_TAXONOMY`, `AGGREGATE_REPORT` | python, pandas, numpy, matplotlib, seaborn |

Nextflow builds each env itself the first time it's needed (`conda.enabled`
in `cambridge.config`/`nextflow.config`'s `local` profile) and caches it —
no manual `conda install`, no shared pre-built env to hand-configure. Only
`FASTQC`/`MULTIQC` are the exception: they run via their own Singularity
container instead (declared directly on those two processes), so Singularity
is needed too, but only for those.

---

## HPC tips (CSD3)

- Set `workDir` in `nextflow.config` to a path outside your home quota (scratch or RDS)
- `process.cache = 'lenient'` prevents spurious cache misses from RDS NFS jitter
- For very large BLAST TSVs, switch `ANNOTATE_BLAST` to `icelake-himem` in `cambridge.config`
- Monitor jobs: `squeue -u <your_username>` and `tail -f logs/nextflow_*.log`
- The pipeline declares the `nf-schema` plugin (`nextflow.config`), used to validate
  `samples.csv`. Nextflow downloads it automatically the first time the pipeline runs,
  so the login node needs outbound internet access on that first run; after that it's
  cached locally (`~/.nextflow/plugins`) and no further downloads are needed.
- Similarly, `FASTQC`/`MULTIQC`'s Singularity images are pulled on first use and cached under
  `NXF_SINGULARITY_CACHEDIR` (or `~/.singularity/cache` if that env var isn't set)
  — set `NXF_SINGULARITY_CACHEDIR` to somewhere outside your home quota
  (e.g. `export NXF_SINGULARITY_CACHEDIR=<path>/hpc-work/singularity_cache` in
  `~/.bashrc`) before the first run.
- Same idea for the per-process conda envs (`config/envs/*.yaml`): set
  `params.conda_cache_dir` in `params.yml` to a persistent path outside your
  home quota before the first run — each env is built there once (a few
  minutes total) and reused after that, on every subsequent run.
