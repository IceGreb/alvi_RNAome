*Manuscript in preparation.* Raw sequencing data and intermediate files will be made publicly available upon publication.

This repository documents the steps used to process RNA-seq data from *Apis mellifera*
and identify transmissible RNA candidates. Prior setup of tools and environments is
required (see _Config_ and _Software requirements_ sections below). For questions or
clarifications, contact Nikolaos Vergoulidis (nv358@cam.ac.uk). The code is provided
for transparency and reproducibility of the analysis.

# Honeybee *Apis mellifera* Transmissible RNAome

The aim of this repository is to document the steps used to process raw paired-end RNA-seq
reads from Royal Jelly and systemic larval tissue samples and identify candidate transmissible
RNA sequences. Part of the process was carried out using nf-core/rnaseq and nf-core/mag. The custom
transRNA identification pipeline is implemented in Nextflow (DSL2) v1.2 and parallelized on
the CSD3 icelake cluster at the University of Cambridge (kernel: Linux 4.18.0-553.125.1.el8_10.x86_64).
Information about software versions and conda environments is found in the `config/` directory
described in the _Config_ section below.

## Overview

The `workflow/` directory contains the Nextflow pipelines used for bioinformatic processing.
The `workflow/nf-core/` subdirectory holds configuration for nf-core/rnaseq (trimming and host
alignment) and nf-core/mag (metagenome-assembled genome removal). The main transRNA identification
pipeline is in `workflow/transRNA_pipeline/`, which also contains all `bin/` scripts and a
detailed step-by-step README.

The `downstream_analysis/` directory contains scripts for post-pipeline analysis and figure
generation.

The `config/` directory contains shared configuration files, environment specifications, and
metadata used across pipelines.

All paths are specified relative to the working directory unless absolute paths are required.
The pipeline was parallelized on CSD3 using SLURM via the `cambridge.config` profile.

## Raw data

### DNA sequencing

For each sequenced DNA sample, two systemic larval DNA extracts were pooled to create a single
DNA library.

### RNA sequencing

* **Mid-RNAseq** — Novogene
* **RISC-seq** — Cambridge Genomics Service, University of Cambridge
* **RNAseq (differential expression)** — Novogene

### Samples

Two sample groups were sequenced for this study:

| Group | Description | Samples |
|-------|-------------|---------|
| **RJ** | Royal Jelly | RJ1, RJ2, RJ3 |
| **ST** | Systemic larval tissue | T1GMN, T4GMN, TA1, TB1, TC1, TD1 |

This reflects the dataset used in the manuscript, not a constraint of the
pipeline: samples, groups, and raw read paths are all defined in
`workflow/transRNA_pipeline/samples.csv`, and the pipeline supports any
number of groups with any labels — see the pipeline README's *Samples*
section for the sample sheet format.

## Code Overview

### nf-core workflows

* `workflow/nf-core/rnaseq/` — nf-core/rnaseq configuration for adapter trimming (TrimGalore)
  and host genome alignment (STAR). STAR-unmapped reads are carried forward as input to the
  transRNA pipeline.
* `workflow/nf-core/mag/` — nf-core/mag configuration for metagenome-assembled genome (MAG)
  identification and removal via BBsplit, reducing contamination from environmental microbiota.

### transRNA pipeline

* `workflow/transRNA_pipeline/main.nf` — Main Nextflow (DSL2) workflow. Orchestrates all 14
  steps from read counting through taxonomy plotting and aggregate reporting. Processes all
  samples in parallel across RJ and ST groups.
* `workflow/transRNA_pipeline/nextflow.config` — Global Nextflow configuration: profiles,
  process labels, work directory, cache mode (`lenient` for RDS NFS jitter), and pipeline
  reports (HTML, timeline, DAG, trace).
* `workflow/transRNA_pipeline/cambridge.config` — CSD3-specific SLURM executor settings,
  partition and resource assignments per process label (`count_only`, `low`, `med`).
* `workflow/transRNA_pipeline/params.yml` — All user-configurable parameters: input data paths,
  output directory, filtering thresholds (`min_occ`, `min_len`), taxonomy priority, plot
  parameters, and skip flags for modular reruns.
* `workflow/transRNA_pipeline/samples.csv` — Sample sheet listing sample IDs and group
  assignments (RJ/ST).
* `workflow/transRNA_pipeline/run.sh` — Launch script for CSD3. Submits the Nextflow master
  job to SLURM; supports `--resume` for restarting after failure.

### Scripts (`workflow/transRNA_pipeline/bin/`)

* `annotate_kraken_lineage.py` — Extracts taxids from Kraken2 output (field 3) and appends
  an 8-rank lineage column (Domain→Species) via taxonkit.
* `filter_kraken_invertebrates.py` — Filters Kraken-annotated TSVs: keeps classified (C) reads
  only; removes viruses, invertebrates, host genera, and reads shorter than `min_len`. Emits a
  per-sample filter stats TSV.
* `annotate_blast_lineage.py` — Appends 8 lineage columns to BLAST tabular output using the
  taxid in column 13.
* `filter_blast_all_conditions.py` — Filters BLAST-annotated TSVs: removes viruses,
  invertebrates, host species, reads with mismatches or gaps, non-full-span alignments, and
  reads shorter than `min_len`. Emits a per-sample filter stats TSV.
* `collapse_stats.py` — Summarises `seqkit rmdup` output into a per-sample collapse stats TSV
  (merged reads, total groups, ge5 groups).
* `ge5_18nt_filtering_pipeline_for_both_mates.py` — Applies the ≥`min_occ` duplicate and
  ≥`min_len` length thresholds to the collapsed read pool to produce the final candidate set.
  Outputs candidate IDs, a weighted histogram, and a weighted IDs TSV.
* `make_reads_posttrim_tab.py` — Aggregates per-sample trimmed-read seqkit stats into a single
  `reads_posttrim_tab.tsv` used for RPM normalisation in the taxonomy parser.
* `04_05_2026_transRNA_taxonomy_parser.py` — Global taxonomy parser (called once across all
  samples). Re-inflates candidate sequences using duplicate weights, RPM-normalises by trimmed
  reads, and produces per-sample and per-group (RJ/ST) taxonomy summary and top-10 TSVs at
  four ranks (Domain, Kingdom, Order, Species).
* `plot_top10_taxa_global_colors.py` — Generates multi-panel broken-axis horizontal bar plots
  of top-10 taxa for RJ and ST groups across all ranks.
* `plot_lengths.py` — Plots weighted length distributions (15–100 nt) of candidate transRNAs
  per sample and as a dataset average.
* `aggregate_report.py` — Aggregates per-step read counts into a single per-sample and
  per-group TSV/HTML report. Metrics include raw → trimmed → STAR → BBsplit → Kraken/BLAST
  classification → candidate selection rates (RPM-normalised).
* `report_to_html.py` — Converts the aggregate TSV report to a styled HTML table.

### Config

* `config/` — Configuration files and conda environment specifications for the workflows and
  downstream analysis scripts.

## Some useful intermediate files

Several large output directories are not included in this repository but are available through
the accompanying [Zenodo](https://doi.org/10.5281/zenodo.XXXXXXX) repository.

Key outputs per run:

* `reports/pipeline_read_counts_report.tsv` and `.html` — Per-sample read counts at every
  depletion step with RPM-normalised transRNA metrics
* `reports/dataset_summary_report.tsv` — Per-group (RJ/ST) averages across all metrics
* `final_transRNAs_taxonomies/` — Per-sample and per-group taxonomy tables at four ranks;
  `fetch_reinflate_report.tsv` summarises re-inflation counts
* `final_transRNAs_fasta_collapsed/` — FASTA files of representative transmissible RNA
  sequences per sample
* `plots/` — Length distribution and taxonomy bar plots
