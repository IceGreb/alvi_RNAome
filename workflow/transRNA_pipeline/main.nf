#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
================================================================================
  transRNA Pipeline  v1.2
  Paired-end RNA-seq → filtering → taxonomy → final transRNA identification

  Samples and groups are defined entirely by the sample sheet (see
  params.sample_sheet / samples.csv: sample,group,fastq_1,fastq_2). Any
  number of groups, with any labels, is supported. group is optional — a
  sample with no group is treated as its own group (see resolveGroup()).

  

  Key design notes:
    • Read counts reported at EVERY step (raw → trimmed → STAR → bbsplit ×3)
    • Collapse runs on the final clean (no_MAGs) reads; counts reported there
    • seqkit rmdup -D writes duplicated.detail.txt directly — no reformatting needed
    • Kraken C-reads → Kraken filter path
      Kraken U-reads → silently passed; handled by the BLAST path
    • BLAST path is independent and processes all reads in the BLAST TSV
    • Candidate set = UNION(Kraken-passing, BLAST-passing) with ≥5 dups & ≥18 nt
    • reads_posttrim_tab.tsv produced fresh from step-2 seqkit stats
    • Taxonomy parser called ONCE GLOBALLY (faster, designed for directory scanning)
    • top10 TSVs produced by taxonomy parser → fed directly to plotter
================================================================================
*/

include { samplesheetToList } from 'plugin/nf-schema'
include { FASTQC            } from './modules/nf-core/fastqc/main.nf'

// ── Helpers ───────────────────────────────────────────────────────────────────

def isBlank(v) {
    return v == null || v == [] || (v instanceof String && v.trim() == "")
}

// group is optional in the schema; a sample with no group becomes its own
// group, so every downstream per-group summary/plot still has something to
// key on without any pipeline logic needing to special-case "no group".
def resolveGroup(meta) {
    return isBlank(meta.group) ? meta.id : meta.group
}

// ============================================================================
//  READ COUNT REPORTS  (seqkit stats -T at every pre-computed step)
// ============================================================================

process COUNT_RAW {
    tag "${meta.id}"
    label 'count_only'
    publishDir "${params.outdir}/reports/01_raw", mode: 'copy'
    input:  tuple val(meta), path(reads)
    output: tuple val(meta), path("${meta.id}_raw_stats.tsv")
    script:
    """
    seqkit stats -T -j ${task.cpus} \
        ${reads} \
        > ${meta.id}_raw_stats.tsv
    """
}

process COUNT_TRIMMED {
    tag "${meta.id}"
    label 'count_only'
    publishDir "${params.outdir}/reports/02_trimmed", mode: 'copy'
    input:  val(meta)
    output: tuple val(meta), path("${meta.id}_trimmed_stats.tsv")
    script:
    """
    seqkit stats -T -j ${task.cpus} \
        ${params.trimmed_dir}/${meta.id}_1_trimmed.fq.gz \
        ${params.trimmed_dir}/${meta.id}_2_trimmed.fq.gz \
        > ${meta.id}_trimmed_stats.tsv
    """
}

process COUNT_STAR_UNMAPPED {
    tag "${meta.id}"
    label 'count_only'
    publishDir "${params.outdir}/reports/03_star", mode: 'copy'
    input:  val(meta)
    output: tuple val(meta),
                  path("${meta.id}_star_unmapped_stats.tsv"),
                  path("${meta.id}_Log.final.out")
    script:
    """
    seqkit stats -T -j ${task.cpus} \
        ${params.star_dir}/${meta.id}/${meta.id}_Unmapped.out.mate1 \
        ${params.star_dir}/${meta.id}/${meta.id}_Unmapped.out.mate2 \
        > ${meta.id}_star_unmapped_stats.tsv

    cp ${params.star_dir}/${meta.id}/Log.final.out \
       ${meta.id}_Log.final.out
    """
}

process COUNT_BBSPLIT_HOST {
    tag "${meta.id}"
    label 'count_only'
    publishDir "${params.outdir}/reports/04a_bbsplit_host", mode: 'copy'
    input:  val(meta)
    output: tuple val(meta), path("${meta.id}_bbsplit_host_stats.tsv")
    script:
    """
    seqkit stats -T -j ${task.cpus} \
        ${params.bbsplit_host_dir}/Filtered_${meta.id}_2MM_bbsplit_clean1.fq \
        ${params.bbsplit_host_dir}/Filtered_${meta.id}_2MM_bbsplit_clean2.fq \
        > ${meta.id}_bbsplit_host_stats.tsv
    """
}

process COUNT_BBSPLIT_VIRUS {
    tag "${meta.id}"
    label 'count_only'
    publishDir "${params.outdir}/reports/04b_bbsplit_virus", mode: 'copy'
    input:  val(meta)
    output: tuple val(meta), path("${meta.id}_bbsplit_virus_stats.tsv")
    script:
    """
    seqkit stats -T -j ${task.cpus} \
        ${params.bbsplit_virus_dir}/Clean_${meta.id}_2MM_bbsplit_1.fq \
        ${params.bbsplit_virus_dir}/Clean_${meta.id}_2MM_bbsplit_2.fq \
        > ${meta.id}_bbsplit_virus_stats.tsv
    """
}

process COUNT_BBSPLIT_MAGS {
    tag "${meta.id}"
    label 'count_only'
    publishDir "${params.outdir}/reports/05_no_mags", mode: 'copy'
    input:  val(meta)
    output: tuple val(meta), path("${meta.id}_noMAGs_stats.tsv")
    script:
    """
    seqkit stats -T -j ${task.cpus} \
        ${params.bbsplit_mags_dir}/no_MAGs_${meta.id}_2MM_clean1.fq \
        ${params.bbsplit_mags_dir}/no_MAGs_${meta.id}_2MM_clean2.fq \
        > ${meta.id}_noMAGs_stats.tsv
    """
}

// ============================================================================
//  READS POSTTRIM TABLE
//  Produced fresh from the trimmed-read seqkit stats. Written to reports/ so
//  it is available as input to the taxonomy parser.
// ============================================================================

process MAKE_READS_POSTTRIM_TAB {
    label 'count_only'
    publishDir "${params.outdir}/reports", mode: 'copy'
    input:  path(trimmed_stats_files)   // all *_trimmed_stats.tsv collected
    output: path("reads_posttrim_tab.tsv")
    script:
    """
    python3 ${projectDir}/bin/make_reads_posttrim_tab.py
    """
}

// ============================================================================
//  COLLAPSE READS
//  seqkit rmdup -s -D on the clean no_MAGs reads.
//  -D writes duplicated.detail.txt directly in the exact format needed by
//  the taxonomy parser and ge5 pipeline — no reformatting step required.
// ============================================================================

// =============================================================================
//  COLLAPSE_READS
//
//  New logic (replacing no_MAGs-based collapse):
//    1. Extract passing seq IDs from BLAST filtered TSVs (col 1, per mate)
//    2. Extract passing seq IDs from Kraken filtered TSV (col 2, mate 1 only)
//    3. seqkit grep each ID set against the corresponding trimmed FASTQ
//    4. cat blast_mate1 + blast_mate2 + kraken_mate1 → merged.fq
//    5. seqkit rmdup -s -D on merged.fq
//
//  This means only BLAST/Kraken-validated reads enter the duplicate pool,
//  so any sequence reaching ≥ min_occ duplicates is already a candidate.
//  No further cross-referencing is needed in SELECT_CANDIDATES.
// =============================================================================

process COLLAPSE_READS {
    tag "${meta.id}"
    label 'med'
    publishDir "${params.outdir}/collapsed", mode: 'copy'
    input:
    tuple val(meta),
          path(blast1_tsv),
          path(blast2_tsv),
          path(kraken_tsv)
    output:
    tuple val(meta),
          path("${meta.id}_merged_duplicated.detail.txt"),
          path("${meta.id}_merged_collapsed_clean.fq"),
          path("${meta.id}_merged.fq"),
          path("${meta.id}_collapse_stats.tsv")
    script:
    """
    # ── Extract passing IDs from BLAST filtered TSVs (col 1, no header) ──────
    awk -F'\t' '{print \$1}' ${blast1_tsv} | sort -u > blast1_ids.txt
    awk -F'\t' '{print \$1}' ${blast2_tsv} | sort -u > blast2_ids.txt

    # ── Extract passing IDs from Kraken filtered TSV (col 2, no header) ──────
    awk -F'\t' '{print \$2}' ${kraken_tsv} | sort -u > kraken_ids.txt

    # ── Fetch sequences from trimmed reads ────────────────────────────────────
    seqkit grep -j ${task.cpus} -f blast1_ids.txt \
        ${params.trimmed_dir}/${meta.id}_1_trimmed.fq.gz > blast1_passing.fq
    seqkit grep -j ${task.cpus} -f blast2_ids.txt \
        ${params.trimmed_dir}/${meta.id}_2_trimmed.fq.gz > blast2_passing.fq
    seqkit grep -j ${task.cpus} -f kraken_ids.txt \
        ${params.trimmed_dir}/${meta.id}_1_trimmed.fq.gz > kraken_passing.fq

    # ── Merge and collapse ────────────────────────────────────────────────────
    cat blast1_passing.fq blast2_passing.fq kraken_passing.fq \
        > ${meta.id}_merged.fq

    seqkit rmdup -s \
        -j ${task.cpus} \
        -o ${meta.id}_merged_collapsed_clean.fq \
        -d /dev/null \
        -D ${meta.id}_merged_duplicated.detail.txt \
        ${meta.id}_merged.fq

    # ── Collapse summary ──────────────────────────────────────────────────────
    MERGED_READS=\$(awk 'END{print NR/4}' ${meta.id}_merged.fq)
    python3 ${projectDir}/bin/collapse_stats.py ${meta.id} ${params.min_occ} \${MERGED_READS}
    """
}

// ============================================================================
//  KRAKEN ANNOTATION + FILTERING
//
//  annotate_kraken_lineage.py:
//    Extracts taxid from field 3 "Name (taxid N)", runs taxonkit, appends
//    lineage as column 6: Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species
//    Both C and U rows are annotated and passed through unchanged.
//
//  filter_kraken_invertebrates.py:
//    • KEEPS only C (classified) rows
//    • U rows silently skipped — they are NOT lost; handled by BLAST path
//    • Removes: viruses (Domain), invertebrates (Phylum), host genera
//    • Removes: either mate length < min_len nt
// ============================================================================

process ANNOTATE_KRAKEN {
    tag "${meta.id}"
    label 'med'
    //publishDir "${params.outdir}/kraken_annotated", mode: 'copy'
    input:  val(meta)
    output: tuple val(meta), path("${meta.id}_kraken_annotated.tsv")
    script:
    """
    python3 ${projectDir}/bin/annotate_kraken_lineage.py \
        --input       ${params.kraken_dir}/${meta.id}_2MM_CLEAN_noMAGs_kraken_output_005_.txt \
        --output      ${meta.id}_kraken_annotated.tsv \
    """
}

process FILTER_KRAKEN {
    tag "${meta.id}"
    label 'low'
    publishDir "${params.outdir}/kraken_filtered", mode: 'copy'
    input:
    tuple val(meta), path(annotated_tsv)
    output:
    tuple val(meta),
          path("${meta.id}_kraken_invertebrates_filtered.tsv"),
          path("${meta.id}_kraken_filter_stats.tsv")             // ← new
    script:
    """
    python3 ${projectDir}/bin/filter_kraken_invertebrates.py \
        --min-len   ${params.min_len} \
        --stats-out ${meta.id}_kraken_filter_stats.tsv \
        ${params.invertebrate_phyla} \
        ${annotated_tsv} \
        ${meta.id}_kraken_invertebrates_filtered.tsv
    """
}

// ============================================================================
//  BLAST ANNOTATION + FILTERING  (mate 1 and mate 2 independently)
//
//  annotate_blast_lineage.py:
//    taxid from col 13 (may be semicolon-delimited; first used).
//    Appends 8 columns at end: Domain Kingdom Phylum Class Order Family Genus Species
//
//  filter_blast_all_conditions.py:
//    Removes: viruses (Domain), invertebrate phyla (Phylum), host (Species)
//    Removes: mismatch≠0, gapopen≠0, not-full-span, length<min_len
// ============================================================================

process ANNOTATE_BLAST {
    tag "${meta.id}_${mate}"
    label 'med'
    //publishDir "${params.outdir}/blast_annotated", mode: 'copy'
    input:  tuple val(meta), val(mate)
    output: tuple val(meta), val(mate), path("${meta.id}_${mate}_blast_annotated.tsv")
    script:
    """
    python3 ${projectDir}/bin/annotate_blast_lineage.py \
        --input       ${params.blast_dir}/${meta.id}_${mate}.tsv \
        --output      ${meta.id}_${mate}_blast_annotated.tsv \
        --taxid-col   13
    """
}

process FILTER_BLAST {
    tag "${meta.id}_${mate}"
    label 'low'
    publishDir "${params.outdir}/blast_filtered", mode: 'copy'
    input:
    tuple val(meta), val(mate), path(annotated_tsv)
    output:
    tuple val(meta), val(mate),
          path("${meta.id}_${mate}_blast_all_lengths_filtered.tsv"),
          path("${meta.id}_${mate}_blast_filter_stats.tsv")      // ← new
    script:
    """
    python3 ${projectDir}/bin/filter_blast_all_conditions.py \
        --min-len   ${params.min_len} \
        --stats-out ${meta.id}_${mate}_blast_filter_stats.tsv \
        ${params.invertebrate_phyla} \
        ${annotated_tsv} \
        ${meta.id}_${mate}_blast_all_lengths_filtered.tsv
    """
}

// ============================================================================
//  CANDIDATE SELECTION
//  Since COLLAPSE_READS now collapses only BLAST/Kraken-validated reads,
//  any sequence reaching ≥ min_occ duplicates is already a candidate.
//  The ge5 script no longer cross-references BLAST/Kraken TSVs — it simply
//  applies the ≥ min_occ threshold and reads lengths from the collapsed FASTQ.
// ============================================================================

process SELECT_CANDIDATES {
    tag "${meta.id}"
    label 'low'
    publishDir "${params.outdir}/candidates", mode: 'copy'
    input:
    tuple val(meta),
          path(dup_detail),
          path(collapsed_clean_fq),
          path(collapse_stats)
    output:
    tuple val(meta),
          path("${meta.id}_ge5_detected_ids.txt"),
          path("${meta.id}_ge5_detected.hist"),
          path("${meta.id}_ge5_detected_weighted_ids.tsv")
    script:
    """
    python3 ${projectDir}/bin/ge5_18nt_filtering_pipeline_for_both_mates.py \
        --mode       filtered \
        --dup-file   ${dup_detail} \
        --fq-file    ${collapsed_clean_fq} \
        --sample     ${meta.id} \
        --output-dir . \
        --min-occ    ${params.min_occ} \
        --min-len    ${params.min_len}
    """
}

// ============================================================================
//  EXTRACT FINAL FASTA
//  Unique read IDs from the ge5 detected ID files → seqkit grep against
//  the clean no_MAGs mate-1 reads → FASTA.
//  Re-inflation is done in the taxonomy step using dup weights.
// ============================================================================

process EXTRACT_FASTA {
    tag "${meta.id}"
    label 'low'
    publishDir "${params.outdir}/final_transRNAs_fasta_collapsed", mode: 'copy'
    input:
    tuple val(meta),
          path(ids),
          path(hist),
          path(wids),
          path(merged_fq)
    output:
    tuple val(meta),
          path("${meta.id}_final_transRNAs.fasta"),
          path("${meta.id}_fasta_fetch_report.txt")
    script:
    """
    EXPECTED=\$(wc -l < ${ids})

    # Grep from the same merged fq used for collapse.
    # seqkit rmdup picks representative IDs from this file,
    # so all candidate IDs are guaranteed to be present here.
    seqkit grep -f ${ids} \
        -j ${task.cpus} \
        ${merged_fq} \
        | seqkit fq2fa \
        > ${meta.id}_final_transRNAs.fasta

    FETCHED=\$(grep -c '^>' ${meta.id}_final_transRNAs.fasta || echo 0)

    echo "Sample: ${meta.id}"                           > ${meta.id}_fasta_fetch_report.txt
    echo "Expected IDs: \${EXPECTED}"                >> ${meta.id}_fasta_fetch_report.txt
    echo "Sequences fetched: \${FETCHED}"            >> ${meta.id}_fasta_fetch_report.txt
    if [ "\${FETCHED}" -eq "\${EXPECTED}" ]; then
        echo "Status: OK - all sequences fetched"     >> ${meta.id}_fasta_fetch_report.txt
    else
        MISSING=\$(( EXPECTED - FETCHED ))
        echo "Status: WARNING - \${MISSING} IDs not found in merged fq" >> ${meta.id}_fasta_fetch_report.txt
    fi
    cat ${meta.id}_fasta_fetch_report.txt
    """
}

// ============================================================================
//  LENGTH HISTOGRAMS — publish per-sample hists and collect for global plot
// ============================================================================

process PUBLISH_HISTS {
    tag "${meta.id}"
    label 'count_only'
    publishDir "${params.outdir}/final_transRNAs_length_hists", mode: 'copy'
    input:
    tuple val(meta),
          path(ids), path(hist), path(wids)
    output:
    path("${meta.id}_ge5_detected.hist")
    script:
    """
    ls ${hist}
    """
}

process PLOT_LENGTH_DIST {
    label 'count_only'
    publishDir "${params.outdir}/plots", mode: 'copy'
    input:  path(hist_files)
            path(sample_sheet)
    output: path("*.png")          // ← collect all PNGs (both plots)
    script:
    """
    python3 ${projectDir}/bin/plot_lengths.py . --samplesheet ${sample_sheet}
    """
}

// ============================================================================
//  TAXONOMY PARSER — called ONCE GLOBALLY
//
//  Inputs collected from all samples:
//    --ids-dir    : folder with all *_ge5_detected_ids.txt files
//    --blast-dir  : folder with all *_blast_all_lengths_filtered.tsv files
//    --kraken-dir : folder with all *_kraken_invertebrates_filtered.tsv files
//    --dup-dir    : folder with all no_MAGs_*_duplicated.detail.txt files
//    --reads-table: reads_posttrim_tab.tsv (fresh from this run)
//    --samplesheet: resolved (post-group-default) sample,group table
//
//  Outputs per sample:  {sample}_{blast,kraken,combined}_{Rank}_summary.tsv
//                       {sample}_{blast,kraken,combined}_{Rank}_top10.tsv
//  Outputs per group:   {group}_{blast,kraken,combined}_{Rank}_summary.tsv
//                       {group}_{blast,kraken,combined}_{Rank}_top10.tsv
//                       (skipped for a group with only one sample — its
//                       per-sample output already covers it; see script)
//  Report:              fetch_reinflate_report.tsv
// ============================================================================

process TAXONOMY_PARSER {
    label 'med'
    publishDir "${params.outdir}/final_transRNAs_taxonomies", mode: 'copy'
    input:
    // ids_files: all *_ge5_detected_ids.txt from all samples
    path(ids_files)
    // kraken_files: all *_kraken_invertebrates_filtered.tsv
    path(kraken_files)
    // blast_files: all *_blast_all_lengths_filtered.tsv
    path(blast_files)
    // dup_files: all no_MAGs_*_duplicated.detail.txt
    path(dup_files)
    // reads table produced by MAKE_READS_POSTTRIM_TAB
    path(reads_table)
    // resolved sample sheet, so the parser can look up each sample's group
    path(sample_sheet)
    output:
    path("*_top10.tsv"),           emit: top10_tsvs
    path("*_summary.tsv"),         emit: summary_tsvs
    path("fetch_reinflate_report.tsv"), emit: report
    script:
    // All files are staged into the work dir by Nextflow.
    // Point all --*-dir args to "." so the parser scans the work dir.
    """
    python3 ${projectDir}/bin/04_05_2026_transRNA_taxonomy_parser.py \
        --ids-dir     . \
        --blast-dir   . \
        --kraken-dir  . \
        --dup-dir     . \
        --reads-table ${reads_table} \
        --samplesheet ${sample_sheet} \
        --outdir      . \
        --priority    ${params.priority}
    """
}

// ============================================================================
//  TAXONOMY PLOT — top10 TSVs → multi-panel broken-axis plot (one column per group)
// ============================================================================

process PLOT_TAXONOMY {
    label 'count_only'
    publishDir "${params.outdir}/plots", mode: 'copy'
    input:  path(top10_tsvs)
            path(sample_sheet)
    output: path("*.png")
    script:
    """
    python3 ${projectDir}/bin/plot_top10_taxa_global_colors.py \
        --top10_dir . \
        --outdir    . \
        --threshold ${params.tax_plot_threshold} \
        --max-normals ${params.tax_plot_max_normals} \
        --xbreak    ${params.tax_plot_xbreak} \
        --samplesheet ${sample_sheet}
    """
}

// ============================================================================
//  AGGREGATE READ-COUNT REPORT
// ============================================================================

process AGGREGATE_REPORT {
    label 'count_only'
    publishDir "${params.outdir}/reports", mode: 'copy'
    input:
    path(raw_stats)
    path(trimmed_stats)
    path(star_stats)
    path(star_logs)              // Log.final.out files for STAR mapped %
    path(host_stats)
    path(virus_stats)
    path(nomags_stats)
    path(collapse_stats)
    path(all_filter_stats_ch)   // filter stats: contains classified counts pre-filter
    path(wids_files)             // weighted_ids: for Total transRNAs after all filters
    path(sample_sheet)
    output:
    path("pipeline_read_counts_report.tsv")
    path("pipeline_read_counts_report.html")
    path("dataset_summary_report.tsv")
    path("virus_exclusion_report.tsv")
    script:
    """
    python3 ${projectDir}/bin/aggregate_report.py \
        --output          pipeline_read_counts_report.tsv \
        --summary-output  dataset_summary_report.tsv \
        --virus-output    virus_exclusion_report.tsv \
        --samplesheet     ${sample_sheet}

    python3 ${projectDir}/bin/report_to_html.py \
        pipeline_read_counts_report.tsv \
        pipeline_read_counts_report.html \
        ${sample_sheet}
    """
}

// ============================================================================
//  WORKFLOW
// ============================================================================

workflow {

    // ── Read + validate the sample sheet (nf-schema) ──────────────────────────
    ch_samplesheet = Channel.fromList(
            samplesheetToList(params.sample_sheet, "${projectDir}/assets/schema_input.json")
        )
        .map { meta, fastq_1, fastq_2 ->
            meta.group = resolveGroup(meta)
            tuple(meta, [fastq_1, fastq_2])
        }

    samples_full_ch = ch_samplesheet                                // (meta, [fastq_1, fastq_2])
    samples_ch      = samples_full_ch.map { meta, reads -> meta }   // meta only

        resolved_samplesheet_ch = Channel.value("sample,group")
        .concat(samples_ch.map { meta -> "${meta.id},${meta.group}" })
        .collectFile(name: "resolved_samplesheet.csv", newLine: true)
        .first()

    // ── Raw-read QC (FastQC, vendored nf-core module) ─────────────────────────
    // First thing done to the raw reads, before anything else. Independent
    // of skip_count_reports: this is QC, not a read-count accounting step.
    if (!params.skip_fastqc) {
        FASTQC(samples_full_ch)
    }

    // ── Read count reports at every pre-computed step ─────────────────────────
    // Skippable via params.skip_count_reports.
    // When skipped, AGGREGATE_REPORT is also disabled (no stats to aggregate).
    // If taxonomy is still needed, reads_posttrim_tab.tsv is loaded from disk.
    if (!params.skip_count_reports) {
        raw_ch     = COUNT_RAW(samples_full_ch)
        trimmed_ch = COUNT_TRIMMED(samples_ch)
        star_ch    = COUNT_STAR_UNMAPPED(samples_ch)
        host_ch    = COUNT_BBSPLIT_HOST(samples_ch)
        virus_ch   = COUNT_BBSPLIT_VIRUS(samples_ch)
        nomags_ch  = COUNT_BBSPLIT_MAGS(samples_ch)

        reads_table_ch = MAKE_READS_POSTTRIM_TAB(
            trimmed_ch.map { it[1] }.collect()
        )
    } else if (!params.skip_taxonomy) {
        // Load the pre-existing reads_posttrim_tab.tsv for the taxonomy parser
        reads_table_ch = Channel.fromPath(
            "${params.outdir}/reports/reads_posttrim_tab.tsv",
            checkIfExists: true
        )
    }

    // ── Kraken annotation + filtering ─────────────────────────────────────────
    kraken_ann_ch          = ANNOTATE_KRAKEN(samples_ch)
    kraken_filtered_raw_ch = FILTER_KRAKEN(kraken_ann_ch)
    kraken_filtered_ch     = kraken_filtered_raw_ch.map { meta, tsv, stats -> tuple(meta, tsv) }
    // raw shape: (meta, filtered_tsv, filter_stats)
    // downstream shape: (meta, filtered_tsv)

    // ── BLAST annotation + filtering (mate 1 and 2 in parallel) ──────────────
    blast_mates_ch = samples_ch.flatMap { meta ->
        [ tuple(meta, "1"), tuple(meta, "2") ]
    }
    blast_ann_ch      = ANNOTATE_BLAST(blast_mates_ch)
    blast_filtered_ch = FILTER_BLAST(blast_ann_ch)

    blast1_ch = blast_filtered_ch.filter { it[1] == "1" }
                    .map { meta, mate, tsv, stats -> tuple(meta, tsv) }
    blast2_ch = blast_filtered_ch.filter { it[1] == "2" }
                    .map { meta, mate, tsv, stats -> tuple(meta, tsv) }

    // Also collect all filter stats for AGGREGATE_REPORT (virus report).
    all_filter_stats_ch = kraken_filtered_raw_ch
        .map { meta, tsv, stats -> stats }
        .mix(blast_filtered_ch.map { meta, mate, tsv, stats -> stats })
        .collect()

    // ── Collapse (fetch from trimmed reads, cat, seqkit rmdup) ────────────────
    // Skippable via params.skip_collapse — loads pre-existing outputs from
    // outdir/collapsed/ so that SELECT_CANDIDATES and downstream can still run.
    // Annotation + filtering above still runs when skip_collapse is true because
    // those outputs are needed by the taxonomy parser and aggregate report.
    if (!params.skip_collapse) {
        collapse_input_ch = blast1_ch
            .join(blast2_ch,          by: 0)
            .join(kraken_filtered_ch, by: 0)
        collapse_ch = COLLAPSE_READS(collapse_input_ch)
    } else {
        collapse_ch = samples_ch.map { meta ->
            tuple(
                meta,
                file("${params.outdir}/collapsed/${meta.id}_merged_duplicated.detail.txt",       checkIfExists: true),
                file("${params.outdir}/collapsed/${meta.id}_merged_collapsed_clean.fq",           checkIfExists: true),
                file("${params.outdir}/collapsed/${meta.id}_merged.fq",                           checkIfExists: true),
                file("${params.outdir}/collapsed/${meta.id}_collapse_stats.tsv",                  checkIfExists: true)
            )
        }
    }
    // shape: (meta, dup_detail, collapsed_clean_fq, merged_fq, collapse_stats)

    // ── Candidate selection ───────────────────────────────────────────────────
    // No longer needs BLAST/Kraken dirs — pool is already filtered.
    candidates_input_ch = collapse_ch
        .map { meta, dup_detail, collapsed_clean_fq, merged_fq, collapse_stats ->
               tuple(meta, dup_detail, collapsed_clean_fq, collapse_stats) }

    candidates_ch = SELECT_CANDIDATES(candidates_input_ch)
    // shape: (meta, ids, hist, wids)

    // ── Extract FASTA ─────────────────────────────────────────────────────────
    extract_input_ch = candidates_ch
        .join(collapse_ch.map { meta, dup_detail, collapsed_clean_fq, merged_fq, collapse_stats ->
            tuple(meta, collapsed_clean_fq) }, by: 0)
        .map { meta, ids, hist, wids, collapsed_clean_fq ->
            tuple(meta, ids, hist, wids, collapsed_clean_fq) }

    if (!params.skip_extract_fasta) {
        EXTRACT_FASTA(extract_input_ch)
    }

    // ── Length histograms + dataset plot ─────────────────────────────────────
    PUBLISH_HISTS(candidates_ch)

    if (!params.skip_length_plot) {
        all_hists_ch = PUBLISH_HISTS.out.collect()
        PLOT_LENGTH_DIST(all_hists_ch, resolved_samplesheet_ch)
    }

    // ── Taxonomy parser (once globally) ──────────────────────────────────────
    all_ids_ch = candidates_ch
        .map    { meta, ids, hist, wids -> ids }
        .collect()

    // Collect all kraken filtered files
    all_kraken_ch = kraken_filtered_ch
        .map    { meta, tsv -> tsv }
        .collect()

    // Collect all blast filtered files (both mates)
    all_blast_ch = blast_filtered_ch
        .map    { meta, mate, tsv, stats -> tsv }
        .collect()

    // Collect all dup detail files (one per sample in new flow)
    all_dups_ch = collapse_ch
        .map    { meta, dup_detail, collapsed_clean_fq, merged_fq, collapse_stats -> dup_detail }
        .collect()

    if (!params.skip_taxonomy) {
        taxonomy_ch = TAXONOMY_PARSER(
            all_ids_ch,
            all_kraken_ch,
            all_blast_ch,
            all_dups_ch,
            reads_table_ch,
            resolved_samplesheet_ch
        )

        // Pick out per-group combined top10 files (e.g. "RJ_blast_Domain_top10.tsv")
        // from per-sample ones (e.g. "RJ1_blast_Domain_top10.tsv"), for any
        // number/labels of groups. Filenames always end in a fixed
        // "_{tool}_{rank}_top10.tsv" suffix from a known small vocabulary, so
        // the label is extracted exactly (not by prefix) and checked against
        // the resolved group set — this also correctly excludes a sample
        // whose own ID happens to start with a real group's name followed by
        // "_" (e.g. sample "RJ_1" in group "RJ"), which a prefix check would
        // wrongly sweep in as if it were RJ's group-level file.
        // (A group with only one sample has no combined file at all — the
        // taxonomy parser skips it since the per-sample file already covers
        // it; see 04_05_2026_transRNA_taxonomy_parser.py.)
        group_labels_ch = samples_ch
            .map    { meta -> meta.group.toString() }
            .unique()
            .collect()

        top10_label_re = ~/^(.+)_(blast|kraken|combined)_(Domain|Kingdom|Order|Species)_top10\.tsv$/

        top10_for_plot_ch = taxonomy_ch.top10_tsvs
            .flatten()
            .combine(group_labels_ch)
            .filter  { f, groups ->
                def m = f.name =~ top10_label_re
                m.matches() && groups.contains(m.group(1))
            }
            .map     { f, groups -> f }
            .collect()

        PLOT_TAXONOMY(top10_for_plot_ch, resolved_samplesheet_ch)
    }

    // ── Aggregate read count report ───────────────────────────────────────────
    // Skipped automatically when skip_count_reports is true (no stats to read).
    // Can also be skipped independently via skip_report.
    if (!params.skip_report && !params.skip_count_reports) {
        all_wids_ch = candidates_ch
            .map { meta, ids, hist, wids -> wids }
            .collect()

        AGGREGATE_REPORT(
            raw_ch    .map { it[1] }.collect(),
            trimmed_ch.map { it[1] }.collect(),
            star_ch   .map { it[1] }.collect(),
            star_ch   .map { it[2] }.collect(),
            host_ch   .map { it[1] }.collect(),
            virus_ch  .map { it[1] }.collect(),
            nomags_ch .map { it[1] }.collect(),
            collapse_ch.map { it[4] }.collect(),
            all_filter_stats_ch,
            all_wids_ch,
            resolved_samplesheet_ch
        )
    }
}
