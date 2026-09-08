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
//  RAW-READ QC (FastQC)
//  Skippable via params.skip_fastqc.
// ============================================================================

process FASTQC {
    tag "${meta.id}"
    label 'count_only'
    container 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/reports/00_fastqc_raw", mode: 'copy'
    input:  tuple val(meta), path(reads)
    output: tuple val(meta), path("*_fastqc.html"), path("*_fastqc.zip")
    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}

// Aggregates one step's FastQC zips into one combined report. Called once
// per QC'd step (report_name picks the output subdir), so it stays a single
// process instead of being duplicated per step.
process MULTIQC {
    tag "${report_name}"
    label 'count_only'
    container 'https://depot.galaxyproject.org/singularity/multiqc:1.27--pyhdfd78af_0'
    publishDir "${params.outdir}/reports/${report_name}", mode: 'copy'
    input:  val(report_name)
            path(fastqc_zips)
    output: path("multiqc_report.html")
            path("multiqc_data")
    script:
    """
    multiqc .
    """
}

// ============================================================================
//  READ COUNT REPORTS  (seqkit stats -T at every pre-computed step)
// ============================================================================

process COUNT_RAW {
    tag "${meta.id}"
    label 'count_only'
    conda "${moduleDir}/config/envs/seqkit.yaml"
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

// Runs with TrimGalore's own defaults (quality/stringency/length shown
// explicitly below); trim_extra_args covers anything else (e.g. explicit
// adapters). --cores is deliberately not tied to task.cpus: TrimGalore's
// own docs note actual thread usage runs to ~4x --cores (see cambridge.config).
process TRIMGALORE {
    tag "${meta.id}"
    label 'med'
    conda "${moduleDir}/config/envs/trim_galore.yaml"
    publishDir "${params.outdir}/trimmed",            mode: 'copy', pattern: '*_trimmed.fq.gz'
    publishDir "${params.outdir}/reports/02_trimmed", mode: 'copy', pattern: '*.{txt,html,zip}'
    input:  tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.id}_1_trimmed.fq.gz"), path("${meta.id}_2_trimmed.fq.gz"), emit: reads
    path("*trimming_report.txt"), emit: reports
    path("*_fastqc.zip"),         emit: fastqc_zip
    path("*_fastqc.html"),        emit: fastqc_html
    script:
    """
    trim_galore --paired \
        --quality ${params.trim_quality} \
        --stringency ${params.trim_stringency} \
        --length ${params.trim_length} \
        --fastqc \
        --cores ${params.trim_cores} \
        ${params.trim_extra_args} \
        ${reads}

    mv *_val_1.f*q.gz ${meta.id}_1_trimmed.fq.gz
    mv *_val_2.f*q.gz ${meta.id}_2_trimmed.fq.gz
    """
}

process COUNT_TRIMMED {
    tag "${meta.id}"
    label 'count_only'
    conda "${moduleDir}/config/envs/seqkit.yaml"
    publishDir "${params.outdir}/reports/02_trimmed", mode: 'copy'
    input:  tuple val(meta), path(reads)
    output: tuple val(meta), path("${meta.id}_trimmed_stats.tsv")
    script:
    """
    seqkit stats -T -j ${task.cpus} \
        ${reads} \
        > ${meta.id}_trimmed_stats.tsv
    """
}

// Structural flags only (genomeDir, readFilesCommand, outSAMtype,
// outReadsUnmapped are needed for this pipeline's own data flow); alignment
// sensitivity stays at STAR's own defaults. star_extra_args covers anything
// else (e.g. project-specific filtering — see config/params/alvi_rnaome.yml).
process STAR {
    tag "${meta.id}"
    label 'high'
    conda "${moduleDir}/config/envs/star.yaml"
    publishDir "${params.outdir}/star", mode: 'copy', pattern: '*.bam'
    input:  tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.id}_Unmapped.out.mate1"), path("${meta.id}_Unmapped.out.mate2"), emit: unmapped
    tuple val(meta), path("${meta.id}_Log.final.out"),                                             emit: log_final
    path("${meta.id}_Aligned.sortedByCoord.out.bam"),                                               emit: bam
    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${params.star_genome_dir} \
        --readFilesIn ${reads} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outFileNamePrefix ${meta.id}_ \
        ${params.star_extra_args}
    """
}

process COUNT_STAR_UNMAPPED {
    tag "${meta.id}"
    label 'count_only'
    conda "${moduleDir}/config/envs/seqkit.yaml"
    publishDir "${params.outdir}/reports/03_star", mode: 'copy'
    input:  tuple val(meta), path(unmapped_reads), path(log_final, stageAs: 'input_log_final.out')
    output: tuple val(meta),
                  path("${meta.id}_star_unmapped_stats.tsv"),
                  path("${meta.id}_Log.final.out")
    script:
    """
    seqkit stats -T -j ${task.cpus} \
        ${unmapped_reads} \
        > ${meta.id}_star_unmapped_stats.tsv

    cp ${log_final} ${meta.id}_Log.final.out
    """
}

// ============================================================================
//  BBSPLIT  (pipeline step 4)
//  Decontaminates STAR's unmapped reads against host (bee) + human + viral
//  references.
//
//  Only structural flags (in1/in2, ref, basename, outu1/outu2, refstats) are
//  hard-coded here -- alignment sensitivity is left at BBSplit's own
//  defaults, with project-specific tuning passed through bbsplit_extra_args
//  (see config/params/alvi_rnaome.yml).
//
//  Per-reference match % for the report comes straight from BBSPLIT's own
//  refstats.txt (see AGGREGATE_REPORT) -- no separate COUNT_* process needed.
//  COUNT_BBSPLIT_MAGS below still reads from the separate externally-computed
//  params.bbsplit_mags_dir directory -- not yet wired to this process's own
//  output. That rewiring is a deliberate follow-up step, not done here.
// ============================================================================

process BBSPLIT {
    tag "${meta.id}"
    label 'high'
    conda "${moduleDir}/config/envs/bbsplit.yaml"
    publishDir "${params.outdir}/bbsplit", mode: 'copy',
        pattern: params.bbsplit_keep_matched ? "*" : "*_unmatched_{1,2}.fq,*_refstats.txt"
    input:  tuple val(meta), path(unmapped_reads)   // STAR.out.unmapped: [mate1, mate2]
    output:
    tuple val(meta), path("${meta.id}_bbsplit_unmatched_1.fq"), path("${meta.id}_bbsplit_unmatched_2.fq"), emit: unmatched
    tuple val(meta), path("${meta.id}_bbsplit_refstats.txt"),                                              emit: refstats
    script:
    """
    bbsplit.sh in1=${unmapped_reads[0]} in2=${unmapped_reads[1]} \
        ref=${params.apis_mellifera_genome_fasta},${params.viral_genomes_fasta},${params.human_genome_fasta} \
        basename=${meta.id}_bbsplit_%.fq \
        outu1=${meta.id}_bbsplit_unmatched_1.fq outu2=${meta.id}_bbsplit_unmatched_2.fq \
        refstats=${meta.id}_bbsplit_refstats.txt \
        ${params.bbsplit_extra_args}
    """
}

process COUNT_BBSPLIT_MAGS {
    tag "${meta.id}"
    label 'count_only'
    conda "${moduleDir}/config/envs/seqkit.yaml"
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
    conda "${moduleDir}/config/envs/python_analysis.yaml"
    publishDir "${params.outdir}/reports", mode: 'copy'
    input:  path(trimmed_stats_files)   // all *_trimmed_stats.tsv collected
    output: path("reads_posttrim_tab.tsv")
    script:
    """
    python3 ${moduleDir}/bin/make_reads_posttrim_tab.py
    """
}

// ============================================================================
//  COLLAPSE READS  —  two selectable variants, chosen by params.collapse_mode
//
//  Both produce the identical output shape (dup_detail, collapsed_clean_fq,
//  merged_fq, collapse_stats), so everything downstream (SELECT_CANDIDATES
//  onward) is unaffected by which one ran.
//
//  "pre_taxonomy"  (default, COLLAPSE_READS_PRE_TAXONOMY) — dedup runs right
//    after BBSplit, before Kraken/BLAST ever see a read. The standard,
//    efficient order for any dataset processed through this pipeline from
//    here on: classify each unique sequence once, not once per duplicate.
//
//  "post_taxonomy" (COLLAPSE_READS_POST_TAXONOMY) — dedup runs AFTER
//    Kraken/BLAST, by extracting the IDs that already passed annotation +
//    filtering and fetching those specific reads back out of the trimmed
//    pool. This is the pipeline's original logic, kept unchanged — needed
//    for the SLT/RJ dataset specifically, whose Kraken/BLAST results were
//    already generated (long before this pipeline existed) on the full,
//    non-deduplicated read set. Set via config/params/alvi_rnaome.yml.
// ============================================================================

process COLLAPSE_READS_PRE_TAXONOMY {
    tag "${meta.id}"
    label 'med'
    conda "${moduleDir}/config/envs/collapse.yaml"
    publishDir "${params.outdir}/collapsed", mode: 'copy'
    input:
    tuple val(meta), path(bbsplit_unmatched_1), path(bbsplit_unmatched_2)
    output:
    tuple val(meta),
          path("${meta.id}_merged_duplicated.detail.txt"),
          path("${meta.id}_merged_collapsed_clean.fq"),
          path("${meta.id}_merged.fq"),
          path("${meta.id}_collapse_stats.tsv")
    script:
    """
    # ── Merge both mates of BBSplit's clean/unmatched output ──────────────────
    # Mate 1 and mate 2 of the same fragment share the same base read ID
    # (Illumina's "<id> 1:N:..." / "<id> 2:N:..." convention -- the mate
    # distinguisher lives in the description field, which every downstream
    # ID-based lookup here strips). Tag /1 and /2 onto each mate's IDs before
    # merging so they can never collide under one shared bare ID once
    # collapsed -- verified previously to cause ~99% ID collision between
    # mates when omitted.
    awk 'NR%4==1{sub(/^\\S+/, "&/1")}1' ${bbsplit_unmatched_1} > mate1_tagged.fq
    awk 'NR%4==1{sub(/^\\S+/, "&/2")}1' ${bbsplit_unmatched_2} > mate2_tagged.fq
    cat mate1_tagged.fq mate2_tagged.fq > ${meta.id}_merged.fq

    seqkit rmdup -s \
        -j ${task.cpus} \
        -o ${meta.id}_merged_collapsed_clean.fq \
        -d /dev/null \
        -D ${meta.id}_merged_duplicated.detail.txt \
        ${meta.id}_merged.fq

    # ── Collapse summary ──────────────────────────────────────────────────────
    MERGED_READS=\$(awk 'END{print NR/4}' ${meta.id}_merged.fq)
    python3 ${moduleDir}/bin/collapse_stats.py ${meta.id} ${params.min_occ} \${MERGED_READS}
    """
}

process COLLAPSE_READS_POST_TAXONOMY {
    tag "${meta.id}"
    label 'med'
    conda "${moduleDir}/config/envs/collapse.yaml"
    publishDir "${params.outdir}/collapsed", mode: 'copy'
    input:
    tuple val(meta),
          path(blast1_tsv),
          path(blast2_tsv),
          path(kraken_tsv),
          path(trimmed_reads)
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

    # ── Fetch sequences from trimmed reads, then tag /1 or /2 on the header ID ──
    # Same mate-collision issue as COLLAPSE_READS_PRE_TAXONOMY above: mate 1
    # and mate 2 share a bare base ID, so untagged IDs from blast1/kraken
    # (both mate-1-sourced) and blast2 (mate-2-sourced) can collide once
    # merged. Kraken hits are always mate-1-scoped in this pipeline (fetched
    # from trimmed_reads[0] only), so they're tagged /1 like blast1.
    seqkit grep -j ${task.cpus} -f blast1_ids.txt \
        ${trimmed_reads[0]} | awk 'NR%4==1{sub(/^\\S+/, "&/1")}1' > blast1_passing.fq
    seqkit grep -j ${task.cpus} -f blast2_ids.txt \
        ${trimmed_reads[1]} | awk 'NR%4==1{sub(/^\\S+/, "&/2")}1' > blast2_passing.fq
    seqkit grep -j ${task.cpus} -f kraken_ids.txt \
        ${trimmed_reads[0]} | awk 'NR%4==1{sub(/^\\S+/, "&/1")}1' > kraken_passing.fq

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
    python3 ${moduleDir}/bin/collapse_stats.py ${meta.id} ${params.min_occ} \${MERGED_READS}
    """
}

// ============================================================================
//  KRAKEN2  —  real process, only runs under collapse_mode=pre_taxonomy
//              and skip_kraken=false (see the workflow block).
//
//  Query is EXTRACT_FASTA's candidate FASTA — the whole point of the
//  pre_taxonomy dedup order is running the expensive classification tools
//  on the smallest possible (deduplicated, ≥min_occ) set, not on every raw
//  read. Structural flags only (--db, --output); tuning via
//  kraken2_extra_args (see config/params/alvi_rnaome.yml for the pattern).
// ============================================================================

process KRAKEN2 {
    tag "${meta.id}"
    label 'high'
    conda "${moduleDir}/config/envs/kraken2.yaml"
    publishDir "${params.outdir}/kraken_raw", mode: 'copy'
    input:  tuple val(meta), path(candidate_fasta)
    output: tuple val(meta), path("${meta.id}_kraken_output.txt")
    script:
    """
    kraken2 --db ${params.kraken2_db} \
        --output ${meta.id}_kraken_output.txt \
        --threads ${task.cpus} \
        ${params.kraken2_extra_args} \
        ${candidate_fasta}
    """
}

// ============================================================================
//  BLASTN  —  real process, only runs under collapse_mode=pre_taxonomy and
//             skip_blast=false (see the workflow block).
//
//  Query is FILTER_KRAKEN's unclassified-ID subset of the candidate FASTA
//  (seqkit grep'd out as this process's own first step) — the real
//  "U rows -> BLAST" handoff, wired for real instead of an unenforced
//  assumption about an external run. 14-column outfmt (no qlen/slen)
//  matches the historical convention ANNOTATE_BLAST's --taxid-col 13
//  already assumes.
// ============================================================================

process BLASTN {
    tag "${meta.id}"
    label 'high'
    conda "${moduleDir}/config/envs/blast.yaml"
    publishDir "${params.outdir}/blast_raw", mode: 'copy'
    input:
    tuple val(meta), path(candidate_fasta), path(unclassified_ids)
    output: tuple val(meta), path("${meta.id}_blast_top_hits.tsv")
    script:
    """
    seqkit grep -j ${task.cpus} -f ${unclassified_ids} \
        ${candidate_fasta} > ${meta.id}_unclassified.fasta

    blastn -query ${meta.id}_unclassified.fasta \
        -db ${params.blast_nt_db} \
        -task blastn-short \
        -out ${meta.id}_blast_top_hits.tsv \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames" \
        -num_threads ${task.cpus} \
        ${params.blastn_extra_args}
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
//    • U rows silently skipped from the kept output — but now also written to
//      a real unclassified-IDs output, which is BLASTN's actual query list
//      under collapse_mode=pre_taxonomy (previously just an assumption about
//      how an external BLAST run was set up)
//    • Removes: viruses (Domain), invertebrates (Phylum), host genera
//    • Removes: either mate length < min_len nt
// ============================================================================

process ANNOTATE_KRAKEN {
    tag "${meta.id}"
    label 'med'
    conda "${moduleDir}/config/envs/taxonkit.yaml"
    //publishDir "${params.outdir}/kraken_annotated", mode: 'copy'
    // Source-agnostic: kraken_raw is either a real KRAKEN2 process output
    // (collapse_mode=pre_taxonomy, skip_kraken=false) or a file() lookup
    // into params.kraken_dir (post_taxonomy, or pre_taxonomy+skip_kraken) —
    // see the workflow block for which.
    input:  tuple val(meta), path(kraken_raw)
    output: tuple val(meta), path("${meta.id}_kraken_annotated.tsv")
    script:
    """
    python3 ${moduleDir}/bin/annotate_kraken_lineage.py \
        --input       ${kraken_raw} \
        --output      ${meta.id}_kraken_annotated.tsv \
    """
}

process FILTER_KRAKEN {
    tag "${meta.id}"
    label 'low'
    conda "${moduleDir}/config/envs/python_core.yaml"
    publishDir "${params.outdir}/kraken_filtered", mode: 'copy'
    input:
    tuple val(meta), path(annotated_tsv)
    output:
    tuple val(meta),
          path("${meta.id}_kraken_invertebrates_filtered.tsv"),
          path("${meta.id}_kraken_filter_stats.tsv"),
          path("${meta.id}_kraken_unclassified_ids.txt")        // ← new: real
          // BLAST query list under collapse_mode=pre_taxonomy (unused, but
          // still produced, under post_taxonomy — harmless small file).
    script:
    """
    python3 ${moduleDir}/bin/filter_kraken_invertebrates.py \
        --min-len          ${params.min_len} \
        --stats-out        ${meta.id}_kraken_filter_stats.tsv \
        --unclassified-out ${meta.id}_kraken_unclassified_ids.txt \
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
    conda "${moduleDir}/config/envs/taxonkit.yaml"
    //publishDir "${params.outdir}/blast_annotated", mode: 'copy'
    // Source-agnostic, same pattern as ANNOTATE_KRAKEN: blast_raw is either a
    // real BLASTN process output (pre_taxonomy, skip_blast=false — "mate" is
    // just "candidates" there, no real mate split) or a file() lookup into
    // params.blast_dir (post_taxonomy, or pre_taxonomy+skip_blast).
    input:  tuple val(meta), val(mate), path(blast_raw)
    output: tuple val(meta), val(mate), path("${meta.id}_${mate}_blast_annotated.tsv")
    script:
    """
    python3 ${moduleDir}/bin/annotate_blast_lineage.py \
        --input       ${blast_raw} \
        --output      ${meta.id}_${mate}_blast_annotated.tsv \
        --taxid-col   13
    """
}

process FILTER_BLAST {
    tag "${meta.id}_${mate}"
    label 'low'
    conda "${moduleDir}/config/envs/python_core.yaml"
    publishDir "${params.outdir}/blast_filtered", mode: 'copy'
    input:
    tuple val(meta), val(mate), path(annotated_tsv)
    output:
    tuple val(meta), val(mate),
          path("${meta.id}_${mate}_blast_all_lengths_filtered.tsv"),
          path("${meta.id}_${mate}_blast_filter_stats.tsv")      // ← new
    script:
    """
    python3 ${moduleDir}/bin/filter_blast_all_conditions.py \
        --min-len   ${params.min_len} \
        --stats-out ${meta.id}_${mate}_blast_filter_stats.tsv \
        ${params.invertebrate_phyla} \
        ${annotated_tsv} \
        ${meta.id}_${mate}_blast_all_lengths_filtered.tsv
    """
}

// ============================================================================
//  CANDIDATE SELECTION
//  Never cross-references BLAST/Kraken TSVs itself — just applies the
//  ≥ min_occ duplicate-count threshold and reads lengths from the collapsed
//  FASTQ. Why that's sufficient depends on collapse_mode (see the COLLAPSE
//  READS comment block above): under "post_taxonomy", the collapsed pool is
//  already BLAST/Kraken-validated, so ≥ min_occ duplicates alone makes a
//  candidate; under "pre_taxonomy", candidates are selected before taxonomy
//  even runs, and Kraken/BLAST classify only this already-deduplicated set.
// ============================================================================

process SELECT_CANDIDATES {
    tag "${meta.id}"
    label 'low'
    conda "${moduleDir}/config/envs/python_core.yaml"
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
    python3 ${moduleDir}/bin/ge5_18nt_filtering_pipeline_for_both_mates.py \
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
    conda "${moduleDir}/config/envs/seqkit.yaml"
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
    conda "${moduleDir}/config/envs/python_analysis.yaml"
    publishDir "${params.outdir}/plots", mode: 'copy'
    input:  path(hist_files)
            path(sample_sheet)
    output: path("*.png")          // ← collect all PNGs (both plots)
    script:
    """
    python3 ${moduleDir}/bin/plot_lengths.py . --samplesheet ${sample_sheet}
    """
}

// ============================================================================
//  TAXONOMY PARSER — called ONCE GLOBALLY
//
//  Inputs collected from all samples:
//    --ids-dir    : folder with all *_ge5_detected_ids.txt files
//    --blast-dir  : folder with all *_blast_all_lengths_filtered.tsv files
//    --kraken-dir : folder with all *_kraken_invertebrates_filtered.tsv files
//    --dup-dir    : folder with all *_merged_duplicated.detail.txt files
//                   (one per sample, from COLLAPSE_READS_{PRE,POST}_TAXONOMY;
//                   IDs inside are /1 /2 mate-suffixed)
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
    conda "${moduleDir}/config/envs/python_analysis.yaml"
    publishDir "${params.outdir}/final_transRNAs_taxonomies", mode: 'copy'
    input:
    // ids_files: all *_ge5_detected_ids.txt from all samples
    path(ids_files)
    // kraken_files: all *_kraken_invertebrates_filtered.tsv
    path(kraken_files)
    // blast_files: all *_blast_all_lengths_filtered.tsv
    path(blast_files)
    // dup_files: all *_merged_duplicated.detail.txt
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
    python3 ${moduleDir}/bin/04_05_2026_transRNA_taxonomy_parser.py \
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
    conda "${moduleDir}/config/envs/python_analysis.yaml"
    publishDir "${params.outdir}/plots", mode: 'copy'
    input:  path(top10_tsvs)
            path(sample_sheet)
    output: path("*.png")
    script:
    """
    python3 ${moduleDir}/bin/plot_top10_taxa_global_colors.py \
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
    conda "${moduleDir}/config/envs/python_analysis.yaml"
    publishDir "${params.outdir}/reports", mode: 'copy'
    input:
    path(raw_stats)
    path(trimmed_stats)
    path(star_stats)
    path(star_logs)              // Log.final.out files for STAR mapped %
    path(bbsplit_refstats)       // raw BBSPLIT refstats.txt, one per sample
    path(nomags_stats)
    path(collapse_stats)
    path(all_filter_stats_ch)   // filter stats: contains classified counts pre-filter
    path(wids_files)             // weighted_ids: for Total transRNAs after all filters
    path(sample_sheet)
    output:
    path("pipeline_read_counts_report.tsv"),   emit: tsv
    path("pipeline_read_counts_report.html"),  emit: html
    path("dataset_summary_report.tsv"),        emit: dataset_summary
    path("virus_exclusion_report.tsv"),        emit: virus_exclusion
    script:
    """
    python3 ${moduleDir}/bin/aggregate_report.py \
        --output          pipeline_read_counts_report.tsv \
        --summary-output  dataset_summary_report.tsv \
        --virus-output    virus_exclusion_report.tsv \
        --samplesheet     ${sample_sheet}

    python3 ${moduleDir}/bin/report_to_html.py \
        pipeline_read_counts_report.tsv \
        pipeline_read_counts_report.html \
        ${sample_sheet}
    """
}

// ============================================================================
//  PLOT REPORT SUMMARY — two-panel cascade + classification-metrics heatmap,
//  built from AGGREGATE_REPORT's own TSV. Group colouring/legend comes from
//  the real sample sheet (general-purpose across any dataset).
// ============================================================================

process PLOT_REPORT_SUMMARY {
    label 'count_only'
    conda "${moduleDir}/config/envs/python_analysis.yaml"
    publishDir "${params.outdir}/reports", mode: 'copy'
    input:
    path(report_tsv)
    path(sample_sheet)
    output:
    path("*_pipeline_report_summary.png")
    script:
    """
    python3 ${moduleDir}/bin/plot_report_summary.py \
        ${report_tsv} \
        --samplesheet ${sample_sheet}
    """
}

// ============================================================================
//  WORKFLOW
// ============================================================================

workflow {

    // ── Read + validate the sample sheet (nf-schema) ──────────────────────────
    ch_samplesheet = Channel.fromList(
            samplesheetToList(params.sample_sheet, "${moduleDir}/assets/schema_input.json")
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

    // ── Raw-read QC (FastQC) ───────────────────────────────────────────────────
    // First thing done to the raw reads, before anything else. Independent
    // of skip_count_reports: this is QC, not a read-count accounting step.
    if (!params.skip_fastqc) {
        FASTQC(samples_full_ch)
        MULTIQC('00_fastqc_raw', FASTQC.out.map { meta, html, zip -> zip }.flatten().collect())
    }

    // ── Adapter/quality trimming (TrimGalore) ─────────────────────────────────
    // Skippable via skip_trimming — falls back to pre-trimmed reads already
    // in params.trimmed_dir, matching {sample}_{1,2}_trimmed.fq.gz.
    if (!params.skip_trimming) {
        TRIMGALORE(samples_full_ch)
        trimmed_reads_ch = TRIMGALORE.out.reads.map { meta, r1, r2 -> tuple(meta, [r1, r2]) }

        if (!params.skip_fastqc) {
            MULTIQC('02_trimmed', TRIMGALORE.out.fastqc_zip.collect())
        }
    } else {
        trimmed_reads_ch = samples_ch.map { meta ->
            tuple(meta, [
                file("${params.trimmed_dir}/${meta.id}_1_trimmed.fq.gz", checkIfExists: true),
                file("${params.trimmed_dir}/${meta.id}_2_trimmed.fq.gz", checkIfExists: true)
            ])
        }
    }

    // ── Host-genome alignment (STAR) ──────────────────────────────────────────
    // Skippable via skip_star — falls back to pre-computed unmapped reads
    // already in params.star_dir, matching {sample}/{sample}_Unmapped.out.mate{1,2}.
    if (!params.skip_star) {
        STAR(trimmed_reads_ch)
        star_out_ch = STAR.out.unmapped
            .join(STAR.out.log_final, by: 0)
            .map { meta, r1, r2, log -> tuple(meta, [r1, r2], log) }

        if (!params.skip_fastqc) {
            MULTIQC('03_star', STAR.out.log_final.map { meta, log -> log }.collect())
        }
    } else {
        star_out_ch = samples_ch.map { meta ->
            tuple(meta,
                [
                    file("${params.star_dir}/${meta.id}/${meta.id}_Unmapped.out.mate1", checkIfExists: true),
                    file("${params.star_dir}/${meta.id}/${meta.id}_Unmapped.out.mate2", checkIfExists: true)
                ],
                file("${params.star_dir}/${meta.id}/Log.final.out", checkIfExists: true)
            )
        }
    }

    // ── Decontamination (BBSplit) — pipeline step 4 ───────────────────────────
    // Skippable via skip_bbsplit — falls back to pre-computed unmatched reads
    // already in params.bbsplit_dir, matching
    // {sample}/{sample}_bbsplit_unmatched_{1,2}.fq.
    // Consumed downstream by COLLAPSE_READS_PRE_TAXONOMY (see the COLLAPSE
    // READS comment block below) when params.collapse_mode == "pre_taxonomy".
    if (!params.skip_bbsplit) {
        bbsplit_in_ch = star_out_ch.map { meta, reads, log -> tuple(meta, reads) }
        BBSPLIT(bbsplit_in_ch)
        bbsplit_out_ch      = BBSPLIT.out.unmatched
        bbsplit_refstats_ch = BBSPLIT.out.refstats.map { it[1] }.collect()
    } else {
        bbsplit_out_ch = samples_ch.map { meta ->
            tuple(meta,
                file("${params.bbsplit_dir}/${meta.id}/${meta.id}_bbsplit_unmatched_1.fq", checkIfExists: true),
                file("${params.bbsplit_dir}/${meta.id}/${meta.id}_bbsplit_unmatched_2.fq", checkIfExists: true)
            )
        }
        bbsplit_refstats_ch = samples_ch.map { meta ->
            file("${params.bbsplit_dir}/${meta.id}/${meta.id}_bbsplit_refstats.txt", checkIfExists: true)
        }.collect()
    }

    // ── Read count reports at every pre-computed step ─────────────────────────
    // Skippable via params.skip_count_reports.
    // When skipped, AGGREGATE_REPORT is also disabled (no stats to aggregate).
    // If taxonomy is still needed, reads_posttrim_tab.tsv is loaded from disk.
    if (!params.skip_count_reports) {
        raw_ch     = COUNT_RAW(samples_full_ch)
        trimmed_ch = COUNT_TRIMMED(trimmed_reads_ch)
        star_ch    = COUNT_STAR_UNMAPPED(star_out_ch)
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

    // ── Kraken + BLAST, post_taxonomy mode only ───────────────────────────────
    // Legacy datasets (collapse_mode=post_taxonomy) already have real
    // Kraken2/BLAST results, computed long before this pipeline existed, on
    // the full non-deduplicated pool — read from params.kraken_dir/blast_dir
    // unconditionally here, same as always. COLLAPSE_READS_POST_TAXONOMY
    // (below) needs kraken_filtered_ch/blast1_ch/blast2_ch as real inputs, so
    // this has to happen before collapse in this mode — structurally the
    // opposite of pre_taxonomy mode, where KRAKEN2/BLASTN run for real, but
    // only after EXTRACT_FASTA produces something to classify (see the
    // matching block right after EXTRACT_FASTA, below).
    if (params.collapse_mode == "post_taxonomy") {
        kraken_raw_ch = samples_ch.map { meta ->
            tuple(meta, file("${params.kraken_dir}/${meta.id}_2MM_CLEAN_noMAGs_kraken_output_005_.txt", checkIfExists: true))
        }
        kraken_ann_ch          = ANNOTATE_KRAKEN(kraken_raw_ch)
        kraken_filtered_raw_ch = FILTER_KRAKEN(kraken_ann_ch)
        kraken_filtered_ch     = kraken_filtered_raw_ch.map { meta, tsv, stats, unclassified -> tuple(meta, tsv) }
        // raw shape: (meta, filtered_tsv, filter_stats, unclassified_ids)
        // downstream shape: (meta, filtered_tsv)

        blast_mates_ch = samples_ch.flatMap { meta ->
            [ tuple(meta, "1"), tuple(meta, "2") ]
        }.map { meta, mate ->
            tuple(meta, mate, file("${params.blast_dir}/${meta.id}_${mate}.tsv", checkIfExists: true))
        }
        blast_ann_ch      = ANNOTATE_BLAST(blast_mates_ch)
        blast_filtered_ch = FILTER_BLAST(blast_ann_ch)

        blast1_ch = blast_filtered_ch.filter { it[1] == "1" }
                        .map { meta, mate, tsv, stats -> tuple(meta, tsv) }
        blast2_ch = blast_filtered_ch.filter { it[1] == "2" }
                        .map { meta, mate, tsv, stats -> tuple(meta, tsv) }

        // Also collect all filter stats for AGGREGATE_REPORT (virus report).
        all_filter_stats_ch = kraken_filtered_raw_ch
            .map { meta, tsv, stats, unclassified -> stats }
            .mix(blast_filtered_ch.map { meta, mate, tsv, stats -> stats })
            .collect()
    }

    // ── Collapse (dedup) ───────────────────────────────────────────────────────
    // Skippable via params.skip_collapse — loads pre-existing outputs from
    // outdir/collapsed/ so that SELECT_CANDIDATES and downstream can still run.
    // Annotation + filtering above still runs when skip_collapse is true because
    // those outputs are needed by the taxonomy parser and aggregate report.
    //
    // Which collapse process actually runs is chosen by params.collapse_mode
    // (see the COLLAPSE READS comment block above) — both branches converge
    // to the same collapse_ch shape, so nothing past this point cares which
    // one ran.
    if (!params.skip_collapse) {
        if (params.collapse_mode == "pre_taxonomy") {
            collapse_ch = COLLAPSE_READS_PRE_TAXONOMY(bbsplit_out_ch)
        } else {
            collapse_input_ch = blast1_ch
                .join(blast2_ch,          by: 0)
                .join(kraken_filtered_ch, by: 0)
                .join(trimmed_reads_ch,   by: 0)
            collapse_ch = COLLAPSE_READS_POST_TAXONOMY(collapse_input_ch)
        }
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

    // fasta_ch: (meta, candidate_fasta) — only actually needed downstream
    // under collapse_mode=pre_taxonomy (KRAKEN2's query, below); still built
    // either way since EXTRACT_FASTA's own skip-fallback needs the same shape.
    if (!params.skip_extract_fasta) {
        fasta_ch = EXTRACT_FASTA(extract_input_ch).map { meta, fasta, report -> tuple(meta, fasta) }
    } else {
        fasta_ch = samples_ch.map { meta ->
            tuple(meta, file("${params.outdir}/final_transRNAs_fasta_collapsed/${meta.id}_final_transRNAs.fasta", checkIfExists: true))
        }
    }

    // ── Kraken + BLAST, pre_taxonomy mode only ────────────────────────────────
    // Mirrors the post_taxonomy block far above, but the query is
    // EXTRACT_FASTA's candidate set instead of the pre-computed external
    // dirs — the actual point of the pre_taxonomy dedup order: classify the
    // smallest possible (deduplicated, >=min_occ) set once, not every raw
    // read. skip_kraken/skip_blast are independent — each still falls back
    // to a pre-computed directory if set, same fallback shape as every other
    // skip_* flag in this pipeline.
    if (params.collapse_mode == "pre_taxonomy") {
        if (!params.skip_kraken) {
            kraken_raw_ch = KRAKEN2(fasta_ch)
        } else {
            kraken_raw_ch = samples_ch.map { meta ->
                tuple(meta, file("${params.kraken_dir}/${meta.id}_2MM_CLEAN_noMAGs_kraken_output_005_.txt", checkIfExists: true))
            }
        }
        kraken_ann_ch          = ANNOTATE_KRAKEN(kraken_raw_ch)
        kraken_filtered_raw_ch = FILTER_KRAKEN(kraken_ann_ch)
        kraken_filtered_ch     = kraken_filtered_raw_ch.map { meta, tsv, stats, unclassified -> tuple(meta, tsv) }
        kraken_unclassified_ch = kraken_filtered_raw_ch.map { meta, tsv, stats, unclassified -> tuple(meta, unclassified) }

        if (!params.skip_blast) {
            // BLASTN's own query = the unclassified subset of fasta_ch,
            // extracted inside the process itself (see BLASTN's script).
            blastn_in_ch = fasta_ch.join(kraken_unclassified_ch, by: 0)
            blast_raw_ch = BLASTN(blastn_in_ch).map { meta, tsv -> tuple(meta, "candidates", tsv) }
        } else {
            blast_raw_ch = samples_ch.map { meta ->
                tuple(meta, "candidates", file("${params.blast_dir}/${meta.id}_candidates.tsv", checkIfExists: true))
            }
        }
        blast_ann_ch      = ANNOTATE_BLAST(blast_raw_ch)
        blast_filtered_ch = FILTER_BLAST(blast_ann_ch)
        // No real mate-1/mate-2 split under pre_taxonomy — "candidates" is a
        // single pool tag, not a real mate. blast1_ch/blast2_ch aren't
        // needed here: they only ever feed COLLAPSE_READS_POST_TAXONOMY,
        // which doesn't run in this mode.

        all_filter_stats_ch = kraken_filtered_raw_ch
            .map { meta, tsv, stats, unclassified -> stats }
            .mix(blast_filtered_ch.map { meta, mate, tsv, stats -> stats })
            .collect()
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
            bbsplit_refstats_ch,
            nomags_ch .map { it[1] }.collect(),
            collapse_ch.map { it[4] }.collect(),
            all_filter_stats_ch,
            all_wids_ch,
            resolved_samplesheet_ch
        )

        PLOT_REPORT_SUMMARY(
            AGGREGATE_REPORT.out.tsv,
            resolved_samplesheet_ch
        )
    }
}
