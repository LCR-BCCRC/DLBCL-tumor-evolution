##### PYTHON MODULES #####

import oncopipe as op
import os
import pandas as pd

##### SETUP VARIABLES #####

# Raw data
RAW_BIOPSY      = "../private_data/clinical_md_biopsy_FINAL.xlsx"
RAW_FISH        = "../private_data/Assay_FISH_consensus.tsv"
RAW_NANOSTRING  = "../private_data/Assay_Nanostring.xlsx"
RAW_MTP         = "../private_data/DLBCL_Trio_DS_Mastertable.xlsx"


# Metadata tables
MD_SEQUENCING   = "data/metadata/metadata_sequencing.tsv"
MD_PATIENT      = "data/metadata/metadata_clinical_patient.tsv"
MD_BIOPSY       = "data/metadata/metadata_clinical_biopsy.tsv"
MD_FISH_NS      = "data/metadata/metadata_fish_nanostring.tsv"
MD_QC           = "data/qc/gambl_qc_combined.tsv"
MD_SURV         = "../private_data/DLBCL_ITT_LH_Aug2022.xlsx"

# GAMBL
GAMBL_BASE      = "/projects/rmorin/projects/gambl-repos/gambl-crushton-canary/"
GAMBL_LYMPHGEN  = GAMBL_BASE + "versioned_results/LymphGen/GAMBL_all_the_things.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv"
GAMBL_METADATA  = GAMBL_BASE + "data/metadata/gambl_samples_available.tsv"
GAMBL_QC        = GAMBL_BASE + "results/gambl/qc-1.0/99-outputs/{seq_type}.qc_metrics.tsv"

# LySeq
LYSEQ_BASE      = "../LySeq_Validation/"
LYSEQ_LYMPHGEN  = LYSEQ_BASE + "results/lymphgen-1.0/99-outputs/LyseqST_genome_callable.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv"
LYSEQ_QC        = LYSEQ_BASE + "results/QC/04-mergedMetrics/capture--grch37.hs_metrics.txt"
LYSEQ_BED       = LYSEQ_BASE + "data/bed/lyseqst_targets.grch37.bed"

# MAF data
GENOME_MAF      = "data/maf/genomes_augmented.grch37.maf"
MAF_BASE        = "data/maf/"
GENOME_CALL_MAF = MAF_BASE + "genomes_callable.maf"
LYSEQ_MAF       = MAF_BASE + "LySeqST_only_read_threshold.maf"
COMBINED_MAF    = MAF_BASE + "LySeqST_genomes_merged.maf"
COMB_CALL_MAF   = MAF_BASE + "LySeqST_genomes_merged_callable.maf"
PHYCLONE_MAF    = "data/phyclone/phyclone_mafs.tsv"


##### CREATE METADATA #####

rule metadata_sequencing:
    output:
        metadata = MD_SEQUENCING
    script:
        "src/create_sequencing_metadata.R"

rule metadata_clinical:
    input:
        md = MD_SEQUENCING,
        raw = RAW_BIOPSY
    output:
        patient = MD_PATIENT,
        biopsy = MD_BIOPSY
    script:
        "src/finalize_clinical_md.R"

rule metadata_fish_nanostring:
    input:
        RAW_MTP,
        RAW_NANOSTRING,
        RAW_FISH
    output:
        MD_FISH_NS
    script:
        "src/prepare_fish_nanostring_data.R"

##### PLOTS #####

# Swimmer plot

rule swimmer_plot:
    input:
        MD_SEQUENCING,
        MD_BIOPSY,
        MD_PATIENT
    output:
        "figures/swimmer_plot.pdf",
        "figures/swimmer_plot_alt_layout.pdf",
        "figures/swimmer_plot_breaks.pdf",
        "figures/swimmer_plot_breaks_alt_layout.pdf"
    script:
        "src/swimmer_plot.R"
        

# FISH and NanoString plots
rule fish_ns_plots:
    input:
        MD_FISH_NS
    output:
        "figures/coo_alluvial_score.pdf",
        "figures/dhitsig_alluvial_score.pdf",
        "figures/fish_pairs.pdf",
        "data/fish_nanostring/fish_summary.tsv",
        "data/fish_nanostring/ns_summary.tsv", 
        "data/fish_nanostring/multi_FISH_summary.tsv"
    script:
        "src/analyze_fish_ns.R"

# Assemble a merged maf of genome variant calls with min read support of >= 3 reads
rule genome_maf:
    input:
        MD_SEQUENCING
    output:
        GENOME_MAF
    script:
        "src/genome_maf.R"

# Plot QC metrics
rule qc_data:
    input:
        MD_SEQUENCING, 
        GENOME_MAF
    output:
        MD_QC,  
        "figures/qc_boxplots.pdf" 
    script: 
        "src/get_qc_data.R"
        
# Assemble SV results
rule wgs_svs:
    input:
        MD_SEQUENCING
    output:
        "data/svs/summarize_all_svs.tsv",
        "data/svs/summarize_multi_FISH_pos.tsv",
        "data/svs/summarize_multi_FISH_pos_relapse_timing.tsv"
    script:
        "src/get_sv_results.R"

# Assemble MiXCR IG rearrangement calls
rule assemble_mixcr:
    input:
        MD_SEQUENCING
    output:
        "data/ig_rearrangement/lightchain_alluvial.tsv",
        "data/ig_rearrangement/heavychain_alluvial.tsv",
        "data/ig_rearrangement/all_mixcr_results.tsv"
    script:
        "src/assemble_mixcr_results.R"

rule plot_mixcr:
    input:
        rules.assemble_mixcr.output
    output:
        "data/ig_rearrangement/summary_ig_concordance.tsv",
        "figures/mixcr_lightchain_alluvial.pdf",
        "figures/mixcr_heavychain_alluvial.pdf"
    script:
        "src/plot_mixcr.R"

# Compare LymphGen calls 
rule compare_lymphgen:
    input:
        MD_SEQUENCING,
        GAMBL_LYMPHGEN,
        LYSEQ_LYMPHGEN
    output:
        "data/lymphgen/lymphgen_results.tsv",
        "data/lymphgen/lymphgen_summary.tsv",
        "figures/lymphgen_comparison.pdf"
    script:
        "src/compare_lymphgen.R"

# Survival Analyses
rule survival_plots: 
    input: 
       MD_SURV 
    output: 
        "data/survival/BMT_summary.tsv", 
        "data/survival/response_summary.tsv", 
        "data/survival/surv_wide.tsv", 
        "data/survival/cox_models.tsv", 
        "figures/surv_forestplot.pdf", 
        "figures/bmt_barplot.pdf", 
        "figures/response_barplot.pdf", 
        "figures/survival_HSCT_faceted.pdf", 
        "figures/survival_prog_faceted.pdf", 
        "figures/survival_faceted.pdf", 
    script: 
        "src/survival_analyses.R"

# Shared mutation analysis

rule shared_mutations: 
    input: 
        MD_SEQUENCING,
        MD_PATIENT,
        MD_QC, 
        GENOME_MAF
    output: 
        "data/shared_mutations/unique_mutations_per_patient.tsv",
        "data/shared_mutations/mutation_tally_by_group.tsv", 
        "data/shared_mutations/lm_mutations_vs_relapse_coverage_categorical.tsv", 
        "data/shared_mutations/lm_mutations_vs_relapse_coverage_continuous.tsv", 
        "data/shared_mutations/mutation_tally_by_group_wilcox.tsv", 
        "data/shared_mutations/mutation_percent_by_group_wilcox.tsv", 
        "figures/percent_unique_boxplot.pdf", 
        "figures/count_unique_boxplot.pdf", 
        "figures/percent_unique_boxplot_wide.pdf", 
        "figures/count_unique_linear.pdf", 
        "figures/count_unique_linear_coverage.pdf", 
        "figures/total_mutations_vs_age.pdf", 
        "figures/unique_mutations_vs_age.pdf"
    script: 
        "src/analyze_shared_mutations.R"

# PhyClone        
rule pyclone_maf: 
    input: 
        MD_SEQUENCING, 
        LYSEQ_LYMPHGEN
    output: 
        "data/phyclone/phyclone_mafs.tsv"
    script: 
        "src/all_phyclone_plots.R"
        
# LySeq Validation

rule lyseq_qc: 
    input: 
        MD_SEQUENCING, 
        LYSEQ_QC
    output: 
        "data/qc/lst_qc_coverage.tsv", 
        "figures/lst_qc.pdf"
    script: 
        "src/LySeq_qc.R"

rule lyseq_validation: 
    input: 
        MD_SEQUENCING, 
        MD_QC, 
        LYSEQ_QC, 
        GENOME_MAF, 
        LYSEQ_BED
    output: 
        "data/lyseq/genome_sensitivity_lyseq_lm.tsv", 
        "figures/lyseq_vaf_violin_sensitivity.pdf", 
        "figures/genome_sensitivity_lyseq.pdf"
    script: 
        "src/lyseq_mutations2.R"
        
# Trunk Mutations
rule trunk_mutations: 
    input: 
        MD_SEQUENCING, 
        rules.pyclone_maf.output[0], 
        rules.shared_mutations.output[0]
    output: 
        "data/phyclone/mutations_divergent_patients.tsv", 
        "data/phyclone/tabulate_trunk.tsv", 
        "data/phyclone/tabulate_constrained.tsv", 
        "data/phyclone/tabulate_constrained_pt.tsv"
    script: 
        "src/trunk_mutations.R"
        
# Summarize everything
rule summarize_all: 
    input: 
        MD_SEQUENCING, 
        MD_FISH_NS, 
        MD_PATIENT
    output: 
        "data/metadata/summarize_all_assays.tsv", 
        "data/metadata/all_assays_results.tsv"
    script: 
        "src/summarize_all.R"
        
# Assay UpSet Plot
rule upset_plot: 
    input: 
        "data/metadata/all_assays_results.tsv"
    output: 
        "figures/assay_upset_plot.pdf"
    script: 
        "src/summarize_cohort_plots.R"

##### TARGETS #####
rule all:
    input:
        # Metadata
        MD_SEQUENCING,
        MD_PATIENT,
        MD_BIOPSY,
        MD_FISH_NS,
        # Genome maf file (unfiltered from GAMBL)
        GENOME_MAF,
        # QC plots
        rules.qc_data.output,
        # Swimmer plot
        rules.swimmer_plot.output,
        # FISH Nanostring
        rules.fish_ns_plots.output,
        # MiXCR Plots
        rules.plot_mixcr.output,
        # LymphGen
        rules.compare_lymphgen.output,
        # Shared mutations
        rules.shared_mutations.output,
        # PhyClone
        rules.pyclone_maf.output, 
        # Survival
        rules.survival_plots.output,
        # LySeq Validation
        rules.lyseq_validation.output,
        rules.lyseq_qc.output,
        # Trunk mutations
        rules.trunk_mutations.output,
        # Summary of everything
        rules.summarize_all.output, 
        rules.upset_plot.output
