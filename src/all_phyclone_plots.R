##### Generate publication-ready clonal evolution figures #####

source("src/libs.R")
source("src/phyclone_plot_function.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>%
  select(-matches("RNAseq")) %>%
  rename_with(~ str_remove(.x, "DNAseq_")) %>%
  rename("seq_type" = "type") %>%
  filter(tissue_status == "tumour")

biopsy_timing_all <- trios_md %>%
  select(
    patient_id,
    Tumor_Sample_Barcode = sample_id,
    time_since_diagnosis_years,
    relapse_timing
  )

trios_lymphgen <- read_tsv(
  "../LySeq_Validation/results/lymphgen-1.0/99-outputs/LyseqST_genome_callable.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv"
)

results_dir <- "../LySeq_Validation/results/pyclone_vi-1.0_indels/99-outputs/genome--grch37/"
maf_dir <- "../LySeq_Validation/results/pyclone_vi-1.0_indels/00-inputs/maf/genome--grch37/"
plot_dir <- "figures/phyclone_plots_lymphgen_lst_callable_indels/"

pyclone_mafs <- lapply(
  unique(trios_md$patient_id),
  # "03-23488",
  generate_plots,
  maf_dir = maf_dir,
  results_dir = results_dir,
  lymphgen_data = trios_lymphgen,
  biopsy_timing = biopsy_timing_all,
  plot_dir = plot_dir
  # feature_set = "all"
)

pyclone_mafs2 <- bind_rows(pyclone_mafs) %>%
  select(-matches("label")) %>%
  mutate(across(everything(), ~ str_replace_all(.x, "\n", "-")))

write_tsv(pyclone_mafs2, "data/phyclone/phyclone_mafs.tsv")

# 03-23488
# Has many coding mutations
# Specify good ones to plot here

pt_03_23488 <- generate_plots(
  patient = "03-23488",
  maf_dir = maf_dir,
  results_dir = results_dir,
  lymphgen_data = trios_lymphgen,
  biopsy_timing = biopsy_timing_all,
  plot_dir = plot_dir,
  custom_genes = c("FOXO1", "SGK1", "RB1", "CREBBP", "MYC", "KLHL6", "BTG2", "BCL6")
)

pt_13_26835 <- generate_plots(
  patient = "13-26835",
  maf_dir = maf_dir,
  results_dir = results_dir,
  lymphgen_data = trios_lymphgen,
  biopsy_timing = biopsy_timing_all,
  plot_dir = plot_dir,
  custom_genes = c("KLF2", "SOCS1", "OSBPL10", "ACTB", "HLA-A", "HLA-B")
)



# pt_14_20552 <- generate_plots(
#   patient = "14-20552",
#   maf_dir = maf_dir,
#   results_dir = results_dir,
#   lymphgen_data = trios_lymphgen,
#   biopsy_timing = biopsy_timing_all,
#   plot_dir = plot_dir
#   # custom_genes = c("KLF2", "SOCS1", "OSBPL10", "ACTB", "HLA-A", "HLA-B")
# )

# pt_07_32561 <- generate_plots(
#   patient = "07-32561",
#   maf_dir = maf_dir,
#   results_dir = results_dir,
#   lymphgen_data = trios_lymphgen,
#   biopsy_timing = biopsy_timing_all,
#   plot_dir = plot_dir
#   # custom_genes = c("KLF2", "SOCS1", "OSBPL10", "ACTB", "HLA-A", "HLA-B")
# )
