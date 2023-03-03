
library(tidyverse)
library(RefManageR)

# Sequencing summary
cohort_summary <- read_tsv("data/metadata/cohort_summary.tsv")
num_patients_seq <- length(cohort_summary$patient_id)
num_per_seq_type <- cohort_summary %>%
  count(seq_type) %>%
  pull(n, name = seq_type)

all_summary <- read_tsv("data/metadata/summarize_all_assays.tsv") %>%
  mutate(denovo_percent = round(is_denovo / lowgrade_data_available * 100, digits = 1)) %>%
  mutate(transformed_percent = 100 - denovo_percent)

seq_md <- read_tsv("data/metadata/metadata_sequencing.tsv")
lst_count <- seq_md %>%
  filter(
    tissue_status == "tumour",
    !lyseq_library_id %in% c("failed", "pending", NA)
  ) %>%
  group_by(patient_id) %>%
  summarize(lyseq_assays = sum(n() > 1)) %>%
  ungroup() %>%
  count(lyseq_assays) %>%
  pull(n)
names(lst_count) <- c("one_lst", "two_or_more_lst")

# QC
qc <- read_tsv("data/qc/gambl_qc_combined.tsv")
mean_coverage <- qc %>%
  filter(metric == "MeanVariantDepth") %>%
  group_by(seq_type) %>%
  summarize(mean_depth = mean(value, na.rm = TRUE), num_samples = n()) %>%
  ungroup() %>%
  column_to_rownames("seq_type")

lst_qc <- read_tsv("data/qc/lst_qc_coverage.tsv")
lst_coverage <- mean(lst_qc$MEAN_TARGET_COVERAGE)


# FISH

fish_summary <- read_tsv("data/fish_nanostring/fish_summary.tsv")

fish_md <- read_tsv("data/metadata/metadata_fish_nanostring.tsv")
ever_pos <- fish_md %>%
  select(patient_id, biopsy_id, matches("ba_consensus")) %>%
  pivot_longer(
    matches("ba_consensus"),
    names_to = "assay",
    values_to = "result"
  ) %>%
  select(-biopsy_id) %>%
    group_by(patient_id) %>%
    slice_max(result == "POS", n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(assay = toupper(str_remove(assay, "_ba_consensus"))) %>%
    group_by(assay) %>% 
    summarize(ever_positive = sum(result == "POS")) %>% 
    ungroup()

total_discordant <- fish_summary %>%
  add_row(
    assay = "BCL2",
    discordant = "Discordant",
    n = 0
  ) %>%
  left_join(ever_pos) %>% 
  group_by(assay) %>%
  mutate(total_assay = min(total_assay, na.rm = TRUE)) %>%
  filter(discordant == "Discordant") %>%
  summarize(total_discordant = sum(n), total_assay, ever_positive) %>%
  mutate(percent_discordant = round(total_discordant / total_assay * 100), 
  percent_discordant_pos = round(total_discordant / ever_positive * 100)) %>%
  ungroup() %>%
  distinct() %>%
  column_to_rownames("assay")

# Structural Variants

sv_summary <- read_tsv("data/svs/summarize_multi_FISH_pos.tsv")

sv_list <- sv_summary %>%
  column_to_rownames("FISH")

sv_sensitivity <- read_tsv("data/svs/all_svs.tsv")
sv_sensitivity <- sv_sensitivity %>%
  filter(
    FISH_result == "POS",
    seq_type == "genome"
  ) %>%
  group_by(FISH, tumour_sample_id) %>%
  slice_max(!is.na(partner), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(FISH) %>%
  summarize(
    bp_found = sum(!is.na(partner)),
    total_pos = n()
  ) %>%
  mutate(sensitivity = round((bp_found / total_pos) * 100, digits = 1)) %>%
  column_to_rownames("FISH")


# IG Rearrangements

mixcr <- read_tsv("data/ig_rearrangement/summary_ig_concordance.tsv")

mixcr <- mixcr %>%
  pivot_wider(
    names_from = discordant,
    values_from = n,
    values_fill = 0
  ) %>%
  rename(discordant = `TRUE`, concordant = `FALSE`) %>%
  mutate(total = discordant + concordant) %>%
  mutate(item = str_c(str_replace(tolower(relapse_timing), " ", "_"), IG_gene, sep = "_"))

mixcr_list <- mixcr %>%
  column_to_rownames("item") %>%
  split(., seq(nrow(.)))
mixcr_list <- setNames(mixcr_list, mixcr$item)

# COO

d90 <- read_tsv("data/fish_nanostring/ns_summary.tsv")

d90_assays <- d90$total_assay[1]

d90 <- d90 %>%
  filter(!is.na(discordant), time_point == "Diagnosis") %>%
  select(-total_time_point, -percent, -total_assay, -time_point) %>%
  pivot_wider(
    names_from = "discordant",
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    total_assayed = Discordant + Concordant,
    percent_discordant = round(Discordant / total_assayed * 100, digits = 1),
    percent_concordant = round(Concordant / total_assayed * 100, digits = 1)
  )

dhitsig_list <- d90 %>%
  filter(assay == "dhitsig") %>%
  column_to_rownames("relapse_timing")

coo_list <- d90 %>%
  filter(assay == "dlbcl") %>%
  column_to_rownames("relapse_timing")


# LymphGen

summary_lymphgen <- read_tsv("data/lymphgen/lymphgen_summary.tsv")

lg_table <- table(summary_lymphgen$relapse_timing, summary_lymphgen$discordant_type)

# Constrained Mutations

constrained <- read_tsv("data/phyclone/tabulate_constrained_pt.tsv") %>%
  filter(num_constrained_genes >= 2)
divergent <- read_tsv("data/phyclone/mutations_divergent_patients.tsv")
divergent_count <- divergent %>%
  select(patient_id) %>%
  left_join(select(summary_lymphgen, patient_id, relapse_timing)) %>%
  distinct() %>%
  count(relapse_timing) %>%
  pull(n, name = relapse_timing)


# Sequencing Outcomes Overlap
outcomes <- read_tsv("data/survival/surv_wide.tsv") %>%
  pull(ID)

overlap <- cohort_summary %>%
filter(patient_id %in% outcomes) %>%
count(relapse_timing) %>%
pull(n, relapse_timing)
overlap = c(overlap, "total" = sum(unname(overlap)))
