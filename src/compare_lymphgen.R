##### Compare LymphGen classifications #####

source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>%
  filter(tissue_status == "tumour")

lymphgen_raw <-
  read_tsv("/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/results/all_the_things/lymphgen-1.0/99-outputs/GAMBL_all_the_things.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv")

lymphgen_capture_raw <-
  read_tsv("../LySeq_Validation/results/lymphgen-1.0/99-outputs/LyseqST_genome_callable.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv")

lymphgen_raw <- lymphgen_raw %>%
  filter(!Sample.Name %in% lymphgen_capture_raw$Sample.Name) %>%
  bind_rows(lymphgen_capture_raw)

lymphgen <- trios_md %>%
  left_join(lymphgen_raw, by = c("DNAseq_sample_id" = "Sample.Name")) %>%
  filter(!is.na(Subtype.Prediction)) %>%
  tidy_lymphgen(
    lymphgen_column_in = "Subtype.Prediction",
    lymphgen_column_out = "LymphGen",
    relevel = TRUE
  ) %>%
  mutate(
    relapse = case_when(
      time_since_diagnosis_years == 0 ~ "Diagnosis",
      time_since_diagnosis_years > 0 ~ "Relapse",
      time_since_diagnosis_years < 0 ~ "Pre-Diagnosis"
    ),
    relapse_timing = factor(relapse_timing, levels = names(relapse_colours))
  ) %>%
  group_by(patient_id) %>%
  mutate(num_tumours = n()) %>%
  ungroup() %>%
  group_by(patient_id, LymphGen) %>%
  mutate(discordant = n() < num_tumours) %>%
  ungroup() %>%
  select(patient_id,
    sample_id = DNAseq_sample_id,
    relapse_timing,
    relapse,
    time_since_diagnosis_years,
    DNAseq_type,
    LymphGen,
    discordant,
    any_of(colnames(lymphgen_raw))
  )

summarize_lg <- lymphgen %>%
  filter(time_since_diagnosis_years >= 0) %>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  slice_min(time_since_diagnosis_years, n = 2, with_ties = FALSE) %>%
  ungroup() %>%
  select(patient_id, relapse, Subtype.Prediction, discordant, relapse_timing) %>%
  mutate(all_subclasses = str_replace_all(Subtype.Prediction, "/", "|")) %>%
  pivot_wider(
    names_from = relapse,
    values_from = c(Subtype.Prediction, all_subclasses),
    names_glue = "{relapse}_{.value}"
  ) %>%
  mutate(discordant_type = case_when(
    discordant &
      (
        Relapse_Subtype.Prediction == "Other" |
          Diagnosis_Subtype.Prediction == "Other"
      ) ~ "Partial",
    discordant &
      !str_detect(
        Diagnosis_Subtype.Prediction,
        Relapse_all_subclasses
      ) &
      !str_detect(
        Relapse_Subtype.Prediction,
        Diagnosis_all_subclasses
      ) ~ "Frank",
    discordant ~ "Partial",
    !discordant & Relapse_Subtype.Prediction == "Other" ~ "Unclassified",
    !discordant ~ "Concordant"
  )) %>%
  select(-matches("all_subclasses"))

write_tsv(summarize_lg, "data/lymphgen/lymphgen_summary.tsv")

lymphgen <- lymphgen %>%
  left_join(select(summarize_lg, patient_id, discordant_type))

write_tsv(lymphgen, "data/lymphgen/lymphgen_results.tsv")

lg_alluvial <- lymphgen %>%
  filter(time_since_diagnosis_years >= 0) %>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  slice_min(time_since_diagnosis_years, n = 2, with_ties = FALSE) %>%
  ungroup()

lg_alluvial_plot <- lg_alluvial %>%
  ggplot(aes(
    x = relapse,
    stratum = LymphGen,
    alluvium = patient_id,
    fill = LymphGen
  )) +
  geom_flow(
    stat = "alluvium",
    lode.guidance = "frontback",
    aes(alpha = discordant_type == "Frank")
  ) +
  geom_stratum(colour = NA) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  scale_fill_manual(values = get_gambl_colours(classification = "LymphGen")[levels(lg_alluvial$LymphGen)[levels(lg_alluvial$LymphGen) %in% unique(lg_alluvial$LymphGen)]]) +
  facet_grid2(~relapse_timing, strip = relapse_strips) +
  xlab("") +
  guides(alpha = "none")

lg_alluvial_plot

ggsave("figures/lymphgen_comparison.pdf", height = 4, width = 10)
saveRDS(lg_alluvial_plot, "figures/lymphgen_comparison.RDS")
