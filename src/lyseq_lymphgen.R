source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv")

metadata <- trios_md %>%
  filter(tissue_status == "tumour", DNAseq_type == "genome") %>%
  group_by(patient_id) %>%
  filter(n() > 1) %>%
  slice_min(time_since_diagnosis_years, n=2) %>%
  mutate(relapse = ifelse(time_since_diagnosis_years == 0, "Diagnosis", "Relapse")) %>%
  ungroup() %>%
  select(patient_id, sample_id = DNAseq_sample_id, relapse_timing, relapse, time_since_diagnosis_years)

genome_only <- read_tsv("../../gambl-repos/gambl-mutect2-lhilton/versioned_results/LymphGen/GAMBL_all_the_things.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv") %>%
  filter(Sample.Name %in% metadata$sample_id) %>%
  rename(sample_id = Sample.Name) %>%
  mutate(mutations = "genome_only") %>%
  tidy_lymphgen(lymphgen_column_in = "Subtype.Prediction",
                lymphgen_column_out = "LymphGen",
                relevel = TRUE)


lyseq_only <- read_tsv("../LySeq_Validation/results/lymphgen-1.0/99-outputs/LySeqST_only.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv") %>%
  mutate(mutations = "lyseq_only") %>%
  rename(sample_id = Sample.Name) %>%
  left_join(metadata) %>%
  filter(!is.na(patient_id)) %>%
  tidy_lymphgen(lymphgen_column_in = "Subtype.Prediction",
                lymphgen_column_out = "LymphGen",
                relevel = TRUE)

genome_callable <- read_tsv("../LySeq_Validation/results/lymphgen-1.0/99-outputs/genome_callable.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv")%>%
  mutate(mutations = "genome_callable") %>%
  rename(sample_id = Sample.Name) %>%
  left_join(metadata) %>%
  filter(!is.na(patient_id)) %>%
  tidy_lymphgen(lymphgen_column_in = "Subtype.Prediction",
                lymphgen_column_out = "LymphGen",
                relevel = TRUE)

combined <- read_tsv("../LySeq_Validation/results/lymphgen-1.0/99-outputs/LyseqST_genome_callable.lymphgen_calls.with_cnvs.with_sv.with_A53.tsv") %>%
  mutate(mutations = "combined") %>%
  rename(sample_id = Sample.Name) %>%
  left_join(metadata) %>%
  filter(!is.na(patient_id)) %>%
  tidy_lymphgen(lymphgen_column_in = "Subtype.Prediction",
                lymphgen_column_out = "LymphGen",
                relevel = TRUE)

all_lymphgen <- genome_only %>%
  filter(sample_id %in% combined$sample_id) %>%
  mutate(LymphGen = factor(LymphGen, levels = levels(combined$LymphGen))) %>%
  bind_rows(lyseq_only) %>%
  bind_rows(combined) %>%
  bind_rows(genome_callable) %>%
  group_by(sample_id, LymphGen) %>%
  mutate(discordant = n() < 4) %>%
  ungroup() %>%
  mutate(mutations = factor(
    mutations,
    levels = c("genome_only", "combined", "lyseq_only", "genome_callable"),
    labels = c("Genome", "Combined", "LySeqST", "Genome Callable")
  ))

lg_alluvial_plot <- all_lymphgen %>%
  filter(mutations %in% c("Genome", "Combined")) %>%
  group_by(sample_id) %>% filter(n() == 2) %>%
  ungroup() %>%
  group_by(sample_id, LymphGen) %>%
  mutate(discordant = n() < 2) %>%
  ungroup() %>%
  ggplot(aes(
    x = mutations,
    stratum = LymphGen,
    alluvium = sample_id,
    fill = LymphGen
  )) +
  geom_flow(stat = "alluvium",
            lode.guidance = "frontback",
            aes(alpha = discordant)) +
  geom_stratum(colour = NA) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  scale_fill_manual(values = get_gambl_colours(classification = "LymphGen")[levels(all_lymphgen$LymphGen)[levels(all_lymphgen$LymphGen) %in% unique(all_lymphgen$LymphGen)]]) +
  xlab("Source of Variant Calls") +
  guides(alpha = "none")

lg_alluvial_plot

pivot_lymphgen_features <- function(df){
  df_out <- df %>%
    select(
      sample_id,
      LymphGen,
      matches("[.]Features")
    ) %>%
    pivot_longer(
      matches("[.]Features"),
      names_to = "feature_class",
      values_to = "feature"
    ) %>%
    separate(
      feature,
      into = as.character(seq_len(max(str_count(.$feature, ","), na.rm = TRUE) + 1)),
      sep = ","
    ) %>%
    pivot_longer(
      matches(as.character(1:50)),
      names_to = "junk",
      values_to = "feature"
    ) %>%
    select(-junk) %>%
    drop_na(feature) %>%
    mutate(feature_class = str_remove(feature_class, "[.]Features"))
  return(df_out)
}

lyseq_features <- pivot_lymphgen_features(lyseq_only) %>%
  mutate(lyseq_only = TRUE) %>%
  rename(LymphGen_LySeq = LymphGen)
genome_features <- pivot_lymphgen_features(genome_only) %>%
  mutate(genome_only = TRUE) %>%
  rename(LymphGen_genome = LymphGen)
combined_features <- pivot_lymphgen_features(combined) %>%
  mutate(combined = TRUE) %>%
  rename(LymphGen_combined = LymphGen)
callable_features <- pivot_lymphgen_features(genome_callable) %>%
  mutate(genome_callable = TRUE) %>%
  rename(LymphGen_callable = LymphGen)


all_features <- lyseq_features %>%
  full_join(genome_features) %>%
  full_join(combined_features) %>%
  full_join(callable_features) %>%
  mutate(across(c(lyseq_only, genome_only, combined, genome_callable), ~replace_na(.x, FALSE))) %>%
  group_by(sample_id) %>%
  mutate(across(matches("LymphGen"), ~ .x[!is.na(.x)][1])) %>%
  ungroup() %>%
  filter(!is.na(LymphGen_LySeq))

all_features %>%
  group_by(feature) %>%
  count(feature_class, genome_only, combined) %>%
  mutate(
    total = sum(n),
    percent = round(n/sum(n)*100)
  ) %>%
  ungroup() %>%
  View()









