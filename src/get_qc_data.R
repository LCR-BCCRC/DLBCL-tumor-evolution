source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>%
  filter(tissue_status == "tumour") %>%
  rename(
    sample_id = DNAseq_sample_id,
    seq_type = DNAseq_type
  )

qc_data <- read_tsv("/projects/rmorin/projects/gambl-repos/gambl-lhilton/results/gambl/qc-1.0/99-outputs/genome.qc_metrics.tsv") %>%
  select(sample_id = UID, AverageBaseQuality:ProportionCoverage30x) %>%
  filter(sample_id %in% trios_md$sample_id)

qc_exomes <- read_tsv("/projects/rmorin/projects/gambl-repos/gambl-lhilton/results/gambl/qc-1.0/99-outputs/capture.qc_metrics.tsv") %>%
  select(sample_id = UID, AverageBaseQuality:ProportionCoverage30x) %>%
  filter(
    sample_id %in% trios_md$sample_id,
    !sample_id %in% qc_data$sample_id
  )

trios_maf <- read_tsv("data/maf/genomes_augmented.grch37.maf")

depth_from_maf <- trios_maf %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(coding_depth = ifelse(Variant_Classification %in% GAMBLR:::coding_class, t_depth, NA)) %>%
  summarize(
    MeanVariantDepth = mean(t_depth),
    MeanCodingDepth = mean(coding_depth, na.rm = TRUE),
    NumberOfVariants = n()
  )


trios_qc <- trios_md %>%
  left_join(bind_rows(qc_data, qc_exomes)) %>%
  left_join(depth_from_maf, by = c("sample_id" = "Tumor_Sample_Barcode")) %>%
  mutate(AltCovEstimate = (TotalUniquelyMapped - TotalDuplicatedreads) * AverageReadLength / 3101788170) %>%
  select(
    sample_id,
    patient_id,
    seq_type,
    relapse_timing,
    MeanVariantDepth,
    MeanCorrectedCoverage,
    ProportionCoverage10x,
    ProportionCoverage30x,
    AltCovEstimate,
    NumberOfVariants
  ) %>%
  pivot_longer(matches("Mean|Prop|Alt"), names_to = "metric", values_to = "value") %>%
  mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours)))

write_tsv(trios_qc, "data/qc/gambl_qc_combined.tsv")

compare_depths <- trios_qc %>%
  filter(seq_type == "genome", str_detect(metric, "Mean|Alt")) %>%
  group_by(metric) %>%
  rstatix::pairwise_wilcox_test(value ~ relapse_timing) %>%
  add_y_position(step.increase = 0.15)

strip_labs <- c("Mean Variant Depth", "Mean Corrected Coverage", "Raw Coverage Estimate")
names(strip_labs) <- c("MeanVariantDepth", "MeanCorrectedCoverage", "AltCovEstimate")


qc_boxplot <-
  trios_qc %>%
  filter(seq_type == "genome", str_detect(metric, "Mean|Alt")) %>%
  ggplot(
    aes(
      x = relapse_timing,
      y = value
    )
  ) +
  geom_boxplot(aes(fill = relapse_timing)) +
  geom_quasirandom() +
  geom_hline(yintercept = 25) +
  facet_wrap(~metric, scales = "fixed", labeller = labeller(metric = strip_labs)) +
  stat_pvalue_manual(filter(compare_depths, str_detect(metric, "Mean|Alt")), tip.length = 0.01) +
  scale_fill_manual(values = relapse_colours) +
  scale_x_discrete(labels = relapse_timing_labels) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125), limits = c(0, 150)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")

ggsave("figures/qc_boxplots.pdf", qc_boxplot, height = 6, width = 12)
saveRDS(qc_boxplot, "figures/qc_boxplots.RDS")

trios_qc %>%
  filter(metric == "MeanCorrectedCoverage") %>%
  filter(NumberOfVariants < 50000) %>%
  ggplot(aes(x = value, y = NumberOfVariants)) +
  geom_point(alpha = 0.8) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black"
  ) +
  stat_cor() +
  xlab("MeanCorrectedCoverage") +
  facet_wrap(~seq_type, scales = "free")

ggsave("figures/coverage_vs_mutations.pdf", height = 5, width = 9)

# shared_vars <- read_tsv("data/shared_mutations/unique_mutations_per_patient.tsv")

# shared_vars_vs_depth <- trios_qc %>%
#   filter(!str_detect(sample_id, "normal")) %>%
#   left_join(shared_vars, by = c("patient_id", "relapse_timing")) %>%
#   mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours)))

# broom::tidy(lm(unique_muts_percent ~ value + relapse_timing, data = filter(shared_vars_vs_depth, metric.x == "MeanCorrectedCoverage"))) %>%
#   mutate(term = str_replace(term, "value", "MeanCorrectedCoverage"))

# broom::tidy(lm(unique_mutations ~ value + relapse_timing, data = filter(shared_vars_vs_depth, metric.x == "MeanCorrectedCoverage"))) %>%
#   mutate(term = str_replace(term, "value", "MeanCorrectedCoverage"))

# ggplot(shared_vars_vs_depth, aes(x = value, y = unique_muts_percent)) +
#   geom_point() +
#   facet_wrap(~metric.x)
