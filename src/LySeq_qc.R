source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv")
trios_samples <- trios_md %>%
  filter(!lyseq_library_id %in% c("pending", "failed", NA)) %>%
  pull(DNAseq_sample_id)

qc <- read_tsv("../LySeq_Validation/results/QC/04-mergedMetrics/capture--grch37.hs_metrics.txt") %>%
  filter(sampleID %in% trios_samples) %>%
  left_join(select(trios_md, sampleID = DNAseq_sample_id, relapse_timing)) %>%
  mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours))) %>%
  select(sampleID, MEAN_TARGET_COVERAGE, PCT_TARGET_BASES_50X, PCT_TARGET_BASES_100X, relapse_timing) %>%
  group_by(relapse_timing) %>%
  mutate(MEAN_COVERAGE_GROUP = mean(MEAN_TARGET_COVERAGE)) %>%
  ungroup()

write_tsv(qc, "data/qc/lst_qc_coverage.tsv")


compare_depths <- qc %>%
  rstatix::pairwise_wilcox_test(MEAN_TARGET_COVERAGE ~ relapse_timing) %>%
  add_y_position(step.increase = 0.5) %>%
  mutate(y.position = log10(y.position))

lst_qc_plot <- qc %>%
  ggplot(aes(x = relapse_timing, y = MEAN_TARGET_COVERAGE)) +
  geom_boxplot(aes(fill = relapse_timing)) +
  geom_quasirandom() +
  geom_hline(yintercept = 100) +
  scale_y_log10(breaks = c(10, 100, 1000, 10000)) +
  stat_pvalue_manual(compare_depths, tip.length = 0.005) +
  scale_fill_manual(values = relapse_colours, name = "") +
  scale_x_discrete(labels = relapse_timing_labels) +
  xlab("") +
  ylab("Mean Target Coverage (Corected)")

ggsave("figures/lst_qc.pdf", lst_qc_plot, height = 5, width = 6)
saveRDS(lst_qc_plot, "figures/lst_qc.RDS")
