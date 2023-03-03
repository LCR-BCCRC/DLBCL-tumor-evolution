##### Shared mutation analysis

source("src/libs.R")
library("lemon")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>%
  filter(
    tissue_status == "tumour",
    !is.na(DNAseq_sample_id),
    pathology != "FL",
    time_since_diagnosis_years >= 0
  ) %>%
  group_by(patient_id, DNAseq_type) %>%
  slice_min(dtbx, n = 2) %>%
  filter(n() == 2) %>%
  mutate(time_point = ifelse(
    time_since_diagnosis_years == 0,
    "Diagnosis",
    "Relapse"
  )) %>%
  ungroup() %>%
  filter(!patient_id %in% c(
    "01-16433" # One time point is an extreme outlier for number of mutations
  )) %>%
  rename(seq_type = DNAseq_type)

genome_cov <- read_tsv("data/qc/gambl_qc_combined.tsv") %>%
  filter(
    metric == "MeanCorrectedCoverage",
    sample_id %in% trios_md$DNAseq_sample_id
  ) %>%
  select(Tumor_Sample_Barcode = sample_id, genome_coverage = value)

pairs_table <- select(trios_md, patient_id, DNAseq_sample_id)
pairs_table <- left_join(
  pairs_table,
  pairs_table,
  by = "patient_id"
) %>%
  filter(DNAseq_sample_id.x != DNAseq_sample_id.y) %>%
  rename(sample_id = DNAseq_sample_id.x, partner_sample_id = DNAseq_sample_id.y) %>%
  select(-patient_id) %>%
  left_join(genome_cov, by = c("sample_id" = "Tumor_Sample_Barcode")) %>%
  left_join(genome_cov, by = c("partner_sample_id" = "Tumor_Sample_Barcode"), suffix = c("", "_partner")) %>%
  select(-partner_sample_id)


# genome_maf <- read_tsv("data/maf/genomes_augmented.grch37.maf")
genome_maf <- read_tsv("../LySeq_Validation/results/slms_3-1.0_vcf2maf-1.3/level_3/capture--grch37/LySeqST_genomes_merged_with_0_reads.maf")

maf_dt <- genome_maf %>%
  as.data.table()

setkey(maf_dt, Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position)

full_cn <- get_sample_cn_segments() %>%
  filter(ID %in% trios_md$DNAseq_sample_id) %>%
  rename(
    Tumor_Sample_Barcode = ID,
    Chromosome = chrom,
    Start_Position = start,
    End_Position = end
  ) %>%
  as.data.table()

setkey(full_cn, Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position)

maf_cn <- foverlaps(maf_dt, full_cn) %>%
  select(-c(Start_Position, End_Position)) %>%
  select(
    Start_Position = i.Start_Position,
    End_Position = i.End_Position,
    everything()
  ) %>%
  select(
    all_of(colnames(genome_maf)),
    LOH_flag,
    CN
  ) %>%
  distinct()

maf_md <- trios_md %>%
  select(patient_id,
    seq_type,
    Tumor_Sample_Barcode = DNAseq_sample_id,
    relapse_timing,
    time_point
  ) %>%
  left_join(maf_cn) %>%
  left_join(pairs_table, by = c("Tumor_Sample_Barcode" = "sample_id")) %>%
  filter(!is.na(Start_Position))

vaf_df <- maf_md %>%
  mutate(VAF = ifelse(
    t_alt_count >= 3,
    t_alt_count / t_depth,
    0
  )) %>%
  select(
    patient_id,
    seq_type,
    relapse_timing,
    time_point,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    Hugo_Symbol,
    HGVSp_Short,
    Variant_Classification,
    Variant_Type,
    t_depth,
    VAF,
    CN,
    LOH_flag,
    matches("genome_coverage")
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = time_point,
    values_from = c(t_depth, VAF, CN, LOH_flag, genome_coverage, genome_coverage_partner),
    names_glue = "{time_point}_{.value}",
    id_cols =
    ) %>%
  mutate(callable = ifelse(
    Diagnosis_t_depth >= 10 & Relapse_t_depth >= 10,
    TRUE,
    FALSE
  )) %>%
  mutate(shared = case_when(
    Diagnosis_VAF > 0 & Relapse_VAF > 0 & callable ~ 1,
    callable ~ 0
  )) %>%
  filter(!is.na(Diagnosis_t_depth) & !is.na(Relapse_t_depth)) %>%
  mutate(coding = Variant_Classification %in% GAMBLR:::coding_class)


dropout <- vaf_df %>%
  group_by(patient_id) %>%
  summarize(
    total_mutations = n(),
    total_callable = sum(callable == TRUE, na.rm = TRUE),
    Diagnosis_genome_coverage,
    Relapse_genome_coverage,
    seq_type
  ) %>%
  ungroup() %>%
  distinct() %>%
  mutate(dropout = round(100 * (1 - (total_callable / total_mutations)), digits = 1))

dropout %>%
  pivot_longer(matches("genome_coverage"),
    names_to = "time_point",
    values_to = "genome_coverage"
  ) %>%
  mutate(time_point = str_remove(time_point, "_genome_coverage")) %>%
  ggplot(aes(x = genome_coverage, y = dropout)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(p.accuracy = 0.001) +
  facet_wrap(~time_point) +
  geom_hline(yintercept = 30) +
  ylim(0, 40)

vaf_df <- vaf_df %>%
  filter(callable) %>%
  filter(!patient_id %in% dropout[dropout$dropout > 30, ]$patient_id)

# Make some oncoplots of what was kept vs. dropped
maf_callable <- maf_md %>%
  left_join(select(
    vaf_df,
    patient_id,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    callable
  )) %>%
  filter(
    !is.na(callable),
    t_alt_count > 2
  )

trios_md$sample_id <- trios_md$DNAseq_sample_id
trios_md$Tumor_Sample_Barcode <- trios_md$DNAseq_sample_id
maf_callable <- maftools::read.maf(maf_callable, clinicalData = trios_md)

pdf("figures/oncoplot_callable_genes.pdf", height = 7, width = 9)
prettyOncoplot(
  maf_callable,
  these_samples_metadata = trios_md,
  metadataColumns = "relapse_timing",
  minMutationPercent = 5,
  genes = unique(GAMBLR::lymphoma_genes_comprehensive$Gene),
  custom_colours = list("relapse_timing" = relapse_colours),
  splitColumnName = "relapse_timing"
)
dev.off()

maf_dropped <- maf_md %>%
  left_join(select(
    vaf_df,
    patient_id,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    callable
  )) %>%
  filter(
    is.na(callable),
    t_alt_count > 2
  )

maf_dropped_obj <- maftools::read.maf(maf_dropped, clinicalData = trios_md)

# What's been dropped from at least two patients?
dropped_genes <- maf_dropped %>%
  filter(Variant_Classification %in% GAMBLR:::coding_class) %>%
  select(Hugo_Symbol, patient_id) %>%
  distinct() %>%
  count(Hugo_Symbol) %>%
  filter(n > 1) %>%
  pull(Hugo_Symbol)

pdf("figures/oncoplot_dropped_genes.pdf", height = 7, width = 9)
prettyOncoplot(
  maf_dropped_obj,
  these_samples_metadata = trios_md,
  metadataColumns = "relapse_timing",
  genes = unique(GAMBLR::lymphoma_genes_comprehensive$Gene),
  minMutationPercent = 1,
  custom_colours = list("relapse_timing" = relapse_colours),
  splitColumnName = "relapse_timing"
)
dev.off()

tabulate_shared <- vaf_df %>%
  select(-matches("CN|LOH")) %>%
  pivot_longer(
    matches("VAF|coverage"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  extract(metric, into = c("time_point", "metric"), "(Relapse|Diagnosis)_(.*+)") %>%
  pivot_wider(
    names_from = "metric",
    values_from = "value"
  ) %>%
  filter(VAF > 0, !is.na(VAF)) %>%
  group_by(patient_id, time_point, seq_type) %>%
  mutate(
    total.muts_all = n(),
    shared.muts_all = sum(shared, na.rm = TRUE),
    total.muts_coding = sum(coding, na.rm = TRUE),
    shared.muts_coding = sum(shared & coding, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(patient_id, relapse_timing, time_point, seq_type, matches("muts|coverage")) %>%
  distinct() %>%
  pivot_longer(
    cols = matches("muts"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  separate(metric, into = c("tally", "metric"), sep = "_") %>%
  pivot_wider(
    names_from = tally,
    values_from = value
  ) %>%
  mutate(unique_muts_percent = 100 - round(shared.muts / total.muts * 100, digits = 1)) %>%
  mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours))) %>%
  mutate(metric = ifelse(metric == "coding", "Coding", "All")) %>%
  mutate(unique_mutations = total.muts - shared.muts) %>%
  filter(!(seq_type == "capture" & metric == "All"))

write_tsv(
  tabulate_shared,
  "data/shared_mutations/unique_mutations_per_patient.tsv"
)

# Make MAFs and Oncoplots of unique vs. shared variants
maf_annotated_shared <- maf_callable@data %>%
  group_by(
    patient_id,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2
  ) %>%
  mutate(shared = case_when(
    n() == 2 ~ "Shared",
    TRUE ~ time_point
  )) %>%
  slice_max(time_point == "Diagnosis", n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Tumor_Sample_Barcode = patient_id) %>%
  mutate(hot_spot = case_when(
    Hugo_Symbol == "EZH2" & str_detect(HGVSp_Short, "p.Y646") ~ TRUE,
    Hugo_Symbol == "MYD88" & HGVSp_Short == "p.L265P" ~ TRUE,
    Hugo_Symbol == "CD79B" & str_detect(HGVSp_Short, "p.Y197") ~ TRUE,
    Hugo_Symbol %in% c("NOTCH1", "NOTCH2") & Exon_Number == "34/34" & str_detect(Variant_Classification, ("Nonsense|Frame_Shift")) ~ TRUE,
    Hugo_Symbol == "CREBBP" & Start_Position > GAMBLR::hotspot_regions_grch37["CREBBP", "start"] & End_Position < GAMBLR::hotspot_regions_grch37["CREBBP", "end"] ~ TRUE,
    Hugo_Symbol == "MEF2B" & str_detect(HGVSp_Short, "D83") ~ TRUE,
    Hugo_Symbol == "STAT6" & str_detect(HGVSp_Short, "D419") ~ TRUE
  ))

patient_metadata <- tabulate_shared %>%
  filter(metric == "Coding") %>%
  select(
    Tumor_Sample_Barcode = patient_id,
    relapse_timing,
    time_point,
    unique_muts_percent
  ) %>%
  pivot_wider(
    names_from = time_point,
    values_from = unique_muts_percent,
    names_glue = "{time_point}_{.value}"
  ) %>%
  rowwise() %>%
  mutate(
    mean_unique = mean(c(
      Relapse_unique_muts_percent,
      Diagnosis_unique_muts_percent
    ))
  ) %>%
  ungroup() %>%
  mutate(
    relapse_timing = factor(relapse_timing, levels = names(relapse_colours))
  ) %>%
  mutate(
    sample_id = Tumor_Sample_Barcode
  ) %>%
  rename(
    `Diagnosis % Unique` = Diagnosis_unique_muts_percent,
    `Relapse % Unique` = Relapse_unique_muts_percent,
    `Mean % Unique` = mean_unique,
    `Relapse Timing` = relapse_timing
  )

load_mafs <- function(mut_cat, maf_list) {
  maf_filt <- maf_annotated_shared %>%
    filter(shared == mut_cat)
  maf_obj <- maftools::read.maf(
    maf_filt,
    clinicalData = patient_metadata
  )
  maf_list[[mut_cat]] <- maf_obj
  return(maf_list)
}

maf_list <- list()
maf_list <- unlist(lapply(c("Shared", "Diagnosis", "Relapse"), load_mafs, maf_list = maf_list))

# Define the most frequently mutated genes among shared variants

plotgenes <-
  c(
    "ACTB", "ARID1A", "BCL2", "BCL6", "BTG1",
    "BTG2", "CD70", "CD79B", "CREBBP", "DTX1", "EP300", "ETS1",
    "HLA-B", "IRF8", "KLF2", "KLHL14", "KMT2D", "MEF2B",
    "MYD88", "NOTCH2", "OSBPL10", "PIM1", "SETD1B", "SGK1",
    "SOCS1", "STAT6", "TBL1XR1", "TNFAIP3", "TNFRSF14", "TP53",
    "HLA-A", "DUSP2", "ETV6", "ITPKB", "NOTCH1"
  )

col_fun <- colorRamp2(c(0, 50, 100), c("blue", "white", "red"))

# Run this once to get the gene order
oncoplot_shared <- prettyOncoplot(
  maf_list$Shared,
  these_samples_metadata = patient_metadata,
  metadataColumns = "Relapse Timing",
  # highlightHotspots = TRUE,
  genes = plotgenes,
  removeNonMutated = FALSE
)

# Make gene order vector usable by next plot functions
roworder <- c(plotgenes)[row_order(oncoplot_shared)]

oncoplot_fun <- function(mut_cat) {
  outfile <- paste0("figures/oncoplot_", mut_cat, ".pdf")
  pdf(outfile, height = 6, width = 7)
  col_fun <- colorRamp2(c(0, 50, 100), c("blue", "white", "red"))
  prettyOncoplot(
    maf_list[[mut_cat]],
    these_samples_metadata = with(
      patient_metadata,
      patient_metadata[order(`Relapse Timing`), ]
    ),
    metadataColumns = "Relapse Timing",
    numericMetadataColumns = c(
      "Diagnosis % Unique",
      "Relapse % Unique",
      "Mean % Unique"
    ),
    keepSampleOrder = TRUE,
    sortByColumns = "Mean % Unique",
    highlightHotspots = TRUE,
    genes = roworder,
    keepGeneOrder = TRUE,
    splitColumnName = "Relapse Timing",
    removeNonMutated = FALSE,
    custom_colours = list(
      "Relapse Timing" = relapse_colours,
      "Diagnosis % Unique" = col_fun,
      "Relapse % Unique" = col_fun,
      "Mean % Unique" = col_fun
    ),
    hideTopBarplot = FALSE,
    tally_all_mutations = TRUE,
    legendFontSize = 8,
    fontSizeGene = 8,
    metadataBarFontsize = 8,
    metadataBarHeight = 2
  )
  dev.off()
}

lapply(c("Shared", "Diagnosis", "Relapse"), oncoplot_fun)



summary_means <- tabulate_shared %>%
  group_by(relapse_timing, metric, time_point) %>%
  summarize(across(where(is.numeric),
    ~ round(mean(.x, na.rm = TRUE), digits = 2),
    .names = "mean_{.col}"
  ))

write_tsv(
  summary_means,
  "data/shared_mutations/mutation_tally_by_group.tsv"
)

big_lm <- tabulate_shared %>%
  rename(mutation_type = metric) %>%
  select(-total.muts) %>%
  pivot_longer(matches("mut") & where(is.numeric),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(mutation_type, time_point, metric) %>%
  nest() %>%
  mutate(model = map(
    data, ~ broom::tidy(
      lm(value ~ relapse_timing + genome_coverage + genome_coverage_partner,
        data = .
      )
    )
  )) %>%
  select(-data) %>%
  unnest(c(model)) %>%
  mutate(p.format = p_format(p.value, accuracy = 0.001))

write_tsv(big_lm, "data/shared_mutations/lm_mutations_vs_relapse_coverage_categorical.tsv")

big_lm_cont <- tabulate_shared %>%
  left_join(select(
    trios_md,
    patient_id,
    time_point,
    time_since_diagnosis_years
  )) %>%
  group_by(patient_id) %>%
  mutate(time_since_diagnosis_years = max(time_since_diagnosis_years)) %>%
  ungroup() %>%
  rename(mutation_type = metric) %>%
  select(-total.muts) %>%
  pivot_longer(matches("mut") & where(is.numeric),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(mutation_type, time_point, metric) %>%
  nest() %>%
  mutate(model = map(
    data, ~ broom::tidy(
      lm(value ~ time_since_diagnosis_years + genome_coverage + genome_coverage_partner,
        data = .
      )
    )
  )) %>%
  select(-data) %>%
  unnest(c(model)) %>%
  mutate(p.format = p_format(p.value, accuracy = 0.001))

write_tsv(big_lm_cont, "data/shared_mutations/lm_mutations_vs_relapse_coverage_continuous.tsv")

summary_stats_counts <- tabulate_shared %>%
  group_by(metric, time_point) %>%
  rstatix::pairwise_wilcox_test(unique_mutations ~ relapse_timing, p.adjust.method = "BH") %>%
  rstatix::add_y_position(fun = "max", scales = "free_y") %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "****", "***", p.adj.signif)) %>%
  ungroup()

write_tsv(
  summary_stats_counts,
  "data/shared_mutations/mutation_tally_by_group_wilcox.tsv"
)

count_unique_boxplot <- tabulate_shared %>%
  ggplot(aes(x = relapse_timing, y = unique_mutations)) +
  geom_boxplot(aes(fill = relapse_timing), outlier.shape = NA) +
  geom_quasirandom() +
  scale_fill_manual(values = relapse_colours) +
  scale_x_discrete(labels = relapse_timing_labels) +
  facet_grid(metric ~ time_point, scales = "free_y") +
  stat_pvalue_manual(summary_stats_counts, label = "p.adj.signif") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Unique Variants")

count_unique_boxplot

ggsave("figures/count_unique_boxplot.pdf",
  count_unique_boxplot,
  height = 6, width = 6
)
saveRDS(
  count_unique_boxplot,
  "figures/count_unique_boxplot.RDS"
)

summary_stats_percent <- tabulate_shared %>%
  group_by(metric, time_point) %>%
  rstatix::pairwise_wilcox_test(unique_muts_percent ~ relapse_timing, p.adjust.method = "BH") %>%
  rstatix::add_y_position(fun = "max") %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "****", "***", p.adj.signif)) %>%
  ungroup()

write_tsv(
  summary_stats_percent,
  "data/shared_mutations/mutation_percent_by_group_wilcox.tsv"
)

percent_unique_boxplot <- tabulate_shared %>%
  ggplot(aes(x = relapse_timing, y = unique_muts_percent)) +
  geom_boxplot(aes(fill = relapse_timing), outlier.shape = NA) +
  geom_quasirandom() +
  scale_fill_manual(values = relapse_colours) +
  scale_x_discrete(labels = relapse_timing_labels) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  facet_grid(metric ~ time_point, scales = "free_y") +
  stat_pvalue_manual(summary_stats_percent, label = "p.adj.signif", tip.length = 0) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Unique Variants (%)")

percent_unique_boxplot

ggsave("figures/percent_unique_boxplot.pdf",
  percent_unique_boxplot,
  height = 6, width = 6
)
saveRDS(
  percent_unique_boxplot,
  "figures/percent_unique_boxplot.RDS"
)

metric_labs <- c("All Variants", "Coding Only")
names(metric_labs) <- c("All", "Coding")

percent_unique_boxplot_wide <- percent_unique_boxplot +
  facet_grid(~ metric + time_point,
    labeller = labeller(
      metric = metric_labs
    )
  )


percent_unique_boxplot_wide

ggsave("figures/percent_unique_boxplot_wide.pdf",
  percent_unique_boxplot_wide,
  height = 4, width = 12
)
saveRDS(
  percent_unique_boxplot_wide,
  "figures/percent_unique_boxplot_wide.RDS"
)

shading_df <- data.frame(
  xmin = c(0, 0.75, 2),
  xmax = c(0.75, 2, max(trios_md$time_since_diagnosis_years) + 1),
  ymin = 0,
  ymax = 105,
  relapse_timing = names(relapse_colours)
)

shading_df2 <- bind_rows(
  shading_df,
  shading_df
) %>%
  mutate(metric = c(rep.int("All", 3), rep.int("Coding", 3))) %>%
  mutate(ymax = c(rep.int(25000, 3), rep.int(300, 3)))

count_unique_linear <- tabulate_shared %>%
  left_join(select(
    trios_md,
    patient_id,
    time_point,
    time_since_diagnosis_years
  )) %>%
  group_by(patient_id) %>%
  mutate(time_since_diagnosis_years = max(time_since_diagnosis_years)) %>%
  ungroup() %>%
  ggplot(aes(
    x = time_since_diagnosis_years,
    y = total.muts - shared.muts
  )) +
  geom_rect(
    data = shading_df2,
    aes(
      xmin = xmin,
      ymin = ymin,
      xmax = xmax,
      ymax = ymax,
      fill = relapse_timing
    ),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black"
  ) +
  facet_grid(metric ~ time_point, scales = "free_y") +
  stat_cor(p.accuracy = 0.001) +
  scale_fill_manual(
    values = relapse_colours,
    name = "Relapse Timing"
  ) +
  xlab("Time Between Biopsies (Years)") +
  ylab("Unique Variants") +
  theme(legend.position = "bottom")

count_unique_linear

ggsave("figures/count_unique_linear.pdf",
  count_unique_linear,
  height = 8, width = 8
)
saveRDS(
  count_unique_linear,
  "figures/count_unique_linear.RDS"
)

count_unique_linear_coverage <- tabulate_shared %>%
  filter(seq_type == "genome") %>%
  ggplot(aes(
    x = genome_coverage,
    y = total.muts - shared.muts
  )) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black"
  ) +
  facet_grid(metric ~ time_point, scales = "free") +
  stat_cor(p.accuracy = 0.001) +
  xlab("Mean Corrected Genome Coverage") +
  ylab("Unique Variants") +
  theme(legend.position = "bottom")

count_unique_linear_coverage

ggsave("figures/count_unique_linear_coverage.pdf",
  count_unique_linear_coverage,
  height = 8, width = 8
)
saveRDS(
  count_unique_linear_coverage,
  "figures/count_unique_linear_coverage.RDS"
)

percent_unique_linear <- tabulate_shared %>%
  left_join(select(
    trios_md,
    patient_id,
    time_point,
    time_since_diagnosis_years
  )) %>%
  group_by(patient_id) %>%
  mutate(time_since_diagnosis_years = max(time_since_diagnosis_years)) %>%
  ungroup() %>%
  ggplot(aes(
    x = time_since_diagnosis_years,
    y = unique_muts_percent
  )) +
  geom_rect(
    data = shading_df,
    aes(
      xmin = xmin,
      ymin = ymin,
      xmax = xmax,
      ymax = ymax,
      fill = relapse_timing
    ),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black"
  ) +
  facet_grid(~ time_point + metric, scales = "free_y") +
  stat_cor(label.x = 5, label.y = 15, p.accuracy = 0.001) +
  scale_fill_manual(
    values = relapse_colours,
    name = "Relapse Timing"
  ) +
  xlab("Time Between Biopsies (Years)") +
  ylab("Unique Variants (%)") +
  theme(legend.position = "right") +
  ylim(0, 105)

percent_unique_linear

ggsave("figures/percent_unique_linear.pdf",
  percent_unique_linear,
  height = 4, width = 14
)
saveRDS(
  percent_unique_linear,
  "figures/percent_unique_linear.RDS"
)

lowgrade <- read_tsv("data/metadata/metadata_lowgrade.tsv") %>%
  select(patient_id, is_denovo)

colours <- get_gambl_colours(classification = "blood")[c("Red", "Blue")]
names(colours) <- c("FALSE", "TRUE")

pct_unique_linear_lowgrade <- tabulate_shared %>%
  left_join(select(
    trios_md,
    patient_id,
    time_point,
    time_since_diagnosis_years
  )) %>%
  left_join(lowgrade) %>%
  group_by(patient_id) %>%
  mutate(time_since_diagnosis_years = max(time_since_diagnosis_years)) %>%
  ungroup() %>%
  ggplot(aes(
    x = time_since_diagnosis_years,
    y = unique_muts_percent,
    colour = is_denovo
  )) +
  # geom_rect(
  #   data = shading_df,
  #   aes(
  #     xmin = xmin,
  #     ymin = ymin,
  #     xmax = xmax,
  #     ymax = ymax,
  #     fill = relapse_timing
  #   ),
  #   inherit.aes = FALSE,
  #   alpha = 0.5
  # ) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE,
    # colour = "black"
  ) +
  facet_grid(~ time_point + metric, scales = "free_y") +
  stat_cor(label.x = 5, label.y = c(10, 20), p.accuracy = 0.001) +
  scale_fill_manual(
    values = relapse_colours,
    name = "Relapse Timing"
  ) +
  xlab("Time Between Biopsies (Years)") +
  ylab("Unique Variants (%)") +
  theme(legend.position = "bottom") +
  scale_colour_manual(name = expression(paste(italic("de novo"), " DLBCL")), values = colours) +
  ylim(0, 105) +
  guides(
    colour = guide_legend(
      override.aes = aes(label = "")
    )
  )
pct_unique_linear_lowgrade

ggsave("figures/percent_unique_boxplot_lowgrade.pdf",
  pct_unique_linear_lowgrade,
  height = 4, width = 12
)
saveRDS(
  pct_unique_linear_lowgrade,
  "figures/percent_unique_boxplot_lowgrade.RDS"
)

lm_lowgrade <- tabulate_shared %>%
  left_join(select(
    trios_md,
    patient_id,
    time_point,
    time_since_diagnosis_years
  )) %>%
  group_by(patient_id) %>%
  mutate(time_since_diagnosis_years = max(time_since_diagnosis_years)) %>%
  ungroup() %>%
  rename(mutation_type = metric) %>%
  select(-total.muts) %>%
  pivot_longer(matches("mut") & where(is.numeric),
    names_to = "metric",
    values_to = "value"
  ) %>%
  left_join(lowgrade) %>%
  group_by(mutation_type, time_point, metric) %>%
  nest() %>%
  mutate(model = map(
    data, ~ broom::tidy(
      lm(value ~ relapse_timing + is_denovo,
        data = .
      )
    )
  )) %>%
  select(-data) %>%
  unnest(c(model)) %>%
  mutate(p.format = p_format(p.value, accuracy = 0.001))

write_tsv(lm_lowgrade, "data/shared_mutations/lm_relapse_timing_lowgrade.tsv")

percent_unique_linear_coverage <- tabulate_shared %>%
  ggplot(aes(
    x = genome_coverage,
    y = unique_muts_percent
  )) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black"
  ) +
  facet_grid(metric ~ time_point, scales = "free_y") +
  stat_cor(p.accuracy = 0.001) +
  xlab("Mean Corrected Genome Coverage") +
  ylab("Unique Variants (%)") +
  theme(legend.position = "bottom")

percent_unique_linear_coverage

ggsave("figures/percent_unique_linear_coverage.pdf",
  percent_unique_linear_coverage,
  height = 8, width = 8
)
saveRDS(
  percent_unique_linear_coverage,
  "figures/percent_unique_linear_coverage.RDS"
)


# Plot shared vs. total mutations

shared_vs_total <- tabulate_shared %>%
  select(patient_id, relapse_timing, time_point, total.muts, shared.muts, mutation_type = metric) %>%
  filter(patient_id %in% unique(tabulate_shared[tabulate_shared$total.muts < 50000 & tabulate_shared$metric == "All", ]$patient_id)) %>%
  ggplot(aes(
    x = total.muts,
    y = shared.muts,
    colour = relapse_timing
  )) +
  geom_abline(slope = 1, intercept = 0, colour = "darkgrey", linetype = 2) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  stat_regline_equation(show.legend = FALSE) +
  facet_wrap(time_point ~ mutation_type,
    scales = "free",
    nrow = 1,
    labeller = labeller(
      mutation_type = metric_labs
    )
  ) +
  scale_colour_manual(values = relapse_colours, name = "Relapse Timing") +
  theme(legend.position = "bottom") +
  xlab("Total Mutations") +
  ylab("Shared Mutations")

ggsave("figures/shared_vs_total.pdf", height = 4, width = 13)
saveRDS(shared_vs_total, "figures/shared_vs_total.RDS")


unique_vs_total <- tabulate_shared %>%
  select(patient_id, relapse_timing, time_point, total.muts, unique_mutations, mutation_type = metric) %>%
  filter(patient_id %in% unique(tabulate_shared[tabulate_shared$total.muts < 50000 & tabulate_shared$metric == "All", ]$patient_id)) %>%
  ggplot(aes(
    x = total.muts,
    y = unique_mutations,
    colour = relapse_timing
  )) +
  # geom_abline(slope = -1, intercept = 0, colour = "darkgrey", linetype = 2) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  stat_regline_equation(show.legend = FALSE) +
  facet_wrap(time_point ~ mutation_type,
    scales = "free",
    nrow = 1,
    labeller = labeller(
      mutation_type = metric_labs
    )
  ) +
  scale_colour_manual(values = relapse_colours, name = "Relapse Timing") +
  theme(legend.position = "bottom") +
  xlab("Total Mutations") +
  ylab("Unique Mutations")

ggsave("figures/unique_vs_total.pdf", height = 4, width = 13)
saveRDS(unique_vs_total, "figures/unique_vs_total.RDS")


##### Only CN=2 and no LOH #####

tabulate_shared_CN <- vaf_df %>%
  filter(Diagnosis_CN == 2 & Relapse_CN == 2 &
    Diagnosis_LOH_flag == 0 & Relapse_LOH_flag == 0) %>%
  select(-matches("CN|LOH")) %>%
  pivot_longer(
    matches("VAF"),
    names_to = "time_point",
    values_to = "VAF"
  ) %>%
  mutate(time_point = str_remove(time_point, "_VAF")) %>%
  filter(VAF > 0, !is.na(VAF)) %>%
  group_by(patient_id, time_point) %>%
  mutate(
    total.muts_all = n(),
    shared.muts_all = sum(shared, na.rm = TRUE),
    total.muts_coding = sum(coding, na.rm = TRUE),
    shared.muts_coding = sum(shared & coding, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(patient_id, relapse_timing, time_point, matches("muts")) %>%
  distinct() %>%
  pivot_longer(
    cols = matches("muts"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  separate(metric, into = c("tally", "metric"), sep = "_") %>%
  pivot_wider(
    names_from = tally,
    values_from = value
  ) %>%
  mutate(unique_muts_percent = 100 - round(shared.muts / total.muts * 100, digits = 1)) %>%
  mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours))) %>%
  mutate(metric = ifelse(metric == "coding", "Coding", "All"))

summary_means_CN <- tabulate_shared_CN %>%
  group_by(relapse_timing, metric, time_point) %>%
  summarize(mean_percent_unique = round(mean(unique_muts_percent, na.rm = TRUE), digits = 1)) %>%
  select(relapse_timing, metric, time_point, mean_percent_unique)

summary_stats_CN <- tabulate_shared_CN %>%
  group_by(metric, time_point) %>%
  rstatix::pairwise_wilcox_test(unique_muts_percent ~ relapse_timing, p.adjust.method = "BH") %>%
  rstatix::add_y_position(fun = "max", scales = "free_y") %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "****", "***", p.adj.signif)) %>%
  ungroup()

diploid_mutations_boxplot <- tabulate_shared_CN %>%
  ggplot(aes(x = relapse_timing, y = unique_muts_percent)) +
  geom_boxplot(aes(fill = relapse_timing), outlier.shape = NA) +
  geom_quasirandom() +
  scale_fill_manual(values = relapse_colours) +
  scale_x_discrete(labels = relapse_timing_labels) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  facet_grid(metric ~ time_point) +
  stat_pvalue_manual(summary_stats_CN, label = "p.adj.signif") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Unique Variants (%)") +
  ggtitle("Diploid Loci Only")

ggsave("figures/diploid_mutations_boxplot.pdf",
  diploid_mutations_boxplot,
  height = 6, width = 6
)


##### Only SNPs #####

tabulate_shared_SNP <- vaf_df %>%
  filter(Variant_Type == "SNP") %>%
  select(-matches("CN|LOH")) %>%
  pivot_longer(
    matches("VAF"),
    names_to = "time_point",
    values_to = "VAF"
  ) %>%
  mutate(time_point = str_remove(time_point, "_VAF")) %>%
  filter(VAF > 0, !is.na(VAF)) %>%
  group_by(patient_id, time_point) %>%
  mutate(
    total.muts_all = n(),
    shared.muts_all = sum(shared, na.rm = TRUE),
    total.muts_coding = sum(coding, na.rm = TRUE),
    shared.muts_coding = sum(shared & coding, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(patient_id, relapse_timing, time_point, matches("muts")) %>%
  distinct() %>%
  pivot_longer(
    cols = matches("muts"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  separate(metric, into = c("tally", "metric"), sep = "_") %>%
  pivot_wider(
    names_from = tally,
    values_from = value
  ) %>%
  mutate(unique_muts_percent = 100 - round(shared.muts / total.muts * 100, digits = 1)) %>%
  mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours))) %>%
  mutate(metric = ifelse(metric == "coding", "Coding", "All"))

summary_means_SNP <- tabulate_shared_SNP %>%
  group_by(relapse_timing, metric, time_point) %>%
  summarize(mean_percent_unique = round(mean(unique_muts_percent, na.rm = TRUE), digits = 1)) %>%
  select(relapse_timing, metric, time_point, mean_percent_unique)

summary_stats_SNP <- tabulate_shared_SNP %>%
  group_by(metric, time_point) %>%
  rstatix::pairwise_wilcox_test(unique_muts_percent ~ relapse_timing) %>%
  rstatix::add_y_position(fun = "max") %>%
  ungroup()

snp_only_boxplot <- tabulate_shared_SNP %>%
  ggplot(aes(x = relapse_timing, y = unique_muts_percent)) +
  geom_boxplot(aes(fill = relapse_timing), outlier.shape = NA) +
  geom_quasirandom() +
  scale_fill_manual(values = relapse_colours) +
  scale_x_discrete(labels = relapse_timing_labels) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  facet_grid(metric ~ time_point) +
  stat_pvalue_manual(summary_stats_SNP, label = "p.adj.signif") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Unique Variants (%)") +
  ggtitle("SNVs Only")

ggsave("figures/snv_mutations_boxplot.pdf",
  snp_only_boxplot,
  height = 6, width = 6
)


# Ensure the correlation is not with patient age

pt_age <- read_tsv("data/metadata/metadata_clinical_patient.tsv") %>%
  select(patient_id, age) %>%
  distinct()

tabulate_shared_age <- tabulate_shared %>%
  left_join(pt_age)

tabulate_shared_age %>%
  filter(
    patient_id != "LY_RELY_121",
    time_point == "Diagnosis"
  ) %>%
  ggplot(aes(
    x = age,
    y = total.muts
  )) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black"
  ) +
  stat_cor(p.accuracy = 0.001) +
  facet_grid(metric ~ 1, scales = "free") +
  ylab("Total Variants") +
  xlab("Age at Diagnosis") +
  ggtitle("Total Variants vs. \nAge at Diagnosis") +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_blank()
  )

ggsave("figures/total_mutations_vs_age.pdf",
  height = 6, width = 4
)

tabulate_shared_age %>%
  ggplot(aes(
    x = age,
    y = unique_mutations
  )) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black"
  ) +
  stat_cor(p.accuracy = 0.001) +
  facet_grid(metric ~ time_point, scales = "free_y") +
  ylab("Unique Variants") +
  xlab("Age at Diagnosis") +
  ggtitle("Unique Variants vs. Age at Diagnosis")

ggsave("figures/unique_mutations_vs_age.pdf",
  height = 6, width = 6
)
