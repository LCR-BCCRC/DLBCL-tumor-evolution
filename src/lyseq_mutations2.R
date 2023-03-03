source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>%
  filter(
    tissue_status == "tumour",
    !is.na(DNAseq_sample_id),
    pathology != "FL"
  )

relapse_groups <- trios_md %>%
  select(
    Tumor_Sample_Barcode = DNAseq_sample_id,
    relapse_timing
  )


lyseq_maf <- read_tsv("../LySeq_Validation/results/slms_3-1.0_vcf2maf-1.3/level_3/capture--grch37/LySeqST_merged.maf") %>%
  filter(Tumor_Sample_Barcode %in% trios_md$DNAseq_sample_id) %>%
  mutate(t_alt_count = ifelse(t_alt_count > 2, t_alt_count, 0))

lyseq_cov <- read_tsv("../LySeq_Validation/results/QC/04-mergedMetrics/capture--grch37.hs_metrics.txt") %>%
  select(
    Tumor_Sample_Barcode = sampleID,
    lyseq_coverage = MEAN_TARGET_COVERAGE
  )

genome_cov <- read_tsv("data/qc/gambl_qc_combined.tsv") %>%
  filter(sample_id %in% trios_md$DNAseq_sample_id, metric == "MeanCorrectedCoverage") %>%
  select(Tumor_Sample_Barcode = sample_id, genome_coverage = value)

get_battenberg_purity <- function(sample_id) {
  # sample_id <- unique(lyseq_maf$Tumor_Sample_Barcode)[2]
  batt_dir <- "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/results/gambl/battenberg-1.1/99-outputs/txt/genome--grch37/"
  batt_file <- dir(batt_dir, pattern = sample_id, full.names = TRUE)
  if (!str_detect(batt_file[1], "cellularity") | length(batt_file) == 0) {
    return()
  } else {
    batt_file <- read_tsv(batt_file[1])
    purity <- batt_file[1]$cellularity
    df <- data.frame(Tumor_Sample_Barcode = sample_id, purity = purity)
    return(df)
  }
}

all_purities <- lapply(unique(lyseq_maf$Tumor_Sample_Barcode), get_battenberg_purity)
all_purities <- bind_rows(all_purities)

genome_maf <- read_tsv("data/maf/genomes_augmented.grch37.maf")

genome_vafs <- genome_maf %>%
  select(
    Hugo_Symbol,
    HGVSp_Short,
    Variant_Classification,
    Variant_Type,
    Chromosome,
    Start_Position,
    End_Position,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    t_depth_genome = t_depth,
    t_alt_count_genome = t_alt_count,
    t_ref_count_genome = t_ref_count,
    Tumor_Sample_Barcode
  ) %>%
  mutate(
    t_alt_count_genome = ifelse(t_alt_count_genome < 3, 0, t_alt_count_genome),
    VAF_genome = t_alt_count_genome / t_depth_genome,
    HGVSp_Short = ifelse(HGVSp_Short == "", NA, HGVSp_Short)
  )

load_aug_genomes_fn <- function(sample_id, spec_maf) {
  genome_aug_dir <- "../LySeq_Validation/results/slms_3-1.0_vcf2maf-1.3/level_3/genome--grch37/genome_augmented"
  maf_file <- dir(genome_aug_dir, pattern = sample_id, full.names = TRUE)
  maf <- read_tsv(maf_file, col_types = spec_maf)
  return(maf)
}

load_aug_genomes <- lapply(
  unique(lyseq_maf$Tumor_Sample_Barcode),
  load_aug_genomes_fn,
  spec_maf = spec(lyseq_maf)
)

aug_genomes <- bind_rows(load_aug_genomes) %>%
  select(
    Hugo_Symbol,
    HGVSp_Short,
    Variant_Type,
    Variant_Classification,
    Chromosome,
    Start_Position,
    End_Position,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    t_depth_genome_aug = t_depth,
    t_alt_count_genome_aug = t_alt_count,
    t_ref_count_genome_aug = t_ref_count,
    n_depth_genome = n_depth,
    n_ref_count_genome = n_ref_count,
    n_alt_count_genome = n_alt_count,
    Tumor_Sample_Barcode
  ) %>%
  mutate(
    VAF_genome_aug = t_alt_count_genome_aug / t_depth_genome_aug
  )

capspace <- read_tsv("../LySeq_Validation/data/bed/lyseqst_targets.grch37.bed",
  col_names = c("Chromosome", "Start_Position", "End_Position", "Name")
) %>%
  mutate(
    # Start_Postition = Start_Position - 100,
    # End_Position = End_Position + 100,
    in_capspace = TRUE
  ) %>%
  select(-Name)
capspace <- data.table(capspace)
setkey(capspace, Chromosome, Start_Position, End_Position)

lyseq_vs_genomes <- lyseq_maf %>%
  select(
    Hugo_Symbol,
    HGVSp_Short,
    Variant_Type,
    Variant_Classification,
    Chromosome,
    Start_Position,
    End_Position,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    t_depth_lyseq = t_depth,
    t_alt_count_lyseq = t_alt_count,
    t_ref_count_lyseq = t_ref_count,
    n_depth_lyseq = n_depth,
    n_ref_count_lyseq = n_ref_count,
    n_alt_count_lyseq = n_alt_count,
    Tumor_Sample_Barcode
  ) %>%
  mutate(VAF_lyseq = t_alt_count_lyseq / t_depth_lyseq) %>%
  full_join(genome_vafs) %>%
  full_join(aug_genomes) %>%
  left_join(all_purities) %>%
  left_join(genome_cov) %>%
  left_join(lyseq_cov) %>%
  mutate(across(matches("VAF"), ~ replace_na(.x, 0)))

lyseq_vs_genomes <- data.table(lyseq_vs_genomes)
setkey(lyseq_vs_genomes, Chromosome, Start_Position, End_Position)

lyseq_vs_genomes <- foverlaps(
  lyseq_vs_genomes,
  capspace
) %>%
  select(-Start_Position, -End_Position) %>%
  rename_with(~ str_replace(.x, "^i[.]", "")) %>%
  filter(!is.na(Hugo_Symbol)) %>%
  select(Hugo_Symbol, HGVSp_Short, Chromosome, Start_Position, End_Position, everything()) %>%
  mutate(t_alt_count_genome_theoretical = round(VAF_lyseq * t_depth_genome_aug)) %>%
  mutate(genome_LOD = t_alt_count_genome_theoretical > 3) %>%
  mutate(
    genome_unaugmented = VAF_genome > 0,
    genome_augmented = VAF_genome_aug > 0,
    in_capspace = replace_na(in_capspace, FALSE)
  ) %>%
  mutate(lyseq = case_when(
    Tumor_Sample_Barcode %in% unique(lyseq_maf$Tumor_Sample_Barcode) & VAF_lyseq > 0 ~ TRUE,
    Tumor_Sample_Barcode %in% unique(lyseq_maf$Tumor_Sample_Barcode) & VAF_lyseq == 0 ~ FALSE
  ))

lyseq_vs_genomes %>%
  filter(in_capspace) %>%
  count(
    lyseq,
    genome_unaugmented,
    genome_augmented
  )
lyseq_vs_genomes %>%
  filter(in_capspace) %>%
  filter(!lyseq, genome_unaugmented)



# igv_bed <- lyseq_vs_genomes %>%
#   # mutate(name = str_c(
#   #   Tumor_Sample_Barcode,
#   #   Tumor_Seq_Allele1,
#   #   Tumor_Seq_Allele2,
#   #   ifelse(!is.na(HGVSp_Short), str_remove(HGVSp_Short, "p[.]"), "synon"),
#   #   sep = "_"
#   # )) %>%
#   mutate(
#     name = Tumor_Sample_Barcode,
#     score = 1000,
#     strand = ".",
#     Start_Position = Start_Position - 1,
#     thickStart = Start_Position - 1,
#     thickEnd = End_Position
#   ) %>%
#   mutate(itemRGB = case_when(
#     genome_LOD & !genome_unaugmented & genome_augmented ~ "178,34,34", # firebrick
#     genome_LOD & genome_unaugmented & genome_augmented ~  "0,139,0", # green
#     TRUE ~ "255, 215, 0" # gold
#   )) %>%
#   select(
#     chrom = Chromosome,
#     chromStart = Start_Position,
#     chromEnd = End_Position,
#     name,
#     score,
#     strand,
#     thickStart,
#     thickEnd,
#     itemRGB
#   )
#
# write_tsv(igv_bed, "../LySeq_Validation/check_in_igv.bed", col_names = FALSE)\
#
#
sensitivity <- lyseq_vs_genomes %>%
  filter(VAF_lyseq >= 0.4 * purity) %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(
    in_genome = sum(genome_unaugmented),
    in_lyseq = n(),
    sensitivity = in_genome / in_lyseq
  ) %>%
  select(
    Tumor_Sample_Barcode,
    genome_coverage,
    in_lyseq,
    sensitivity,
    purity
  ) %>%
  distinct() %>%
  filter(in_lyseq >= 10) %>%
  left_join(relapse_groups) %>%
  mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours)))

sens_scatter <- sensitivity %>%
  left_join(select(trios_md, Tumor_Sample_Barcode = DNAseq_sample_id, DNAseq_preservation)) %>%
  filter(genome_coverage < 80) %>%
  ggplot(aes(x = genome_coverage, y = sensitivity)) +
  geom_point(aes(colour = relapse_timing, shape = DNAseq_preservation)) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  scale_colour_manual(
    values = relapse_colours,
    name = "Relapse timing"
  ) +
  stat_cor(label.x = 30, label.y = 0.4) +
  ylab("Sensitivity") +
  xlab("Genome Coverage") +
  ggtitle("Sensitivity for capture-detected\nclonal mutations in genome")

ggsave(
  "figures/genome_sensitivity_lyseq.pdf",
  sens_scatter,
  height = 5, width = 7
)
saveRDS(
  sens_scatter,
  "figures/genome_sensitivity_lyseq.RDS"
)

sense_lm <- broom::tidy(
  lm(
    sensitivity ~ genome_coverage + relapse_timing,
    data = filter(sensitivity, genome_coverage < 75)
  )
) %>%
  mutate(p.format = p_format(p.value))

write_tsv(sense_lm, "data/lyseq/genome_sensitivity_lyseq_lm.tsv")

by_gene <- lyseq_vs_genomes %>%
  filter(
    in_capspace,
    Variant_Classification %in% GAMBLR:::coding_class,
    lyseq | genome_unaugmented
  ) %>%
  group_by(Hugo_Symbol) %>%
  mutate(total_muts_this_gene = n()) %>%
  ungroup() %>%
  group_by(Hugo_Symbol, lyseq, genome_unaugmented) %>%
  mutate(
    num_muts = n(),
    percent = round(num_muts / total_muts_this_gene * 100, digits = 1),
    mean_VAF_lyseq = mean(VAF_lyseq, na.rm = TRUE),
    mean_VAF_genome = mean(VAF_genome, na.rm = TRUE)
  ) %>%
  select(
    Hugo_Symbol,
    lyseq,
    genome_unaugmented,
    percent,
    num_muts,
    total_muts_this_gene,
    mean_VAF_lyseq,
    mean_VAF_genome
  ) %>%
  distinct()

scatter <- lyseq_vs_genomes %>%
  filter(lyseq | genome_unaugmented) %>%
  filter(in_capspace) %>%
  mutate(group = case_when(
    lyseq & genome_unaugmented ~ "Both",
    lyseq ~ "LySeq Only",
    genome_unaugmented ~ "Genome Only"
  )) %>%
  ggplot(aes(
    x = VAF_lyseq,
    y = VAF_genome,
    colour = group
  )) +
  geom_point(alpha = 0.2) +
  xlim(-0.1, 0.9) +
  ylim(-0.1, 0.9) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  xlab("VAF LySeq") +
  ylab("VAF Genome")

ggExtra::ggMarginal(scatter,
  type = "boxplot",
  groupColour = TRUE,
  groupFill = TRUE
)

ggsave(
  "figures/genome_sensitivity_lyseq_marginal.pdf",
  sens_scatter,
  height = 5, width = 7
)
saveRDS(
  sens_scatter,
  "figures/genome_sensitivity_lyseq_marginal.RDS"
)

lyseq_vaf_stats <- lyseq_vs_genomes %>%
  filter(lyseq) %>%
  filter(in_capspace) %>%
  wilcox_test(VAF_lyseq ~ genome_unaugmented) %>%
  add_y_position() %>%
  mutate(p.format = p_format(p))

lyseq_vaf_boxplot <- lyseq_vs_genomes %>%
  filter(lyseq) %>%
  filter(in_capspace) %>%
  ggplot(aes(x = genome_unaugmented, y = VAF_lyseq)) +
  geom_violin(aes(fill = genome_unaugmented)) +
  geom_boxplot(fill = "white", width = 0.2, outlier.colour = NA) +
  stat_pvalue_manual(lyseq_vaf_stats, label = "p.format") +
  ylab("LySeq VAF") +
  xlab("Detected in Genome") +
  theme(legend.position = "none") +
  geom_text(
    data = data.frame(
      x = c(0.75, 1.75),
      y = c(0.4, 0.6),
      label = paste0("n=", c(lyseq_vaf_stats$n1, lyseq_vaf_stats$n2))
    ),
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE
  )

ggsave("figures/lyseq_vaf_violin_sensitivity.pdf",
  lyseq_vaf_boxplot,
  height = 5, width = 5
)

saveRDS(
  lyseq_vaf_boxplot,
  "figures/lyseq_vaf_violin_sensitivity.RDS"
)

# summarize_sensitivity <- lyseq_vs_genomes %>%
#   mutate(detectable = genome_augmented & genome_LOD,
#          missing = !genome_unaugmented & genome_augmented & genome_LOD) %>%
#   group_by(Tumor_Sample_Barcode) %>%
#   mutate(
#     total_detectable = sum(detectable),
#     total_missing = sum(missing),
#     total_lyseq = n(),
#   ) %>%
#   filter(detectable) %>%
#   group_by(Tumor_Sample_Barcode, missing) %>%
#   mutate(mean_depth = mean(t_depth_genome_aug),
#          mean_vaf = mean(VAF_lyseq)) %>%
#   ungroup() %>%
#   # mutate(sensitivity = 1-(total_missing / total_detectable)) %>%
#   mutate(sensitivity = 1-(total_missing / total_lyseq)) %>%
#   select(Tumor_Sample_Barcode, matches("detectable|missing|sensitivity|mean|coverage|total"), purity) %>%
#   distinct()
#
# scatter <- lyseq_vs_genomes %>%
#   filter(Variant_Classification %in% GAMBLR:::coding_class) %>%
#   mutate(`Variant Status Genome` = case_when(
#     !genome_unaugmented & genome_augmented & genome_LOD ~ "detectable missing",
#     genome_unaugmented & genome_augmented & genome_LOD ~ "detectable",
#     TRUE ~ "below LOD"
#   )) %>%
#   ggplot(aes(x = VAF_lyseq, y = VAF_genome, colour = `Variant Status Genome`)) +
#   geom_point() +
#   theme(legend.position = "bottom")
#
# ggExtra::ggMarginal(scatter, type = "density", groupColour = TRUE, groupFill = TRUE)
#
# lyseq_vs_genomes %>%
#   filter(Variant_Classification %in% GAMBLR:::coding_class) %>%
#   mutate(detectable = genome_augmented & genome_LOD,
#          detectable_missing = !genome_unaugmented & genome_augmented & genome_LOD) %>%
#   group_by(Tumor_Sample_Barcode) %>%
#   summarize(
#     detectable_muts = sum(detectable),
#     detectable_missing = sum(detectable_missing)
#   ) %>%
#   mutate(sensitivity = 1-(detectable_missing / detectable_muts)) %>% View()
#
# missing_coding <- lyseq_vs_genomes %>%
#   filter(Variant_Classification %in% GAMBLR:::coding_class) %>%
#   filter(!genome_unaugmented & genome_augmented & genome_LOD)
#
# lyseq_vs_genomes %>%
#   mutate(detectable = genome_augmented & genome_LOD,
#          missing = !genome_unaugmented & genome_augmented & genome_LOD) %>%
#   filter(detectable) %>%
#   ggplot(aes(x = missing, y = t_depth_genome_aug)) +
#   geom_boxplot() +
#   geom_quasirandom() +
#   ylim(c(0, 200)) +
#   stat_compare_means(method = "wilcox.test")
#
# lyseq_vs_genomes %>%
#   mutate(detectable = genome_augmented & genome_LOD,
#          missing = !genome_unaugmented & genome_augmented & genome_LOD) %>%
#   filter(detectable) %>%
#   ggplot(aes(x = missing, y = VAF_genome_aug)) +
#   geom_boxplot() +
#   geom_quasirandom() +
#   stat_compare_means(method = "wilcox.test")
#
# lyseq_vs_genomes %>%
#   mutate(detectable = genome_augmented & genome_LOD,
#          missing = !genome_unaugmented & genome_augmented & genome_LOD) %>%
#   filter(detectable) %>%
#   ggplot(aes(x = t_depth_genome_aug, y = VAF_genome_aug, colour = missing)) +
#   geom_point() +
#   xlim(0, 200)
#
# lyseq_vs_genomes %>%
#   mutate(detectable = genome_augmented & genome_LOD,
#          missing = !genome_unaugmented & genome_augmented & genome_LOD) %>%
#   filter(detectable & missing) %>%
#   View()
#
