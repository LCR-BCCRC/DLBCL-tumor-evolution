# Analyze trunk mutations

source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv")

shared_muts <- read_tsv("data/shared_mutations/unique_mutations_per_patient.tsv")

phyclone_maf <- read_tsv("data/phyclone/phyclone_mafs.tsv")

trios_md %>%
  filter(tissue_status == "tumour") %>%
  left_join(shared_muts, by = c("patient_id", "relapse_timing", "relapse" = "time_point")) %>%
  filter(!is.na(unique_muts_percent), metric == "All") %>%
  select(Tumor_Sample_Barcode = DNAseq_sample_id, metric, matches("muts"), relapse_timing) %>%
  ggplot(aes(x = relapse_timing, y = unique_muts_percent)) +
  geom_boxplot() +
  geom_quasirandom() +
  geom_hline(yintercept = 25)

# Based on this plot, use a threshold of 25% unique mutations
# to explore shared clonal mutations

divergent <- trios_md %>%
  filter(tissue_status == "tumour") %>%
  left_join(shared_muts, by = c("patient_id", "relapse_timing", "relapse" = "time_point")) %>%
  group_by(patient_id) %>%
  filter(min(unique_muts_percent) > 25, metric == "All") %>%
  ungroup()

divergent_maf <- phyclone_maf %>%
  filter(Tumor_Sample_Barcode %in% divergent$DNAseq_sample_id) %>%
  filter(Hugo_Symbol %in% unique(c(GAMBLR::lymphoma_genes_comprehensive$Gene, GAMBLR::grch37_ashm_regions$gene))) %>%
  group_by(patient_id, clone_id, Chromosome, Start_Position, Tumor_Seq_Allele2) %>%
  mutate(num_high_ccf = sum(ccf >= 0.1)) %>%
  mutate(shared = case_when(
    num_high_ccf == n() ~ "clonal",
    num_high_ccf > 1 ~ "branch_clonal",
    num_high_ccf == 1 ~ "unique"
  )) %>%
  ungroup() %>%
  select(-num_high_ccf) %>%
  mutate(Hugo_Symbol = ifelse(!is.na(hot_spot), str_c(Hugo_Symbol, "*"), Hugo_Symbol)) %>%
  GAMBLR::tidy_lymphgen(lymphgen_column_out = "LymphGen")

divergent_mutations <- divergent_maf %>%
  select(
    patient_id,
    Tumor_Sample_Barcode,
    Hugo_Symbol,
    clone_id,
    ccf,
    HGVSp_Short,
    vaf,
    source,
    shared,
    Subtype.Prediction,
    LymphGen
  ) %>%
  left_join(select(divergent, Tumor_Sample_Barcode = DNAseq_sample_id, unique_muts_percent, relapse, time_since_diagnosis_years))

write_tsv(divergent_mutations, "data/phyclone/mutations_divergent_patients.tsv")

# lg_feats <- GAMBLR::lymphoma_genes_comprehensive
# lg_feats <- lg_feats %>%
#   filter(LymphGen) %>%
#   pull(Gene) %>%
#   unique()
# lg_feats <- c(lg_feats, paste0(lg_feats, "*"))

lg_feats <- read_tsv("data/lymphgen//lg_features.tsv")

hotspot_feats <- divergent_maf %>%
  filter(str_detect(Hugo_Symbol, "[*]")) %>%
  count(Hugo_Symbol) %>%
  pull(Hugo_Symbol)

lg_feats <- lg_feats %>%
  mutate(Hugo_Symbol = str_c(Hugo_Symbol, "*")) %>%
  bind_rows(lg_feats) %>%
  filter(Hugo_Symbol %in% c(lg_feats$Hugo_Symbol, hotspot_feats))

lg_feat_pairs <- data.frame(t(combn(lg_feats$Hugo_Symbol, 2, simplify = TRUE)))
colnames(lg_feat_pairs) <- c("gene1", "gene2")
lg_feat_pairs <- lg_feat_pairs %>%
  compare_divergent() <- divergent_mutations %>%
  filter(shared == c("unique"), vaf > 0) %>%
  filter(time_since_diagnosis_years >= 0) %>%
  group_by(patient_id, Hugo_Symbol) %>%
  slice_min(time_since_diagnosis_years, n = 2) %>%
  ungroup()

t1_muts <- compare_divergent %>%
  filter(relapse == "Diagnosis", Hugo_Symbol %in% lg_feats) %>%
  select(patient_id, Hugo_Symbol) %>%
  distinct() %>%
  rename(T1 = Hugo_Symbol)

t2_muts <- compare_divergent %>%
  filter(relapse == "Relapse", Hugo_Symbol %in% lg_feats) %>%
  select(patient_id, Hugo_Symbol) %>%
  distinct() %>%
  rename(T2 = Hugo_Symbol)

t1_t2_pairs <- t1_muts %>%
  left_join(t2_muts) %>%
  count(T1, T2) %>%
  arrange(desc(n))

t1_mat <- compare_divergent %>%
  filter(relapse == "Diagnosis", Hugo_Symbol %in% lg_feats) %>%
  select(patient_id, Hugo_Symbol) %>%
  mutate(T1 = 1) %>%
  distinct() %>%
  full_join(distinct(select(divergent, patient_id))) %>%
  group_by(Hugo_Symbol) %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = T1,
    values_fill = 0
  ) %>%
  pivot_longer(
    -patient_id,
    names_to = "T1_mut",
    values_to = "T1_mutated"
  ) %>%
  filter(!is.na(T1_mut))

t2_mat <- compare_divergent %>%
  filter(relapse == "Relapse", Hugo_Symbol %in% lg_feats) %>%
  select(patient_id, Hugo_Symbol) %>%
  mutate(T2 = 1) %>%
  distinct() %>%
  group_by(Hugo_Symbol) %>%
  full_join(distinct(select(divergent, patient_id))) %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = T2,
    values_fill = 0
  ) %>%
  pivot_longer(
    -patient_id,
    names_to = "T2_mut",
    values_to = "T2_mutated"
  ) %>%
  filter(!is.na(T2_mut))

test_mat_t12 <- left_join(t1_mat, t2_mat) %>%
  group_by(T1_mut, T2_mut) %>%
  filter(sum(T1_mutated == 1 & T2_mutated == 1) >= 2) %>%
  ungroup()

test_mat_summary <- test_mat_t12 %>%
  group_by(T1_mut, T2_mut) %>%
  count(T1_mutated, T2_mutated)

fish_test_t12 <- test_mat_t12 %>%
  group_by(T1_mut, T2_mut) %>%
  summarize(table = list(table(T1_mutated, T2_mutated))) %>%
  mutate(test = map(table, fisher.test)) %>%
  mutate(tidy = map(test, broom::tidy)) %>%
  unnest(tidy) %>%
  select(-table, -test) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  filter(p.value < 0.05, q.value < 0.1, estimate > 1) %>%
  arrange(q.value)

clonal_muts <- divergent_mutations %>%
  filter(shared == c("clonal")) %>%
  select(patient_id, Hugo_Symbol) %>%
  distinct() %>%
  rename(clonal = Hugo_Symbol)

unique_muts <- divergent_mutations %>%
  filter(shared == c("unique")) %>%
  select(patient_id, Hugo_Symbol) %>%
  distinct() %>%
  rename(unique = Hugo_Symbol)

unique_mat <- divergent_mutations %>%
  filter(shared == "unique") %>%
  select(patient_id, Hugo_Symbol) %>%
  mutate(unique = 1) %>%
  distinct() %>%
  group_by(Hugo_Symbol) %>%
  filter(n() >= 3) %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = unique,
    values_fill = 0
  ) %>%
  pivot_longer(
    -patient_id,
    names_to = "unique_mut",
    values_to = "unique_mutated"
  )

clonal_mat <- divergent_mutations %>%
  filter(shared == "clonal") %>%
  select(patient_id, Hugo_Symbol) %>%
  mutate(clonal = 1) %>%
  distinct() %>%
  group_by(Hugo_Symbol) %>%
  filter(n() >= 3) %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = clonal,
    values_fill = 0
  ) %>%
  pivot_longer(
    -patient_id,
    names_to = "clonal_mut",
    values_to = "clonal_mutated"
  )

test_mat <- left_join(unique_mat, clonal_mat) %>%
  group_by(unique_mut, clonal_mut) %>%
  filter(sum(unique_mutated == 1 & clonal_mutated == 1) >= 3) %>%
  ungroup() %>%
  group_by(clonal_mut, unique_mut) %>%
  summarize(table = list(table(clonal_mutated, unique_mutated))) %>%
  mutate(test = map(table, fisher.test)) %>%
  mutate(tidy = map(test, broom::tidy)) %>%
  unnest(tidy) %>%
  select(-table, -test)


# With these mutation pairs, generate a matrix

tabulate_trunk <- divergent_maf %>%
  select(
    patient_id,
    Hugo_Symbol,
    shared,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    Start_Position,
    End_Position
  ) %>%
  distinct() %>%
  drop_na(shared) %>%
  group_by(Hugo_Symbol) %>%
  count(shared) %>%
  mutate(
    total_mutations = sum(n),
    percent = round(n / sum(n) * 100, digits = 1)
  ) %>%
  ungroup()

tabulate_trunk %>%
  pivot_wider(
    names_from = shared,
    values_from = c(n, percent),
    names_glue = "{shared}_{.value}",
    values_fill = 0
  ) %>%
  write_tsv("data/phyclone/tabulate_trunk.tsv")

lg_feats <- read_tsv("data/lymphgen//lg_features.tsv")
lg_feats <- lg_feats %>%
  mutate(Hugo_Symbol = str_c(Hugo_Symbol, "*")) %>%
  bind_rows(lg_feats)

lg_feat_pairs <- mutate(lg_feats, lymphgen = 1) %>%
  left_join(mutate(lg_feats, lymphgen = 1), by = "lymphgen") %>%
  mutate(pair_number = row_number()) %>%
  select(matches("Hugo_Symbol"), pair_number) %>%
  pivot_longer(
    matches("Hugo_Symbol"),
    names_to = "which_item",
    values_to = "Hugo_Symbol"
  ) %>%
  select(-which_item)

lg_colours <- GAMBLR::get_gambl_colours("lymphgen")
lg_colours <- lg_colours[c("EZB", "ST2", "BN2", "MCD")]

hotspot_feats <- divergent_maf %>%
  filter(str_detect(Hugo_Symbol, "[*]")) %>%
  count(Hugo_Symbol) %>%
  pull(Hugo_Symbol)

lg_feat_clonal <- lg_feats %>%
  filter(!Hugo_Symbol %in% c("EZH2", "NOTCH2")) %>%
  select(-feature) %>%
  filter(!Hugo_Symbol %in% c("BCL6", "BCL2")) %>%
  filter(featureClass != "N1") %>%
  mutate(featureClass = factor(
    featureClass,
    levels = names(lg_colours)
  )) %>%
  left_join(tabulate_trunk) %>%
  filter(!is.na(shared)) %>%
  distinct() %>%
  mutate(shared = factor(
    shared,
    levels = c(
      "unique",
      "branch_clonal",
      "clonal"
    ),
    labels = c(
      "Unique",
      "Branch Clonal",
      "Clonal"
    )
  )) %>%
  ggplot(aes(
    x = Hugo_Symbol,
    y = percent,
    fill = featureClass,
    alpha = shared
  )) +
  geom_col() +
  geom_text(aes(
    x = Hugo_Symbol,
    y = 90,
    label = total_mutations
  ),
  inherit.aes = FALSE,
  show.legend = FALSE
  ) +
  scale_fill_manual(
    values = lg_colours
  ) +
  scale_alpha_manual(
    values = c(
      "Clonal" = 1,
      "Branch Clonal" = 0.7,
      "Unique" = 0.4
    ), name = ""
  ) +
  facet_wrap(~featureClass, scales = "free_y", nrow = 1) +
  coord_flip() +
  theme(legend.position = "bottom") +
  guides(fill = "none", label = "none") +
  xlab("LymphGen Feature") +
  ylab("Percent Clonal Variants")

ggsave("figures/clonal_features.pdf", height = 5, width = 9)
saveRDS(lg_feat_clonal, "figures/clonal_features.RDS")

tabulate_constrained <- divergent_maf %>%
  filter(shared == "unique") %>%
  filter(vaf > 0.1) %>%
  select(patient_id, Hugo_Symbol, Tumor_Sample_Barcode, LymphGen) %>%
  distinct() %>%
  group_by(patient_id, Hugo_Symbol) %>%
  mutate(tumours_mutated = n()) %>%
  filter(tumours_mutated > 1) %>%
  ungroup() %>%
  group_by(patient_id) %>%
  mutate(LymphGen = paste0(unique(LymphGen), collapse = "--")) %>%
  ungroup() %>%
  select(-Tumor_Sample_Barcode) %>%
  distinct() %>%
  group_by(Hugo_Symbol) %>%
  mutate(number_of_patients = n()) %>%
  ungroup()

write_tsv(tabulate_constrained, "data/phyclone/tabulate_constrained.tsv")

lg_feat_constrained <- lg_feats %>%
  mutate(Hugo_Symbol = str_c(Hugo_Symbol, "*")) %>%
  bind_rows(lg_feats) %>%
  filter(Hugo_Symbol %in% c(lg_feats$Hugo_Symbol, hotspot_feats)) %>%
  filter(!Hugo_Symbol %in% c("EZH2", "NOTCH2", "BCL2")) %>%
  select(-feature) %>%
  filter(featureClass != "N1") %>%
  mutate(featureClass = factor(
    featureClass,
    levels = names(lg_colours)
  )) %>%
  left_join(tabulate_constrained) %>%
  filter(!is.na(number_of_patients)) %>%
  mutate(class_informing = case_when(
    str_detect(LymphGen, as.character(featureClass)) ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  distinct() %>%
  group_by(Hugo_Symbol, class_informing, featureClass) %>%
  mutate(num_pts_informing = n()) %>%
  ggplot(aes(
    x = Hugo_Symbol,
    y = as.integer(num_pts_informing),
    fill = featureClass,
    alpha = class_informing
  )) +
  geom_col() +
  # geom_text(aes(
  #   x = Hugo_Symbol,
  #   y = 90,
  #   label = total_mutations
  # ),
  # inherit.aes = FALSE,
  # show.legend = FALSE
  # ) +
  scale_fill_manual(
    values = lg_colours
  ) +
  scale_alpha_manual(
    values = c(
      "TRUE" = 1,
      "FALSE" = 0.7
    ), name = "Class Informing"
  ) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  facet_wrap(~featureClass, scales = "free", nrow = 1) +
  coord_flip() +
  theme(legend.position = "bottom") +
  guides(fill = "none", label = "none") +
  xlab("LymphGen Feature") +
  ylab("Number of Patients")

ggsave("figures/constrained_features.pdf", height = 5, width = 9)
saveRDS(lg_feat_constrained, "figures/constrained_features.RDS")

tabulate_constrained_pt <- divergent_maf %>%
  filter(shared == "unique") %>%
  filter(vaf > 0.1) %>%
  select(patient_id, Hugo_Symbol, Tumor_Sample_Barcode, LymphGen) %>%
  distinct() %>%
  group_by(patient_id, Hugo_Symbol) %>%
  mutate(tumours_mutated = n()) %>%
  # filter(tumours_mutated > 1) %>%
  ungroup() %>%
  group_by(patient_id) %>%
  mutate(LymphGen = paste0(unique(LymphGen), collapse = "--")) %>%
  ungroup() %>%
  select(-Tumor_Sample_Barcode) %>%
  distinct() %>%
  group_by(patient_id, LymphGen) %>%
  summarize(num_constrained_genes = sum(tumours_mutated > 1)) %>%
  ungroup()

tabulate_constrained_pt %>%
  left_join(read_tsv("data/lymphgen//lymphgen_summary.tsv")) %>%
  ggplot(aes(x = num_constrained_genes)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~discordant_type, ncol = 1)

write_tsv(tabulate_constrained_pt, "data/phyclone/tabulate_constrained_pt.tsv")
