source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>%
  filter(tissue_status == "tumour",
         DNAseq_type == "genome", 
         !is.na(DNAseq_sample_id), 
         pathology != "FL") %>% 
  group_by(patient_id, DNAseq_type) %>% 
  slice_min(dtbx, n=2) %>% 
  filter(n() == 2) %>% 
  mutate(time_point = ifelse(
    time_since_diagnosis_years == 0, 
    "Diagnosis", 
    "Relapse"
  )) %>% 
  ungroup() %>% 
  filter(!patient_id %in% c(
    "01-16433" # One time point is an extreme outlier for number of mutations
  ))

lyseq_maf <- read_tsv("../LySeq_Validation/results/slms_3-1.0_vcf2maf-1.3/level_3/capture--grch37/LySeqST_merged.maf")

lyseq_md <- trios_md %>%
  filter(DNAseq_sample_id %in% unique(lyseq_maf$Tumor_Sample_Barcode), 
         tissue_status == "tumour",
         DNAseq_type == "genome", 
         !is.na(DNAseq_sample_id), 
         pathology != "FL") %>% 
  group_by(patient_id, DNAseq_type) %>% 
  slice_min(dtbx, n=2) %>% 
  filter(n() == 2) %>% 
  mutate(time_point = ifelse(
    time_since_diagnosis_years == 0, 
    "Diagnosis", 
    "Relapse"
  )) %>% 
  ungroup() %>% 
  filter(!patient_id %in% c(
    "01-16433" # One time point is an extreme outlier for number of mutations
  ))


maf_md <- lyseq_md %>% 
  select(patient_id, 
         Tumor_Sample_Barcode = DNAseq_sample_id, 
         relapse_timing, 
         time_point) %>% 
  left_join(lyseq_maf) %>% 
  filter(!is.na(Start_Position))


vaf_df <- maf_md %>% 
  mutate(VAF = ifelse(
    t_alt_count >= 3, 
    t_alt_count / t_depth, 
    0
  )) %>% 
  select(
    patient_id, 
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
    t_alt_count, 
    t_depth,
    VAF
  ) %>% 
  distinct() %>% 
  pivot_wider(
    names_from = time_point, 
    values_from = c(VAF, t_alt_count, t_depth), 
    names_glue = "{time_point}_{.value}"
  ) %>% 
  mutate(shared = ifelse(
    Diagnosis_VAF > 0 & Relapse_VAF > 0, 
    1, 
    0
  )) %>% 
  filter(!is.na(shared)) %>% 
  mutate(coding = Variant_Classification %in% GAMBLR:::coding_class)

genomes_vaf_df <- read_tsv("../private_data/genome_mutations_shared.maf")

lyseq_vs_genomes_vafs <- vaf_df %>% 
  left_join(
    genomes_vaf_df, 
    by = c(
      "patient_id", 
      "relapse_timing", 
      "Chromosome", 
      "Start_Position", 
      "End_Position", 
      "Reference_Allele", 
      "Tumor_Seq_Allele1", 
      "Tumor_Seq_Allele2", 
      "Hugo_Symbol", 
      "HGVSp_Short", 
      "Variant_Classification", 
      "Variant_Type"
    ), 
  suffix = c("_LST", "_genome"))

tabulate_shared <- vaf_df %>% 
  pivot_longer(
    matches("VAF"), 
    names_to = "time_point", 
    values_to = "VAF"
  ) %>% 
  mutate(time_point = str_remove(time_point, "_VAF")) %>% 
  filter(VAF > 0, !is.na(VAF)) %>% 
  group_by(patient_id, time_point) %>% 
  mutate(total.muts = n(), 
         shared.muts = sum(shared, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(patient_id, relapse_timing, time_point, matches("muts")) %>% 
  distinct() %>% 
  mutate(unique_muts_percent = 100 - round(shared.muts / total.muts * 100, digits = 1)) %>%
  mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours))) %>% 
  mutate(unique_mutations = total.muts - shared.muts) %>% 
  mutate(metric = "LST") %>% 
  filter(total.muts >= 20)

genomes_shared <- read_tsv("data/shared_mutations/unique_mutations_per_patient.tsv")

all_shared <- genomes_shared %>% 
  filter(patient_id %in% tabulate_shared$patient_id) %>% 
  bind_rows(tabulate_shared)

all_shared %>% 
  filter(metric %in% c("Coding", "LST")) %>% 
  select(patient_id, time_point, metric, unique_muts_percent, relapse_timing) %>% 
  pivot_wider(names_from = metric, values_from = unique_muts_percent) %>% 
  ggplot(aes(x = Coding, y = LST)) + 
  geom_point(alpha = 0.5, aes(colour = relapse_timing)) + 
  geom_smooth(method = "lm", se = FALSE, colour = "black") + 
  stat_cor() + 
  geom_abline(slope = 1, intercept = 0) + 
  ggtitle("Unique mutations per tumour (%)")

all_shared %>% 
  filter(metric %in% c("Coding", "LST")) %>% 
  select(patient_id, time_point, metric, unique_muts_percent, relapse_timing) %>% 
  pivot_wider(names_from = metric, values_from = unique_muts_percent) %>%
  filter(LST == 0, Coding > 70)


lyseq_vs_genomes_vafs %>% 
  group_by(patient_id, relapse_timing) %>% 
  mutate(across(matches("coverage"), ~max(.x, na.rm = TRUE))) %>% 
  filter(shared_LST == 1, is.na(shared_genome)) %>%  
  summarize(missing = n(), across(matches("coverage$"), max)) %>% 
  View()  

lyseq_vs_genomes_vafs %>% count(shared_LST, shared_genome)
