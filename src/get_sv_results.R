##### Generate SV summary table #####

source("src/libs.R")

# Load stable trios metadata
trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>% 
  filter(tissue_status == "tumour") %>% 
  select(-matches("RNAseq")) %>% 
  rename_with(~str_remove(.x, "DNAseq_")) %>% 
  rename(seq_type = type)

# Load combined manta/gridss SV data from gambl svar_master
svs <- get_combined_sv(sample_ids = trios_md$sample_id, 
                       oncogenes = c("MYC", "BCL2", "BCL6"), 
                       min_vaf = 0.05)

# Load breakpoints identified by custom capture in the LLMPP/CCSRI project
llmpp_capture <- read_tsv("/projects/dscott_prj/CCSRI_1500/capture/results/annotate_sv_vcf/myc_annotated_unique_with_bcl2_bcl6_partners.tsv")

llmpp_trios_overlap_myc <- trios_md %>% 
  dplyr::select(sample_id, biopsy_id, patient_id) %>% 
  left_join(dplyr::select(llmpp_capture, -sample_id), by = c("biopsy_id" = "surgno", "patient_id")) %>%
  filter(str_detect(seq_type, "capture")) %>% 
  rename_with( ~toupper(str_replace(., "_target", "_A")), matches("_target")) %>% 
  rename_with( ~toupper(str_replace(., "_partner", "_B")), matches("_partner")) %>% 
  dplyr::select(any_of(colnames(svs)), 
                SCORE = score, 
                VAF = vaf, 
                DP = dp, 
                tumour_sample_id = sample_id, 
                myc_ba = llmpp_myc_ba) %>% 
  filter(!is.na(CHROM_A))

# Since the capture data are hg38-aligned, I will export them and lift them over
# to grch37 to match the genomes
write_tsv(llmpp_trios_overlap_myc, "data/svs/capture_myc_sv_trios.hg38.bedpe")

# Load the grch37 bedpe file from capture containing only breaks not identified
# in the genomes. These were manually examined in IGV for their veracity. 
llmpp_cap_37 <- read_tsv("data/svs/capture_myc_sv_trios.grch37.bedpe") %>% 
  rename(VAF_tumour = VAF) %>% 
  mutate(across(matches("CHROM"), as.character))

svs <- bind_rows(svs, llmpp_cap_37)

# Add additional partner regions to capture two non-recurrent MYC partners
grch37_partners_addl <- data.frame(
  chrom = c("15", "8"), 
  start = c(72890335, 8237664), 
  end = c(72900336, 8247666), 
  gene = c("ARIH", "PRAG1"), 
  entrez = c(0, 0)
)

grch37_partners <- bind_rows(
  grch37_partners, 
  grch37_partners_addl
)

BCL2_tx <- svs %>% 
  rename(SOMATIC_SCORE = SCORE) %>%
  annotate_sv() %>% 
  filter(str_detect(fusion, "IG[HKL]-BCL2"))

MYC_tx <- svs %>% 
  rename(SOMATIC_SCORE = SCORE) %>% 
  annotate_sv() %>% 
  filter(str_detect(fusion, "^MYC-|-MYC$"), 
         !str_detect(fusion, "^NA-|-NA$"))

BCL6_tx <- svs %>% 
  rename(SOMATIC_SCORE = SCORE) %>% 
  annotate_sv() %>% 
  filter(str_detect(fusion, "^BCL6-|-BCL6$"), 
         !str_detect(fusion, "^NA-|-NA$|FOXP1"))

myc_svs_merged <- trios_md %>%
  filter(tissue_status == "tumour") %>%
  mutate(normal_sample_id = str_c(patient_id, "_normal")) %>% 
  select(patient_id, tumour_sample_id = sample_id, normal_sample_id, myc_ba, seq_type, relapse_timing) %>% 
  full_join(MYC_tx) %>% 
  mutate(FISH = "myc_ba") %>% 
  rename(FISH_result = myc_ba)

bcl2_svs_merged <- trios_md %>%
  filter(tissue_status == "tumour") %>%
  mutate(normal_sample_id = str_c(patient_id, "_normal")) %>% 
  select(patient_id, tumour_sample_id = sample_id, normal_sample_id, bcl2_ba, seq_type, relapse_timing) %>% 
  full_join(BCL2_tx) %>% 
  mutate(FISH = "bcl2_ba") %>% 
  rename(FISH_result = bcl2_ba)

bcl6_svs_merged <- trios_md %>%
  filter(tissue_status == "tumour") %>%
  mutate(normal_sample_id = str_c(patient_id, "_normal")) %>% 
  select(patient_id, tumour_sample_id = sample_id, normal_sample_id, bcl6_ba, seq_type, relapse_timing) %>% 
  full_join(BCL6_tx) %>% 
  mutate(FISH = "bcl6_ba") %>% 
  rename(FISH_result = bcl6_ba)

all_svs <- bind_rows(
  bcl2_svs_merged, 
  bcl6_svs_merged, 
  myc_svs_merged
)

write_tsv(all_svs, "data/svs/all_svs.tsv")

##### Generate some summary values #####

multi_ba <- trios_md %>% 
  filter(tissue_status == "tumour") %>% 
  select(patient_id, matches("_ba$")) %>% 
  pivot_longer(matches("_ba"), 
               names_to = "FISH", 
               values_to = "FISH_result") %>% 
  group_by(patient_id, FISH) %>% 
  mutate(num_FISH_pos = sum(FISH_result == "POS", na.rm = TRUE), 
         num_tumours = n()) %>% 
  ungroup() %>% 
  select(-FISH_result) %>% 
  distinct()


summ_table <- all_svs %>%
  left_join(multi_ba) %>%
  group_by(
    patient_id,
    relapse_timing,
    seq_type,
    num_tumours,
    num_FISH_pos,
    FISH,
    gene,
    partner,
    chrom1,
    start1,
    end1,
    chrom2,
    start2,
    end2,
    strand1,
    strand2
  ) %>%
  filter(!(tumour_sample_id == "18-19313_tumorA" & FISH == "bcl6_ba" & gene == "MYC")) %>% 
  summarize(
    num_tumours_this_breakpoint = sum(!is.na(gene)),
    sample_ids = paste0(tumour_sample_id, collapse = ", "),
    FISH_results = paste0(FISH_result, collapse = ", ")
  ) %>%
  ungroup() %>% 
  group_by(patient_id, FISH, gene, partner) %>% 
  slice_max(num_tumours_this_breakpoint, n=1, with_ties = FALSE) %>% 
  ungroup() %>% 
  group_by(patient_id, FISH) %>% 
  mutate(num_unique_breakpoints = sum(!is.na(partner))) %>% 
  ungroup() %>%
  select(
    patient_id,
    relapse_timing,
    sample_ids,
    seq_type,
    num_tumours,
    num_FISH_pos,
    num_tumours_this_breakpoint,
    num_unique_breakpoints,
    FISH,
    FISH_results,
    everything()
  ) 

write_tsv(summ_table, "data/svs/summarize_all_svs.tsv")
         
summ_table2 <- summ_table %>%
  filter(num_unique_breakpoints > 1 | num_tumours_this_breakpoint > 1) %>%
  filter(num_unique_breakpoints > 0) %>%
  select(patient_id, FISH, num_unique_breakpoints, num_FISH_pos) %>%
  distinct() %>%
  group_by(FISH) %>%
  mutate(multi_FISH_pos = sum(num_FISH_pos > 0)) %>% 
  count(multi_FISH_pos, num_unique_breakpoints) %>%
  ungroup() %>% 
  mutate(discordant = ifelse(num_unique_breakpoints > 1, "Discordant", "Concordant")) %>% 
  select(-num_unique_breakpoints) %>% 
  pivot_wider(names_from = discordant, values_from = n, values_fill = 0)

write_tsv(summ_table2, "data/svs/summarize_multi_FISH_pos.tsv")

summ_table3 <- summ_table %>%
  filter(num_unique_breakpoints > 1 | num_tumours_this_breakpoint > 1) %>%
  filter(num_unique_breakpoints > 0) %>%
  select(patient_id, FISH, num_unique_breakpoints, num_FISH_pos, relapse_timing) %>%
  distinct() %>%
  group_by(FISH, relapse_timing) %>%
  mutate(multi_FISH_pos = sum(num_FISH_pos > 0)) %>%
  count(multi_FISH_pos, num_unique_breakpoints) %>%
  ungroup() %>%
  mutate(discordant = ifelse(num_unique_breakpoints > 1, "Discordant", "Concordant")) %>%
  select(-num_unique_breakpoints) %>%
  pivot_wider(names_from = discordant, values_from = n, values_fill = 0)

write_tsv(summ_table3, "data/svs/summarize_multi_FISH_pos_relapse_timing.tsv")

