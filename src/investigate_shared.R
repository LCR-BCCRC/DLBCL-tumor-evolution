source("src/libs.R")

md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>% 
  select(-matches("RNAseq")) %>% 
  rename_with(~ str_remove(.x, "DNAseq_")) %>% 
  rename(seq_type = type) %>% 
  filter(tissue_status == "tumour")

gambl_md <- get_gambl_metadata(seq_type_filter = c("genome", "capture")) %>% 
  filter(sample_id %in% md$sample_id)
tumourA <- get_ssm_by_sample(this_sample_id = "00-15201_tumorA", 
                             these_samples_metadata = gambl_md, 
                             min_read_support = 0)

tumourA_min <- tumourA %>% 
  select(Hugo_Symbol, 
         Chromosome, 
         Start_Position, 
         End_Position, 
         HGVSp_Short,
         Variant_Classification,
         t_depth, 
         t_alt_count)

tumourB <- get_ssm_by_sample(this_sample_id = "00-15201_tumorB", 
                             these_samples_metadata = gambl_md, 
                             min_read_support = 0) 

tumourB_min <- tumourB %>% 
  select(Hugo_Symbol, 
         Chromosome, 
         Start_Position, 
         End_Position, 
         HGVSp_Short, 
         Variant_Classification,
         t_depth, 
         t_alt_count)

compare <- left_join(tumourA_min, tumourB_min, 
                     by = c(
                       "Hugo_Symbol", 
                       "Chromosome", 
                       "Start_Position", 
                       "End_Position"
                     ), 
                     suffix = c("_A", "_B"))

compare %>% filter(t_alt_count_A > 0 & t_alt_count_B > 0) %>% View()
