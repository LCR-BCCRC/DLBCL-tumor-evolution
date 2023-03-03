##### Creates a stable metadata table for the DLBCL Trios Manuscript #####

source("src/libs.R")
source("src/swap_biopsy_id.R")

gambl_md <- get_gambl_metadata(
  seq_type_filter = c("genome", "capture"),
  tissue_status_filter = c("tumour", "normal")
) %>%
  swap_biopsy_id()

gambl_normals <- gambl_md %>%
  filter(tissue_status == "normal") %>%
  pull(patient_id)

gambl_md <- gambl_md %>%
  filter(patient_id %in% gambl_normals)

# Patients with too low coverage to include here:
low_coverage <- c(
  "00-23442",
  "01-12047",
  "01-23942",
  "17-23711",
  "76-10146",
  "94-25764",
  "08-10249", # Very low TC and mutation yield at T1
  "08-26563", # Very low TC and mutation yield at T2
  "15-21665", # Very low TC and mutation yield at T2
  "18-14625", # Very low TC and mutation yield at T2
  "09-32452" # 0 mutations at T2
)

pts_exclude <- c(pts_exclude, low_coverage)

trios_md <- gambl_md %>%
  filter(patient_id %in%
    pull(filter(gambl_md, cohort == "DLBCL_LSARP_Trios"), patient_id)) %>%
  filter(!patient_id %in% pts_exclude) %>%
  filter(!sample_id == "14-27873T") # Redundant with exome for this biopsy


multi_dlbcl <- trios_md %>%
  arrange(patient_id, time_since_diagnosis_years) %>%
  group_by(patient_id, seq_type) %>%
  mutate(num_tumours = sum(tissue_status == "tumour")) %>%
  mutate(num_normals = sum(tissue_status == "normal")) %>%
  mutate(num_dlbcl = sum(pathology %in% c("DLBCL", "COMFL", "HGBL", "BLL"))) %>%
  filter(
    num_tumours > 1,
    (seq_type == "genome" & num_normals == 1) | seq_type == "capture"
  ) %>%
  filter(tissue_status == "tumour") %>%
  mutate(path_list = paste0(pathology, collapse = ",")) %>%
  mutate(times = paste0(round(time_since_diagnosis_years, digits = 1), collapse = ",")) %>%
  mutate(time_points = paste0(time_point, collapse = ",")) %>%
  slice_max(time_since_diagnosis_years) %>%
  ungroup() %>%
  select(patient_id, seq_type, num_dlbcl, num_tumours, num_normals, relapse_timing, time_points, path_list, times) %>%
  filter(num_dlbcl > 1) %>%
  distinct()

write_tsv(multi_dlbcl, "data/metadata/cohort_summary.tsv")

count(multi_dlbcl, relapse_timing, seq_type)


trios_md_dlbcl <- trios_md %>%
  filter(patient_id %in% multi_dlbcl$patient_id) %>%
  # mutate(time_since_diagnosis_years =
  #          as.numeric(dtbx - min(dtbx, na.rm = TRUE)) / 365) %>%
  select(-relapse_timing, -time_since_diagnosis_years)

index_dlbcl <- trios_md_dlbcl %>%
  filter(tissue_status == "tumour") %>%
  filter(pathology %in% c("DLBCL", "COMFL", "HGBL")) %>%
  group_by(patient_id) %>%
  slice_min(dtbx, with_ties = FALSE, n = 1) %>%
  ungroup() %>%
  select(patient_id, index_dtbx = dtbx)

time_since_diagnosis <- trios_md_dlbcl %>%
  select(patient_id, biopsy_id, dtbx) %>%
  left_join(index_dlbcl) %>%
  mutate(
    time_since_diagnosis_years =
      round(as.numeric(dtbx - index_dtbx) / 365, digits = 2)
  ) %>%
  # Fix time since diagnosis for biopsy with missing dtbx
  mutate(time_since_diagnosis_years = case_when(
    biopsy_id == "Bx_31842" ~ 1.15,
    biopsy_id == "Bx_31840" ~ 0,
    TRUE ~ time_since_diagnosis_years
  ))

relapse_timing <- time_since_diagnosis %>%
  filter(time_since_diagnosis_years > 0) %>%
  group_by(patient_id) %>%
  slice_min(time_since_diagnosis_years, n = 1, with_ties = FALSE) %>%
  mutate(relapse_timing = case_when(
    time_since_diagnosis_years == 0 ~ "Diagnostic",
    time_since_diagnosis_years < (10 / 12) ~ "Primary Refractory",
    time_since_diagnosis_years < 2 ~ "Early Relapse",
    time_since_diagnosis_years >= 2 ~ "Late Relapse"
  )) %>%
  ungroup() %>%
  select(patient_id, relapse_timing) %>%
  distinct()

trios_mrna <- get_gambl_metadata(seq_type_filter = "mrna") %>%
  swap_biopsy_id() %>%
  filter(
    patient_id %in% trios_md_dlbcl$patient_id,
    protocol == "Ribodepletion"
  ) %>%
  select(patient_id,
    biopsy_id,
    RNAseq_sample_id = sample_id,
    RNAseq_protocol = protocol,
    RNAseq_preservation = ffpe_or_frozen,
    RNAseq_source_type = seq_source_type
  )

trios_lst <- read_tsv("../LySeq_Validation/data/metadata/LST_metadata.tsv") %>%
  filter(
    tissue_status == "tumour",
  ) %>%
  select(patient_id, biopsy_id, lyseq_library_id) %>%
  mutate(lyseq_library_id = ifelse(
    lyseq_library_id == "failed",
    NA,
    lyseq_library_id
  )) %>%
  swap_biopsy_id()


trios_md_min <- trios_md_dlbcl %>%
  left_join(time_since_diagnosis) %>%
  left_join(relapse_timing) %>%
  select(
    patient_id,
    biopsy_id,
    pathology,
    tissue_status,
    time_since_diagnosis_years,
    relapse_timing,
    time_point,
    matches("_ba$|_cn$"),
    matches("DLBCL90"),
    DNAseq_sample_id = sample_id,
    DNAseq_type = seq_type,
    DNAseq_preservation = ffpe_or_frozen,
    DNAseq_source_type = seq_source_type
  ) %>%
  left_join(trios_mrna) %>%
  left_join(trios_lst) %>%
  mutate(relapse = case_when(
    time_since_diagnosis_years == 0 ~ "Diagnosis",
    time_since_diagnosis_years > 0 ~ "Relapse",
    time_since_diagnosis_years < 0 ~ "Prior Indolent"
  )) %>%
  distinct()

write_tsv(trios_md_min, "data/metadata/metadata_sequencing.tsv")

# Redo multi_dlbcl to account for updated relapse_timing values from above
multi_dlbcl <- trios_md_min %>%
  arrange(patient_id, time_since_diagnosis_years) %>%
  group_by(patient_id, DNAseq_type) %>%
  mutate(num_tumours = sum(tissue_status == "tumour")) %>%
  mutate(num_normals = sum(tissue_status == "normal")) %>%
  mutate(num_RNAseq = sum(!is.na(RNAseq_sample_id))) %>%
  mutate(num_dlbcl = sum(pathology %in% c("DLBCL", "COMFL", "HGBL", "BLL"))) %>%
  filter(
    num_tumours > 1,
    (DNAseq_type == "genome" & num_normals == 1) | DNAseq_type == "capture"
  ) %>%
  filter(tissue_status == "tumour") %>%
  mutate(path_list = paste0(pathology, collapse = ",")) %>%
  mutate(times = paste0(round(time_since_diagnosis_years, digits = 1), collapse = ",")) %>%
  mutate(time_points = paste0(time_point, collapse = ",")) %>%
  slice_max(time_since_diagnosis_years) %>%
  ungroup() %>%
  select(patient_id, seq_type = DNAseq_type, num_dlbcl, num_tumours, num_normals, num_RNAseq, relapse_timing, time_points, path_list, times) %>%
  filter(num_dlbcl > 1) %>%
  distinct()

write_tsv(multi_dlbcl, "data/metadata/cohort_summary.tsv")


trios_md_min %>%
  filter(tissue_status == "tumour") %>%
  mutate(lyseq_library_id = case_when(
    is.na(lyseq_library_id) ~ "ND",
    lyseq_library_id %in% c("failed") ~ lyseq_library_id,
    TRUE ~ "complete_or_pending"
  )) %>%
  group_by(patient_id) %>%
  mutate(total_tumours = n()) %>%
  count(total_tumours, lyseq_library_id) %>%
  filter(total_tumours == n, lyseq_library_id == "complete_or_pending")
