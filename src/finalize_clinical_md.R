##### Finalize clinical metadata. #####
# Following extensive clean up of redundant biopsies and
# reviewing treatment courses and pathology, this script
# loads and finalizes the clinical metadata.

source("src/libs.R")
source("src/swap_biopsy_id.R")

seq_md <- read_tsv("data/metadata/metadata_sequencing.tsv")

clinical_md_final <- read_excel("../private_data/clinical_md_biopsy_FINAL.xlsx") %>%
  filter(is.na(EXCLUDE)) %>%
  filter(patient_id %in% seq_md$patient_id) %>%
  filter(!patient_id %in% pts_exclude) %>%
  swap_biopsy_id()

# Anonymize dates
clinical_md_final <- clinical_md_final %>%
  mutate(DTDX_anon = as.POSIXct("2000-01-01")) %>%
  mutate(time_diff = as.numeric(as.Date(DTDX_pt) - as.Date(DTDX_anon))) %>%
  mutate(
    DTDX_pt_anon = as.Date(DTDX_pt) - time_diff,
    DTDX_bx_anon = as.Date(DTDX_bx) - time_diff,
    DTBX_anon = as.Date(DTBX) - time_diff,
    DTFUP_anon = as.Date(DTFUP) - time_diff,
    DTBMT_anon = as.Date(DTBMT) - time_diff
  )

# Confirm this anonymization preserved time differences
check_dates <- clinical_md_final %>%
  mutate(
    time_since_diag = DTBX - DTDX_pt,
    time_since_diag_anon = DTBX_anon - DTDX_pt_anon
  ) %>%
  filter(time_since_diag != time_since_diag_anon)

write_tsv(clinical_md_final, "../private_data/anonymizing_dates.tsv")

clinical_md_final <- clinical_md_final %>%
  select(-c(DTDX_pt, DTDX_bx, DTBX, DTFUP, DTBMT, time_diff, DTDX_anon)) %>%
  rename_with(~ str_remove(.x, "_anon"))

lowgrade <- read_tsv("data/metadata/metadata_lowgrade.tsv")

clinical_md_pt <- clinical_md_final %>%
  select(patient_id,
    age,
    sex,
    DTDX = DTDX_pt,
    DIAG,
    DTFUP,
    CODE_OS,
    OS_YEARS,
    CAUSE_DEATH,
    BMT,
    DTBMT
  ) %>%
  # mutate(DTFUP = as.POSIXct(as.Date(as.numeric(DTFUP), origin = "1899-12-30"), format="%Y-%m-%d")) %>%
  mutate(DIAG = case_when(
    DIAG == "DLBC" ~ "DLBCL",
    str_detect(DIAG, "FOLL") ~ "FL",
    DIAG == "COM" ~ "COMFL",
    TRUE ~ DIAG
  )) %>%
  left_join(lowgrade) %>%
  distinct()

write_tsv(clinical_md_pt, "data/metadata/metadata_clinical_patient.tsv")

colnames_biopsy <- c(
  "patient_id",
  "DTBX",
  "biopsy_id",
  "biopsy_site",
  "DTDX",
  "TX",
  "CT/PET",
  "NOTE",
  "sample_id",
  "pathology"
)

write_tsv(
  select(rename(clinical_md_final, DTDX = DTDX_bx), all_of(colnames_biopsy)),
  "data/metadata/metadata_clinical_biopsy.tsv"
)
