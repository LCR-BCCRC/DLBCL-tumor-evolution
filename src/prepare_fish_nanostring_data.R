##### Prepare FISH and NanoString data for all tumours with and without sequencing #####

source("src/libs.R")
source("src/swap_biopsy_id.R")

# Exclude two biopsies that are either pre-DLBCL or concurrent with diagnostic
excluded_biopsies <- c("Bx_1634", "Bx_13206")
# Exclude two patients with multiple diagnostic biopsies
excluded_patients <- c("07-25994", "07-10924", pts_exclude)

# Load mastertable
trios_all <- read_excel("../private_data/DLBCL_Trio_DS_Mastertable.xlsx") %>%
  swap_biopsy_id(biopsy_col = "Biopsy_ID", patient_col = "Patient_ID") %>%
  filter(
    is.na(EXCLUDED),
    !Biopsy_ID %in% excluded_biopsies,
    !Patient_ID %in% excluded_patients
  ) %>%
  rename_with(tolower) %>%
  rename_with(~ str_replace_all(.x, " ", "_")) %>%
  select(
    patient_id,
    biopsy_id,
    source,
    dtbx
  ) %>%
  group_by(patient_id) %>%
  filter(n() > 1) %>%
  ungroup()

nms <- names(read_excel("../private_data/Assay_Nanostring.xlsx"))
ct <- ifelse(grepl("Duplicate_include", nms), "text", "guess")

# Load and de-duplicate NanoString data
nanostring <- read_excel("../private_data/Assay_Nanostring.xlsx",
  col_types = ct
) %>%
  rename_with(tolower) %>%
  filter(
    str_detect(code_set, "DLBCL90"),
    duplicate_include %in% c("YES", NA, "yes"),
    dlbcl_call != "fail"
  ) %>%
  select(
    patient_id,
    code_set,
    biopsy_id,
    dlbcl_call,
    dlbcl_score,
    dhitsig_call,
    dhitsig_score = dhitisg_score
  ) %>%
  swap_biopsy_id() %>%
  mutate(dhitsig_call = case_when(
    dhitsig_call == "POS" ~ "DHITsig+",
    dhitsig_call == "NEG" ~ "DHITsig-",
    dhitsig_call == "UNCLASS" ~ "DHITsig-IND"
  )) %>%
  group_by(patient_id, biopsy_id) %>%
slice_min(code_set, n = 1, with_ties = FALSE) %>%
  ungroup()

# Load de-duplicated FISH data
fish_all <- read_tsv("../private_data/Assay_FISH_consensus.tsv") %>%
  rename_with(tolower) %>%
  select(-matches("IGH")) %>%
  swap_biopsy_id()

# Combine metadata with FISH and NanoString and calculate relapse timing
trios_fish_ns <- trios_all %>%
  left_join(fish_all) %>%
  left_join(nanostring) %>%
  group_by(patient_id) %>%
  group_by(patient_id) %>%
  mutate(
    time_since_diagnosis_years =
      as.numeric(as.Date(dtbx) - min(as.Date(dtbx), na.rm = TRUE)) / 365
  )

relapse_timing <- trios_fish_ns %>%
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

trios_fish_ns <- trios_fish_ns %>%
  left_join(relapse_timing) %>%
  group_by(patient_id) %>%
  slice_min(dtbx, n = 2) %>%
  filter(n() == 2) %>%
  mutate(time_point = ifelse(dtbx == min(dtbx), "Diagnosis", "Relapse")) %>%
  ungroup() %>%
  select(-dtbx)

write_tsv(trios_fish_ns, "data/metadata/metadata_fish_nanostring.tsv")
