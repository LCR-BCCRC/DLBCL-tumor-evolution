source("src/libs.R")

seq_md <- read_tsv("data/metadata/metadata_sequencing.tsv")

fish_ns_md <- read_tsv("data/metadata/metadata_fish_nanostring.tsv")

all_pts <- unique(c(seq_md$patient_id, fish_ns_md$patient_id))

full_path <- read_excel("../private_data/0_Internal_Biopsies_local.xlsx") %>%
  filter(Patient_ID %in% all_pts) %>%
  filter(!PATHa %in% c("REACT", NA, "NSQ", "NEG", "??", "???", "CARCINOMA", "ALH")) %>%
  group_by(Patient_ID) %>%
  mutate(is_diagnostic = ifelse(DTBX == min(DTBX), TRUE, FALSE)) %>%
  ungroup()

diagnostic_entity <- full_path %>%
  filter(is_diagnostic) %>%
  select(Patient_ID, PATHa) %>%
  distinct()

diagnostic_lowgrade <- diagnostic_entity %>%
  filter(!PATHa %in% c("DL?", "DLBC", "IBL", "NS2", "NS3", "HDNSLP") |
    Patient_ID == "07-40609") %>%
  rename(
    transformed_from = PATHa,
    patient_id = Patient_ID
  ) %>%
  # Include FL for 07-40609 which had FL in 2007
  mutate(transformed_from = ifelse(patient_id == "07-40609", "FL", transformed_from)) %>%
  group_by(patient_id) %>%
  slice_max(!transformed_from %in% c("NS1", "COM"), n = 1, with_ties = FALSE) %>%
  ungroup()

diagnostic_denovo <- diagnostic_entity %>%
  filter(PATHa %in% c("DL?", "DLBC", "IBL", "NS2", "NS3")) %>%
  # Exclude 07-40609  from de novo because it had FL in 2007
  filter(Patient_ID != "07-40609") %>%
  pull(Patient_ID)

diagnostic_denovo <- c(
  diagnostic_denovo,
  unique(pull(filter(fish_ns_md, str_detect(patient_id, "^CA")), patient_id)),
)

relapse_lowgrade <- full_path %>%
  filter(Patient_ID %in% diagnostic_denovo) %>%
  group_by(Patient_ID) %>%
  filter(DTBX != min(DTBX)) %>%
  slice_max(!PATHa %in% c("NS1", "COM")) %>%
  ungroup() %>%
  filter(str_detect(PATHa, "COM|FOLL|MALT|MZL|NS1")) %>%
  select(
    patient_id = Patient_ID,
    subsequent_lowgrade = PATHa
  ) %>%
  distinct()

lowgrade_raw <- read_tsv("../private_data/lowgrade_summary.tsv")
lowgrade <- lowgrade_raw %>%
  left_join(diagnostic_lowgrade) %>%
  left_join(relapse_lowgrade) %>%
  filter(patient_id %in% all_pts) %>%
  mutate(is_denovo = ifelse(patient_id %in% diagnostic_denovo, TRUE, FALSE)) %>%
  mutate(BM_infiltrates = case_when(
    !str_detect(comment, "Possible|MGUS") & str_detect(comment, "BM at A|BM at diag") ~ "at_diagnosis",
    !str_detect(comment, "Possible|MGUS") & str_detect(comment, "BM at relapse") ~ "at_relapse"
  )) %>%
  mutate(transformed_from = case_when(
    str_detect(comment, "tMZL") | str_detect(transformed_from, "MZL") ~ "MZL",
    str_detect(comment, "tFL") & !is_denovo ~ "FL",
    str_detect(transformed_from, "FOLL|COMFL") ~ "FL",
    TRUE ~ transformed_from
  )) %>%
  mutate(subsequent_lowgrade = case_when(
    str_detect(subsequent_lowgrade, "MZL") ~ "MZL",
    str_detect(subsequent_lowgrade, "FOLL|COMFL") ~ "FL",
    TRUE ~ subsequent_lowgrade
  )) %>%
  select(patient_id, is_denovo, transformed_from, subsequent_lowgrade, BM_infiltrates) %>%
  distinct() %>%
  group_by(patient_id) %>%
  slice_max(!transformed_from %in% c("NS1", "COM"), n = 1, with_ties = FALSE) %>%
  ungroup()


write_tsv(lowgrade, "data/metadata/metadata_lowgrade.tsv")

# Import table with David's edits and write to final tsv file

lowgrade_manual <- read_excel("../private_data/metadata_lowgrade_DWS.xlsx") %>%
  select(-`DS comment`)

write_tsv(lowgrade_manual, "data/metadata/metadata_lowgrade.tsv", na = "")
