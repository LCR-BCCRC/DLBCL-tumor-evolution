source("src/libs.R")

seq_md <- read_tsv("data/metadata/metadata_sequencing.tsv")

seq_md <- seq_md %>%
    filter(tissue_status == "tumour") %>%
    mutate(DLBCL90_dhitsig_call = case_when(
        DLBCL90_dhitsig_call == "POS" ~ "DHITsig+",
        DLBCL90_dhitsig_call == "NEG" ~ "DHITsig-",
        DLBCL90_dhitsig_call == "UNCLASS" ~ "DHITsig-IND",
        TRUE ~ DLBCL90_dhitsig_call
    )) %>%
    rename_with(~ str_replace(.x, "dhitsig", "DZsig"))

fish_ns <- read_tsv("data/metadata/metadata_fish_nanostring.tsv") %>%
    rename_with(~ str_replace(.x, "dhitsig", "DZsig")) %>%
    filter(
        !patient_id %in% pts_exclude
    ) %>%
    select(
        patient_id,
        biopsy_id,
        matches("_call|_score|_ba|_cn"),
        relapse_timing,
        time_since_diagnosis_years,
    ) %>%
    rename_with(~ str_remove(., "_consensus"), everything()) %>%
    rename_with(~ str_c("DLBCL90_", .), matches("_score|_call")) %>%
    select(any_of(colnames(seq_md))) %>%
    mutate(pathology = "DLBCL")

lowgrade <- read_tsv("data/metadata/metadata_lowgrade.tsv")

col_order <- c(
    "patient_id",
    "biopsy_id",
    "pathology",
    "tissue_status",
    "time_since_diagnosis_years",
    "DNAseq_sample_id",
    "DNAseq_type",
    "DNAseq_preservation",
    "DNAseq_source_type",
    "RNAseq_sample_id",
    "RNAseq_protocol",
    "RNAseq_preservation",
    "RNAseq_source_type",
    "lyseq_library_id",
    "relapse",
    "myc_ba",
    "myc_cn",
    "bcl2_ba",
    "bcl2_cn",
    "bcl6_ba",
    "bcl6_cn",
    "DLBCL90_dlbcl_call",
    "DLBCL90_dlbcl_score",
    "DLBCL90_DZsig_call",
    "DLBCL90_DZsig_score",
    "multi_FISH",
    "is_denovo",
    "transformed_from",
    "subsequent_lowgrade",
    "BM_infiltrates",
    "relapse_timing",
    "multi_d90",
    "multi_RNAseq",
    "multi_LySeq"
)
all_md <- seq_md %>%
    full_join(fish_ns, by = c(colnames(seq_md)[colnames(seq_md) %in% colnames(fish_ns)])) %>%
    left_join(lowgrade) %>%
    mutate(across(matches("_score"), ~ round(., digits = 1))) %>%
    group_by(patient_id, biopsy_id) %>%
    slice_max(!is.na(DNAseq_sample_id)) %>%
    ungroup() %>%
    group_by(patient_id) %>%
    mutate(multi_FISH = ifelse(
        (
            sum(!is.na(myc_ba)) > 1 |
                sum(!is.na(bcl2_ba)) > 1 |
                sum(!is.na(bcl6_ba)) > 1
        ),
        TRUE, FALSE
    )) %>%
    mutate(multi_d90 = ifelse(sum(!is.na(DLBCL90_DZsig_call)) > 1, TRUE, FALSE)) %>%
    mutate(multi_RNAseq = ifelse(sum(!is.na(RNAseq_protocol)) > 1, TRUE, FALSE)) %>%
    mutate(multi_LySeq = ifelse(sum(!is.na(lyseq_library_id)) > 1, TRUE, FALSE)) %>%
    ungroup() %>%
    select(all_of(col_order)) %>%
    filter(multi_d90 | multi_FISH | !is.na(DNAseq_sample_id))

summarize_all <- all_md %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    summarize(
        total_patients = n(),
        multi_d90 = sum(multi_d90),
        multi_RNAseq = sum(multi_RNAseq),
        multi_FISH = sum(multi_FISH, na.rm = TRUE),
        multi_LySeq = sum(multi_LySeq, na.rm = TRUE),
        lowgrade_data_available = sum(!is.na(is_denovo)),
        is_denovo = sum(is_denovo, na.rm = TRUE)
    )

write_tsv(summarize_all, "data/metadata/summarize_all_assays.tsv")

all_md %>%
    select(
        -tissue_status,
    ) %>%
    write_tsv("data/metadata/all_assays_results.tsv")
