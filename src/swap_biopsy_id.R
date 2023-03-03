swap_biopsy_id <- function(df, biopsy_col = "biopsy_id", patient_col = "patient_id") {
    biopsy_key <- "../private_data/0_Internal_Biopsies_local.xlsx"
    cctg_biopsy_key <- "../private_data/CCTG_biopsy_key.tsv"
    biopsy_key <- read_excel(biopsy_key) %>%
        select(Patient_ID, Biopsy_ID, Biopsy_pub) %>%
        bind_rows(mutate(read_tsv(cctg_biopsy_key), Biopsy_pub = as.character(Biopsy_pub)))
    join_key <- c("Patient_ID", "Biopsy_ID")
    names(join_key) <- c(patient_col, biopsy_col)
    df2 <- df %>%
        left_join(biopsy_key, by = join_key) %>%
        mutate(., !!biopsy_col := ifelse(!is.na(Biopsy_pub), Biopsy_pub, get(biopsy_col))) %>%
        select(all_of(colnames(df)))
}
