# Explore mixture modelling to identify divergent "relapses"

library(tidyverse)
library(mclust)

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>%
    filter(
        tissue_status == "tumour",
        !is.na(DNAseq_sample_id),
        pathology != "FL",
        time_since_diagnosis_years >= 0
    ) %>%
    group_by(patient_id, DNAseq_type) %>%
    slice_min(dtbx, n = 2) %>%
    filter(n() == 2) %>%
    mutate(time_point = ifelse(
        time_since_diagnosis_years == 0,
        "Diagnosis",
        "Relapse"
    )) %>%
    ungroup() %>%
    filter(!patient_id %in% c(
        "01-16433" # One time point is an extreme outlier for number of mutations
    )) %>%
    rename(seq_type = DNAseq_type)

shared_vars <- read_tsv("data/shared_mutations/unique_mutations_per_patient.tsv")

shared_vars_pt <- shared_vars %>%
    select(-genome_coverage_partner) %>%
    pivot_wider(
        names_from = c(time_point, metric),
        values_from = c(genome_coverage, total.muts:unique_mutations),
        names_glue = "{time_point}_{metric}_{.value}"
    )


shared_vars_scatter <- shared_vars_pt %>%
    ggplot(aes(x = Diagnosis_Coding_total.muts, y = Diagnosis_Coding_shared.muts, colour = relapse_timing)) +
    geom_point()

shared_vars_scatter

mclust_model_fun <- function(data_subset) {
    model_data <- shared_vars_pt %>%
        filter(seq_type == "genome") %>%
        select(
            patient_id,
            relapse_timing,
            matches(data_subset),
            -matches("percent")
        ) %>%
        left_join(select(trios_md, patient_id, time_point, time_since_diagnosis_years)) %>%
        filter(time_since_diagnosis_years > 0) %>%
        drop_na() %>%
        select(patient_id, relapse_timing, time_point, everything()) %>%
        rename_with(~ str_remove(.x, data_subset))

    clustered <- Mclust(model_data[, c(4:8)], 2)

    outdir <- "figures/mixture_models/"

    pdf(paste0(outdir, data_subset, "mclust2.pdf"), height = 10, width = 10)
    plot(clustered, what = "classification", main = "Mclust Classification")
    dev.off()

    model_data$mclust2 <- clustered$classification
    colnames(model_data)[9] <- paste0(data_subset, "mclust2")

    clustered3 <- Mclust(model_data[, c(4:8)], 3)
    pdf(paste0(outdir, data_subset, "mclust3.pdf"), height = 10, width = 10)
    plot(clustered3, what = "classification", main = "Mclust Classification")
    dev.off()

    model_data$mclust3 <- clustered3$classification
    colnames(model_data)[10] <- paste0(data_subset, "mclust3")

    model_data <- select(
        model_data,
        patient_id,
        relapse_timing,
        time_point,
        matches("mclust")
    )

    return(model_data)
}

all_mods <- lapply(
    c("Diagnosis_All_", "Diagnosis_Coding_", "Relapse_All_", "Relapse_Coding_"),
    mclust_model_fun
)

all_mods <- full_join(all_mods[[1]], all_mods[[2]]) %>%
    full_join(all_mods[[3]]) %>%
    full_join(all_mods[[4]])
