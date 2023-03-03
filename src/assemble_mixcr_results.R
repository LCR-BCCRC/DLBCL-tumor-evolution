##### Get and tidy MiXCR results #####
##### MUST BE RUN ON A GPHOST #####

source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv") %>%
  filter(
    tissue_status == "tumour",
    !is.na(RNAseq_sample_id)
  ) %>%
  select(-matches("DNAseq")) %>%
  rename_with(~ str_remove(.x, "RNAseq_")) %>%
  group_by(patient_id, biopsy_id) %>%
  slice_max(str_length(sample_id), n = 1, with_ties = FALSE) %>%
  ungroup()

mixcr_dir <- "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/results/gambl/mixcr-1.2/99-outputs/txt/mrna/"

mixcr_trios <- unlist(lapply(
  trios_md[trios_md$tissue_status == "tumour", ]$sample_id,
  function(x) paste0(mixcr_dir, x)
))

colspec <- cols(
  .default = col_character(),
  cloneId = col_double(),
  cloneCount = col_double(),
  cloneFraction = col_double(),
  bestVHitScore = col_double(),
  bestDHitScore = col_double(),
  bestJHitScore = col_double(),
  bestCHitScore = col_double()
)

mixcr_IGH <- lapply(
  trios_md$sample_id,
  function(x) {
    read_tsv(dir(
      mixcr_dir,
      pattern = paste0(
        "mixcr.",
        x,
        ".clonotypes.IGH.igblast.txt"
      ),
      full.names = TRUE
    ), col_types = colspec) %>%
      mutate(sample_id = x)
  }
)

mixcr_IGH <- bind_rows(mixcr_IGH)

mixcr_IGH_filt <- mixcr_IGH %>%
  group_by(sample_id, bestVGene) %>%
  mutate(
    cloneCount = sum(cloneCount),
    cloneFraction = sum(cloneFraction)
  ) %>%
  ungroup() %>%
  filter(
    cloneCount >= 100,
    cloneFraction > 0.1
  ) %>%
  mutate(
    igblastnVAllele = str_remove(str_remove(igblastnVAllele, "[*].*"), ",.*"),
    igblastnVAllele = case_when(
      !str_detect(igblastnVAllele, "IGH") ~ bestVGene,
      TRUE ~ igblastnVAllele
    )
  ) %>%
  select(
    sample_id,
    matches("clone"),
    igblastnVAllele,
    matches("Gene"),
    matches("Score"),
    -matches("Seq")
  ) %>%
  group_by(sample_id) %>%
  slice_max(cloneCount, n = 1, with_ties = FALSE) %>%
  ungroup()

mixcr_IGK <- lapply(
  trios_md$sample_id,
  function(x) {
    read_tsv(dir(
      mixcr_dir,
      pattern = paste0(
        "mixcr.",
        x,
        ".clonotypes.IGK.igblast.txt"
      ),
      full.names = TRUE
    ), col_types = colspec) %>%
      mutate(sample_id = x)
  }
)

mixcr_IGK <- bind_rows(mixcr_IGK)

mixcr_IGK_filt <- mixcr_IGK %>%
  group_by(sample_id, bestVGene) %>%
  mutate(
    cloneCount = sum(cloneCount),
    cloneFraction = sum(cloneFraction)
  ) %>%
  ungroup() %>%
  filter(
    cloneCount >= 100,
    cloneFraction > 0.1
  ) %>%
  mutate(
    igblastnVAllele = str_remove(str_remove(igblastnVAllele, "[*].*"), ",.*"),
    igblastnVAllele = case_when(
      !str_detect(igblastnVAllele, "IGK") ~ bestVGene,
      TRUE ~ igblastnVAllele
    )
  ) %>%
  select(
    sample_id,
    matches("clone"),
    igblastnVAllele,
    matches("Gene"),
    matches("Score"),
    -matches("Seq")
  ) %>%
  group_by(sample_id) %>%
  slice_max(cloneCount, n = 1, with_ties = FALSE) %>%
  ungroup()

mixcr_IGL <- lapply(
  trios_md$sample_id,
  function(x) {
    read_tsv(dir(
      mixcr_dir,
      pattern = paste0(
        "mixcr.",
        x,
        ".clonotypes.IGL.igblast.txt"
      ),
      full.names = TRUE
    ), col_types = colspec) %>%
      mutate(sample_id = x)
  }
)

mixcr_IGL <- bind_rows(mixcr_IGL)

mixcr_IGL_filt <- mixcr_IGL %>%
  group_by(sample_id, bestVGene) %>%
  mutate(
    cloneCount = sum(cloneCount),
    cloneFraction = sum(cloneFraction)
  ) %>%
  ungroup() %>%
  filter(
    cloneCount >= 50,
    cloneFraction > 0.1
  ) %>%
  mutate(
    igblastnVAllele = str_remove(str_remove(igblastnVAllele, "[*].*"), ",.*"),
    igblastnVAllele = case_when(
      !str_detect(igblastnVAllele, "IGL") ~ bestVGene,
      TRUE ~ igblastnVAllele
    )
  ) %>%
  select(
    sample_id,
    matches("clone"),
    igblastnVAllele,
    matches("Gene"),
    matches("Score"),
    -matches("Seq")
  ) %>%
  group_by(sample_id) %>%
  slice_max(cloneCount, n = 1, with_ties = FALSE) %>%
  ungroup()

mixcr_lightchain <- bind_rows(
  mixcr_IGL_filt,
  mixcr_IGK_filt
) %>%
  group_by(sample_id) %>%
  slice_max(cloneCount, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(
    select(trios_md,
      sample_id,
      patient_id,
      time_point,
      time_since_diag = time_since_diagnosis_years,
      relapse_timing
    )
  )


lightchain_alluvial <-
  mixcr_lightchain %>%
  group_by(patient_id) %>%
  filter(n() > 1) %>%
  slice_min(time_since_diag, n = 2, with_ties = FALSE) %>%
  mutate(time_point = case_when(
    time_since_diag == min(time_since_diag) ~ "Diagnosis",
    TRUE ~ "Relapse"
  )) %>%
  mutate(relapse_timing = factor(relapse_timing,
    levels = c(
      "Primary Refractory",
      "Early Relapse",
      "Late Relapse"
    )
  )) %>%
  ungroup() %>%
  group_by(patient_id, igblastnVAllele) %>%
  mutate(discordant = n() < 2) %>%
  ungroup() %>%
  mutate(`V Gene` = str_extract(igblastnVAllele, "IG[KL]V[[:digit:]]"))

write_tsv(lightchain_alluvial, "data/ig_rearrangement/lightchain_alluvial.tsv")


mixcr_heavychain <- mixcr_IGH_filt %>%
  group_by(sample_id) %>%
  slice_max(cloneCount, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(
    select(trios_md,
      sample_id,
      patient_id,
      time_point,
      time_since_diag = time_since_diagnosis_years,
      relapse_timing
    )
  )

heavychain_alluvial <-
  mixcr_heavychain %>%
  group_by(patient_id) %>%
  filter(n() > 1) %>%
  slice_min(time_since_diag, n = 2, with_ties = FALSE) %>%
  mutate(time_point = case_when(
    time_since_diag == min(time_since_diag) ~ "Diagnosis",
    TRUE ~ "Relapse"
  )) %>%
  mutate(relapse_timing = factor(relapse_timing,
    levels = c(
      "Primary Refractory",
      "Early Relapse",
      "Late Relapse"
    )
  )) %>%
  ungroup() %>%
  group_by(patient_id, igblastnVAllele) %>%
  mutate(discordant = n() < 2) %>%
  ungroup() %>%
  mutate(`V Gene` = str_extract(igblastnVAllele, "IGHV[[:digit:]]"))

write_tsv(heavychain_alluvial, "data/ig_rearrangement/heavychain_alluvial.tsv")

alluvial_all <- bind_rows(
  mutate(lightchain_alluvial, locus = "IGK/L"),
  mutate(heavychain_alluvial, locus = "IGH")
)

both_discordant <- alluvial_all %>%
  select(patient_id, locus, discordant) %>%
  distinct() %>%
  pivot_wider(names_from = "locus", values_from = "discordant") %>%
  filter(IGH & `IGK/L`) %>%
  pull(patient_id)

bind_rows(mixcr_heavychain, mixcr_lightchain) %>%
  write_tsv("data/ig_rearrangement/all_mixcr_results.tsv")
