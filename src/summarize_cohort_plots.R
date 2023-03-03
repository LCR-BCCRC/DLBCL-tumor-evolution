# Create UpSet plots to show what analyses were completed on which samples

source("src/libs.R")

all_results <- read_tsv("data/metadata/all_assays_results.tsv")

results_matrix <- all_results %>%
    select(
        patient_id,
        DNAseq_type,
        multi_FISH,
        multi_d90,
        multi_RNAseq,
        multi_LySeq
    ) %>%
    distinct() %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(-patient_id) %>%
    drop_na() %>%
    mutate(name = ifelse(name == "DNAseq_type", value, name)) %>%
    mutate(value = ifelse(value == FALSE, 0, 1)) %>%
    pivot_wider(names_from = name, values_from = value, values_fill = 0) %>%
    rename(Name = patient_id) %>%
    mutate(across(!matches("Name"), as.integer)) %>%
    as.data.frame() %>%
    rowwise() %>%
    filter(sum(c_across(genome:capture)) > 0) %>%
    ungroup() %>%
    select(
        Name,
        `Whole Genome` = genome,
        `Whole Exome` = capture,
        `LySeq Capture` = multi_LySeq,
        `FISH` = multi_FISH,
        NanoString = multi_d90,
        RNAseq = multi_RNAseq
    )

comb_mat <- make_comb_mat(results_matrix, mode = "distinct")
comb_names <- names(comb_size(comb_mat))
has_genome <- comb_names[str_detect(comb_names, "^1")]
has_exome <- comb_names[str_detect(comb_names, "^01")]
has_neither <- comb_names[!comb_names %in% c(has_genome, has_exome)]
colours <- c(
    rep("red", length(has_genome)),
    rep("blue", length(has_exome)),
    rep("black", length(has_neither))
)
names(colours) <- c(has_genome, has_exome, has_neither)

pdf("figures/assay_upset_plot.pdf", height = 5, width = 7)
UpSet(
    comb_mat,
    set_order = c(colnames(results_matrix)[2:7]),
    comb_order = rev(order(comb_names)),
    comb_col = colours[names(comb_degree(comb_mat))],
    top_annotation = upset_top_annotation(
        comb_mat,
        add_numbers = TRUE,
        numbers_rot = 0
    ),
    right_annotation = upset_right_annotation(comb_mat, add_numbers = TRUE)
)
dev.off()
