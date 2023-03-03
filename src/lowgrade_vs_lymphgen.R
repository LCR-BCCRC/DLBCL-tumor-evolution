source("src/libs.R")

lg_summary <- read_tsv("data/lymphgen/lymphgen_summary.tsv")

indolent_ever <- read_tsv("data/metadata/metadata_lowgrade.tsv")

lg_matrix <- lg_summary %>%
    left_join(indolent_ever) %>%
    group_by(is_denovo) %>%
    pivot_longer(c(Diagnosis_Subtype.Prediction, Relapse_Subtype.Prediction)) %>%
    separate(value, into = paste0("value", 1:5), sep = "/") %>%
    pivot_longer(matches("value"), names_to = "which_val", values_to = "class") %>%
    drop_na(class) %>%
    select(-name, -which_val) %>%
    distinct() %>%
    mutate(is_class = 1) %>%
    pivot_wider(names_from = class, values_from = is_class, values_fill = 0)

lg_colours <- GAMBLR::get_gambl_colours("lymphgen")
lg_colours <- lg_colours[names(lg_colours) %in% colnames(lg_matrix)]

model <- lg_matrix %>%
    select(patient_id, is_denovo, EZB:MCD, -Other) %>%
    pivot_longer(c(EZB:MCD)) %>%
    group_by(name) %>%
    summarize(table = list(table(is_denovo, value))) %>%
    mutate(
        model = purrr::map(table, fisher.test),
        tidy = purrr::map(model, ~ broom::tidy(.x, exponentiate = TRUE))
    ) %>%
    unnest(cols = c(tidy)) %>%
    select(-model, -table, -method, -alternative)

percent_classified <- lg_matrix %>%
    pivot_longer(c(EZB:MCD),
        names_to = "class",
        values_to = "is_class"
    ) %>%
    filter(class != "Other") %>%
    group_by(class, is_denovo) %>%
    summarize(percent_class = round(sum(is_class) / n() * 100, digits = 1)) %>%
    ungroup()


bar <- percent_classified %>%
    ggplot(aes(x = class, y = percent_class, fill = class, alpha = is_denovo)) +
    geom_col(position = "dodge", width = 0.5, colour = "black") +
    xlab("") +
    ylab("% Classified") +
    coord_flip() +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = lg_colours) +
    scale_alpha_manual(
        values = c("TRUE" = 1, "FALSE" = 0.6),
        name = "de novo"
    ) +
    theme(axis.text.y = element_blank(), legend.position = "bottom") +
    guides(fill = "none")


forest <- model %>%
    ggplot(aes(
        x = name,
        y = log(estimate)
    )) +
    geom_point(
        size = 2,
        aes(shape = p.value < 0.05),
        position = position_dodge(width = 0.75)
    ) +
    geom_hline(
        yintercept = 0,
        lty = 2
    ) +
    coord_flip() +
    geom_errorbar(aes(
        ymin = log(conf.low),
        ymax = log(conf.high),
        width = 0.2
    ), position = position_dodge(width = 0.75)) +
    ylab("ln(Odds Ratio)") +
    xlab("LymphGen Class") +
    guides(colour = "none") +
    cowplot::theme_cowplot() +
    theme(legend.position = c(0.03, 0.9))


plots <- plot_grid(
    forest,
    bar + theme(legend.position = c(0.4, 0.7)),
    rel_widths = c(1, 0.8),
    nrow = 1,
    align = "hv",
    axis = "b"
)
plots

ggsave("figures/lg_vs_low_grade.pdf", height = 5, width = 6)
saveRDS(plots, "figures/lg_vs_low_grade.RDS")

lg_vs_lowgrade_entity <- lg_matrix %>%
    mutate(lowgrade_entity = case_when(
        !is.na(transformed_from) ~ transformed_from,
        !is.na(subsequent_lowgrade) ~ subsequent_lowgrade,
        !is.na(BM_infiltrates) ~ str_c("BM infiltrates\n", str_replace(BM_infiltrates, "_", " ")),
        TRUE ~ "None"
    )) %>%
    mutate(lowgrade_entity = str_replace(lowgrade_entity, "COMFL", "FL")) %>%
    mutate(is_denovo = ifelse(is_denovo, "de novo", "Transformed")) %>%
    group_by(lowgrade_entity, is_denovo) %>%
    summarize(
        across(EZB:MCD, sum)
    ) %>%
    ungroup() %>%
    pivot_longer(
        EZB:MCD,
        names_to = "LymphGen",
        values_to = "Count"
    ) %>%
    mutate(LymphGen = factor(
        LymphGen,
        levels = c("EZB", "ST2", "BN2", "N1", "MCD", "A53", "Other")
    )) %>%
    ggplot(aes(x = lowgrade_entity, y = Count, fill = LymphGen)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~is_denovo) +
    scale_fill_manual(values = lg_colours) +
    xlab("") +
    theme(legend.position = "bottom") +
    coord_flip()

ggsave("figures/lg_vs_lowgrade_entity.pdf", height = 4, width = 6)
saveRDS(lg_vs_lowgrade_entity, "figures/lg_vs_lowgrade_entity.RDS")
