# Survival analyses

source("src/libs.R")
library(survival)
library(survminer)
library(table1)
library(kableExtra)

clindata_raw <- read_excel("../private_data/DLBCL_ITT_LH_Aug2022_DWS.xlsx") %>%
    filter(!str_detect(TXPROG, "ALLO"))
    
refined_treatment <- read_tsv("../private_data/refined_salvage.tsv")

clindata_tidy <- data.frame(
    ID = clindata_raw$res_id,
    TTP = clindata_raw$Timetoprogressiony,
    Age_Dx = clindata_raw$AGE,
    Stage_Dx = clindata_raw$STAGE,
    Stage_Relapse = clindata_raw$Stageatrelapse,
    IPI_Dx = clindata_raw$IPI,
    IPI_Relapse = clindata_raw$IPIatrelapse,
    BMT = clindata_raw$BMT,
    ResponseToSalvage = clindata_raw$ResponseToFristRxatrelapse,
    Sex = clindata_raw$SEX,
    DTDX = clindata_raw$DTDX,
    DTPROG = clindata_raw$DTPROG
) %>%
left_join(refined_treatment) %>%
    # mutate(time_to_relapse = as.numeric(DTPROG - DTDX) / 365) %>%
    mutate(time_to_relapse = TTP) %>%
    mutate(relapse_timing = case_when(
        time_to_relapse < 0.75 ~ "Primary Refractory",
        time_to_relapse < 2 ~ "Early Relapse",
        time_to_relapse >= 2 ~ "Late Relapse"
    )) %>%
    mutate(relapse_timing = factor(
        relapse_timing,
        levels = names(relapse_colours)
    )) %>%
    mutate(across(matches("Stage"), ~ case_when(
        .x %in% c(1, 2) ~ "I-II",
        .x %in% c(3, 4) ~ "III-IV"
    ))) %>%
    mutate(across(matches("Stage"), ~ factor(
        .x,
        levels = c("I-II", "III-IV")
    ))) %>%
    mutate(across(matches("IPI"), ~ case_when(
        .x == 0 ~ "0",
        .x %in% c(1, 2) ~ "1-2",
        .x %in% c(3:5) ~ "3-5"
    ))) %>%
    mutate(across(matches("IPI"), ~ factor(
        .x,
        levels = c("0", "1-2", "3-5")
    ))) %>%
    mutate(Salvage_Tx = case_when(
        str_detect(Salvage_regimen, "GDP") ~ "GDP +/- R",
        TRUE ~ "Other"
    )) %>%
    mutate(Salvage_Tx = factor(
        Salvage_Tx,
        levels = c("GDP +/- R", "Other")
    )) %>%
    mutate(BMT = factor(
        BMT,
        levels = c("Y", "N"),
        labels = c("Yes", "No")
    )) %>%
    mutate(ResponseToSalvage = factor(
        ResponseToSalvage,
        levels = c("COMP", "PART", "PROG"),
        labels = c("CR", "PR", "Progression")
    )) %>%
    mutate(Sex = factor(
        Sex,
        levels = c("M", "F"),
        labels = c("Male", "Female")
    )) %>%
    mutate(ORR = case_when(
        ResponseToSalvage %in% c("CR", "PR") ~ "Responder",
        ResponseToSalvage == "Progression" ~ "Non-Responder"
    )) %>%
    select(-DTDX, -DTPROG)

pvals <- clindata_tidy %>%
    mutate(Age_Dx = kruskal.test(Age_Dx ~ relapse_timing)$p.value) %>%
    mutate(across(
        c(Sex, Stage_Dx, IPI_Dx, Stage_Relapse, IPI_Relapse, Salvage_Tx, ResponseToSalvage, BMT),
        # ~ chisq.test(table(.x, relapse_timing), simulate.p.value = TRUE)$p.value
        ~ fisher.test(table(.x, relapse_timing))$p.value
    )) %>%
    select(
        Age_Dx,
        Sex,
        Stage_Dx,
        Stage_Relapse,
        IPI_Dx,
        IPI_Relapse,
        Salvage_Tx,
        ResponseToSalvage,
        BMT
    ) %>%
    distinct() %>%
    pivot_longer(everything(),
        names_to = "Category",
        values_to = "pvalue"
    ) %>%
    mutate(`qvalue` = p.adjust(`pvalue`, method = "BH")) %>%
    mutate(across(where(is.numeric), ~ format.pval(.x, digits = 3, eps = 0.0001))) %>%
    column_to_rownames("Category")

label(clindata_tidy$BMT) <- "Proceed to HSCT"
label(clindata_tidy$ResponseToSalvage) <- "Response to Salvage Therapy"
label(clindata_tidy$Salvage_Tx) <- "Salvage Chemotherapy"
label(clindata_tidy$ORR) <- "Overall Response Rate"
label(clindata_tidy$time_to_relapse) <- "Time to Relapse"
label(clindata_tidy$IPI_Relapse) <- "IPI at Relapse"
label(clindata_tidy$IPI_Dx) <- "IPI at Diagnosis"
label(clindata_tidy$Age_Dx) <- "Age at Diagnosis"
label(clindata_tidy$Stage_Relapse) <- "Stage at Relapse"
label(clindata_tidy$Stage_Dx) <- "Stage at Diagnosis"



pvalue <- function(x, ...) {
    # Hacky method to iterate over variables and add adjusted p-values to table
    item <<- item + 1
    c("", pvals[item - 1, "qvalue"])
}

item <- 1
char_table <- table1(
    ~ Age_Dx + Sex + Stage_Dx + IPI_Dx + Stage_Relapse + IPI_Relapse + Salvage_Tx + ResponseToSalvage + BMT + time_to_relapse | relapse_timing,
    data = clindata_tidy,
    overall = NULL,
    extra.col = list(`Adjusted P-value` = pvalue)
)

saveRDS(char_table, "figures/surv_characteristic_table.RDS")

char_kable <- t1kable(char_table)


percent_BMT <- clindata_tidy %>%
    group_by(relapse_timing) %>%
    count(BMT) %>%
    mutate(percent = round(n / sum(n) * 100, digits = 2)) %>%
    filter(BMT == "Yes")

write_tsv(percent_BMT, "data/survival/BMT_summary.tsv")

bmt_table <- table(
    BMT = clindata_tidy$BMT,
    relapse_timing = clindata_tidy$relapse_timing
)

bmt_pairwise <- pairwise_fisher_test(bmt_table, p.adjust.method = "BH") %>%
    left_join(
        select(percent_BMT, relapse_timing, percent),
        by = c("group2" = "relapse_timing")
    ) %>%
    group_by(group2) %>%
    mutate(y.position = case_when(
        row_number(group2) != min(row_number(group2)) ~ percent + 6 * (1 + row_number(group2) - min(row_number(group2))),
        TRUE ~ percent + 6
    )) %>%
    ungroup() %>%
    mutate(p.adj.signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))

bmt_barplot <- percent_BMT %>%
    ggplot(aes(
        x = relapse_timing,
        y = percent,
        fill = relapse_timing
    )) +
    geom_col(colour = "darkgrey") +
    geom_text(
        aes(
            label = str_c("N=", n),
            x = relapse_timing,
            y = percent
        ),
        position = position_stack(reverse = FALSE, vjust = 0.5),
        colour = "black",
        inherit.aes = FALSE
    ) +
    stat_pvalue_manual(
        bmt_pairwise,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        label = "p.adj.signif",
        inherit.aes = FALSE
    ) +
    scale_fill_manual(
        values = relapse_colours,
        name = "Relapse Timing",
        guide = guide_none()
    ) +
    xlab("") +
    ylab("Received HSCT (%)") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

bmt_barplot

ggsave("figures/bmt_barplot.pdf", height = 4, width = 4)
saveRDS(bmt_barplot, "figures/bmt_barplot.RDS")

percent_response <- clindata_tidy %>%
    filter(!is.na(ResponseToSalvage)) %>%
    group_by(relapse_timing) %>%
    count(ResponseToSalvage) %>%
    mutate(percent = round(n / sum(n) * 100, digits = 2))

write_tsv(percent_response, "data/survival/response_summary.tsv")

response_table <- table(
    ResponseToSalvage = ifelse(
        clindata_tidy$ResponseToSalvage %in% c("CR", "PR"),
        "Response",
        as.character(clindata_tidy$ResponseToSalvage)
    ),
    relapse_timing = clindata_tidy$relapse_timing
)

response_pairwise <- pairwise_fisher_test(response_table, p.adjust.method = "BH") %>%
    left_join(
        select(filter(
            percent_response,
            ResponseToSalvage == "Progression"
        ), relapse_timing, percent),
        by = c("group2" = "relapse_timing")
    ) %>%
    group_by(group2) %>%
    mutate(y.position = case_when(
        row_number(group2) != min(row_number(group2)) ~ (100 - percent) + 6 * (1 + row_number(group2) - min(row_number(group2))),
        TRUE ~ (100 - percent) + 5
    )) %>%
    ungroup() %>%
    mutate(p.adj.signif = ifelse(p.adj.signif == "****", "***", p.adj.signif))

response_plot <- percent_response %>%
    filter(ResponseToSalvage %in% c("PR", "CR")) %>%
    mutate(ResponseToSalvage = fct_rev(ResponseToSalvage)) %>%
    ggplot(aes(
        x = relapse_timing,
        y = percent,
        fill = relapse_timing,
        alpha = ResponseToSalvage
    )) +
    geom_col(colour = "darkgrey") +
    geom_text(
        aes(
            label = str_c("N=", n),
            x = relapse_timing,
            y = percent
        ),
        position = position_stack(reverse = FALSE, vjust = 0.5),
        colour = "black",
        inherit.aes = FALSE
    ) +
    stat_pvalue_manual(
        response_pairwise,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        label = "p.adj.signif",
        inherit.aes = FALSE
    ) +
    scale_fill_manual(
        values = relapse_colours,
        name = "Relapse Timing"
    ) +
    scale_alpha_manual(
        name = "Response",
        values = c(0.6, 1),
        labels = c("Partial", "Complete")
    ) +
    xlab("") +
    ylab("Response (%)") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

response_plot

ggsave("figures/response_barplot.pdf", response_plot, height = 4, width = 6)
saveRDS(response_plot, "figures/response_barplot.RDS")



surv_wide <- clindata_raw %>%
    select(
        ID = res_id,
        BMT,
        ResponseToSalvage = ResponseToFristRxatrelapse,
        OSSCT_Y = OSpostSCT,
        OSSCT_CODE = CODE_OS,
        PFSSCT_Y = PFSPostSCT,
        PFSSCT_CODE = CODE_ProgPostBMT,
        OSPROG_Y = OSPostProg,
        OSPROG_CODE = CODE_OS,
        PFSPROG_Y = PFSpostProg,
        PFSPROG_CODE = CODEPROG2,
        TreatmentPrimary = TXPRIMARY,
        TreatmentProgression = TXPROG
    ) %>%
    left_join(select(clindata_tidy, ID, relapse_timing, Age_Dx, IPI_Relapse))

write_tsv(surv_wide, "data/survival/surv_wide.tsv")



surv_tidy <- surv_wide %>%
    # mutate(across(matches("_Y"), as.character)) %>%
    pivot_longer(
        matches("PFS|OS", ignore.case = FALSE),
        names_to = "metric",
        values_to = "value"
    ) %>%
    separate(
        metric,
        into = c("metric", "parameter"),
        sep = "_"
    ) %>%
    pivot_wider(
        names_from = parameter,
        values_from = value
    ) %>%
    mutate(from = factor(
        ifelse(str_detect(metric, "SCT"), "HSCT", "Progression"),
        levels = c("Progression", "HSCT")
    )) %>%
    mutate(metric = factor(
        str_remove(metric, "SCT|PROG"),
        levels = c("PFS", "OS")
    )) %>%
    mutate(relapse_timing = factor(
        relapse_timing,
        levels = names(relapse_colours)
    ))

survtidy_cut10 <- surv_tidy %>%
    mutate(CODE = ifelse(
        Y > 10,
        0,
        CODE
    )) %>%
    mutate(Y = ifelse(
        Y > 10,
        10,
        Y
    ))

surv_wide_cut <- survtidy_cut10 %>%
    pivot_wider(
        names_from = c("metric", "from"),
        values_from = c("Y", "CODE"),
        names_glue = c("{metric}{from}_{.value}")
    )

responses <- colnames(surv_wide_cut)[str_detect(colnames(surv_wide_cut), "_CODE")]

survplot_fun <- function(response) {
    formula <- surv_fit(as.formula(paste0(
        "Surv(", str_replace(response, "_CODE", "_Y"), ", ", response, ") ~ relapse_timing"
    )), data = surv_wide_cut)
    since <- ifelse(str_detect(response, "SCT"), "ASCT", "Progression")
    outcome <- ifelse(str_detect(response, "PFS"), "PFS", "OS")
    title <- paste0(outcome, " from time of ", since)
    plot <- ggsurvplot(
        formula,
        surv_wide_cut,
        palette = palette(relapse_colours),
        risk.table = TRUE,
        legend.title = "",
        legend.labs = c("REFR", "ER", "LR"),
        xlab = "Time (y)",
        title = title,
        break.time.by = 2
    )
    plot_grid(
        plot$plot +
            theme(
                legend.position = "none",
                plot.margin = unit(c(0, 0, 0, 0), "cm")
            ),
        plot$table +
            theme(
                plot.margin = unit(c(0, 0, 0, 0), "cm"),
                plot.title = element_text(size = 14)
            ) +
            xlab(""),
        nrow = 2,
        rel_heights = c(1, 0.5),
        align = "hv"
    )
}
surv_with_tables <- plot_grid(
    survplot_fun(responses[4]),
    survplot_fun(responses[3]),
    survplot_fun(responses[2]),
    survplot_fun(responses[1]),
    nrow = 2,
    ncol = 2,
    labels = c("C", "D", "E", "F"),
    label_size = 22
)

ggsave("figures/surv_with_tables.pdf", surv_with_tables, height = 10, width = 10)
saveRDS(surv_with_tables, "figures/surv_with_tables.RDS")

survdata_HSCT <- data.frame(filter(
    survtidy_cut10,
    from == "HSCT"
))

fit_HSCT <- survfit(
    Surv(Y, CODE) ~ relapse_timing,
    data = survdata_HSCT
)

surv_HSCT <- ggsurvplot(
    fit_HSCT,
    data = survdata_HSCT,
    facet.by = "metric",
    # legend.labs = names(relapse_colours),
    # palette = (relapse_colours),
    pval = TRUE,
    risk.table = TRUE,
    short.panel.labs = TRUE,
    legend.title = "Relapse Timing",
    pval.coord = c(6, 75),
    fun = "pct"
) +
    theme_cowplot() +
    scale_colour_manual(values = relapse_colours) +
    theme(
        legend.position = "right"
    ) +
    xlab("Time since HSCT (years)")

surv_HSCT

ggsave("figures/survival_HSCT_faceted.pdf", height = 4, width = 8)
saveRDS(surv_HSCT, "figures/survival_HSCT_faceted.RDS")

survdata_prog <- data.frame(filter(
    survtidy_cut10,
    from == "Progression"
)) %>%
    mutate(metric = factor(metric, levels = c("PFS", "OS")))

fit_prog <- survfit(
    Surv(Y, CODE) ~ relapse_timing,
    data = survdata_prog
)

surv_prog <- ggsurvplot(
    fit_prog,
    data = survdata_prog,
    facet.by = "metric",
    pval = TRUE,
    risk.table = TRUE,
    short.panel.labs = TRUE,
    legend.title = "Relapse Timing",
    pval.coord = c(6, 75),
    cumevents = TRUE,
    combine = TRUE,
    fun = "pct"
) +
    theme_cowplot() +
    scale_colour_manual(values = relapse_colours) +
    theme(
        legend.position = "right"
    ) +
    xlab("Time since progression (years)")

surv_prog

ggsave("figures/survival_prog_faceted.pdf", height = 4, width = 8)
saveRDS(surv_prog, "figures/survival_prog_faceted.RDS")

fit_facet <- survfit(
    Surv(Y, CODE) ~ relapse_timing + metric + from,
    data = survtidy_cut10
)


surv_facet <- ggsurvplot_facet(
    fit_facet,
    data = surv_tidy,
    facet.by = c("from", "metric"),
    # legend.labs = names(relapse_colours),
    # palette = (relapse_colours),
    pval = FALSE,
    risk.table = FALSE,
    short.panel.labs = TRUE,
    fun = "pct",
    xlim = c(0, 10)
) +
    theme_cowplot() +
    scale_colour_manual(values = relapse_colours) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank()
    ) +
    xlab("Time (years)")

surv_facet

ggsave("figures/survival_faceted.pdf", surv_facet, height = 8, width = 8)
saveRDS(surv_facet, "figures/survival_faceted.RDS")

survmodel <- function(df) {
    coxph(Surv(df$Y, df$CODE) ~ df$relapse_timing)
}

test <- surv_tidy %>%
    mutate(relapse_timing = fct_rev(relapse_timing)) %>%
    group_by(metric, from) %>%
    nest() %>%
    mutate(
        cox = map(data, survmodel),
        cox_tidy = map(cox, broom::tidy, exponentiate = TRUE, conf.int = TRUE)
    ) %>%
    select(-data, -cox) %>%
    unnest(cols = c(cox_tidy)) %>%
    mutate(term = str_replace(term, "df\\$relapse_timing", "")) %>%
    # mutate(p.value = round(p.value, digits = 4)) %>%
    mutate(p.value.format = format.pval(
        round(p.value, digits = 3),
        digits = 2,
        eps = 0.001
    )) %>%
    mutate(p.term = str_c(
        "HR = ", round(estimate, digits = 2),
        " (95% CI ", round(conf.low, digits = 2), " - ", round(conf.high, digits = 2),
        "), P = ", p.value.format
    ))

write_tsv(test, "data/survival/cox_models.tsv")

forestplot <- test %>%
    ggplot(aes(
        x = estimate,
        y = term,
        colour = term
    )) +
    geom_segment(
        aes(
            x = conf.low,
            xend = conf.high,
            y = term,
            yend = term
        ),
        inherit.aes = FALSE
    ) +
    geom_point(
        aes(size = -log10(p.value))
    ) +
    geom_text(
        aes(
            label = p.value.format,
            x = 9,
            y = term
        ),
        inherit.aes = FALSE
    ) +
    geom_vline(
        xintercept = 1,
        linetype = "dashed"
    ) +
    facet_grid(
        from + metric ~ 1
    ) +
    scale_size(
        breaks = c(-log10(0.05), 2, 3, 4),
        labels = c("0.05", "0.01", "0.001", "0.0001"),
        name = "P value\nvs. Late Relapse"
    ) +
    scale_colour_manual(
        values = relapse_colours[1:2]
    ) +
    scale_x_continuous(
        expand = expansion(add = c(0, 0.25)),
        trans = "log",
        breaks = c(0.5, 1, 2, 4, 6, 8)
    ) +
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    ) +
    guides(colour = "none", size = "none") +
    xlab("Hazard Ratio") +
    ylab("")

ggsave("figures/surv_forestplot.pdf", forestplot, height = 6, width = 6)

survmodel_adj <- function(df) {
    coxph(Surv(df$Y, df$CODE) ~ df$relapse_timing + df$IPI_Relapse + df$Age_Dx)
}

test_adj <- surv_tidy %>%
    mutate(relapse_timing = fct_rev(relapse_timing)) %>%
    group_by(metric, from) %>%
    nest() %>%
    mutate(
        cox = map(data, survmodel_adj),
        cox_tidy = map(cox, broom::tidy, exponentiate = TRUE, conf.int = TRUE)
    ) %>%
    select(-data, -cox) %>%
    unnest(cols = c(cox_tidy)) %>%
    mutate(term = str_replace(term, "df\\$relapse_timing", "")) %>%
    mutate(term = str_replace(term, "df\\$", "")) %>%
    mutate(term = str_replace(term, "_", " ")) %>%
    # mutate(p.value = round(p.value, digits = 3)) %>%
    mutate(p.value.format = format.pval(
        round(p.value, digits = 3),
        digits = 2,
        eps = 0.001
    )) %>%
    mutate(p.term = str_c(
        "HR = ", round(estimate, digits = 2),
        " (95% CI ", round(conf.low, digits = 2), " - ", round(conf.high, digits = 2),
        "), P = ", p.value.format
    ))

write_tsv(test_adj, "data/survival/cox_models_adjusted.tsv")

forestplot_adj <- test_adj %>%
    mutate(term = str_replace(term, "Age Dx", "Age at Diagnosis")) %>%
    mutate(term = factor(term, levels = c("Age at Diagnosis", rev(unique(test_adj$term))))) %>%
    ggplot(aes(
        x = estimate,
        y = term,
        colour = term
    )) +
    geom_segment(
        aes(
            x = conf.low,
            xend = conf.high,
            y = term,
            yend = term
        ),
        inherit.aes = FALSE
    ) +
    geom_point(
        aes(size = -log10(p.value))
    ) +
    geom_text(
        aes(
            label = p.value.format,
            x = 40,
            y = term
        ),
        inherit.aes = FALSE
    ) +
    geom_vline(
        xintercept = 1,
        linetype = "dashed"
    ) +
    facet_grid(
        from + metric ~ 1
    ) +
    scale_size(
        breaks = c(-log10(0.5), 2, 3, 4),
        labels = c("0.05", "0.01", "0.001", "0.0001"),
        name = "P value\nvs. Late Relapse"
    ) +
    scale_colour_manual(
        values = relapse_colours[1:2]
    ) +
    scale_x_continuous(
        trans = "log",
        breaks = c(0.5, 1, 5, 10, 20),
        expand = expansion(add = c(0, 0.5))
    ) +
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    ) +
    guides(colour = "none", size = "none") +
    xlab("Hazard Ratio") +
    ylab("")

ggsave("figures/surv_forestplot_adjusted.pdf", forestplot_adj, height = 6, width = 7)


survmodel_cont <- function(df) {
    coxph(Surv(df$Y, df$CODE) ~ df$time_to_relapse)
}

test_cont <- surv_tidy %>%
    left_join(select(clindata_tidy, ID, time_to_relapse)) %>%
    group_by(metric, from) %>%
    nest() %>%
    mutate(
        cox = map(data, survmodel_cont),
        cox_tidy = map(cox, broom::tidy, exponentiate = TRUE, conf.int = TRUE)
    ) %>%
    select(-data, -cox) %>%
    unnest(cols = c(cox_tidy)) %>%
    mutate(term = "Time to Relapse") %>%
    # mutate(p.value = round(p.value, digits = 4)) %>%
    mutate(p.value.format = format.pval(
        round(p.value, digits = 3),
        digits = 2,
        eps = 0.001
    )) %>%
    mutate(p.term = str_c(
        "HR = ", round(estimate, digits = 2),
        " (95% CI ", round(conf.low, digits = 2), " - ", round(conf.high, digits = 2),
        "), P = ", p.value.format
    ))

write_tsv(test_cont, "data/survival/cox_models_cont.tsv")

surv_5yr <- surv_tidy %>%
    mutate(relapse_timing = fct_rev(relapse_timing)) %>%
    group_by(metric, from, relapse_timing) %>%
    nest() %>%
    mutate(
        summary = map(data, ~ summary(survfit(Surv(Y, CODE) ~ 1, data = .), times = 5)),
        survival = map(summary, ~ .$surv),
        conf_low = map(summary, ~ .$lower),
        conf_high = map(summary, ~ .$upper)
    ) %>%
    select(-data, -summary) %>%
    unnest(cols = c(survival, conf_low, conf_high)) %>%
    ungroup()

write_tsv(surv_5yr, "data/survival/five_year_survival.tsv")

surv_2yr <- surv_tidy %>%
    mutate(relapse_timing = fct_rev(relapse_timing)) %>%
    group_by(metric, from, relapse_timing) %>%
    nest() %>%
    mutate(
        summary = map(data, ~ summary(survfit(Surv(Y, CODE) ~ 1, data = .), times = 2)),
        survival = map(summary, ~ .$surv),
        conf_low = map(summary, ~ .$lower),
        conf_high = map(summary, ~ .$upper)
    ) %>%
    select(-data, -summary) %>%
    unnest(cols = c(survival, conf_low, conf_high)) %>%
    ungroup()

write_tsv(surv_2yr, "data/survival/two_year_survival.tsv")

lrmodel <- function(df) {
    pairwise_survdiff(Surv(Y, CODE) ~ relapse_timing, data = df)
}

test_lr <- surv_tidy %>%
    mutate(relapse_timing = fct_rev(relapse_timing)) %>%
    group_by(metric, from) %>%
    nest() %>%
    mutate(
        lr = map(data, lrmodel),
        lr_tidy = map(lr, broom::tidy, conf.int = TRUE)
    ) %>%
    select(-data, -lr) %>%
    unnest(cols = c(lr_tidy)) %>%
    mutate(p.value.format = format.pval(
        p.value,
        digits = 3,
        eps = 0.001
    ))

write_tsv(test_lr, "data/survival/logrank_tests.tsv")

test_lr_2yr <- surv_tidy %>%
    mutate(CODE = ifelse(Y > 2, 0, CODE)) %>%
    mutate(relapse_timing = fct_rev(relapse_timing)) %>%
    group_by(metric, from) %>%
    nest() %>%
    mutate(
        lr = map(data, lrmodel),
        lr_tidy = map(lr, broom::tidy, conf.int = TRUE)
    ) %>%
    select(-data, -lr) %>%
    unnest(cols = c(lr_tidy)) %>%
    mutate(p.value.format = format.pval(
        p.value,
        digits = 3,
        eps = 0.001
    ))

write_tsv(test_lr_2yr, "data/survival/logrank_tests_2yr.tsv")



# Make a two-panel figure for ASH 2022 abstract
# lymphgen_comp <- readRDS("figures/lymphgen_comparison.RDS")


# surv_pfs_prog <- survtidy_cut10 %>%
#     filter(metric == "PFS",
#            from == "Progression") %>%
#     data.frame()

# fit_pfs_prog <- survfit(
#     Surv(Y, CODE) ~ relapse_timing,
#     data = surv_pfs_prog
# )


# plot_pfs_prog <- ggsurvplot(
#     fit_pfs_prog,
#     data = surv_pfs_prog,
#     legend.labs = names(relapse_colours),
#     palette = (relapse_colours),
#     pval = TRUE,
#     risk.table = FALSE,
#     short.panel.labs = TRUE,
#     pval.coord = c(6, 75),
#     fun = "pct",
#     xlim = c(0, 10),
#     ggtheme = theme_cowplot(),
#     legend = "right",
#     legend.title = "Relapse Timing",
#     xlab = "Time since progression (years)",
#     title = "Progression-free survival from first progression"
# )

# plot_pfs_prog

# bottom_row = plot_grid(
#     NULL,
#     plot(plot_pfs_prog$plot),
#     NULL,
#     rel_widths = c(0.2, 1, 0.2),
#     nrow = 1,
#     labels = c("", "B", ""),
#     label_size = 18
# )

# top_row = plot_grid(
#     lymphgen_comp +
#         theme(legend.position = "bottom") +
#         ggtitle("LymphGen Classifications at Relapse vs Diagnosis"),
#     labels = c("A"),
#     label_size = 18
# )

# plot_grid(
#     top_row,
#     bottom_row,
#     nrow = 2,
#     rel_heights = c(1, 0.7)
# )

# ggsave("figures/ASH_figure.pdf", height = 8, width = 8)
