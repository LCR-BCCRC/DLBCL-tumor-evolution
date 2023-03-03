##### Plot FISH and NanoString data #####

source("src/libs.R")

fish_ns_pairs <- read_tsv("data/metadata/metadata_fish_nanostring.tsv") %>%
  mutate(relapse_timing = factor(relapse_timing, levels = names(relapse_colours)))

fish_pairs <- fish_ns_pairs %>%
  select(
    -matches("call|score|code|cn")
  ) %>%
  group_by(patient_id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  pivot_longer(matches("_ba"),
    names_to = "assay",
    values_to = "result"
  ) %>%
  group_by(patient_id, assay) %>%
  mutate(discordant = ifelse(
    max(result) == min(result),
    "Concordant",
    "Discordant"
  )) %>%
  ungroup() %>%
  mutate(
    assay = factor(toupper(str_remove(assay, "_ba_consensus")),
      levels = c("MYC", "BCL2", "BCL6")
    ),
    result = factor(result, levels = c("POS", "NEG"))
  )


fish_plot <- fish_pairs %>%
  filter(!is.na(result)) %>%
  group_by(assay, patient_id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  ggplot(aes(
    x = time_point,
    stratum = result,
    alluvium = patient_id,
    fill = result
  )) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(alpha = discordant)) +
  geom_stratum(colour = NA) +
  scale_fill_manual(
    values = get_gambl_colours(classification = "pos_neg")[c("POS", "NEG")],
    name = "BA FISH"
  ) +
  scale_alpha_manual(values = c("Discordant" = 1, "Concordant" = 0.5)) +
  xlab("") +
  facet_grid2(assay ~ relapse_timing, strip = relapse_strips) +
  theme(
    strip.text.y = element_text(face = "italic"),
    legend.position = c(0.05, 0.9),
    legend.direction = "horizontal"
  ) +
  guides(alpha = "none")

ggsave("figures/fish_pairs.pdf", fish_plot, height = 8, width = 8)
saveRDS(fish_plot, "figures/fish_pairs.RDS")

summarize_fish <- fish_pairs %>%
  filter(!is.na(result)) %>%
  group_by(assay, patient_id) %>%
  filter(n() == 2) %>%
  mutate(multi_FISH = TRUE) %>%
  ungroup() %>%
  group_by(assay, time_point, relapse_timing) %>%
  count(discordant) %>%
  mutate(total_time_point = sum(n)) %>%
  mutate(percent = round(n / total_time_point * 100, digits = 2)) %>%
  ungroup() %>%
  group_by(time_point, assay) %>%
  mutate(total_assay = sum(n)) %>%
  ungroup() %>%
  select(-time_point) %>%
  distinct()

write_tsv(summarize_fish, "data/fish_nanostring/fish_summary.tsv")

fish_pairs %>%
  filter(!is.na(result)) %>%
  group_by(assay, patient_id) %>%
  filter(n() == 2) %>%
  mutate(multi_FISH = TRUE) %>%
  ungroup() %>%
  select(patient_id, assay, multi_FISH) %>%
  distinct() %>%
  pivot_wider(names_from = assay, values_from = multi_FISH) %>%
  write_tsv("data/fish_nanostring/multi_FISH_summary.tsv")


ns_pairs <- fish_ns_pairs %>%
  select(-matches("consensus")) %>%
  drop_na(code_set) %>%
  mutate(across(
    matches("dhitsig_"),
    ~ ifelse(dlbcl_call %in% c("GCB", "UNCLASS"), .x, NA)
  )) %>%
  group_by(patient_id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  mutate(across(matches("_score"), as.character)) %>%
  pivot_longer(matches("call|score"),
    names_to = c("assay", "measure"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = "measure",
    values_from = "value"
  ) %>%
  mutate(score = as.numeric(score)) %>%
  group_by(patient_id, assay) %>%
  mutate(discordant = case_when(
    (max(call) == "ABC" & min(call) == "GCB") |
      (max(call) == "GCB" & min(call) == "ABC") ~ "Discordant",
    (max(call) == "DHITsig+" & min(call) != "DHITsig+") |
      (min(call) == "DHITsig+" & max(call) != "DHITsig+") ~ "Discordant",
    !is.na(max(call)) & !is.na(min(call)) ~ "Concordant"
  )) %>%
  ungroup()

summarize_ns <- ns_pairs %>%
  group_by(assay, time_point, relapse_timing) %>%
  count(discordant) %>%
  mutate(total_time_point = sum(n)) %>%
  mutate(percent = round(n / total_time_point * 100, digits = 2)) %>%
  ungroup() %>%
  group_by(time_point, assay) %>%
  mutate(total_assay = sum(n)) %>%
  ungroup()

write_tsv(summarize_ns, "data/fish_nanostring/ns_summary.tsv")


dhitsig_levels <- c("DHITsig-", "DHITsig-IND", "DHITsig+")

dhitsig_score_plot <- ns_pairs %>%
  filter(assay == "dhitsig") %>%
  select(
    patient_id,
    relapse_timing,
    time_point,
    score,
    call,
    discordant
  ) %>%
  pivot_wider(
    names_from = time_point,
    values_from = c(score, call),
    names_glue = "{time_point}_{.value}"
  ) %>%
  mutate(Diagnosis_call = factor(Diagnosis_call, levels = dhitsig_levels)) %>%
  ggplot(aes(
    x = Diagnosis_score,
    y = Relapse_score
  )) +
  geom_hline(yintercept = dhitsig_cuts, colour = "lightgrey") +
  geom_vline(xintercept = dhitsig_cuts, colour = "lightgrey") +
  geom_smooth(
    method = "lm",
    colour = "black",
    se = FALSE
  ) +
  geom_point(aes(
    colour = discordant,
    size = discordant
  ),
  shape = 1
  ) +
  geom_point(aes(colour = Diagnosis_call),
    size = 2
  ) +
  facet_wrap2(~relapse_timing, strip = relapse_strips) +
  scale_colour_manual(
    values = c(coo_colours[dhitsig_levels],
      "Discordant" = "red",
      "Concordant" = "white"
    ),
    name = ""
  ) +
  scale_size_manual(values = c("Discordant" = 4, "Concordant" = 0)) +
  stat_cor(p.accuracy = 0.001) +
  xlab("Diagnostic Score") +
  ylab("Relapse Score") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  theme(legend.position = "none")

dhitsig_alluvial_plot <- ns_pairs %>%
  filter(
    assay == "dhitsig",
    !is.na(call)
  ) %>%
  mutate(call = factor(call, levels = dhitsig_levels)) %>%
  group_by(patient_id) %>%
  filter(n() == 2) %>%
  select(
    patient_id,
    relapse_timing,
    time_point,
    discordant,
    call
  ) %>%
  ggplot(aes(
    x = time_point,
    alluvium = patient_id,
    stratum = call,
    fill = call
  )) +
  geom_stratum(colour = NA) +
  geom_flow(aes(alpha = discordant)) +
  scale_fill_manual(
    values = coo_colours[dhitsig_levels], name = "",
    labels = c(bquote(DZsig^NEG), bquote(DZsig^IND), bquote(DZsig^POS))
  ) +
  scale_alpha_manual(values = c("Discordant" = 1, "Concordant" = 0.5)) +
  xlab("") +
  ylab("Count") +
  facet_wrap2(~relapse_timing, strip = relapse_strips) +
  guides(alpha = "none")

legend <- get_legend(
  dhitsig_alluvial_plot +
    theme(
      legend.position = c(0.35, 0),
      legend.direction = "horizontal",
      legend.margin = margin(0, 0, 40, 0)
    )
)

dhitsig_plot <- plot_grid(dhitsig_alluvial_plot + theme(legend.position = "none"),
  legend,
  dhitsig_score_plot + theme(legend.position = "none"),
  nrow = 3,
  rel_heights = c(1, 0.1, 1),
  labels = c("A", "", "B"),
  align = "v",
  axis = "l"
)

ggsave("figures/dhitsig_alluvial_score.pdf", dhitsig_plot, height = 8, width = 10)
saveRDS(dhitsig_plot, "figures/dhitsig_alluvial_score.RDS")

coo_levels <- c("ABC", "UNCLASS", "GCB")

coo_score_plot <- ns_pairs %>%
  filter(assay == "dlbcl") %>%
  select(
    patient_id,
    relapse_timing,
    time_point,
    score,
    discordant,
    call
  ) %>%
  pivot_wider(
    names_from = time_point,
    values_from = c(score, call),
    names_glue = "{time_point}_{.value}"
  ) %>%
  mutate(Diagnosis_call = factor(Diagnosis_call, levels = coo_levels)) %>%
  ggplot(aes(
    x = Diagnosis_score,
    y = Relapse_score
  )) +
  geom_hline(yintercept = coo_cuts, colour = "lightgrey") +
  geom_vline(xintercept = coo_cuts, colour = "lightgrey") +
  geom_smooth(
    method = "lm",
    colour = "black",
    se = FALSE
  ) +
  geom_point(aes(
    colour = discordant,
    size = discordant
  ),
  shape = 1
  ) +
  geom_point(aes(colour = Diagnosis_call),
    size = 2
  ) +
  facet_wrap2(~relapse_timing, strip = relapse_strips) +
  scale_colour_manual(
    values = c(coo_colours[coo_levels],
      "Discordant" = "red",
      "Concordant" = "white"
    ),
    name = ""
  ) +
  scale_x_continuous(breaks = c(-1500, 0, 1500, 3000), limits = c(-1500, 4000)) +
  scale_y_continuous(breaks = c(-1500, 0, 1500, 3000), limits = c(-1500, 4000)) +
  scale_size_manual(values = c("Discordant" = 4, "Concordant" = 0)) +
  stat_cor(p.accuracy = 0.001, label.y = -1400) +
  xlab("Diagnostic Score") +
  ylab("Relapse Score") +
  theme(legend.position = "none")

saveRDS(coo_score_plot, "figures/coo_score.RDS")

coo_alluvial_plot <- ns_pairs %>%
  filter(
    assay == "dlbcl",
    !is.na(call)
  ) %>%
  mutate(call = factor(call, levels = coo_levels)) %>%
  group_by(patient_id) %>%
  filter(n() == 2) %>%
  select(
    patient_id,
    relapse_timing,
    time_point,
    discordant,
    call
  ) %>%
  ggplot(aes(
    x = time_point,
    alluvium = patient_id,
    stratum = call,
    fill = call
  )) +
  geom_stratum(colour = NA) +
  geom_flow(aes(alpha = discordant)) +
  scale_fill_manual(values = coo_colours[coo_levels], name = "") +
  scale_alpha_manual(values = c("Discordant" = 1, "Concordant" = 0.5)) +
  xlab("") +
  ylab("Count") +
  facet_wrap2(~relapse_timing, strip = relapse_strips) +
  guides(alpha = "none")

saveRDS(coo_alluvial_plot, "figures/coo_alluvial.RDS")

legend <- get_legend(
  coo_alluvial_plot +
    theme(
      legend.position = c(0.35, 0),
      legend.direction = "horizontal",
      legend.margin = margin(0, 0, 40, 0)
    )
)

coo_plot <- plot_grid(coo_alluvial_plot + theme(legend.position = "none"),
  legend,
  coo_score_plot + theme(legend.position = "none") + theme(strip.background = element_blank(), strip.text = element_blank()),
  nrow = 3,
  rel_heights = c(1, 0.1, 1),
  labels = c("A", "", "B"),
  label_size = 18,
  align = "v",
  axis = "l"
)

ggsave("figures/coo_alluvial_score.pdf", coo_plot, height = 8, width = 10)
saveRDS(coo_plot, "figures/coo_alluvial_score.RDS")

ns_pairs %>%
  filter(
    time_point == "Relapse",
    relapse_timing %in% c("Early Relapse", "Late Relapse")
  ) %>%
  mutate(relapse_timing = factor(
    relapse_timing,
    levels = c("Early Relapse", "Late Relapse")
  )) %>%
  with(table(relapse_timing, discordant)) %>%
  fisher.test(.x)

ns_pairs %>%
  filter(
    time_point == "Relapse",
    relapse_timing %in% c("Primary Refractory", "Late Relapse")
  ) %>%
  mutate(relapse_timing = factor(
    relapse_timing,
    levels = c("Primary Refractory", "Late Relapse")
  )) %>%
  with(table(relapse_timing, discordant)) %>%
  fisher.test(.x)

pts_exclude <- c(
  "03-19103", # Composite SLL/cHL at one time point
  "92-38626", # PTLD
  "04-13783", # Cutaneous leg-type DLBCL
  "12-17272", # PCNSL
  "14-19181" # PCNSL
)

ns_pairs %>%
  filter(time_point == "Relapse", !patient_id %in% pts_exclude) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  count(discordant)
