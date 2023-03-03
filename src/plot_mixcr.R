##### Plot MiXCR results #####

source("src/libs.R")
library(ggalluvial)

lightchain_alluvial <- read_tsv("data/ig_rearrangement/lightchain_alluvial.tsv")

lc_plot <- lightchain_alluvial %>%
mutate(relapse_timing = factor(
    relapse_timing,
    levels = names(relapse_colours)
)) %>%
  ggplot(aes(
    x = time_point,
    stratum = `V Gene`,
    alluvium = patient_id,
    fill = `V Gene`
  )) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            aes(alpha = discordant)) +
  geom_stratum(colour = NA) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  facet_wrap2(~ relapse_timing, nrow = 1, strip = relapse_strips) +
  xlab("") +
  guides(alpha = "none")

lc_plot

ggsave("figures/mixcr_lightchain_alluvial.pdf", lc_plot, height = 4, width = 10)
saveRDS(lc_plot, "figures/mixcr_lightchain_alluvial.RDS")

heavychain_alluvial <- read_tsv("data/ig_rearrangement/heavychain_alluvial.tsv")

hc_plot <- heavychain_alluvial %>%
mutate(relapse_timing = factor(
    relapse_timing,
    levels = names(relapse_colours)
)) %>%
  ggplot(aes(
    x = time_point,
    stratum = `V Gene`,
    alluvium = patient_id,
    fill = `V Gene`
  )) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            aes(alpha = discordant)) +
  geom_stratum(colour = NA) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  facet_wrap2(~ relapse_timing, nrow = 1, strip = relapse_strips) +
  scale_y_continuous(breaks = c(0,5,10)) +
  xlab("") +
  guides(alpha = "none")

hc_plot

ggsave("figures/mixcr_heavychain_alluvial.pdf", hc_plot, height = 4, width = 10)
saveRDS(hc_plot, "figures/mixcr_heavychain_alluvial.RDS")

summary_heavychain <- heavychain_alluvial %>%
  select(patient_id, discordant, relapse_timing) %>%
  distinct() %>%
  group_by(relapse_timing) %>%
  count(discordant) %>%
  mutate(IG_gene = "IGH")

summary_lightchain <- lightchain_alluvial %>%
  select(patient_id, discordant, relapse_timing) %>%
  distinct() %>%
  group_by(relapse_timing) %>%
  count(discordant)%>%
  mutate(IG_gene = "IGK/L")

bind_rows(summary_heavychain, summary_lightchain) %>%
  write_tsv("data/ig_rearrangement/summary_ig_concordance.tsv")


