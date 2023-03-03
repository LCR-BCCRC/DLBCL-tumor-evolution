##### Swimmer's Plot #####
# Create a column indicating time of each biopsy 
# relative to the first DLBCL sample with sequencing data

source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv")
clinical_md_biopsy <- read_tsv("data/metadata/metadata_clinical_biopsy.tsv")
clinical_md_patient <- read_tsv("data/metadata/metadata_clinical_patient.tsv")

trios_md_less <- trios_md %>% 
  filter(patient_id %in% clinical_md_patient$patient_id, 
         tissue_status == "tumour") %>% 
  select(-matches("RNAseq")) %>% 
  rename(seq_type = DNAseq_type) %>% 
  rename_with(~str_remove(.x, "DNAseq_")) %>% 
  mutate(biopsy_id = ifelse(
    biopsy_id == "H10-6806_FL", "H10-6806", biopsy_id
  )) %>%
  select(patient_id, biopsy_id, sample_id, seq_type,  
         matches("DLBCL90.*_call"), 
         matches("_ba|_cn"), 
         relapse_timing) %>% 
  full_join(select(clinical_md_biopsy, -sample_id)) %>% 
  rename(DTDX_bx = DTDX) %>% 
  full_join(clinical_md_patient) %>% 
  rename(DTDX_pt = DTDX)


first_dlbcl <- trios_md_less %>% 
  filter(!is.na(sample_id)) %>% 
  filter(pathology %in% c("DLBCL", "COMFL", "HGBL")) %>% 
  group_by(patient_id) %>% 
  slice_min(DTBX) %>% 
  ungroup() %>% 
  select(patient_id, index_DTBX = DTBX)

rel_times <- first_dlbcl %>% 
  left_join(clinical_md_patient) %>% 
  mutate(rel_DTDX = as.numeric(as.Date(DTDX) - as.Date(index_DTBX))/365, 
         rel_DTFUP = as.numeric(as.Date(DTFUP) - as.Date(index_DTBX))/365) %>% 
  select(patient_id, index_DTBX, rel_DTDX, rel_DTFUP)

relapse_timing <- trios_md_less %>% 
  filter(!is.na(sample_id)) %>% 
  select(patient_id, relapse_timing) %>% 
  distinct() %>% 
  mutate(relapse_timing = factor(relapse_timing, 
                                 levels = c("Primary Refractory", "Early Relapse", "Late Relapse"))) %>% 
  select(patient_id, relapse_timing) %>% 
  distinct()

rel_biopsies <- first_dlbcl %>%
  left_join(trios_md_less) %>% 
  select(-relapse_timing) %>% 
  mutate(rel_DTBX = as.numeric(as.Date(DTBX) - as.Date(index_DTBX))/365, 
         rel_DTTX = ifelse(!is.na(TX), rel_DTBX, NA)) %>% 
  mutate(path_summary = case_when(
    pathology %in% c("FL", "NS1") ~ "FL",
    # DLBCL90_dhitsig_call == "POS" & DLBCL90_dlbcl_call == "GCB" ~ "DHITsig+ DLBCL", 
    !is.na(DLBCL90_dlbcl_call) ~ str_c(DLBCL90_dlbcl_call, " DLBCL"), 
    is.na(pathology) ~ "DLBCL-NOS", # temporary for 00-26427 diagnostic
    pathology %in% c("NS2", "NS3", "DLBCL") ~ "DLBCL-NOS",
    pathology == "HDNOS" ~ "HGBL",
    pathology == "PROG" ~ "Progression",
    TRUE ~ pathology
  )) %>% 
  filter(pathology != "MCLD") %>% 
  mutate(sequencing = ifelse(!is.na(sample_id), TRUE, FALSE), 
         treatment = case_when(
           str_detect(TX, "HSCT") ~ "BMT", 
           !is.na(TX) ~ "Other"
         ),
         path_summary = factor(path_summary, 
                               levels = c(
                                 "GCB DLBCL", 
                                 "UNCLASS DLBCL", 
                                 "ABC DLBCL", 
                                 "DLBCL-NOS", 
                                 "HGBL", 
                                 "COMFL",
                                 "cHL", 
                                 "FL", 
                                 "MALT", 
                                 "MZL", 
                                 "SLL", 
                                 "Progression"
                               )))

pt_order <- rel_biopsies %>% 
  filter(sequencing) %>% 
  select(patient_id, rel_DTBX) %>% 
  group_by(patient_id) %>% 
  filter(rel_DTBX > 0) %>%
  slice_min(rel_DTBX) %>% 
  ungroup() %>% 
  mutate(patient_id = fct_reorder(patient_id, rel_DTBX, .desc = TRUE))

plot_data <- rel_times %>% 
  left_join(relapse_timing) %>% 
  left_join(rel_biopsies) %>% 
  filter(patient_id %in% pt_order$patient_id) %>% 
  mutate(patient_id = factor(patient_id, levels = levels(pt_order$patient_id))) 

swimmer_cols <- c(get_gambl_colours()[
  c("ABC", "DLBCL", "FL", "GCB", 
    "CLL", "UNCLASS", "MZL", "THRLBCL", "HGBL", "COMFL", "SCBC")
], "darkgrey")
names(swimmer_cols) <- c(
  "ABC DLBCL", 
  "DLBCL-NOS",
  "FL", 
  "GCB DLBCL", 
  "SLL", 
  "UNCLASS DLBCL",
  "MZL", 
  "MALT", 
  "HGBL", 
  "COMFL", 
  "cHL",
  "Progression"
)

swimmer_cols <- swimmer_cols[levels(plot_data$path_summary)]

extra_seg_pr <- plot_data %>% 
  filter(relapse_timing == "Primary Refractory") %>% 
  filter(rel_DTDX < -5 | rel_DTFUP > 5) %>% 
  mutate(rel_DTFUP = ifelse(rel_DTFUP > 5, rel_DTFUP - 2, NA), 
         rel_DTDX = ifelse(rel_DTDX < -5, rel_DTDX + 4, NA)) %>% 
  select(rel_DTFUP, rel_DTDX, patient_id, relapse_timing) %>% 
  distinct()

pr_plot <- plot_data %>% 
  filter(relapse_timing == "Primary Refractory") %>% 
  ggplot() +
  geom_vline(xintercept = 0, colour = "grey", size = 2) + 
  geom_segment(aes(x = ifelse(rel_DTDX < -5, -1, rel_DTDX), xend = ifelse(rel_DTFUP > 5, 4, rel_DTFUP), y = patient_id, yend = patient_id)) + 
  geom_segment(data = extra_seg_pr, aes(x = -1, xend = rel_DTDX, y = patient_id, yend = patient_id), linetype = "dashed") + 
  geom_segment(data = extra_seg_pr, aes(x = 4, xend = rel_DTFUP, y = patient_id, yend = patient_id), linetype = "dashed") +
  geom_point(aes(x = ifelse(rel_DTFUP > 5, rel_DTFUP - 2, rel_DTFUP), y = patient_id, shape = as.character(CODE_OS)), size = 4) +
  geom_point(aes(x = ifelse(rel_DTBX < -5, rel_DTBX + 4, rel_DTBX), y = patient_id, fill = path_summary, shape = sequencing), size = 5, stroke = 0.5) + 
  geom_point(aes(x = ifelse(rel_DTTX < -5, rel_DTTX + 4, rel_DTTX), y = patient_id, shape = treatment), size = 3) +
  scale_fill_manual(values = swimmer_cols, name = "") +
  scale_shape_manual(values = c("0" = 18, "1" = 4, "TRUE" = 22, "FALSE" = 21, "BMT" = 8, "Other" = 3), 
                     labels = c("Alive", "Deceased", "WGS/WES", "No sequencing", "HSCT", "Other Tx"), 
                     name = "") +
  facet_grid(rows = vars(relapse_timing), scales = "free_y", space = "free") +   
  ylab("") + xlab("") + 
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  theme(legend.position = "none", 
        panel.border = element_blank(), 
        panel.background = element_blank()) 

extra_seg_early <- data.frame(x = 7, xend = 9, y = "94-25764", relapse_timing = factor("Early Relapse", levels = c("Primary Refractory", "Early Relapse", "Late Relapse")))

early_plot <- plot_data %>% 
  filter(relapse_timing == "Early Relapse") %>% 
  ggplot() +
  geom_vline(xintercept = 0, colour = "grey", size = 2) + 
  geom_segment(aes(x = rel_DTDX, xend = ifelse(rel_DTFUP > 7, 7, rel_DTFUP), y = patient_id, yend = patient_id)) + 
  geom_segment(data = extra_seg_early, aes(x = x, xend = xend, y = y, yend = y), linetype = "dashed") + 
  geom_point(aes(x = ifelse(rel_DTFUP > 7, 9, rel_DTFUP), y = patient_id, shape = as.character(CODE_OS)), size = 4) +
  geom_point(aes(x = rel_DTBX, y = patient_id, fill = path_summary, shape = sequencing), size = 5, stroke = 0.5) + 
  geom_point(aes(x = rel_DTTX, y = patient_id, shape = treatment), size = 3) + 
  scale_fill_manual(values = swimmer_cols, name = "") +
  scale_shape_manual(values = c("0" = 18, "1" = 4, "TRUE" = 22, "FALSE" = 21, "BMT" = 8, "Other" = 3), 
                     labels = c("Alive", "Deceased", "WGS/WES", "No sequencing", "HSCT", "Other Tx"), 
                     name = "") +
  facet_grid(rows = vars(relapse_timing), scales = "free_y", space = "free") +   
  ylab("") + xlab("Time (years)") + 
  theme(legend.position = "none", 
        panel.border = element_blank(), 
        panel.background = element_blank()) 

late_plot <- plot_data %>% 
  filter(relapse_timing == "Late Relapse") %>% 
  ggplot() +
  geom_vline(xintercept = 0, colour = "grey", size = 2) + 
  geom_segment(aes(x = rel_DTDX, xend = rel_DTFUP, y = patient_id, yend = patient_id)) + 
  geom_point(aes(x = rel_DTFUP, y = patient_id, shape = as.character(CODE_OS)), size = 4) +
  geom_point(aes(x = rel_DTBX, y = patient_id, fill = path_summary, shape = sequencing), size = 5, stroke = 0.5) + 
  geom_point(aes(x = rel_DTTX, y = patient_id, shape = treatment), size = 3) + 
  scale_fill_manual(values = swimmer_cols, name = "") +
  scale_shape_manual(values = c("0" = 18, "1" = 4, "TRUE" = 22, "FALSE" = 21, "BMT" = 8, "Other" = 3), 
                     labels = c("Alive", "Deceased", "WGS/WES", "No sequencing", "HSCT", "Other Tx"), 
                     name = "") +
  facet_grid(rows = vars(relapse_timing), scales = "free_y", space = "free") + 
  ylab("") + xlab("Time (years)") + 
  guides(fill = guide_legend(override.aes = list(shape = 21)), 
         shape = guide_legend(override.aes = list(fill = "black"))) + 
  theme_cowplot()


path_legend <- cowplot::get_legend(late_plot + theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom") + guides(shape = "none"))
shape_legend <- cowplot::get_legend(late_plot + theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom") + guides(fill = "none"))

early_pr <- cowplot::plot_grid(pr_plot, early_plot, ncol = 1, rel_heights = c(5, 8))
plots <- cowplot::plot_grid(early_pr, late_plot + theme(legend.position = "none"))
legend <- cowplot::plot_grid(path_legend, shape_legend, nrow = 2)

cowplot::plot_grid(plots, legend, nrow = 2, rel_heights = c(10,1.5))  
ggsave("figures/swimmer_plot.pdf", height = 10, width = 10)

early_legend <- plot_grid(early_pr, legend, ncol = 1, rel_heights = c(10, 2))
cowplot::plot_grid(early_legend, late_plot + theme(legend.position = "none"), nrow = 1, rel_widths = c(1,1.2))
ggsave("figures/swimmer_plot_alt_layout.pdf", height = 8, width = 12)




##### Repeat but with broken axes #####

library(ggbreak)

pr_plot_breaks <- plot_data %>% 
  filter(relapse_timing == "Primary Refractory") %>% 
  ggplot() +
  geom_segment(aes(x = rel_DTDX, xend = rel_DTFUP, y = patient_id, yend = patient_id)) + 
  geom_point(aes(x = rel_DTFUP, y = patient_id, shape = as.character(CODE_OS)), size = 4) +
  geom_point(aes(x = rel_DTBX, y = patient_id, fill = path_summary, shape = sequencing), size = 5, stroke = 0.5) + 
  geom_point(aes(x = rel_DTTX, y = patient_id, shape = treatment), size = 3) +
  scale_fill_manual(values = swimmer_cols, name = "") +
  scale_shape_manual(values = c("0" = 18, "1" = 4, "TRUE" = 22, "FALSE" = 21, "BMT" = 8, "Other" = 3), 
                     labels = c("Alive", "Deceased", "WGS/WES", "No sequencing", "HSCT", "Other Tx"), 
                     name = "") +
  facet_grid(rows = vars(relapse_timing), scales = "free_y", space = "free") +   
  ylab("") + xlab("") + 
  theme(legend.position = "none", 
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  scale_x_continuous(breaks = c(-10, -5), limits = c(-10, 8)) + 
  scale_x_break(c(-0.25, -1), scales = 2, ticklabels = c(0, 0.5, 1, 1.5)) + 
  scale_x_break(c(1.9, 2.2), scales = 1, ticklabels = c(2, 4, 6, 8)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

early_plot_breaks <- plot_data %>% 
  filter(relapse_timing == "Early Relapse") %>% 
  ggplot() +
  geom_segment(aes(x = rel_DTDX, xend = rel_DTFUP, y = patient_id, yend = patient_id)) + 
  geom_point(aes(x = rel_DTFUP, y = patient_id, shape = as.character(CODE_OS)), size = 4) +
  geom_point(aes(x = rel_DTBX, y = patient_id, fill = path_summary, shape = sequencing), size = 5, stroke = 0.5) + 
  geom_point(aes(x = rel_DTTX, y = patient_id, shape = treatment), size = 3) + 
  scale_fill_manual(values = swimmer_cols, name = "") +
  scale_shape_manual(values = c("0" = 18, "1" = 4, "TRUE" = 22, "FALSE" = 21, "BMT" = 8, "Other" = 3), 
                     labels = c("Alive", "Deceased", "WGS/WES", "No sequencing", "HSCT", "Other Tx"), 
                     name = "") +
  facet_grid(rows = vars(relapse_timing), scales = "free_y", space = "free") +   
  ylab("") + xlab("Time (years)") + 
  theme(legend.position = "none") + 
  scale_x_continuous(limits = c(-5, 23), breaks = c(-4, -2)) +
  scale_x_break(c(-0.25, -1), scales = 2, ticklabels = c(0, 1, 2, 3, 4, 5)) +
  scale_x_break(c(5, 5.5), scales = 1, ticklabels = c(10, 15, 20)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


late_plot_breaks <- plot_data %>% 
  filter(relapse_timing == "Late Relapse") %>% 
  ggplot() +
  geom_vline(xintercept = 0, colour = "grey", size = 2) + 
  geom_segment(aes(x = rel_DTDX, xend = rel_DTFUP, y = patient_id, yend = patient_id)) + 
  geom_point(aes(x = rel_DTFUP, y = patient_id, shape = as.character(CODE_OS)), size = 4) +
  geom_point(aes(x = rel_DTBX, y = patient_id, fill = path_summary, shape = sequencing), size = 5, stroke = 0.5) + 
  geom_point(aes(x = rel_DTTX, y = patient_id, shape = treatment), size = 3) + 
  scale_fill_manual(values = swimmer_cols, name = "") +
  scale_shape_manual(values = c("0" = 18, "1" = 4, "TRUE" = 22, "FALSE" = 21, "BMT" = 8, "Other" = 3), 
                     labels = c("Alive", "Deceased", "WGS/WES", "No sequencing", "HSCT", "Other Tx"), 
                     name = "") +
  facet_grid(rows = vars(relapse_timing), scales = "free_y", space = "free") + 
  ylab("") + xlab("Time (years)") + 
  guides(fill = guide_legend(override.aes = list(shape = 21)), 
         shape = guide_legend(override.aes = list(fill = "black"))) + 
  theme_cowplot() + 
  xlim(-20,23) + 
  scale_x_cut(-5.25, which = 1, scales = 0.3)


early_pr_breaks <- cowplot::plot_grid(NULL, print(pr_plot_breaks), NULL,  print(early_plot_breaks), ncol = 1, rel_heights = c(-1, 12, -1, 16))
plots_breaks <- cowplot::plot_grid(print(early_pr_breaks), print(late_plot_breaks + theme(legend.position = "none")))

cowplot::plot_grid(plots_breaks, legend, nrow = 2, rel_heights = c(10,1.5))  
ggsave("figures/swimmer_plot_breaks.pdf", height = 10, width = 10)

early_legend <- plot_grid(print(early_pr_breaks), legend, ncol = 1, rel_heights = c(10, 2))
cowplot::plot_grid(early_legend, print(late_plot_breaks + theme(legend.position = "none")), nrow = 1, rel_widths = c(1,1.2))
ggsave("figures/swimmer_plot_breaks_alt_layout.pdf", height = 8, width = 12)