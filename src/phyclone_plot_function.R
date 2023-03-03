generate_plots <- function(patient, # The ID of the patient whose data should be plotted.
                           maf_dir, # The 00-inputs/maf/{seq_type}--{genome_build}/ directory in the PyClone module outputs
                           results_dir, # The 99-outputs/{seq_type}--{genome_build}/ directory in the PyClone module outputs
                           lymphgen_data, # Raw output of the LymphGen module
                           plot_dir, # Where plots should be written
                           biopsy_timing, # Data frame containing patient_id, Tumor_Sample_Barcode to match maf file, relapse_timing, and time_since_diagnosis_years
                           custom_genes # A vector of Hugo Symbols for genes to be labeled on the VAF plot
) {
  # For testing:
  # patient <- "06-11677"
  # maf_dir <-  "../LySeq_Validation/results/pyclone_vi-1.0/00-inputs/maf/genome--grch37/"
  # results_dir <- "../LySeq_Validation/results/pyclone_vi-1.0/99-outputs/genome--grch37/"
  # lymphgen_data <- trios_lymphgen
  # biopsy_timing <- biopsy_timing_all

  message(paste0("Analyzing patient ", patient))
  # Order the biopsies by dtbx
  biopsy_timing <- biopsy_timing %>%
    filter(patient_id == patient) %>%
    filter(!is.na(time_since_diagnosis_years)) %>%
    arrange(time_since_diagnosis_years) %>%
    mutate(sample_label = row_number()) %>%
    select(Tumor_Sample_Barcode, sample_label, relapse_timing)

  # Obtain the results files for this patient
  outfiles <- as.list(dir(results_dir, pattern = patient, full.names = TRUE))
  if (length(outfiles) == 0) {
    return()
  }
  names(outfiles) <- c("clusters", "stats", "tree", "pyclone")

  # Read in the PhyClone cluster file
  clusters <- read_tsv(outfiles$clusters, col_types = cols()) %>%
    filter(clone_id >= 0) %>%
    separate(
      mutation_id,
      into = c("Chromosome", "Start_Position", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"),
      sep = ":"
    ) %>%
    mutate(Start_Position = as.numeric(Start_Position)) %>%
    rename(Tumor_Sample_Barcode = sample_id) %>%
    group_by(Tumor_Sample_Barcode) %>%
    mutate(ccf = ccf / max(ccf)) %>%
    ungroup()

  # Obtain the full paths to the maf for this patient
  patient_mafs <- dir(paste0(maf_dir, patient), pattern = ".maf", full.names = TRUE)

  patient_maf <- lapply(patient_mafs, read_tsv, col_types = cols(
    SOMATIC = col_double(),
    PHENO = col_double(),
    gnomADg_AF = col_double()
  )) %>%
    bind_rows() %>%
    mutate(patient_id = patient) %>%
    left_join(clusters) %>%
    left_join(biopsy_timing) %>%
    filter(!is.na(cluster_id))


  # Get the mutations that are relevant to the LymphGen classification

  lymphgen <- lymphgen_data %>%
    filter(Sample.Name %in% unique(patient_maf$Tumor_Sample_Barcode)) %>%
    GAMBLR::tidy_lymphgen(lymphgen_column_out = "LymphGen", relevel = TRUE)

  lg_classes <- lymphgen$Subtype.Prediction %>% as.list()
  lg_classes <- lapply(lg_classes, function(x) {
    str_split(x, "/")
  }) %>%
    unlist() %>%
    unique()
  lg_feats <- c()
  if (
    sum(lg_classes == "Other") < length(lg_classes) &
      sum(lg_classes == "A53") < length(lg_classes)
  ) {
    lg_classes <- lg_classes %>%
      paste0(., collapse = "|")
    lg_feats <- lymphgen %>%
      select(matches(lg_classes), -matches("Confidence|Count")) %>%
      pivot_longer(
        everything(),
        names_to = "class",
        values_to = "feats"
      ) %>%
      filter(!is.na(feats)) %>%
      pull(feats) %>%
      as.list()
    message(paste0(unlist(lg_feats)))

    lg_feats <- lapply(lg_feats, function(x) {
      str_split(x, ",")
    }) %>%
      unlist() %>%
      str_remove_all(., "_.*") %>%
      str_remove_all(., "L265P") %>%
      str_remove_all(., "Fusion.*") %>%
      str_remove_all(., "Trans.*") %>%
      unique()
  } else {
    message("Using GAMBLR lymphoma genes to label VAF plot. ")
    lg_feats <- GAMBLR::lymphoma_genes$Gene
  }
  if (length(lg_feats) < 10 & missing(custom_genes)) {
    message("Using GAMBLR lymphoma genes to label VAF plot. ")
    lg_feats <- GAMBLR::lymphoma_genes$Gene
  } else if (!missing(custom_genes)) {
    message("Using a custom gene list to label VAF plot. ")
    lg_feats <- custom_genes
  } else {
    message("Using LymphGen feature genes to label VAF plot. ")
  }

  patient_maf <- patient_maf %>%
    mutate(vaf = t_alt_count / t_depth) %>%
    mutate(vaf = ifelse(vaf < 0.05, 0, vaf)) %>%
    mutate(hot_spot = case_when(
      Hugo_Symbol == "EZH2" & str_detect(HGVSp_Short, "p.Y646") ~ TRUE,
      Hugo_Symbol == "MYD88" & HGVSp_Short == "p.L265P" ~ TRUE,
      Hugo_Symbol == "CD79B" & str_detect(HGVSp_Short, "p.Y197") ~ TRUE,
      Hugo_Symbol %in% c("NOTCH1", "NOTCH2") & Exon_Number == "34/34" & str_detect(Variant_Classification, ("Nonsense|Frame_Shift")) ~ TRUE,
      Hugo_Symbol == "CREBBP" & Start_Position > GAMBLR::hotspot_regions_grch37["CREBBP", "start"] & End_Position < GAMBLR::hotspot_regions_grch37["CREBBP", "end"] ~ TRUE,
      Hugo_Symbol == "MEF2B" & str_detect(HGVSp_Short, "D83") ~ TRUE,
      Hugo_Symbol == "STAT6" & str_detect(HGVSp_Short, "D419") ~ TRUE
    )) %>%
    mutate(label = case_when(
      hot_spot ~ str_c(Hugo_Symbol, "*"),
      Variant_Classification %in% GAMBLR:::coding_class ~ Hugo_Symbol
    )) %>%
    group_by(label) %>%
    mutate(label = ifelse(
      Hugo_Symbol %in% lg_feats & sum(vaf) > 0,
      label,
      NA
    )) %>%
    ungroup() %>%
    left_join(select(
      lymphgen,
      Tumor_Sample_Barcode = Sample.Name,
      Subtype.Prediction
    )) %>%
    mutate(Subtype.Prediction = replace_na(Subtype.Prediction, "ND")) %>%
    mutate(
      sample_label = str_c(
        str_remove(sample_label, "-.*"),
        "-",
        Subtype.Prediction
      )
    ) %>%
    mutate(
      sample_label_x = str_c(
        str_remove(sample_label, "-.*"),
        "\n",
        str_replace_all(Subtype.Prediction, "/", "\n")
      )
    )

  # Obtain an edge table from the newick tree

  tree <- ape::read.tree(outfiles$tree)

  if (length(tree$edge[, 1]) > 1) {
    treeplot <- ggtree::ggtree(tree, aes(fill = as.character(node))) +
      ggtree::geom_nodelab(node = "all", geom = "label", colour = "white") +
      # geom_tiplab(geom = "label") +
      coord_cartesian(clip = "off") +
      theme(legend.position = "none")

    g <- ggplot_build(treeplot)
    nodecols <- g$data[[3]][c("label", "node", "fill")]
    # nodecols$fill <- ggsci::pal_d3("category10")(10)[1:length(nodecols$node)]
    nodecols[nodecols$label == "root", ]$fill <- "darkgrey"
    # tipcols <- g$data[[4]][c("label", "fill")]
    # allcols <- bind_rows(nodecols, tipcols)
    tree_cols <- nodecols$fill
    names(tree_cols) <- nodecols$node

    treeplot <- treeplot + scale_fill_manual(values = tree_cols)
  }


  # CCF plots

  # Obtain the clone_ids for clones that are shared between time points

  shared_clones <- patient_maf %>%
    group_by(clone_id) %>%
    mutate(shared = ifelse(min(ccf) > 0.1, TRUE, FALSE)) %>%
    filter(shared) %>%
    pull(clone_id) %>%
    unique()


  if (length(unique(patient_maf$Tumor_Sample_Barcode)) > 2) {
    label_data <- patient_maf %>%
      filter(!is.na(label)) %>%
      group_by(Hugo_Symbol, HGVSp_Short) %>%
      filter(sample_label_x %in% c(min(sample_label_x), max(sample_label_x))) %>%
      slice_max(vaf, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      group_by(label, clone_id) %>%
      slice_max(vaf, n = 1, with_ties = FALSE) %>%
      ungroup()

    vafplot <- patient_maf %>%
      filter(!is.na(sample_label)) %>%
      ggplot(aes(
        x = sample_label_x,
        y = vaf,
        group = str_c(Chromosome, Start_Position, Tumor_Seq_Allele1, Tumor_Seq_Allele2),
        colour = as.character(clone_id)
      )) +
      geom_point(alpha = 0.6) +
      geom_line(aes(alpha = !is.na(label))) +
      geom_label_repel(
        data = filter(label_data, sample_label_x == max(sample_label_x)),
        aes(
          label = label,
          x = sample_label_x,
          y = vaf
        ),
        nudge_x = 0.5,
        direction = "y",
        size = 4,
        hjust = "left"
      ) +
      geom_label_repel(
        data = filter(label_data, sample_label_x == min(sample_label_x)),
        aes(
          label = label,
          x = sample_label_x,
          y = vaf
        ),
        nudge_x = -0.5,
        direction = "y",
        size = 4,
        hjust = "right"
      ) +
      scale_x_discrete(
        expand = expansion(mult = 0.5)
      ) +
      scale_alpha_manual(values = c(0.2, 1)) +
      xlab("") +
      ylab("VAF") +
      theme_cowplot() +
      theme(legend.position = "none")
  }

  if (length(unique(patient_maf$Tumor_Sample_Barcode)) == 2) {
    vafplot_data <- patient_maf %>%
      filter(!is.na(sample_label)) %>%
      select(
        sample_label,
        clone_id,
        label,
        Hugo_Symbol,
        HGVSp_Short,
        vaf
      ) %>%
      distinct(Hugo_Symbol, HGVSp_Short, sample_label, .keep_all = TRUE) %>%
      pivot_wider(
        names_from = sample_label,
        values_from = vaf
      )

    set.seed(123)
    label_data <- vafplot_data %>%
      filter(!is.na(label)) %>%
      group_by(label, clone_id) %>%
      slice_sample(n = 1) %>%
      ungroup()

    vafplot <- vafplot_data %>%
      ggplot(aes(
        x = !!sym(colnames(.)[5]),
        y = !!sym(colnames(.)[6]),
        colour = as.character(clone_id)
      )) +
      geom_point(alpha = 0.6) +
      geom_label_repel(
        data = label_data,
        aes(label = label),
        max.overlaps = Inf,
        size = 4,
        nudge_x = 0.1,
        nudge_y = 0.1,
        show.legend = FALSE
      ) +
      theme_cowplot() +
      theme(legend.position = "none")
  }

  cluster_summary_phy <- patient_maf %>%
    select(clone_id, sample_label_x, ccf) %>%
    distinct() %>%
    select(
      timepoint = sample_label_x,
      clone_id,
      clonal_prev = ccf
    )


  ccf_plot <- cluster_summary_phy %>%
    filter(!is.na(timepoint)) %>%
    ggplot(aes(
      x = timepoint,
      y = clonal_prev,
      group = clone_id,
      colour = as.character(clone_id)
    )) +
    geom_point() +
    geom_line() +
    scale_colour_discrete(name = "Clone ID") +
    xlab("") +
    ylab("CCF") +
    theme_cowplot()

  clone_size <- patient_maf %>%
    filter(vaf > 0.05) %>%
    group_by(sample_label_x) %>%
    count(clone_id) %>%
    mutate(percent = n / sum(n)) %>%
    ungroup() %>%
    filter(clone_id %in% shared_clones) %>%
    ggplot(aes(
      x = sample_label_x,
      y = percent,
      fill = fct_reorder(factor(clone_id), desc(clone_id))
    )) +
    geom_col(width = 0.5) +
    ylab("Fraction of mutations shared") +
    ylim(0, 1) +
    theme_cowplot() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none"
    )

  title <- ggdraw() +
    draw_label(
      unique(patient_maf$relapse_timing),
      x = 0,
      angle = 90,
      size = 16
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )

  if (length(tree$edge[, 1]) > 1) {
    col_vec <- nodecols$fill
    names(col_vec) <- nodecols$label
    cowplot::plot_grid(
      title,
      ccf_plot + scale_colour_manual(values = col_vec) + theme(legend.position = "none"),
      vafplot + scale_colour_manual(values = col_vec),
      clone_size + scale_fill_manual(values = col_vec) + theme(legend.position = "none"),
      treeplot,
      nrow = 1,
      rel_widths = c(0.1, 1, 1.2, 1, 1),
      align = "h",
      axis = "b"
    )
    ggsave(paste0(plot_dir, "/", patient, ".phyclone.pdf"), height = 3.5, width = 14)
  } else {
    col_vec <- pull(unique(ggplot_build(ccf_plot)$data[[1]]["colour"]))
    names(col_vec) <- unique(cluster_summary_phy$clone_id[!is.na(cluster_summary_phy$clone_id)])
    legend <- get_legend(
      ccf_plot +
        scale_colour_manual(values = col_vec, name = "Clone ID")
    )
    cowplot::plot_grid(
      title,
      ccf_plot + theme(legend.position = "none") + scale_colour_manual(values = col_vec),
      vafplot + scale_colour_manual(values = col_vec),
      clone_size + theme(legend.position = "none") + scale_fill_manual(values = col_vec),
      legend,
      nrow = 1,
      rel_widths = c(0.1, 1, 1.2, 1, 0.3),
      align = "h",
      axis = "b"
    )
    ggsave(paste0(plot_dir, "/", patient, ".phyclone.pdf"), height = 3.5, width = 10)
  }

  patient_maf_coding <- patient_maf %>%
    filter(Variant_Classification %in% GAMBLR:::coding_class) %>%
    select(-sample_label)

  return(patient_maf_coding)
}
