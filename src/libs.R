# Load libraries used in pretty much all Rscripts.

# First source .Rprofile so Snakemake an make use of renv
source("renv/activate.R")

# GAMBLR-related
library(GAMBLR)
library(dbplyr)
library(tidyverse)
library(data.table)
library(RMariaDB)
library(DBI)
library(stats)
library(metaviz)

# Other data
library(readxl)


# Visualization
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
library(ggh4x)
library(ggalluvial)


# Parallelization
library(BiocParallel)
register(MulticoreParam(48))

theme_set(theme_cowplot())

# Patients excluded after review for having disease other than DLBCL
pts_exclude <- c(
  "03-19103", # Composite SLL/cHL at one time point
  "92-38626", # PTLD
  "04-13783", # Cutaneous leg-type DLBCL
  "12-17272", # PCNSL
  "14-19181", # PCNSL
  "06-16041",  # Primary testicular with CNS relapse 
  "CAMP0003" # PMBCL
)

# Common colour palette

relapse_colours <- c("#8681BD", "#34A9E0", "#48AF64")
names(relapse_colours) <- c("Primary Refractory", "Early Relapse", "Late Relapse")

relapse_timing_labels <- c("Primary\nRefractory", "Early\nRelapse", "Late\nRelapse")

relapse_strips <- strip_themed(
  background_x = elem_list_rect(fill = relapse_colours)
)

coo_colours <- get_gambl_colours(classification = "COO")
coo_colours <- coo_colours[c(
  "ABC",
  "GCB",
  "UNCLASS",
  "GCB",
  "DHITsig-IND",
  "DHITsig+"
)]
names(coo_colours)[4] <- "DHITsig-"

# NanoString cutoff values for COO and DHITsig classification
dhitsig_cuts <- c(-15.32, -6.33)
coo_cuts <- c(1908, 2440)
