#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
format_base = args[1]
rmarkdown::render("manuscript.Rmd",output_format=paste0("bookdown::",format_base,"_document2"))
