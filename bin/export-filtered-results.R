#!/usr/bin/env Rscript
#
# hreyes June 2020
# export-filtered-results.R
#######################################################################
# Read in a list of filtered hicexp qlf-compared results and export
# to plain text.
########################################################################
# output fixed numbers, not scientific notation
options(scipen = 10)
#
# input options
args = commandArgs(trailingOnly=TRUE)
########################## functions ###################################
# it's dangerous to go alone! take this.
#
############### export significant pairs of interactions. ##############
export_pairs <- function(comparison, outdir) {
  # bring in data frame
  outfile <- sigpairs.list[[comparison]]
  
  # cast pvalue columns as numeric
  outfile$p.value <- as.numeric(outfile$p.value)
  outfile$p.adj <- as.numeric(outfile$p.adj)
  
  # obtain resolution
  hicres <-  unique(outfile$end1 - outfile$start1) + 1
  
  # build name for output text file
  outname <- paste0(outdir, "/", gsub("\\.", "_", comparison), "_", hicres, ".csv")
  
  # export csv
  write.csv(x = outfile, file = outname, quote = FALSE, row.names = FALSE)
}

######################### read in data #################################
sigpairs.list <- readRDS(args[1])
#
########################### export sigpairs ############################
lapply(names(sigpairs.list), export_pairs, outdir = dirname(args[1]))
