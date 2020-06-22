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
export_pairs <- function(comparison) {
  
  write.
  sigpairs.list[[comparison]]
}

############################ read in data #################################
sigpairs.list <- readRDS(args[1])
#
############################# save ouput ##################################
