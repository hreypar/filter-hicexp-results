#!/usr/bin/env Rscript
#
# hreyes Feb 2020
# filter-hicexp-list.R
#######################################################################
# Read in a list of normalized and compared hicexp objects
# and filter the results using specified cutoffs.
#
# The topDirs function helps to filter out less interesting regions and
# retain significant ones.
########################################################################
#
#################### import libraries and set options ##################
suppressMessages(library(multiHiCcompare))
library(optparse)
#library(ggplot2)
#
options(scipen = 10)
######################### Create options ###############################
option_list = list(
  make_option(opt_str = c("-i", "--input"),
              type = "character",
              help = "Rds file with a list of normalized and compared hicexp objects"),
  make_option(opt_str = c("-p", "--pvalue"), 
              type = "numeric", 
              default = 0.01, 
              help = "Adjusted p value cutoff to filter interactions. Default is 0.01"),
  make_option(opt_str = c("-o", "--output"),
              type = "character",
              help = "output file")
)
#
opt <- parse_args(OptionParser(option_list=option_list))
#
if (is.null(opt$input)){
  print_help(OptionParser(option_list=option_list))
  stop("The input file is mandatory.n", call.=FALSE)
}
#
########################## functions ###################################
# it's dangerous to go alone! take this.
#
############### obtain significant pairs of interactions. ##############
obtain_sigpairs <- function(comparison, p) {
  topDirs(hicexp = comparison, p.adj_cutoff = p, return_df = "pairedbed")
}
#
############################ read in data #################################
hicexp.comparison.list <- readRDS(opt$input)
#
################### filter significant interactions #######################
sigpairs.list <- lapply(hicexp.comparison.list, obtain_sigpairs, p = opt$pvalue)
#
names(sigpairs.list) <- gsub("qlf", "sig", names(sigpairs.list))
#
############################# save ouput ##################################
saveRDS(sigpairs.list, file = opt$output)
