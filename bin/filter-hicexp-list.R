#!/usr/bin/env Rscript
#
# hreyes Feb 2020
# filter-hicexp-list.R
#
# Read in a list of normalized and compared hicexp objects
# and filter the results using cutoffs specified by the user.
#
# The topDirs function helps to filter out less interesting regions and
# retain significant ones.
#################### import libraries and set options ####################
library(multiHiCcompare)
library(optparse)
#
#options(scipen = 10)
#
########################## read in data ###################################
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

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)){
  print_help(OptionParser(option_list=option_list))
  stop("The input file is mandatory.n", call.=FALSE)
}

hicexp.comparison.list <- readRDS(opt$input)

################### filter significant interactions #######################

top.qlf.MCF10AT1.MCF10A <- topDirs(qlf.MCF10AT1.MCF10A, return_df = "pairedbed", p.adj_cutoff = 0.05)
top.qlf.MCF10CA1A.MCF10A <- topDirs(qlf.MCF10CA1A.MCF10A, return_df = "pairedbed", p.adj_cutoff = 0.05)



significant.pairs <- topDirs(hicexp = input.hicexp, 
                                    p.adj_cutoff = opt$pvalue, 
                                    logfc_cutoff = opt$logFC, 
                                    D_cutoff = opt$distance,
                                    return_df = "pairedbed")


################  then you save the object and that's that ################
saveRDS(outfile, file = out_bedpe_path)
