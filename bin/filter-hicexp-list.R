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
########################## functions ###################################
# it's dangerous to go alone! take this.
#
# obtain significant pairs of interactions.
obtain_sigpairs <- function(comparison, p) {
  topDirs(hicexp = comparison, p.adj_cutoff = p, return_df = "pairedbed")
}
#
# plot sigpairs by chromosome
plot_sigpairs_chr <- function(sigpairs, p, outdir) {
  
  c <- sigpairs.list[[sigpairs]]$chr1
  c <- replace(c, c=="chr23", "chrX")
  c <- factor(c, levels = c(paste0("chr", seq(1,22,1)), "chrX"))
  
  png(file = paste0(outdir, "/", sigpairs, "-barplot-chrs.png"), 
      height = 10, width = 10, units = "in", res = 300)
  
  barplot(rev(table(c)), horiz = TRUE, las=1, border=F)

  dev.off()
  
}
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
#
opt <- parse_args(OptionParser(option_list=option_list))
#
if (is.null(opt$input)){
  print_help(OptionParser(option_list=option_list))
  stop("The input file is mandatory.n", call.=FALSE)
}
#
hicexp.comparison.list <- readRDS(opt$input)
#
################### filter significant interactions #######################
#
sigpairs.list <- lapply(hicexp.comparison.list, obtain_sigpairs, p = opt$pvalue)
#
names(sigpairs.list) <- gsub("qlf", "sig", names(sigpairs.list))
#
############################# postprocessig ###############################
# all the battery of plots.
a <- sigpairs.list[[1]]

# how many pairs by chromosome (you'll have to make the chromosomes a factor, should use a separate vector)

# median distance and distance distribution

# logFC 

# distance AND logFC 
# distanca AND logFC by CHROMOSOME (gglplot)


# maybe a circos by chromosome? (DIFFERENT MODULE.)

########### THE BETWEEN-COMPARISONS PLOTS require a different script that I can source here.
#(perhaps this should be a different module)
# boxplots comparing logFC and Distance between comparisons. 





################  then you save the object and that's that ################
saveRDS(sigpairs.list, file = opt$output)
