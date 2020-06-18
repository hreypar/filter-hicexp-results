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
library(multiHiCcompare)
library(optparse)
#
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
###################### plot sigpairs by chromosome #####################
plot_chromosomes <- function(sigpairs, r, outdir) {
  # get chromosomes data
  chrs <- sigpairs.list[[sigpairs]]$chr1
  chrs <- replace(chrs, chrs=="chr23", "chrX")
  chrs <- factor(chrs, levels = c(paste0("chr", seq(1,22,1)), "chrX"))
  chrs <- droplevels(chrs)
  # generate file
  png(file = paste0(outdir, "/", sigpairs, "-barplot-chrs-", r, ".png"), 
      height = 10, width = 13, units = "in", res = 300)
  # set margins
  par(mar=c(5,7,5,2))
  # plot
  barplot(rev(table(chrs)), horiz = TRUE, las=1, border=F,
          main = paste("Differentially Interacting Regions", 
                       gsub("\\.", " vs ", gsub("sig.", "", sigpairs))),
          xlab = paste0("Number of DIRs (", r, ")" ), ylab = "Chromosome\n")
  
  dev.off()
}
#
###################### plot distance distribution #########################
plot_distance_distribution <- function(sigpairs, r, u, outdir) {
  # get distance data
  distance <- sigpairs.list[[sigpairs]]$D * r
  # generate file
  png(filename = paste0(outdir, "/", sigpairs, "-hist-distance-", r, u, ".png"),
      height = 10, width = 13, units = "in", res = 300)
  # set margins
  par(mar=c(5,7,5,2))
  # plot
  hist(distance, las=1, col="gray83", 
       main = paste("Distance Between Differentially Interacting Regions\n",
                    gsub("\\.", " vs ", gsub("sig.", "", sigpairs))),
       xlab=paste0("Distance (", u, ") between DIRs" ))
  
  dev.off()
}
#
###################### plot distance vs logFC #############################
plot_distance_vs_logfc <- function(sigpairs, r, u, outdir) {
  # prepare data
  hicpairs <- sigpairs.list[[sigpairs]]
  colnames(hicpairs)[colnames(hicpairs) == "chr1"] <- "chr"
  hicpairs$chr <- factor(hicpairs$chr, levels = c(paste0("chr", seq(1,22,1)), "chrX"))
  hicpairs$chr <- droplevels(hicpairs$chr)
  
  # create captions
  xcaption = paste0("\nDistance (", u, ")")
  ycaption = "logFC\n"
  maincaption = paste("log fold change difference between", 
                      gsub("\\.", " and ", gsub("sig.", "", sigpairs)),
                      "by distance\n")
  
  # plot distance vs logfc
  png(filename = paste0(outdir, "/", sigpairs, "-distance-vs-logfc-", r, u, ".png"),
      height = 10, width = 13, units = "in", res = 300)
  
  plot1 <- ggplot(hicpairs, aes(x=D*r, y=logFC)) + 
    geom_point(alpha=0.75, shape=19) +
    xlab(xcaption) + ylab(ycaption) + ggtitle(maincaption) +
    theme_minimal()
  
  print(plot1)
  dev.off()
  
  # plot distance vs logfc by chromosome
  png(filename = paste0(outdir, "/", sigpairs, "-distance-vs-logfc-by-chromosome", r, u, ".png"),
      height = 10, width = 13, units = "in", res = 300)
  
  plot2 <- ggplot(hicpairs, aes(x=D*r, y=logFC)) + 
    geom_point(alpha=0.75, shape=19) +
    xlab(xcaption) + ylab(ycaption) + ggtitle(maincaption) +
    facet_wrap(~chr)
  
  print(plot2)
  dev.off()
}
#
############################## MAIN CODE ##################################
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
#
######################### descriptive plots ###############################
# resolution of the data
hicres = resolution(hicexp.comparison.list[[1]])
if(hicres >= 1000000) {
  hicres = hicres/1000000
  hicunit = "Mb"
} else {
  hicres = hicres/1000
  hicunit = "kb"
}
#
# Plot DIRs by chromosome
lapply(names(sigpairs.list), plot_chromosomes, r = paste0(hicres, hicunit), outdir = dirname(opt$output))
#
# Plot distance distribution
lapply(names(sigpairs.list), plot_distance_distribution, r = hicres, u = hicunit, outdir= dirname(opt$output))
#
# Plot distance vs logFC
lapply(names(sigpairs.list), plot_distance_vs_logfc, r = hicres, u = hicunit, outdir= dirname(opt$output))


# maybe a circos by chromosome? (DIFFERENT MODULE.)
# or one of those diagrams where the chromosome's contacts are shown

########### THE BETWEEN-COMPARISONS PLOTS require a different script that I could source here.
#(perhaps this should be a different module)
# ggplots comparing logFC and Distance between comparisons. 

