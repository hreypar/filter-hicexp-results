#!/usr/bin/env Rscript
#
# hreyes June 2020
# plot-stats-filtered-results.R
#######################################################################
# Read in a list of filtered hicexp qlf-compared results and plot
# descriptive stats.
#
# Importantly, bear in mind that this plots work because the data
# corresponds to INTRAchromosomal interactions (i.e. chr1 == chr2)
########################################################################
#
#################### import libraries and set options ##################
#suppressMessages(library(multiHiCcompare))
library(ggplot2)
#
options(scipen = 10)
#
# input options
args = commandArgs(trailingOnly=TRUE)
########################## functions ###################################
# it's dangerous to go alone! take this.
#
####################### format chromosome names ########################
format_chromosomes <- function(chrs) {
  chrs <- replace(x = chrs, chrs == "chr23", "chrX")
  chrs <- factor(chrs, levels = c(paste0("chr", seq(1,22,1)), "chrX"))
  chrs <- droplevels(chrs)
  return(chrs)
}
#
####################### map logFC values to factor #####################
logfc_to_category <- function(l) {
    if(l < 0) {
      return("Negative")
    } else if(l > 0) {
    return("Positive") 
    }
}
#
###################### plot sigpairs by chromosome #####################
plot_chromosomes <- function(sigpairs, r, outdir) {
  ####### get chromosomes and logFC data
  chr.fc <- sigpairs.list[[sigpairs]][ , c("chr1", "logFC")]
  
  ####### prepare data
  # remap chr23 to chrX
  chr.fc$chr1 <- format_chromosomes(chr.fc$chr1)
  # remap logFC to a factor
  chr.fc$logFC <- unlist(lapply(chr.fc$logFC, logfc_to_category))
  
  # final data frame
  t <- as.data.frame(table(chr.fc$logFC, chr.fc$chr1))
  colnames(t) <- c("logFC", "Chromosome", "DIRs")
  
  t$logFC <- factor(t$logFC, levels = c("Positive", "Negative"))
  
  ####### prepare variables for plot
  t.main = paste0("Differentially Interacting Regions ",
                  gsub("\\.", " vs ", gsub("sig.", "", sigpairs)), 
                  "\n(Hi-C resolution ", r, ")")

  ####### generate file
  png(file = paste0(outdir, "/", sigpairs, "-barplot-chrs-", r, ".png"),
      height = 9, width = 15, units = "in", res = 300)
  
  ####### plot
  chrs.plot <- ggplot(t, aes(Chromosome, DIRs, fill=logFC)) + 
    geom_bar(stat="identity", col="gray69", position=position_dodge()) + 
    scale_fill_manual(values=c('#58a4b0','#EE1B49')) +
    theme_minimal() + ggtitle(t.main) +
    geom_text(aes(label = DIRs), vjust=-0.25, position = position_dodge(0.9))
    
  print(chrs.plot)
  
  dev.off()
}
#
####################### plot distance distribution #########################
plot_distance_distribution <- function(sigpairs, r, outdir) {
  ####### get chromosomes, distance and logFC data
  distance.fc <- sigpairs.list[[sigpairs]][,c("chr1", "D", "logFC")]
  
  ####### prepare data
  # remap chr23 to chrX
  distance.fc$chr1 <- format_chromosomes(distance.fc$chr1)
  # remap logFC to a factor
  distance.fc$logFC <- unlist(lapply(distance.fc$logFC, logfc_to_category))
  distance.fc$logFC <- factor(distance.fc$logFC, levels = c("Positive", "Negative"))
  
  ####### prepare variables for plot
  t.main = paste0("Distribution of Chromosomal Distance Between Differentially Interacting Regions ",
         gsub("\\.", " vs ", gsub("sig.", "", sigpairs)),
         "\n(Hi-C resolution ", r, ")")
  
  ####### plot1
  # generate file
  # png(filename = paste0(outdir, "/", sigpairs, "-hist-distance-", r, u, ".png"),
  #     height = 10, width = 13, units = "in", res = 300)
  
  plot1 <- ggplot(distance.fc, aes(logFC, D, fill=logFC)) + 
    geom_boxplot() + scale_fill_manual(values=c('#58a4b0','#EE1B49')) +
    theme_minimal() + facet_wrap(~chr1, scales="free") +
    ggtitle(t.main) + ylab("Distance (bins) between DIRs \n") + xlab("\nlog fold-change")
  
  #dev.off()
  
  ####### plot2
  # generate file
  # png(filename = paste0(outdir, "/", sigpairs, "-hist-distance-", r, u, ".png"),
  #     height = 10, width = 13, units = "in", res = 300)
  # 
  
  # PLOT2 is boxplot using chr, D; coloring by logFC (no facet wrap)
  
  #dev.off()
}
# #
# ###################### plot distance vs logFC #############################
# plot_distance_vs_logfc <- function(sigpairs, r, u, outdir) {
#   # prepare data
#   hicpairs <- sigpairs.list[[sigpairs]]
#   colnames(hicpairs)[colnames(hicpairs) == "chr1"] <- "chr"
#   hicpairs$chr <- factor(hicpairs$chr, levels = c(paste0("chr", seq(1,22,1)), "chrX"))
#   hicpairs$chr <- droplevels(hicpairs$chr)
#   
#   # create captions
#   xcaption = paste0("\nDistance (", u, ")")
#   ycaption = "logFC\n"
#   maincaption = paste("log fold change difference between", 
#                       gsub("\\.", " and ", gsub("sig.", "", sigpairs)),
#                       "by distance\n")
#   
#   # plot distance vs logfc
#   png(filename = paste0(outdir, "/", sigpairs, "-distance-vs-logfc-", r, u, ".png"),
#       height = 10, width = 13, units = "in", res = 300)
#   
#   plot1 <- ggplot(hicpairs, aes(x=D*r, y=logFC)) + 
#     geom_point(alpha=0.75, shape=19) +
#     xlab(xcaption) + ylab(ycaption) + ggtitle(maincaption) +
#     theme_minimal()
#   
#   print(plot1)
#   dev.off()
#   
#   # plot distance vs logfc by chromosome
#   png(filename = paste0(outdir, "/", sigpairs, "-distance-vs-logfc-by-chromosome", r, u, ".png"),
#       height = 10, width = 13, units = "in", res = 300)
#   
#   plot2 <- ggplot(hicpairs, aes(x=D*r, y=logFC)) + 
#     geom_point(alpha=0.75, shape=19) +
#     xlab(xcaption) + ylab(ycaption) + ggtitle(maincaption) +
#     facet_wrap(~chr)
#   
#   print(plot2)
#   dev.off()
# }
#
######################### read in data #################################
sigpairs.list <- readRDS(args[1])
#
#################### resolution of the data ############################
hicres = unique(sigpairs.list[[1]]$end1 - sigpairs.list[[1]]$start1) + 1
if(hicres >= 1000000) {
  hicunit = paste0(hicres/1000000, "Mb")
} else {
  hicunit = paste0(hicres/1000, "kb")
}

######################### produce plots #################################
# Plot DIRs by chromosome
lapply(names(sigpairs.list), plot_chromosomes, r = hicunit, outdir = dirname(args[1]))
#
# Plot distance distribution
#lapply(names(sigpairs.list), plot_distance_distribution, r = hicunit, outdir= dirname(args[1]))
#
# Plot distance vs logFC
#lapply(names(sigpairs.list), plot_distance_vs_logfc, r = hicres, u = hicunit, outdir= dirname(opt$output))


# maybe a circos by chromosome? (DIFFERENT MODULE.)
# or one of those diagrams where the chromosome's contacts are shown

########### THE BETWEEN-COMPARISONS PLOTS 
# ggplot comparing significant interactions neg and pos by chr
# ggplots comparing logFC and Distance between comparisons. 

