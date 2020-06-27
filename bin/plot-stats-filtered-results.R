#!/usr/bin/env Rscript
#
# hreyes June 2020
# plot-stats-filtered-results.R
#######################################################################
# Read in a list of filtered hicexp qlf-compared results and plot
# descriptive stats.
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
###################### plot sigpairs by chromosome #####################
plot_chromosomes <- function(sigpairs, r, outdir) {
  # get chromosomes and logFC data
  chr.fc <- sigpairs.list[[sigpairs]][ , c("chr1", "logFC")]
  
  ####### prepare data
  # remap chrX
  chr.fc$chr1 <- replace(chr.fc$chr1, chr.fc$chr1=="chr23", "chrX")
  # convert to factor so it looks ordered in the plot
  chr.fc$chr1 <- factor(chr.fc$chr1, levels = c(paste0("chr", seq(1,22,1)), "chrX"))
  chr.fc$chr1 <- droplevels(chr.fc$chr1)
  
  # remap logFC to a factor
  chr.fc$logFC <- unlist(lapply(chr.fc$logFC, function(x) {
    if(x < 0) { return("Negative")
    } else if(x > 0) { return("Positive") }
  }))
  
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
# #
# ###################### plot distance distribution #########################

# PLOT OF HISTOGRAM OF DISTANCE

# PLOT OF HISTOGRAM OF DISTANCE BY CHROMOSOME

ggplot(test, aes(chr1, D)) + geom_boxplot()


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
#lapply(names(sigpairs.list), plot_distance_distribution, r = hicres, u = hicunit, outdir= dirname(opt$output))
#
# Plot distance vs logFC
#lapply(names(sigpairs.list), plot_distance_vs_logfc, r = hicres, u = hicunit, outdir= dirname(opt$output))


# maybe a circos by chromosome? (DIFFERENT MODULE.)
# or one of those diagrams where the chromosome's contacts are shown

########### THE BETWEEN-COMPARISONS PLOTS 
# ggplot comparing significant interactions neg and pos by chr
# ggplots comparing logFC and Distance between comparisons. 

