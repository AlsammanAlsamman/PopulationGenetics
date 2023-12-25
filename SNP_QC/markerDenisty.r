#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

# What is this system linux or windows
system<-Sys.info()[1]
# if system is windows then change the path separator
if (system=="Windows") {
  pathSeparator<-"\\"
} else {
  pathSeparator<-"/"
}

#######################################################################################


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} #else if (length(args)==1) {
  # default output file
  #args[2] = "out.txt"
#}

if (length(args) < 3) {
  print(paste("Usage: Rscript plot_GWAS_GEMMA.r <gwasfile> <traitName> <ResultFolder> <fdrcutoff> <qthreshold> <chrTable>",sep=""))
  quit()
}


# plot marker denisty across genome

vcfFile<-args[1]

vcfFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData.vcf"


# Documentation using roxygen2
#' Calculate marker density across the genome
#' @param markerInfo a data frame with columns chrom, pos
#' @param window.size the size of the window in base pairs
#' @export
MarkerDenistyByWindow<-function(markerInfo, window.size=1000000)
{
  chroms<-unique(markerInfo$chrom)
  # calculate density of markers in each chromosome
  chromDensity<-data.frame(matrix(ncol=4, nrow=0))
  colnames(chromDensity)<-c("chrom","pos","count","window.num")
  # for each chromosome
  for (i in 1:length(chroms)){
    # select only that chromosome
    chromInfo <- markerInfo[markerInfo$chrom == chroms[i],]
    # calculate number of windows
    num.windows <- ceiling(max(chromInfo$pos)/window.size)
    # create a vector of window positions
    window.pos <- seq(1, max(chromInfo$pos), by=window.size)
    # create a vector of window counts
    window.count <- rep(0, num.windows)
    # loop through each window
    for (j in 1:num.windows){
      # select only markers in that window
      window <- chromInfo[chromInfo$pos >= window.pos[j] & chromInfo$pos < window.pos[j] + window.size,]
      # count number of markers in that window
      window.count[j] <- nrow(window)
    }
    # create a data frame of window positions and counts
    window.df <- data.frame(pos=window.pos, count=window.count, chrom=chroms[i], window.num=1:num.windows)
    # add the data frame 
    chromDensity<-rbind(chromDensity,window.df)
  }
  chromDensity
}

# Documentation using roxygen2
#' Calculate marker density across the genome
#' @param chromDensity a data frame with columns window.num, count, chrom
#' @export
#' @import ggplot2

MarkerDenistyByWindowPlot<-function(chromDensity)
{
  ggplot(chromDensity, aes(x = window.num, y = 1, fill = count)) +
  geom_bar(stat = "identity", position = "fill") +  
  # Create gradient color scale of red, yellow, green
  scale_fill_gradient2(
    low = "#084908",
    mid = "yellow",
    high = "#d30303",
    midpoint = 10,
    na.value = "#d30303",
    space = "Lab",
    limits = c(0, 20),
    labels = c("0", "5", "10", "15", ">=20")
  ) +
  # Set themes
  theme_bw() +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0, "lines"),
    # remove legend title
    legend.title = element_blank()
  ) +  
  # Facet wrap
  facet_wrap(~ chrom, ncol = 1, scales = "free_y", strip.position = "left")
}

# Documentation using roxygen2
#' Calculate marker density across the genome
#' @param markerInfo a data frame with columns chrom, pos
#' @export
MarkerDenistyTotal<-function( markerInfo)
{
  MarkerCount<-markerInfo$chrom %>% table() %>% as.data.frame()
  colnames(MarkerCount)<-c("chrom","count")
  MarkerCount
}

# Documentation using roxygen2
#' Calculate total markers for each chromosome
#' @param eachChromTotal a data frame with columns chrom, count
#' @export
MarkerDenistyTotalPlot<-function(eachChromTotal)
{
  # reverse order of chromosomes
  eachChromTotal$chrom<-factor(eachChromTotal$chrom, levels=rev(levels(eachChromTotal$chrom)))
  # plot total across chromosomes as bar plot
  ggplot(eachChromTotal, aes(y=chrom, x=count, fill=count)) +
  geom_bar(stat="identity")+
  # change color scale from red to green
  scale_fill_gradient2(low="#084908", 
                       mid="yellow", 
                       high="#d30303", 
                       midpoint= mean(eachChromTotal$count),
                       na.value = "#d30303",
                       space = "Lab") +
  # minimal theme
  theme_minimal()+
  # theme(panel.border = element_blank(), 
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(),
  #       axis.text.y = element_blank(), 
  #       axis.ticks.y = element_blank(), 
  #       legend.position="bottom",
  #       # axis.title.x = element_blank(),
  #       axis.title.y = element_blank(),
  #       legend.title = element_blank(),
  #        # legend label size
  #       legend.text = element_text(size=5),
  #       )+
        scale_x_reverse()
}

# Documentation using roxygen2
#' Calculate marker density across the genome
#' @param markerInfo a data frame with columns chrom, pos
#' @param AcrossGenome a boolean indicating whether to calculate across the genome
#' @param ForeachChromosome a boolean indicating whether to calculate for each chromosome
#' @param window.size the size of the window in base pairs
#' @export
MarkerDenisty<-function(markerInfo, AcrossGenome=TRUE, ForeachChromosome=TRUE, window.size=1000000)
{
  MarkerDen<-list()
  chromDensity<-MarkerDenistyByWindow(markerInfo, window.size)
  eachChromTotal<-MarkerDenistyTotal(markerInfo)
  MarkerDen$chromDensity<-chromDensity
  MarkerDen$eachChromTotal<-eachChromTotal
  # description(MarkerDen)
  # description(MarkerDen$chromDensity)<-"Marker density by window"
  # description(MarkerDen$eachChromTotal)<-"Marker density by chromosome"
  MarkerDen
}

MarkerDenistyPlot<-function(MarkerDen)
{
  AcrossGenomePlot<-NULL
  ForeachChromosomePlot<-NULL
  if(!is.null(MarkerDen$eachChromTotal))
    ForeachChromosomePlot<-MarkerDenistyTotalPlot(MarkerDen$eachChromTotal)
  if(!is.null(MarkerDen$chromDensity))
    AcrossGenomePlot<-MarkerDenistyByWindowPlot(MarkerDen$chromDensity)

  if(!is.null(AcrossGenomePlot) && !is.null(ForeachChromosomePlot))
    grid.arrange(ForeachChromosomePlot, AcrossGenomePlot, ncol=2, widths=c(1,4))
  else if(!is.null(AcrossGenomePlot))
    AcrossGenomePlot
  else if(!is.null(ForeachChromosomePlot))
    ForeachChromosomePlot
  else
    NULL
}


snps<-read.table(vcfFile,sep="\t",header=F,blank.lines.skip=TRUE,
                 comment.char = "#")
# get only the first 5 columns
snps<-snps[,1:2]

library(dplyr)
library(ggplot2)
library(gridExtra)

# rename columns
colnames(snps)<-c("chrom","pos")
# rename NA chromosomes to Unknown
snps$chrom[is.na(snps$chrom)]<-"Unknown"
# calculate marker density
MarkerDen<-MarkerDenisty(snps, AcrossGenome=TRUE, ForeachChromosome=TRUE, window.size=50000000)
MarkerDen
# plot marker density
p<-MarkerDenistyPlot(MarkerDen)
ResultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/RESULTS"
ggsave(paste(ResultFolder,"/MarkerDensity.png",sep=""),p, width = 20, height = 20, units = "in", dpi = 300)

