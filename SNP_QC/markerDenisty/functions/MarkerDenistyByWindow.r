#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################


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