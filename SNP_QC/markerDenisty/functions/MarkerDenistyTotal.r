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
#' @export
MarkerDenistyTotal<-function( markerInfo)
{
  MarkerCount<-markerInfo$chrom %>% table() %>% as.data.frame()
  colnames(MarkerCount)<-c("chrom","count")
  MarkerCount
}