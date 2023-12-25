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
