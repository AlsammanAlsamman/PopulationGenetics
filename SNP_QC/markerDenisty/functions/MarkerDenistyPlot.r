#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################


MarkerDenistyPlot<-function(MarkerDen,maxMarkerCount = 50)
{
  AcrossGenomePlot<-NULL
  ForeachChromosomePlot<-NULL
  if(!is.null(MarkerDen$eachChromTotal))
    ForeachChromosomePlot<-MarkerDenistyTotalPlot(MarkerDen$eachChromTotal)
  if(!is.null(MarkerDen$chromDensity))
    AcrossGenomePlot<-MarkerDenistyByWindowPlot(MarkerDen$chromDensity,maxMarkerCount)
  if(!is.null(AcrossGenomePlot) && !is.null(ForeachChromosomePlot))
    grid.arrange(ForeachChromosomePlot, AcrossGenomePlot, ncol=2, widths=c(1,4))
  else if(!is.null(AcrossGenomePlot))
    AcrossGenomePlot
  else if(!is.null(ForeachChromosomePlot))
    ForeachChromosomePlot
  else
    NULL
}
