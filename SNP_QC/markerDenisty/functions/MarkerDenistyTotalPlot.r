#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################



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