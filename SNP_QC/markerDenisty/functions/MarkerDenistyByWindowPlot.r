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
#' @param chromDensity a data frame with columns window.num, count, chrom
#' @export
#' @import ggplot2

MarkerDenistyByWindowPlot<-function(chromDensity,maxMarkerCount = 50)
{
  ggplot(chromDensity, aes(x = window.num, y = 1, fill = count)) +
  geom_bar(stat = "identity", position = "fill", colour="#adaaaa") +  
  # Create gradient color scale of red, yellow, green
  scale_fill_gradient2(
    guide = "colourbar",
    low = "#ffffff",
    na.value = "#360ae5f5",
    space = "Lab",
    limits = c(0, maxMarkerCount),
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
  )+
  facet_wrap(~ chrom, ncol = 1, scales = "free_y", strip.position = "left")
}
