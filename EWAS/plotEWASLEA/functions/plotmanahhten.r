# Documentation using roxygen2
#' @title plotmanahhten for GWAS
#' @description plotmanahhten for GWAS data
#' @param gwasResults a data frame with the GWAS results
#' @param ChromInfo a data frame with the chromosome information generated by getChrInfoFromHap or getChrInfoFromSNPINFO
#' @param logmin a numeric value to set the log10 pvalue threshold for annotation
#' @param log.global a numeric value to set the log10 pvalue global threshold
#' @param TraitName a character string to set the name of the trait
#' @param highLightSignf a logical value to set if the significant SNPs will be highlighted
#' @param annotateSignf a logical value to set if the significant SNPs will be annotated
#' @param annotate.this a character vector with the SNPs to annotate
#' @param annotationCol a character string with the column name to annotate
#' @param dash.this a character vector with the SNPs to dash vertical lines
#' @param xshowtitle a logical value to set if the x-axis title will be shown
#' @param yshowtitle a logical value to set if the y-axis title will be shown
#' @param xfontsize a numeric value to set the x-axis title font size
#' @param yfontsize a numeric value to set the y-axis title font size
#' @param matchChrTable a logical value to set if the chromosome table will be matched with the GWAS results
#' @return a plot
#' @import ggplot2 ggrepel tidyr dplyr
plotmanahhten <- function(gwasResults,
                          fdrcutoff = 0.05,
                          logmin=NA,
                          log.global = 2,
                          TraitName = "TraitName",
                          highLightSignf = T,
                          annotateSignf = F,
                          annotate.this = c(),
                          annotationCol = "rs",
                          dash.this = c(),
                          # show x-axis title
                          xshowtitle = T,
                          # show y-axis title
                          yshowtitle = T,
                          xfontsize = 8, 
                          yfontsize = 8,
                          matchChrTable = T,
                          ColorByCol="chr",
                          colors=c(),
                          fdrSignSNPs=c()){
  # Table containing the chromosome information and thier order
  # Its better to provide this table to save time
  ChrTable<-NULL
  # Table number of snps per chromosome it will be used to order the chromosomes
  ChrTable <- getChrInfoFromSNPGWAS(gwasResults)
  
  fdrColor <- "#F51202"
  fdrDashColor <- "#F51202"

  qthresholdColor <- "#060cbF"
  qthresholdDashColor <- "#060cbF"
  
  # if chr is not a factor
  if(!is.factor(gwasResults$chr))
  {
    gwasResults$chr <- factor(gwasResults$chr, levels = ChrTable$chr)
  } 
  if(is.factor(gwasResults$chr))
  {
    # match the chromosome table with the GWAS results
   print("Chromosomes were sorted according to the chromosome table")
  }
  # log Max
  logmax <- max(-log10(gwasResults$p_wald)) + 2 # max log10 pvalue
  # Minimize for limits proposes
  gwasResults$ps <- gwasResults$ps / 100000
  # number of chromosomes
  nchr <- length(ChrTable$chr)

  # Prepare the dataset for plotting
  don <- gwasResults %>%
    # Compute chromosome size
    group_by(chr) %>%
    summarise(chr_len = max(ps)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwasResults, ., by = c("chr" = "chr")) %>%
    # Add a cumulative position of each SNP
    arrange(chr, ps) %>%
    mutate(Chromosome = ps + tot) %>%
    mutate(is_fdr = ifelse(fdr < fdrcutoff, "yes", "no")) %>%
    # for qthreshold
    mutate(is_qthreshold = ifelse(-log10(p_wald) > log.global, "yes", "no"))

  # Prepare X axis by splitting by chromosome in plot
  axisdf <- don %>%
    group_by(chr) %>%
    summarize(center = (max(Chromosome) + min(Chromosome)) / 2)

  # how much colors are needed
  ncol <- length(unique(don[[ColorByCol]]))
  # get the colors
  mycolors <- getColors(ncol,fixed = T)
  
  # plot object
  p <- ggplot(don, aes(x = Chromosome, y = -log10(p_wald))) +
    # plot point and color by target column
    geom_point(aes(color = get(ColorByCol)),alpha = 0.8, size = 1.3) +
    # add TraitName in the corner
    ggtitle(TraitName) +
    # custom X axis:
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, logmax)) + # remove space between plot area and x axis    
    # Custom the theme:
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(color = "red", size = 8, face = "bold.italic"),
      # do not show x-axis title
      axis.title.x = element_blank(),
      # do not show y-axis title
      axis.title.y = element_blank())
    
  # Draw horizontal line for logmin
  if (!is.na(logmin)) {
    # dashed lines for fdr
    p <- p + geom_hline(yintercept = logmin, linetype = "dashed", color =  fdrDashColor, size = 0.5)
  }  
  # Draw horizontal line for log.global
  if (!is.na(log.global)) {
    # dashed lines for logmin
    p <- p + geom_hline(yintercept = log.global, linetype = "dashed", color = qthresholdDashColor, size = 0.5)
  }
  
  # if xtitle is true
  if (xshowtitle) {
    p <- p + xlab("Chromosome") + theme(axis.title.x = element_text(size = xfontsize, face = "bold.italic"))
  }
  # if ytitle is true
  if (yshowtitle) {
    p <- p + ylab(expression(paste("pvalue", "(-log"[10], ")"))) + theme(axis.title.y = element_text(size = xfontsize, face = "bold.italic", angle = 90))
  }


  # color by fdr and qthreshold
  p <- p + geom_point(data = subset(don, is_fdr == "yes"), color = fdrColor, size = 1.5)
  p <- p + geom_point(data = subset(don, is_qthreshold == "yes" & is_fdr == "no"), color = qthresholdColor, size = 1.5)
  
  p <- p + geom_text_repel(
      data = subset(don, is_fdr == "yes" | is_qthreshold == "yes"),
      aes_string(label = annotationCol),
      nudge_x = .15,
      box.padding = 0.5,
      nudge_y = 1,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.linetype = 3,
      size = 2,
      segment.color = colorTransparancy("#0a6428", 60),
      segment.angle = 20
    )
  
  # if there is loci to dash
  if (length(dash.this) > 0) {
    # add vertical line
    p <- p + geom_vline(xintercept = don[don$rs %in% dash.this, ]$Chromosome, linetype = "dashed", color = colorTransparancy("#000000", 50), size = 0.2)
  }

  # add color legend
  if (length(ColorByCol) > 0 & ColorByCol!="chr") {
    if(length(colors)>0){
      mycolors <- colors
    }
    p <- p + scale_color_manual(values = mycolors)
    p <- p + theme(legend.position = "right")+
      # legend title
      guides(color = guide_legend(title = ColorByCol))+
      # no legend title
      theme(legend.title = element_blank())}

  # if by chr is true
  if (ColorByCol=="chr") {
    mycolors<-rep(c("#cfcf08", "#7d994b"), nchr)
    if(length(colors)>0){
      mycolors <- rep(colors, nchr)
    }
    p <- p +     scale_color_manual(values = mycolors) +
    # remove legend
    theme(legend.position = "none")
  }
  # return plot object
  return(p)
}


# Documentation using roxygen2
#' @title Plot GWAS results in a QQ plot
#' @description Plot GWAS results in a QQ plot
#' @param gwasResults A data frame with the GWAS results
#' @param col Color of the points
#' @param cex Size of the points
#' @param xshowtitle Show x axis title
#' @param yshowtitle Show y axis title
#' @param xfontsize Size of the x axis title
#' @param yfontsize Size of the y axis title
#' @return A ggplot object
#' @export
#' @import ggplot2 dplyr
plotqq <- function(gwasResults, col = "red", coln = 1, colsize = 1, xshowtitle = T, yshowtitle = T,
                   xfontsize = 8, yfontsize = 8) {
  # generate colors
  if (!colsize == 1) {
    col <- getColors(colsize)[coln]
  }
  ci <- 0.95
  nSNPs <- length(gwasResults$p_wald)
  plotdata <- data.frame(
    observed = -log10(sort(gwasResults$p_wald)),
    expected = -log10(ppoints(nSNPs)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))
  )

  qqplot <- ggplot(plotdata, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = colorTransparancy(col, 100), alpha = 0.5) +
    geom_step(color = "#070707", size = 1.1, direction = "vh") +
    geom_segment(
      data = . %>% filter(expected == max(expected)),
      aes(x = 0, xend = expected, y = 0, yend = expected),
      size = 1.25, alpha = 0.5, color = "#ededed", lineend = "round"
    ) +
    # do not show x-axis title
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      # do not show x-axis title
      axis.title.x = element_blank(),
      # do not show y-axis title
      axis.title.y = element_blank()
    )
  # if xtitle is true
  if (xshowtitle) {
    qqplot <- qqplot + xlab(expression(paste("Expected pvalue", "(-log"[10], ")"))) + theme(axis.title.x = element_text(size = xfontsize, face = "bold.italic"))
  }
  if (yshowtitle) {
    qqplot <- qqplot + ylab(expression(paste("Observed pvalue", "(-log"[10], ")"))) + theme(axis.title.y = element_text(size = yfontsize, face = "bold.italic", angle = 90))
  }
  return(qqplot)
}

# Documentation using roxygen2
#' @title Plot GWAS results in a manhatten plots
#' @description Plot GWAS results in a manhatten plots
#' @param GWAS A data frame with the GWAS results
#' @param GWASDetails A data frame with the GWAS details
#' @param ChromInfo A data frame with the chromosome information
#' @param logmin The minimum value for the log scale
#' @param highLightSignf Highlight significant loci
#' @param annotateSignf Annotate significant loci
#' @param annotate.this A vector with the loci to annotate
#' @param annotationCol The column to use for the annotation
#' @param dash.this A vector with the loci to dash
#' @param xfontsize Size of the x axis title
#' @param yfontsize Size of the y axis title
#' @param matchChrTable Match the chromosome table
#' @param mergePlots Merge the plots
#' @param showQQ Show the QQ plot
#' @param showFileNames Show the file names
#' @return A list of ggplot objects
#' @import ggplot2 ggpubr dplyr patchwork
plotmanahhtenList<-function(GWAS,
                            GWASDetails,
                            ChromInfo = NULL,
                            logmin = 3,
                            log.global = 2,
                            highLightSignf = T,
                            annotateSignf = F,
                            annotate.this = c(),
                            annotationCol = "rs",
                            dash.this = c(),
                            xfontsize = 8, 
                            yfontsize = 8,
                            matchChrTable = T, 
                            mergePlots=T,
                            showQQ=T,
                            showFileNames=T,
                            showxtitle=T,
                            showytitle=T)
{

  # its better to handle ChromInfo if it is not provided
  if(is.null(ChromInfo))
  {
    ChromInfo<-getChrInfoFromSNPGWAS(GWAS)
  }


  manhattenPlots <- list()
  plotsInEach<-list()
  gwasFiles<-GWASDetails$files
  gwasFilesCount<-length(gwasFiles)
  if(mergePlots)
  {
      showxtitle <- F 
  }
  for (f in 1:gwasFilesCount)
  {
    # print the current trait
    print(paste0("Data: ", gwasFiles[f]))
    
    # Current trait
    thisFile<- gwasFiles[f]
    # get gwas details for the current trait
    traitGWAS<-GWAS[GWASDetails$trait.rows,]
    traitGWAS<-traitGWAS[traitGWAS$traitData==thisFile,]
    
  
    #-------------------------------------------
    # Manhatten plots and QQ plots for each trait
    #-------------------------------------------
    # create list for gwas
    manhattenPlots[[thisFile]] <- list()
    
    # if last plot show the x title and mergePlots is true
    if (f == gwasFilesCount && mergePlots) {
      showxtitle <- T
    }
      
    
    # Manhatten plot
    manp <- plotmanahhten(traitGWAS,
                            ChromInfo = ChromInfo,
                            dash.this = dash.this, 
                            xshowtitle = showxtitle, 
                            yshowtitle = showytitle,
                            TraitName = ifelse(showFileNames, thisFile, ""),
                            xfontsize = 8, 
                            yfontsize = 8,
                            matchChrTable=F)
      
      # save the plot
      manhattenPlots[[thisFile]]<-manp
      # QQ plot
      if(showQQ)
      {
        qp <- plotqq(traitGWAS,
                    coln = f, 
                    colsize = gwasFilesCount, 
                    xshowtitle = showxtitle, 
                    yshowtitle = showytitle,
                    xfontsize = 8, 
                    yfontsize = 8)
        # merge the plots together horizontally where the first plot takes 2/3 of the width
        manp <- ggarrange(manp, qp, nrow = 1, ncol = 2, widths = c(2, 1))
        # save the plot
        manhattenPlots[[thisFile]]<-manp
      }
  }  

  # merge plots for each trait if mergePlots is true
  if(mergePlots)
  {
    firstPlot<-manhattenPlots[[1]]
    for (i in 2:gwasFilesCount) {
      # vertical merge
      firstPlot <- firstPlot + manhattenPlots[[i]] + theme(legend.position = "none") +  plot_layout(ncol = 1)
    }
    manhattenPlots <- firstPlot
  }
  
  plots<-list()
  plots$manhattenPlots<-manhattenPlots
  plots$qqPlots<-showQQ
  plots$mergePlots<-mergePlots
  plots$plotsInEach<-gwasFilesCount
  return(plots)
}


# saveManhattanPlot<-function(plots,TraitInfoTable,ResultsFolder,width=10,height=3)
# {
#   plots<-plotLists
#   width<-15
#   height<-3
#   plots$manhattenPlots[[1]][[1]]
#   plotNames<-names(plots$manhattenPlots)
#   plots.plots<-plots$manhattenPlots
#   plots.n<-length(plotNames)
#   plots.nn<-plots$plotsInEach
#   if(plots$mergePlots)
#   {
#     for(i in 1:length(plotNames))
#     {
#       traitName<-plotNames[i]
#       print(paste0("Saving plot: ",traitName))
#       # get the plot
#       plot<-plots.plots[[traitName]]
#       # get folder name from TraitInfoTable
#       folderName<-TraitInfoTable[TraitInfoTable$Trait==traitName,"OutFolderName"][1]
#       #get the file name
#       traitNameFile<-gsub(" ","_",traitName)
#       # remove unwanted characters
#       traitNameFile<-gsub("[^[:alnum:]]", "", traitNameFile)
#       # get the file name from TraitInfoTable
#       fileNamePath<-paste0(ResultsFolder,"/",folderName,"/",traitNameFile)
#       # save the plot
#       ggsave(filename = paste0(fileNamePath,"_ManhattanPlot.pdf"), plot = plot, width = width, height = height*plots.nn[[traitName]])
#     }
#   }
#   else
#   {
#     for(i in 1:plots.n)
#     {
#      traitName<-plotNames[i]
#      # loop over the plots for each trait
#      for(p in 1:plots.nn[[traitName]])
#      {
#       dataplot<-plots.plots[[traitName]][[p]]
#       print(paste0("Saving plot: ",traitName," Data ",p))
#       # get folder name from TraitInfoTable
#       folderName<-TraitInfoTable[TraitInfoTable$Trait==traitName,"OutFolderName"][1]
#       #get the file name
#       dataname<-rownames(TraitInfoTable[TraitInfoTable$Trait==traitName,])[p]
#       # remove unwanted characters
#       dataname<-gsub("[^[:alnum:]]", "", dataname)
#       # get the file name
#       fileNamePath<-paste0(ResultsFolder,"/",folderName,"/",dataname)
#       # save the plot
#       ggsave(filename = paste0(fileNamePath,"_ManhattanPlot.pdf"), plot = dataplot, width = width, height = height)
#      }
#     }
#   }
#   print("Plots saved successfully")  
# }



# Documentation using roxygen2
#' match GWAS with chromosome table
#' @Description match GWAS with chromosome table if the GWAS does not have all chromosomes in ChrTable or the start and end positions are not the same add fake data to avoid errors in the plot
#' @param gwasResults a GWAS data
#' @param ChrTable a chromosome table
#' @return a data frame of GWAS data 
matchGWASwithChrTable<-function(gwasResults,ChrTable)
{
  
  # Check if GWAS does not have all chromosomes 
  #in ChrTable or the start and end positions are not the same add fake data
  #to avoid errors in the plot
  
  gwasChrINFO <- getChrInfoFromSNPGWAS(gwasResults)
  # compare the chromosomes
  # chromosomes does exist in ChrTable and not in gwasChrINFO
  chr.not.in.gwas <- ChrTable[!(ChrTable$chr %in% gwasChrINFO$chr), ]
  faken<-1
  faketraitData<-gwasResults$traitData[1]
  # add fake data to gwasResults
  if (nrow(chr.not.in.gwas) > 0) {
    for (i in 1:nrow(chr.not.in.gwas)) {
      # start
      gwasResults <- rbind(gwasResults, data.frame(
        chr = chr.not.in.gwas$chr[i],
        ps = chr.not.in.gwas$start[i],
        effect = 0,
        se = 0,
        traitData = faketraitData,
        pvalue = 1,
        rs = paste("fake",faken,sep="")
      ))
      faken<-faken+1
      # end
      gwasResults <- rbind(gwasResults, data.frame(
        chr = chr.not.in.gwas$chr[i],
        ps = chr.not.in.gwas$end[i],
        effect = 0,
        se = 0,
        traitData = faketraitData,
        pvalue = 1,
        rs = paste("fake",faken,sep="")
      ))
      faken<-faken+1
    }
  }
  # chromosomes in both ChrTable and gwasChrINFO but the start is not the same
  chr.in.gwas <- ChrTable[ChrTable$chr %in% gwasChrINFO$chr, ]
  chr.in.gwas <- chr.in.gwas[!(chr.in.gwas$start == gwasChrINFO$start), ]
  # add fake data to gwasResults
  if (nrow(chr.in.gwas) > 0) {
    for (i in 1:nrow(chr.in.gwas)) {
      # start
      gwasResults <- rbind(gwasResults, data.frame(
        chr = chr.in.gwas$chr[i],
        ps = chr.in.gwas$start[i],
        effect = 0,
        se = 0,
        traitData = faketraitData,
        pvalue = 1,
        rs = paste("fake",faken,sep="")
      ))
      faken<-faken+1
    }
  }
  # chromosomes in both ChrTable and gwasChrINFO but the end is not the same
  chr.in.gwas <- ChrTable[ChrTable$chr %in% gwasChrINFO$chr, ]
  chr.in.gwas <- chr.in.gwas[!(chr.in.gwas$end == gwasChrINFO$end), ]
  # add fake data to gwasResults
  if (nrow(chr.in.gwas) > 0) {
    for (i in 1:nrow(chr.in.gwas)) {
      # end
      gwasResults <- rbind(gwasResults, data.frame(
        chr = chr.in.gwas$chr[i],
        ps = chr.in.gwas$end[i],
        effect = 0,
        se = 0,
        traitData = faketraitData,
        pvalue = 1,
        rs = paste("fake",faken,sep="")
      ))
      faken<-faken+1
    }
  }
  return(gwasResults) 
}