plotLeaR<-function(LEAoutFolder,resultFolder,passportFile,orderFile,bestKnownK,maxK, plotHeight, plotWidth, pathSeparator="\\", labelfontsize = 12)
{
  # get row names from indvOrderFile
  indvOrder<-read.table(orderFile,header=F,sep="\t")
  indvOrder<-indvOrder[,1]

  # read passport file
  passport<-read.csv(passportFile,header=T,sep="\t")
  # rename first column to Taxa
  colnames(passport)[1]<-"Taxa"
  # get passport groups
  passportGroups<-colnames(passport)[2:ncol(passport)]  

  # define color palette with 16 different colors
  colorPalette<-c("#2f4f4f","#8b4513","#006400","#4b0082",
                  "#48d1cc","#ff0000","#ffa500","#ffff00",
                  "#454645","#00fa9a","#0000ff","#ff00ff",
                  "#1e90ff","#f0e68c","#dda0dd","#ff1493")
  # replace NA with unknown
  passport[is.na(passport)]<-"unknown"

  # convert all passport groups columns to factor
  for (i in 2:ncol(passport))
  {
    passport[,i]<-as.factor(passport[,i])
  }
  passport
  indvOrder
  # if some individuals are not in the passport file, add them with unknown passport group
  if (length(setdiff(indvOrder,passport$Taxa))>0)
  {
    print("Some individuals are not in the passport file, they will be added with unknown passport group")
    passport<-rbind(passport,data.frame(Taxa=setdiff(indvOrder,passport$Taxa),matrix(rep("unknown",ncol(passport)-1),nrow=length(setdiff(indvOrder,passport$Taxa)),ncol=ncol(passport)-1)))
  }
  # sort passport file according to indvOrder
  passport<-passport[match(indvOrder,passport$Taxa),]

  ################# Structure plots #################
  # create list for plots
  plotList<-list()

  for (kn in 2:maxK) {

  # read best run for K from bestKs in outFolder 
  bestKs<-read.table(paste(LEAoutFolder, "/","K",kn,".Q",sep=""),header=F,sep="\t")

  # change row names to individual names using indvOrderFile
  rownames(bestKs)<-indvOrder
  # get the max value for each individual
  MaxValues<-apply(bestKs,1,max)
  #calculating the group by max of col and retrieve column name
  bestKs$bestGroup<-colnames(bestKs)[max.col(bestKs)]
  # add column for max value
  bestKs$maxValue<-MaxValues
  # replace V with Q in group column
  bestKs$bestGroup<-sub("V","Q",bestKs$bestGroup)
  # add taxa names according to rownames
  bestKs$Taxa<-rownames(bestKs)
  #re-ashape the table for plotting and keep columns bestGroup and Taxa
  bestKs<-melt(bestKs,id.vars=c("bestGroup","Taxa","maxValue"))
  # convert bestGroup to factor
  bestKs$bestGroup<-as.factor(bestKs$bestGroup)
  # create a numeric column for bestGroup
  bestKs$bestGroupNum<-sub("Q","",bestKs$bestGroup)
  bestKs$bestGroupNum<-as.numeric(bestKs$bestGroupNum)
  # rename variable column to Q
  bestKs$Q<-sub("V"," ",bestKs$variable)
  # convert Q to numeric
  bestKs$Q<-as.numeric(bestKs$Q)
  # convert to factor
  bestKs$Q<-as.factor(bestKs$Q)
  # sort levels of Q
  levels(bestKs$Q)<-1:kn

  # Plot the population structure for K=kn
  structurePlot<-bestKs %>%  # use mutate to order using bestGroupNum and maxValue
    mutate(ordering = -as.numeric(bestGroupNum) + maxValue, Taxa = fct_reorder(Taxa, ordering), .desc = F)%>%
    ggplot(aes(x=Taxa,y=value,fill=Q))+
    geom_bar(stat="identity",position="fill")+ geom_col()+
    # make labels for x axis vertical
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # add upper highlight to bestGroup dacshed line
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")+
    # use color palette
    scale_fill_manual(values=colorPalette)+
    # remove x axis labels (The Taxa names are too long)
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  ################# Plot passport groups ############################
  # We will now calculate the proportion of how each passport group is represented in each Q group
  # table for best group and taxa
  bestGroupCal<-bestKs[,c("bestGroup","Taxa")]
  # remove duplicates
  bestGroupCal<-unique(bestGroupCal)
  # join with passport file
  bestGroupCal<-merge(bestGroupCal,passport,by.x="Taxa",by.y="Taxa")
  # remove Total column
  bestGroupCal<-bestGroupCal[,grep("Total",colnames(bestGroupCal),invert=T)]

  # melt the table to get only Taxa, bestGroup and two columns for passport groups and corresponding values
  bestGroupCal<-melt(bestGroupCal,id.vars=c("bestGroup","Taxa"))
  colnames(bestGroupCal)<-c("bestGroup","Taxa","passportGroup","passportValue")

  # calculate the proportion of each passport group in each Q group
  groupsCountQ<-bestGroupCal %>% group_by(bestGroup,passportGroup,passportValue) %>% summarise(count=n())
  # convert to data frame
  groupsCountQ<-as.data.frame(groupsCountQ)
  # rename BestGroup to Q
  groupsCountQ$Q<-groupsCountQ$bestGroup
  # convert Q to numeric
  groupsCountQ$Q<-as.numeric(sub("Q","",groupsCountQ$Q))
  # convert to factor
  groupsCountQ$Q<-as.factor(groupsCountQ$Q)
  # sort levels of Q
  levels(groupsCountQ$Q)<-1:kn

  # Plotting the results for passport groups representation in each Q group
  # We can use facet_wrap to plot each passport group in a separate plot
  passportPlot<-groupsCountQ %>%
    ggplot(aes(x=passportValue,y=count,fill=Q))+
    geom_bar(stat="identity",position="fill")+ 
    # make labels for x axis vertical
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # use color palette
    scale_fill_manual(values=colorPalette)+
    # add facet wrap
    facet_grid(~passportGroup,scales="free_x")+
    theme(strip.text.x = element_text(size = labelfontsize, face = "bold"),
    strip.text.y = element_text(size = labelfontsize, face = "bold"))+
    # remove x axis title
    theme(axis.title.x = element_blank())

  # join the two plots, structurePlot and passportPlot, with one legend and one legend
  joinedPLot<-ggarrange(structurePlot, passportPlot, 
            ncol = 2, nrow = 1, legend = "right", common.legend = TRUE,
            widths = c(0.7, 0.3))
  # add plot to list
  plotList[[kn]]<-joinedPLot
  }

  # join all plots in plotList
  GrandPlot<-ggarrange(plotlist = plotList, ncol = 1, legend = "right", common.legend = TRUE)
  # save plot
  ggsave(paste(resultFolder,"/","structurePlotGrandPlot.pdf",sep="/"),GrandPlot,width=plotWidth,height=plotHeight)
  # save plot of best Known K
  ggsave(paste(resultFolder,"/","structurePlotBestK.pdf",sep="/"),plotList[[bestKnownK]],width=plotWidth,height=5)

  # create table for all Ks 
  # create a data frame for all Ks
  allKs<-matrix(NA,nrow=length(indvOrder),ncol=maxK-1)
  colnames(allKs)<-paste0("K",2:maxK)
  rownames(allKs)<-indvOrder
  # fill the data frame with bestKs for each K and each individual
  # for  loop for each K
  for (kn in 1:ncol(allKs))
  {
    kn<-colnames(allKs)[kn]
    bestKs<-read.table(paste(LEAoutFolder, "/",kn,".Q",sep=""),header=F,sep="\t")
    # the best K is the one with the highest max value
    bestKs$bestGroup<-colnames(bestKs)[max.col(bestKs)]
    bestKs$bestGroup<-gsub("V","",bestKs$bestGroup)
    # fill the allKs table with the bestKs
    allKs[,kn]<-bestKs$bestGroup
  }

  allKs<-as.data.frame(allKs)
  # merge with passport file
  allKs<-merge(allKs,passport,by.x="row.names",by.y="Taxa")

  # save to results the Ks table with passport groups
  write.table(allKs,paste(resultFolder,"/","allKs_passportGroups.txt",sep="/"),quote=F,row.names=F,col.names=T,sep="\t")

}
