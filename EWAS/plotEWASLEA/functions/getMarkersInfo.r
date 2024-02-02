# Documentation using roxygen2
#' Get Marker information from GWAS data
#' @Description Get Marker information from GWAS data
#' @param GWAS dataframes of GWAS data
#' @return a dataframe of Marker information
#' @export
getMarkersInfoFromGWAS<-function(GWASdata)
{

  # remove the duplicated markers and keep the first one using rs column as key
  GWASdata<-GWASdata[!duplicated(GWASdata$rs),]
  # add snpcode names
  GWASdata$SNPCode<-paste0("SNP",1:nrow(GWASdata))
  rownames(GWASdata)<-GWASdata$rs
  # remove all the columns except SNPCode, rs, chr, ps
  GWASdata<-GWASdata[,c("rs","chr","ps","SNPCode")]
  return(GWASdata)
}

getMarkersInfoFromHap<-function(Hapdata)
{
  # remove the duplicated markers and keep the first one using rs column as key
  Hapdata<-HapData[,c("chrom","pos")]
  colnames(Hapdata)<-c("chr","ps")
  # add snpcode names
  Hapdata$SNPCode<-paste0("SNP",1:nrow(Hapdata))
  Hapdata$rs<-rownames(Hapdata)
  Hapdata<-Hapdata[,c("rs","chr","ps","SNPCode")]
  return(Hapdata)
}
