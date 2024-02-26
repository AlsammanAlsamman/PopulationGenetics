#Documentation using roxygen2
#' @title Get Colors from Colorbrewer
#' @param n the number of colors will be used
#' @return a vector of colors
#' @export
#' @importFrom RColorBrewer
getColorsFromPallete<-function(n,fixed=FALSE){
  # get n colors from colorbrewer
  # n: number of colors  
  # get all colors
  allColors<-brewer.pal.info
  # get all colors with n colors
  allColors<-allColors[allColors$maxcolors>=n,]
  # get the first color
  colorP<-sample(nrow(allColors),1)
  if (fixed){
    colorP<-1
  }
  colorPname<-rownames(allColors)[colorP]
  colorPlength<-allColors$maxcolors[colorP]
  # get the colors
  colors<-brewer.pal(name = colorPname, n = colorPlength)
  colors<-sample(colors, n)
  return(colors)
}

#' @title Get Colors
#' @param n the number of colors will be used
#' @return a vector of colors
#' @export
#' @importFrom RColorBrewer
getColors<-function(n,fixed = FALSE){
  # get n colors
  # n: number of colors
  # if n less than 12 get colors from colorbrewer
  colors<-NULL
  if (n<12){
    colors<-getColorsFromPallete(n,fixed)
  }else{
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if (fixed){
      colors<-col_vector[1:n]
  }else{
    colors<-sample(col_vector, n)
  }
  }
  return(colors)
}
# Documentation using roxygen2
#' @title Add Transparancy to a Colors for ggplot
#' @param color the color will be used
#' @param trans the transparancy percentage will be used
#' @return a vector of colors with transparancy
#' @export 
#' @import RColorBrewer
colorTransparancy <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}







