

LoadLibraries <- function(libs,lib_instruct_file)
{  
  print(libs) 
  libinstruct<-read.table(lib_instruct_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
  # check if the libraries are installed
  for (i in libs) {
      if (!require(i,character.only = TRUE)) {
        # get the repository
        repo<-libinstruct[libinstruct$library==i,"repository"]
        if (repo=="NA" | repo=="CRAN") {
          install.packages(i,dependencies = TRUE)
        } else if (repo=="Bioconductor") {
            if (!require("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            BiocManager::install(i)
        }
    } 
  }
  # loading libraries
  for (i in libs) {
    print(paste("Loading library",i))
    library(i,character.only = TRUE)
    }
}