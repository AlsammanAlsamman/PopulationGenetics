
openPDFfile<-function(file){
  if (.Platform$OS.type == "windows"){
    system(paste("start", file))
  } else if (.Platform$OS.type == "unix"){
    system(paste("xdg-open", file))
  } else {
    system(paste("open", file))
  }
}