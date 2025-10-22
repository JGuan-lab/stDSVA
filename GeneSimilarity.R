
SgIntra <- function(databulk) {
  if(is.null(databulk)){
    stop("error,no data!")
  } 
  databulk <- as.matrix(databulk)
  return(SM(databulk))
  
}



Sg <- function(Spatial.Deconv.Data) {
  ##calculate gene-gene similarity use PCC
  if(is.null(Spatial.Deconv.Data$BatchCorr)) {
    print("No Batch effect correction!")
    if(is.null(Spatial.Deconv.Data$T_combine)){
      stop("no data!")
    }
    data_combine <- Spatial.Deconv.Data$T_combine
    Spatial.Deconv.Data$Sg <- SM(data_combine)
    
  }else if(!is.matrix(Spatial.Deconv.Data$BatchCorr)){
    print("calculate Sg using GraphST embedding data.....")

    Spatial.Deconv.Data$Sg <- SM(Spatial.Deconv.Data$data_combine.gst)
  }else {
    print("calculate Sg using Mnn Correction data.....")
    Spatial.Deconv.Data$Sg <- SM(Spatial.Deconv.Data$BatchCorr)
  }
  return(Spatial.Deconv.Data)
}



# get  similarity matrix
SM <- function(A){
  
  score = 1/(2-cor(t(A)))
  
  return(score)
}
