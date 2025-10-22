
SsInter <- function(databulk,location) {
  if(is.null(databulk)){
    stop("No spatial data")
  }
  databulk <- as.matrix(databulk)
  Ss_expr <- CS(databulk)
  
  if(dim(Ss_expr)[1] != dim(location)[1]){
    nreal <- dim(location)[1]
    npseudo <- dim(Ss_expr)[1] - dim(location)[1]
    Ss_spatial <- Ed(location)
    Ss_expr[-(1:npseudo),-(1:npseudo)] <- (Ss_expr[-(1:npseudo),-(1:npseudo)] + Ss_spatial) / 2
    return(Ss_expr)
    } else {

    Ss_spatial <- Ed(location)
    return((Ss_expr + Ss_spatial) / 2)
  }
  
 
    
  }

  
  



SsIntra <- function(databulk,location) {
  if(is.null(databulk)){
    stop("No spatial data")
  } else if(is.null(location)) {
    message("No location!continue............")
    Ss_spatial <- matrix(0,ncol(databulk),ncol(databulk))
  } else {
    Ss_spatial <- Ed(location)
  }
  databulk <- as.matrix(databulk)
  Ss_expr <- CS(databulk)
  if(dim(Ss_expr)[1] == dim(Ss_spatial)[1]){
    return((Ss_expr + Ss_spatial) / 2)
  } else{
    stop("dim error")
  }
  
  
  
}


Ss <- function(Spatial.Deconv.Data, addSpatial = FALSE) {
  if(is.null(Spatial.Deconv.Data$BatchCorr)) {
    print("No Batch effect correction!")
    data_combine <- as.matrix(Spatial.Deconv.Data$T_combine)
    Ss_expr <- CS(data_combine)
    if(!is.null(Spatial.Deconv.Data$location.pseudo)){
      print("use both pseudo loc and real loc................")
      data_loc_combine <- rbind(Spatial.Deconv.Data$location.pseudo,Spatial.Deconv.Data$location.real)
      Ss_spatial <- Ed((data_loc_combine))
      Spatial.Deconv.Data$Ss <- Ss_spatial + Ss_expr
    } else if(addSpatial) {
      Ss_spatial <- Ed((Spatial.Deconv.Data$location.real))
      add <- (Ss_spatial + Ss_expr[(ncol(Spatial.Deconv.Data$spatial.pseudo) + 1) : ncol(Ss_expr),(ncol(Spatial.Deconv.Data$spatial.pseudo) + 1) : ncol(Ss_expr)]) / 2
      Ss_expr[(ncol(Spatial.Deconv.Data$spatial.pseudo) + 1) : ncol(Ss_expr),(ncol(Spatial.Deconv.Data$spatial.pseudo) + 1) : ncol(Ss_expr)] <- add 
      Spatial.Deconv.Data$Ss <- Ss_expr 
    } else {
      Spatial.Deconv.Data$Ss <- Ss_expr 
    }
    
  }else if(!is.matrix(Spatial.Deconv.Data$BatchCorr)){
    print("calculate Ss using GraphST embedding data(with both PSEUDO AND REAL spatial information).....")
    data_loc_combine <- as.matrix((Spatial.Deconv.Data$data_loc_combine.gst))
    Ss_spatial <- Ed((data_loc_combine))
    data_combine <- Spatial.Deconv.Data$data_combine.gst
    Ss_expr <- CS(data_combine)
    Ss <- Ss_spatial + Ss_expr
    Spatial.Deconv.Data$Ss <- Ss
  }else {
    if(addSpatial) {
      if(is.null(Spatial.Deconv.Data$location.pseudo)){
        message("calculate Ss using Correction data with only REAL spatial information.....")
        
        Ss_spatial <- Ed((Spatial.Deconv.Data$location.real))
        Ss_expr <- CS(Spatial.Deconv.Data$BatchCorr)
        add <- (Ss_spatial + Ss_expr[(ncol(Spatial.Deconv.Data$spatial.pseudo) + 1) : ncol(Ss_expr),(ncol(Spatial.Deconv.Data$spatial.pseudo) + 1) : ncol(Ss_expr)]) / 2
        Ss_expr[(ncol(Spatial.Deconv.Data$spatial.pseudo) + 1) : ncol(Ss_expr),(ncol(Spatial.Deconv.Data$spatial.pseudo) + 1) : ncol(Ss_expr)] <- add 
        Spatial.Deconv.Data$Ss <- Ss_expr 
      }else{
        message("calculate Ss using Correction data with both REAL and PSEUDO spatial information.....")
        data_loc_combine <- rbind(Spatial.Deconv.Data$location.pseudo,Spatial.Deconv.Data$location.real)
        data_loc_combine <- as.matrix(data_loc_combine)
        Ss_spatial <- Ed((data_loc_combine))
        data_combine <- Spatial.Deconv.Data$BatchCorr
        Ss_expr <- CS(data_combine)
        Ss <- (Ss_spatial + Ss_expr) / 2
        Spatial.Deconv.Data$Ss <- Ss
        }
     
    } else{
      print("calculate Ss using Correction data without spatial information.....")
      Ss_expr <- CS(Spatial.Deconv.Data$BatchCorr)
      Spatial.Deconv.Data$Ss <- Ss_expr    
    }

  }
  return(Spatial.Deconv.Data)
      
}



##Euclidean distance
Ed <- function(A){
  distance <- as.matrix(dist(A))
  sm <- 1 / (1 + distance)
  return(sm)
}


CS <- function(A){
  library(lsa)
  return(cosine(A))
}
