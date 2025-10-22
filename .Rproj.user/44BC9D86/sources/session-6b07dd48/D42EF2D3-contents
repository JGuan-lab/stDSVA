
library(GSEABase)
library(GSVA)
library(limma)
library(GSEABase)



gsva_my <- function(Spatial.Deconv.Data,gsva_thre=NULL,kcdf = "Poisson") {
  Spatial.Deconv.Data$gsva_mat.real <- gsva_apply(Spatial.Deconv.Data$spatial.real.count.marker, Spatial.Deconv.Data$markerslist, kcdf_ = kcdf)
  Spatial.Deconv.Data$gsva_mat.pseudo <- gsva_apply(Spatial.Deconv.Data$spatial.pseudo.count.marker, Spatial.Deconv.Data$markerslist, kcdf_ = kcdf)
  return(Spatial.Deconv.Data)
  
  
}
gsva_apply <- function(grid_exp, marker , kcdf_ = "Poisson", maxDiff = TRUE){

  ##use gsva 1.52.3 to analyse
  if("cluster" %in% colnames(marker)){
    go_list <- split(marker$gene, marker$cluster)} 
  else if("CT" %in% colnames(marker)){
    go_list <- split(marker$gene, marker$CT)  
  }
  gsvaPar <- gsvaParam(as.matrix(grid_exp), go_list, maxDiff=maxDiff, kcdf = kcdf_)
  gsva_mat <- gsva(gsvaPar)

  new_row_names <- gsub("\\.", " ", rownames(gsva_mat))

  rownames(gsva_mat) <- new_row_names
  return(gsva_mat)
  
  
}


makeMaskReal <- function(P,threshold_ = 0.02) {
  ## TRUE positions represent the celltype doesn't exists
  return(P <= threshold_)
  
}


makeMask <- function(Spatial.Deconv.Data,gsva_thre=NULL,sd_time = 1) {
  ct_ref <- colnames(Spatial.Deconv.Data$C_ref)
  rownames(Spatial.Deconv.Data$P_train) <- gsub("[^a-zA-Z0-9]",".",rownames(Spatial.Deconv.Data$P_train))
  rownames(Spatial.Deconv.Data$gsva_mat.real) <- gsub("[^a-zA-Z0-9]",".",rownames(Spatial.Deconv.Data$gsva_mat.real))
  mask.True.Pseudo <- makeMaskReal(Spatial.Deconv.Data$P_train[ct_ref,],threshold_ = 0) 
  gsva_es <- Spatial.Deconv.Data$gsva_mat.real
  if(is.null(gsva_thre)){
    gsva_thre = mean(gsva_es) - sd_time * sd(gsva_es)
  }
  mask.Real <- processGSVAMASK(gsva_es[ct_ref,],gsva_thre)
  MASK <- cbind(mask.True.Pseudo,mask.Real)
  Spatial.Deconv.Data$MASK <- MASK
  Spatial.Deconv.Data$MASK[MASK == TRUE] <- 0
  Spatial.Deconv.Data$MASK[MASK == FALSE] <- 1
  return(Spatial.Deconv.Data)
  
}


findBestThre <- function(gsva_es, P_test, threshshold) {
  if(!is.null(P_test)){
  threshsholdBest = 0
  PortionBest = 0
  for (i in threshshold) {
    MASK <- processGSVAMASK(gsva_es, i)
    Portion = getPortion(MASK,P_test)
    if(Portion > PortionBest){
      PortionBest = Portion
      threshsholdBest = i
    }
  }
  
  MASK <- processGSVAMASK(gsva_es, threshsholdBest)
  Portion = getPortion(MASK,P_test)
  return(list(Portion = Portion, threshsholdBest = threshsholdBest))
  }else {
    return(list(Portion = -1, threshsholdBest = 0))  
  }
  
  
}



processGSVAMASK <- function(gsva_es,threshshold) {
  
  return(gsva_es < threshshold)
  
}


getPortion <- function(MASK,P) {
  if(length(rownames(MASK)) !=  length(rownames(P))){
    stop("length(rownames(MASK) !=  length(rownames(P)))")
    return(0)
  }else{
    return(sum((MASK) == (P <= 0.02))/(dim(P)[1] * dim(P)[2]))
  }
  
}


