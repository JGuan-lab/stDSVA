


library(gtools) 
library(sva)
library(scater)
library(reticulate)



source("CDSC.R")
source("function_help.R", encoding = "utf-8")
source("function_help_new.R", encoding = "utf-8")
source("pythonFuc.R")
source("SimulationFuc.R")
source("SampleSimilarity.R")
source("GeneSimilarity.R")
source("GSVAMy.R")
source("BatchEffectCorrection.R")


#' Performs grid search for hyperparameter tuning in spatial transcriptomics deconvolution
#'
#' @description
#' This function performs grid search across multiple hyperparameter combinations
#' to optimize spatial transcriptomics (ST) deconvolution performance.
#' It takes both ST data and reference single-cell (scRNA-seq) data as inputs,
#' optionally aligning coordinates when available,s and returns labeled
#' deconvolution performance metrics for parameter selection.
#'
#' @param scData Reference single-cell dataset.
#' @param meta Metadata corresponding to \code{scData}.
#' @param ST Target spatial transcriptomics dataset for deconvolution.
#' @param windowTrain Training window size (for generate pseudo-spots). Default is 1000.
#' @param lambda1 Vector of candidate values for regularization parameter λ₁.
#' @param lambda2 Vector of candidate values for regularization parameter λ₂.
#' @param lambdaC Vector of candidate values for regularization parameter λC.
#' @param num_iter_max Maximum number of iterations. Default is 3000.
#' @param constraints_l Vector of GSVA-based constraint strength parameters.
#' @param seedd Random seed for reproducibility. Default is 44.
#' @param TerCondition Convergence threshold for optimization. Default is 1e-8.
#' @param plotBatch Logical, whether to visualize batch correction results. Default is TRUE.
#' @param gsva_thre Threshold for GSVA-based filtering. Default is -0.05.
#' @param batch_corr Logical, whether to perform batch correction. Default is TRUE.
#'
#' @return
#' A list or data frame containing deconvolution performance metrics
#' for each parameter combination, which can be used for model selection
#' and hyperparameter tuning.
#'
#' @examples
#' # Example usage:
#' result <- getParas(scData, meta, ST,
#'                    lambda1 = c(0.001, 0.01),
#'                    lambda2 = c(0.1, 1),
#'                    lambdaC = c(1, 10),
#'                    constraints_l = c(0.001, 1, 1000))
#'
#' @export



getParas <- function(scData, meta, ST, windowTrain = 1000,lambda1 = c(0, 0.001, 0.01, 0.1, 1, 10),
                     lambda2 = c(0, 0.001, 0.01, 0.1, 1, 10),
                     lambdaC = c(0, 0.1, 1, 10, 100, 1000),
                     num_iter_max = 3000,
                     constraints_l = c(0.001,1,1000),
                     seedd = 44,
                     TerCondition = 10^-8,plotBatch = T, gsva_thre = NULL, batch_corr = T
){
  comData <- NULL
  comData$Train = list(data=as.matrix(scData),pData=data.frame(meta))
  
  ###1.simulate pseudo spots
  if(is.null(meta$X)){
    comData <- simulateSpotsSampling(comData,cellNum=1:20,n=ncol(ST),UMI = max(colSums(ST)))
  }else{
    comData <- simulateSpotsGridPseudo(comData,window = windowTrain)
  }
  ###2.calculate ref C & find marker genes
  traindata = as.matrix(scData)
  colnames(traindata) <- meta$cellType
  markers <- Find_markerGene_limma(traindata, plotmarker = F, norm1 = "cpm", log2.threshold = 1)
  markerGenes = unique(markers$gene)
  C_ref <- Get_C(traindata,norm = "cpm")[["C_ref"]][markerGenes,]
  ###3.batch correction

  B1 <- comData$T_train
  B2 <- ST
  ###NEED COUNTS DATA
  if(batch_corr){
      gene <- intersect(rownames(B1),rownames(B2))
  B1 <- B1[gene,]
  B2 <- B2[gene,]
  data_combine <- cbind(B1,B2)
  if(plotBatch){
    mnnPlot(data_combine,B1,B2,norm = "None")
  } 
  batch <- c(rep(1, ncol(B1)), rep(2, ncol(B2))) 
  exprCorr <- ComBatSeq_My(data_combine, batch=batch)
  if(plotBatch){
    mnnPlot(exprCorr$counts,B1,B2,norm = "None")
  } 
  B1_corr.marker <- exprCorr$counts[markerGenes,1:ncol(B1)]
  B2_corr.marker <- exprCorr$counts[markerGenes,-(1:ncol(B1))]
  BatchCorr <- calculateCPM(exprCorr$counts)[markerGenes,]
  
  
  }else {
    B1_corr.marker <- B1[markerGenes,]
    B2_corr.marker <- B2[markerGenes,]
    BatchCorr <- cbind(calculateCPM(B1),calculateCPM(B2))[markerGenes,]
    exprCorr <- NULL
    exprCorr$counts <- cbind(B1,B2)
    exprCorr$gamma <- NULL
    
  }

  if(plotBatch){
    mnnPlot(BatchCorr,B1,B2,norm = "None")
  } 
  
  
  location.real <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(ST),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(ST),split="_"),"[",3)))
  location.pseudo <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(comData$T_train),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(comData$T_train),split="_"),"[",3)))
  rownames(location.real) <- colnames(ST)
  rownames(location.pseudo) <- colnames(comData$T_train)
  if(batch_corr & !is.null(comData$Train$X)){
    library(reticulate)
    use_condaenv("D:/softwares/anaconda/envs/gst/python.exe")
    source_python("PasteMy.py")  
    slice_pseudo <- data.frame(t(B1_corr.marker))
    slice_real <- data.frame(t(B2_corr.marker))
    slice_pseudo_coor <- data.frame(location.pseudo)
    slice_real_coor <- data.frame(location.real)
    rownames(slice_pseudo) <- paste(slice_pseudo_coor$x,"x",slice_pseudo_coor$y,sep = "")
    rownames(slice_real) <- paste(slice_real_coor$x,"x",slice_real_coor$y,sep = "")
    slice_pseudo_coor <- as.matrix(slice_pseudo_coor)
    slice_real_coor <- as.matrix(slice_real_coor)
    spatial = py$pasteGetSpatial(slice_pseudo,slice_real,slice_pseudo_coor,slice_real_coor)
    location.pseudo <- spatial[1:nrow(location.pseudo),]
    location.real <- spatial[(nrow(location.pseudo)+1):nrow(spatial),]
    rownames(location.real) <- colnames(ST)
    rownames(location.pseudo) <- colnames(comData$T_train)
  }
  
  ####calculate Ss Sg and gsva matrix
  if(is.null(comData$Train$pData$X)){
    Ss <- SsInter(BatchCorr,location.real)
  }else{
    Ss <- SsInter(BatchCorr,rbind(location.pseudo,location.real))
    
  }
  Sg <- SgIntra(BatchCorr)
  
  set.seed(seedd)  
  n <- ncol(comData$P_train)
  train_indices <- sample(1:n, size = round(0.7 * n))  
  train_data <- comData$P_train[, train_indices]
  test_data <- comData$P_train[,-train_indices]
  
  gsva_mat.real <- gsva_apply(cbind(B1_corr.marker[,-train_indices], B2_corr.marker), markers)
  mask <- matrix(runif(ncol(C_ref)*ncol(BatchCorr)) ,ncol(C_ref),ncol(BatchCorr))
  colnames(mask) <- colnames(BatchCorr)
  rownames(mask) <- colnames(C_ref)
  if(is.null(gsva_thre)) {
    gsva_thre = mean(gsva_mat.real) - sd(gsva_mat.real)
  }
  mask[,colnames(B2)] <- ifelse(gsva_mat.real[,colnames(B2)] > gsva_thre,1,0)
  mask[,colnames(test_data)] <- ifelse(gsva_mat.real[,colnames(test_data)] > gsva_thre,1,0)
  mask[,colnames(train_data)] <- ifelse(comData$P_train[,colnames(train_data)] > 0,1,0)
    ####find parameters
  result1 = NULL
  #------------cyclic to do deconvolution----------
  pb <- txtProgressBar(style = 3)
  star_time <- Sys.time()
  result_para_c = NULL
  result_para_p = NULL
  number_iter = NULL
  para_find = NULL 
  num = 1
  library(scater)
  library(dplyr)
  for(l_i in 1:length(constraints_l)){
    for (dir_i in 1:length(lambda1)){
      for (dir_j in 1:length(lambda2)){
        for (dir_k in 1:length(lambdaC)){
          result_CDSC = CDSC(data_bulk = BatchCorr, data_ref = as.matrix(C_ref), k = dim(C_ref)[2], MASK = mask, l = constraints_l[l_i],  P.Pseudo = train_data,
                                         lambda1 = lambda1[dir_i], lambda2 = lambda2[dir_j], lambdaC = lambdaC[dir_k],
                                         error = TerCondition, seedd = seedd,
                                         Ss = as.matrix(Ss), Sg = as.matrix(Sg), all_number = num_iter_max)
          result_para_c <- result_CDSC[[1]]
          if(!is.null(exprCorr$gamma)){
            result_para_c <-  exp(log(result_para_c) + vec2mat(exprCorr$gamma[markerGenes,][,c("batch1")], ncol(result_para_c)))
          }
          if(all(is.na(result_CDSC[[1]])) || all(is.na(result_CDSC[[2]]))){
            next
          }
          result_para_p <- NORM(result_CDSC[[2]])
          result_para_p_pre <- result_para_p[,1:ncol(comData$P_train)]
          number_iter <- result_CDSC[[3]]
          ct_labels <- Row_label(result_para_c,C_ref)
          rownames(result_para_p_pre) <- ct_labels
          colnames(result_para_c) <- ct_labels
          index_p <- getPearsonRMSE(result_para_p_pre[,colnames(test_data)],test_data)
          index_c <- getPearsonRMSE(result_para_c,C_ref)
          para_find <- cbind(index_p,index_c)
          para_find <- cbind(para_find,number_iter)
          para_find <- cbind(lambdaC[dir_k],para_find)
          para_find <- cbind(lambda2[dir_j],para_find)
          para_find <- cbind(lambda1[dir_i],para_find)
          para_find <- cbind(constraints_l[l_i],para_find)
          colnames(para_find) <- c("l","l1","l2","lC","RMSE_of_Ptrain","PCC_of_Ptrain","RMSE_of_Cref","PCC_of_Cref","num_iter")
          result1 <- rbind(result1,para_find)
          num = num + 1
          print(para_find)
          setTxtProgressBar(pb, num/(length(lambda1)*length(lambda2)*length(lambdaC)*length(constraints_l)))
        }}}}
  
  return(list(C_ref = as.matrix(C_ref), ST = BatchCorr, Ss = as.matrix(Ss),Sg = as.matrix(Sg), mask = mask, paras = result1, P.Pseudo = train_data, gamma = exprCorr$gamma[markerGenes,]))
  
}

###---------------deconvoltuion--------------------------------------------------------
runDeconv_my <- function(dataRDS,constraints_l = 0.001,lambda1 = 0.001, lambda2 = 0.001, lambdaC = 1000, TerCondition = 10^-8, seedd = 44, num_iter_max = 3000){
  result_CDSC = CDSC(data_bulk = dataRDS$ST, data_ref = as.matrix(dataRDS$C_ref), k = dim(dataRDS$C_ref)[2], MASK = dataRDS$mask, l = constraints_l,  
                                 lambda1 = lambda1, lambda2 = lambda2, lambdaC = lambdaC,
                                 error = TerCondition, seedd = seedd,
                                 Ss = dataRDS$Ss, Sg = dataRDS$Sg, all_number = num_iter_max,P.Pseudo = dataRDS$P.Pseudo)
  
  result_para_c <- result_CDSC[[1]]
  result_para_p <- NORM(result_CDSC[[2]])
  number_iter <- result_CDSC[[3]]
  ct_labels <- Row_label(result_para_c,dataRDS$C_ref)
  rownames(result_para_p) <- ct_labels
  colnames(result_para_c) <- ct_labels 
  if(!is.null(dataRDS$gamma)){
    result_para_c <-  exp(log(result_para_c) + vec2mat(dataRDS$gamma[rownames(dataRDS$ST),][,c("batch2")], ncol(result_para_c)))
  }
  return(list(c = result_para_c,p = result_para_p,number_iter = number_iter))
  
}



