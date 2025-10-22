graphSTPasteCorrMy <- function(Spatial.Deconv.Data, sc = TRUE,plot=F,realTest=F,norm=NULL) {
  print("GraphST begins..................")
  library(reticulate)
  use_condaenv("D:/softwares/anaconda/envs/gst/python.exe")
  source_python("H:/LinYifan/spatialDSSC-BE/R/GraphSTMy.py")  
  Spatial.Deconv.Data$comData <- simulateSpots(Spatial.Deconv.Data$comData,window = Spatial.Deconv.Data$comData$windowTest)
  Spatial.Deconv.Data$comData <- simulateSpotsGridPseudo(Spatial.Deconv.Data$comData,window = Spatial.Deconv.Data$comData$windowTrain)


  Spatial.Deconv.Data$spatial.pseudo.count.marker <- Spatial.Deconv.Data$comData$T_train[Spatial.Deconv.Data$markerGenes,]
  Spatial.Deconv.Data$spatial.real.count.marker <- Spatial.Deconv.Data$comData$T_test[Spatial.Deconv.Data$markerGenes,]
  Spatial.Deconv.Data$location.pseudo <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_train),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_train),split="_"),"[",3)))
  Spatial.Deconv.Data$location.real <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_test),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_test),split="_"),"[",3)))
  
  slice_pseudo <- data.frame(t(Spatial.Deconv.Data$spatial.pseudo.count.marker))
  slice_real <- data.frame(t(Spatial.Deconv.Data$spatial.real.count.marker))
  slice_pseudo_coor <- data.frame(Spatial.Deconv.Data$location.pseudo)
  slice_real_coor <- data.frame(Spatial.Deconv.Data$location.real)
  data_combine <- cbind(Spatial.Deconv.Data$spatial.pseudo.count.marker,Spatial.Deconv.Data$spatial.real.count.marker)
  if(plot){
    mnnPlot(data_combine,Spatial.Deconv.Data$spatial.pseudo.count.marker,Spatial.Deconv.Data$spatial.real.count.marker)
  }
  if(is.null(slice_pseudo) || is.null(slice_real) || is.null(slice_pseudo_coor) || is.null(slice_real_coor)) {
    stop("Error! The pseudo data or the real data is not complete!")
  }
  rownames(slice_pseudo) <- paste(slice_pseudo_coor$x,"x",slice_pseudo_coor$y,sep = "")
  rownames(slice_real) <- paste(slice_real_coor$x,"x",slice_real_coor$y,sep = "")
  slice_pseudo_coor <- as.matrix(slice_pseudo_coor)
  slice_real_coor <- as.matrix(slice_real_coor)
  emb = py$graphSTPasteMy(slice_pseudo,slice_real,slice_pseudo_coor,slice_real_coor)
  Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- emb[[1]]
  Spatial.Deconv.Data$BatchCorr$GraphST_embedding[Spatial.Deconv.Data$BatchCorr$GraphST_embedding < 0]  <- 0
  if(norm == "minmax"){
    Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- (emb[[1]] - min(emb[[1]])) / (max(emb[[1]]) - min(emb[[1]]))
    
  }else if(norm == "zero"){
  }else if(norm == "add") {
    Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- emb[[1]] + abs(min(emb[[1]]))
  }
  Spatial.Deconv.Data$BatchCorr$GraphST_spatial <- emb[[2]]
  
  Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- t(Spatial.Deconv.Data$BatchCorr$GraphST_embedding)
  colnames(Spatial.Deconv.Data$BatchCorr$GraphST_embedding) <- colnames(data_combine)
  rownames(Spatial.Deconv.Data$BatchCorr$GraphST_spatial) <- colnames(data_combine)
  
  if(plot){
    mnnPlot( Spatial.Deconv.Data$BatchCorr$GraphST_embedding,Spatial.Deconv.Data$spatial.pseudo.count.marker,Spatial.Deconv.Data$spatial.real.count.marker)
  }
  
  Spatial.Deconv.Data$P_train <- Spatial.Deconv.Data$comData$P_train
  Spatial.Deconv.Data$P_test <- Spatial.Deconv.Data$comData$P_test
  
  Spatial.Deconv.Data$spatial.real <- Spatial.Deconv.Data$comData$T_test
  Spatial.Deconv.Data$spatial.pseudo <- Spatial.Deconv.Data$comData$T_train
  Spatial.Deconv.Data$location.pseudo <- Spatial.Deconv.Data$BatchCorr$GraphST_spatial[colnames(Spatial.Deconv.Data$comData$T_train),]
  Spatial.Deconv.Data$location.real <- Spatial.Deconv.Data$BatchCorr$GraphST_spatial[colnames(Spatial.Deconv.Data$comData$T_test),]
  Spatial.Deconv.Data$T_combine <- Spatial.Deconv.Data$BatchCorr$GraphST_embedding
  Spatial.Deconv.Data$data_loc_combine.gst <- Spatial.Deconv.Data$BatchCorr$GraphST_spatial
  Spatial.Deconv.Data$data_combine.gst <- Spatial.Deconv.Data$BatchCorr$GraphST_embedding
  
  return(Spatial.Deconv.Data)
}

pasteGraphSTCorrMy <- function(Spatial.Deconv.Data, sc = TRUE,plot=F,realTest=F,norm=NULL) {
  print("GraphST begins..................")
  library(reticulate)
  use_condaenv("D:/softwares/anaconda/envs/gst/python.exe")
  source_python("H:/LinYifan/spatialDSSC-BE/R/GraphSTMy.py")  
  
  if(sc) { ## use sinle cell to do graphst
    slice_pseudo <- data.frame(t(Spatial.Deconv.Data$comData$Train$data[Spatial.Deconv.Data$markerGenes,]))
    slice_real <- data.frame(t(Spatial.Deconv.Data$comData$Test$data[Spatial.Deconv.Data$markerGenes,]))
    slice_pseudo_coor <- cbind.data.frame(x = as.numeric(Spatial.Deconv.Data$comData$Train$pData$X), y = as.numeric(Spatial.Deconv.Data$comData$Train$pData$Y))
    slice_real_coor <- cbind.data.frame(x = as.numeric(Spatial.Deconv.Data$comData$Test$pData$X), y = as.numeric(Spatial.Deconv.Data$comData$Test$pData$Y))
    data_combine <- cbind(Spatial.Deconv.Data$comData$Train$data[Spatial.Deconv.Data$markerGenes,],Spatial.Deconv.Data$comData$Test$data[Spatial.Deconv.Data$markerGenes,])
    if(plot){
      mnnPlot(data_combine,Spatial.Deconv.Data$comData$Train$data[Spatial.Deconv.Data$markerGenes,],Spatial.Deconv.Data$comData$Test$data[Spatial.Deconv.Data$markerGenes,])
    }
    
  } else{
    if(realTest){
      ####realTest,Xreal的合成
      Spatial.Deconv.Data$comData <- simulateSpotsGridPseudo(Spatial.Deconv.Data$comData,window = Spatial.Deconv.Data$comData$windowTrain)
      Spatial.Deconv.Data$comData$T_test <- Spatial.Deconv.Data$comData$Test$data
      Spatial.Deconv.Data$comData$P_test <- NULL
    }else{
      Spatial.Deconv.Data$comData <- simulateSpots(Spatial.Deconv.Data$comData,window = Spatial.Deconv.Data$comData$windowTest)
      Spatial.Deconv.Data$comData <- simulateSpotsGridPseudo(Spatial.Deconv.Data$comData,window = Spatial.Deconv.Data$comData$windowTrain)
    }

    Spatial.Deconv.Data$spatial.pseudo.count.marker <- Spatial.Deconv.Data$comData$T_train[Spatial.Deconv.Data$markerGenes,]
    Spatial.Deconv.Data$spatial.real.count.marker <- Spatial.Deconv.Data$comData$T_test[Spatial.Deconv.Data$markerGenes,]
    Spatial.Deconv.Data$location.pseudo <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_train),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_train),split="_"),"[",3)))
    Spatial.Deconv.Data$location.real <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_test),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_test),split="_"),"[",3)))
      
    slice_pseudo <- data.frame(t(Spatial.Deconv.Data$spatial.pseudo.count.marker))
    slice_real <- data.frame(t(Spatial.Deconv.Data$spatial.real.count.marker))
    slice_pseudo_coor <- data.frame(Spatial.Deconv.Data$location.pseudo)
    slice_real_coor <- data.frame(Spatial.Deconv.Data$location.real)
    data_combine <- cbind(Spatial.Deconv.Data$spatial.pseudo.count.marker,Spatial.Deconv.Data$spatial.real.count.marker)
    if(plot){
      mnnPlot(data_combine,Spatial.Deconv.Data$spatial.pseudo.count.marker,Spatial.Deconv.Data$spatial.real.count.marker)
    }
    
  }

  if(is.null(slice_pseudo) || is.null(slice_real) || is.null(slice_pseudo_coor) || is.null(slice_real_coor)) {
    stop("Error! The pseudo data or the real data is not complete!")
  }
  rownames(slice_pseudo) <- paste(slice_pseudo_coor$x,"x",slice_pseudo_coor$y,sep = "")
  rownames(slice_real) <- paste(slice_real_coor$x,"x",slice_real_coor$y,sep = "")
  slice_pseudo_coor <- as.matrix(slice_pseudo_coor)
  slice_real_coor <- as.matrix(slice_real_coor)
  emb = py$pasteGraphSTMy(slice_pseudo,slice_real,slice_pseudo_coor,slice_real_coor)
  Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- emb[[1]]
  Spatial.Deconv.Data$BatchCorr$GraphST_embedding[Spatial.Deconv.Data$BatchCorr$GraphST_embedding < 0]  <- 0
  if(norm == "minmax"){
    Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- (emb[[1]] - min(emb[[1]])) / (max(emb[[1]]) - min(emb[[1]]))
    
  }else if(norm == "zero"){
  }else if(norm == "add") {
    Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- emb[[1]] + abs(min(emb[[1]]))
  }
  Spatial.Deconv.Data$BatchCorr$GraphST_spatial <- emb[[2]]
  
  Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- t(Spatial.Deconv.Data$BatchCorr$GraphST_embedding)
  colnames(Spatial.Deconv.Data$BatchCorr$GraphST_embedding) <- colnames(data_combine)
  rownames(Spatial.Deconv.Data$BatchCorr$GraphST_spatial) <- colnames(data_combine)
  
  if(sc){
    ###change to original cell ids
    cellids <- colnames(data_combine)
    pdata <- rbind(Spatial.Deconv.Data$comData$Train$pData,Spatial.Deconv.Data$comData$Test$pData)
    locName <- paste(pdata$X,pdata$Y,sep = "x")
    Spatial.Deconv.Data$BatchCorr$GraphST_spatial <- as.matrix(Spatial.Deconv.Data$BatchCorr$GraphST_spatial)
    colnames(Spatial.Deconv.Data$BatchCorr$GraphST_embedding) <- colnames(data_combine)
    rownames(Spatial.Deconv.Data$BatchCorr$GraphST_spatial) <- colnames(data_combine)
    if(plot){
    mnnPlot(Spatial.Deconv.Data$BatchCorr$GraphST_embedding,Spatial.Deconv.Data$comData$Train$data[Spatial.Deconv.Data$markerGenes,],Spatial.Deconv.Data$comData$Test$data[Spatial.Deconv.Data$markerGenes,])
    }
  }else {
    if(plot){
      mnnPlot( Spatial.Deconv.Data$BatchCorr$GraphST_embedding,Spatial.Deconv.Data$spatial.pseudo.count.marker,Spatial.Deconv.Data$spatial.real.count.marker)
    }
    
    Spatial.Deconv.Data$P_train <- Spatial.Deconv.Data$comData$P_train
    Spatial.Deconv.Data$P_test <- Spatial.Deconv.Data$comData$P_test
    
    Spatial.Deconv.Data$spatial.real <- Spatial.Deconv.Data$comData$T_test
    Spatial.Deconv.Data$spatial.pseudo <- Spatial.Deconv.Data$comData$T_train
    Spatial.Deconv.Data$location.pseudo <- Spatial.Deconv.Data$BatchCorr$GraphST_spatial[colnames(Spatial.Deconv.Data$comData$T_train),]
    Spatial.Deconv.Data$location.real <- Spatial.Deconv.Data$BatchCorr$GraphST_spatial[colnames(Spatial.Deconv.Data$comData$T_test),]
    Spatial.Deconv.Data$T_combine <- Spatial.Deconv.Data$BatchCorr$GraphST_embedding
    Spatial.Deconv.Data$data_loc_combine.gst <- Spatial.Deconv.Data$BatchCorr$GraphST_spatial
    Spatial.Deconv.Data$data_combine.gst <- Spatial.Deconv.Data$BatchCorr$GraphST_embedding
      }
  # Spatial.Deconv.Data$BatchCorr$GraphST_spatial <- as.data.frame(Spatial.Deconv.Data$BatchCorr$GraphST_spatial)

  return(Spatial.Deconv.Data)  
  
}


pasteAlignment <- function(Spatial.Deconv.Data) {
  ####correct the coordination using paste and gst env
  print("paste alignment begins..................")
  use_condaenv("D:/softwares/anaconda/envs/gst/python.exe")
  source_python("H:/LinYifan/spatialDSSC-BE/R/PasteMy.py")  
  slice_pseudo <- data.frame(t(Spatial.Deconv.Data$spatial.pseudo.count.marker))
  slice_real <- data.frame(t(Spatial.Deconv.Data$spatial.real.count.marker))
  slice_pseudo_coor <- data.frame(Spatial.Deconv.Data$location.pseudo)
  slice_real_coor <- data.frame(Spatial.Deconv.Data$location.real)
  rownames(slice_pseudo) <- paste(slice_pseudo_coor$x,"x",slice_pseudo_coor$y,sep = "")
  rownames(slice_real) <- paste(slice_real_coor$x,"x",slice_real_coor$y,sep = "")
  slice_pseudo_coor <- as.matrix(slice_pseudo_coor)
  slice_real_coor <- as.matrix(slice_real_coor)
  spatial = py$pasteGetSpatial(slice_pseudo,slice_real,slice_pseudo_coor,slice_real_coor)
  Spatial.Deconv.Data$location.pseudo.original <- Spatial.Deconv.Data$location.pseudo
  Spatial.Deconv.Data$location.real.original <- Spatial.Deconv.Data$location.real
  Spatial.Deconv.Data$location.pseudo <- spatial[1:nrow(Spatial.Deconv.Data$location.pseudo),]
  Spatial.Deconv.Data$location.real <- spatial[(nrow(Spatial.Deconv.Data$location.pseudo)+1):nrow(spatial),]
  rownames(Spatial.Deconv.Data$location.real) <- colnames(Spatial.Deconv.Data$comData$T_test)
  rownames(Spatial.Deconv.Data$location.pseudo) <- colnames(Spatial.Deconv.Data$comData$T_train)
  
  return(Spatial.Deconv.Data)
}


graphSTCorrection <- function(Spatial.Deconv.Data) {
  ###correct expression batch effect for spatial data using GraphST
  ###
  
  library(reticulate)
  use_condaenv("D:/softwares/anaconda/envs/gst/python.exe")
  source_python("H:/LinYifan/spatialDSSC-BE/R/GraphSTMy.py")
  anndata <- pasteAlignment(Spatial.Deconv.Data)
  emb <- graphSTMy(anndata)###emb[[1]]: expr after correction emb[[2]]: spatial coor after correction
  Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- emb[[1]]
  Spatial.Deconv.Data$BatchCorr$GraphST_spatial <- emb[[2]]
  
  Spatial.Deconv.Data$BatchCorr$GraphST_embedding <- t(Spatial.Deconv.Data$BatchCorr$GraphST_embedding)
  
  rownames(Spatial.Deconv.Data$BatchCorr$GraphST_spatial) <- colnames(Spatial.Deconv.Data$BatchCorr$GraphST_spatial)
  colnames(Spatial.Deconv.Data$BatchCorr$GraphST_spatial) <- c("x","y")
  return(Spatial.Deconv.Data)
  
}

mnnPlot <- function(data_combine,B1,B2,norm = "None") {
  library(umap)
  library(ggplot2)
  library(scater)
  if(norm %in%  c("cpm","CPM")){
    data_combine <- calculateCPM(data_combine)
  }
  set.seed(123)
  umap <- umap(t(data_combine),method='naive',n_neighbors = 10)
  df1 <- data.frame(umap$layout)
  df1$label <- c(rep("Pseudo",ncol(B1)), rep("Real",ncol(B2)))
  colnames(df1) <- c('X','Y','label') 
  p <- ggplot(df1, aes(x=X, y=Y, colour=label)) + geom_point(size=4)+
    xlab(NULL)+ 
    ylab(NULL)
  p <- p + theme(  panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.title=element_blank(), #图例标签隐藏
                   panel.border = element_blank(),
                   axis.line.x = element_line(color="black", size = 0.5),
                   axis.line.y = element_line(color="black", size = 0.5),
                   panel.background = element_blank())
  p <- p +  scale_colour_manual(values=c("#D58F25","#882424","#7DA92E","#5A7EB3"))
  print(p)
  
  
}
mnnCombatCorrection <- function(Spatial.Deconv.Data,norm_plot = "None", cos.norm.in = TRUE, cos.norm.out = FALSE,window = 10, plot = FALSE, onlySim = F, norm_corr = "zero",useCombat=F) {
  library(batchelor)
  B1 <- Spatial.Deconv.Data$comData$Train$data
  B2 <- Spatial.Deconv.Data$comData$Test$data
  B1_cpm <- calculateCPM(B1)
  B1_log <- log2(B1_cpm + 1)
  
  interGenes <- intersect(rownames(B1),rownames(B2))
  B1_marker <- B1_log[interGenes,]
  B2_cpm <- calculateCPM(B2)
  B2_log <- log2(B2_cpm + 1)
  B2_marker <- B2_log[interGenes,]
  
  B1 <- as.matrix(B1_marker)
  B2 <- as.matrix(B2_marker)
  data_combine <- cbind(B1,B2)
  if(plot){
    mnnPlot(data_combine,B1,B2,norm = norm_plot)
  }
  
  
  
  m.out <- mnnCorrect(B1, B2, cos.norm.in = cos.norm.in,cos.norm.out = cos.norm.out)
  exprCorr <- assay(m.out, "corrected")
  if(norm_corr == "minmax"){
    exprCorr <- (exprCorr - min(exprCorr)) / (max(exprCorr) - min(exprCorr))
    
  }else if(norm_corr == "zero"){
    exprCorr[exprCorr < 0] = 0
  }else if(norm_corr == "add") {
    exprCorr <- exprCorr + abs(min(exprCorr))
  }
  
  if(plot){
    mnnPlot(exprCorr,B1,B2,norm = norm_plot)
  }
  ###calculate c
  comData_ <- Spatial.Deconv.Data$comData
  comData_$Train$data <- exprCorr[,1:ncol(B1)]
  comData_$Test$data <- exprCorr[,(ncol(B1) + 1):ncol(exprCorr)]
  getcTrain <- comData_$Train$data
  getcTest <- comData_$Test$data
  colnames(getcTrain) <- comData_$Train$pData$cellType
  colnames(getcTest) <- comData_$Test$pData$cellType
  if(!onlySim){
    Spatial.Deconv.Data$C_ref <- as.matrix(Get_C(getcTrain,norm = "None")[['C_ref']])[Spatial.Deconv.Data$markerGenes,]
    Spatial.Deconv.Data$C <- as.matrix(Get_C(getcTest, norm = "None")[['C_ref']])[Spatial.Deconv.Data$markerGenes,]
    
  }
  ###simulate
  Spatial.Deconv.Data$comData <- simulateSpots(Spatial.Deconv.Data$comData,window=window)
  Spatial.Deconv.Data$comData <- simulateSpotsSampling(Spatial.Deconv.Data$comData)
  ###cpm & log
  B1 <- as.matrix(log2(calculateCPM(Spatial.Deconv.Data$comData$T_train)+1)[Spatial.Deconv.Data$markerGenes,])
  B2 <- as.matrix(log2(calculateCPM(Spatial.Deconv.Data$comData$T_test)+1)[Spatial.Deconv.Data$markerGenes,])
  
  ###combat
  
  
  message("Use combat.............")
  library(sva)
  library(scater)
  set.seed(44)
  

  data_combine <- cbind(B1,B2)
  
  mnnPlot(data_combine,B1,B2,norm = norm_plot)
  
  batch <- c(rep(1, ncol(B1)), rep(2, ncol(B2))) 
  exprCorr <- ComBat(data_combine, batch=batch)
  Spatial.Deconv.Data$BatchCorr <- exprCorr
  exprCorr[exprCorr < 0] <- 0
  mnnPlot(exprCorr,B1,B2,norm = norm_plot)
 
  Spatial.Deconv.Data$T_combine <- exprCorr
  
  Spatial.Deconv.Data$spatial.pseudo <- Spatial.Deconv.Data$comData$T_train[Spatial.Deconv.Data$markerGenes,]
  Spatial.Deconv.Data$spatial.real <- Spatial.Deconv.Data$comData$T_test[Spatial.Deconv.Data$markerGenes,]
  Spatial.Deconv.Data$location.real <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_test),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_test),split="_"),"[",3)))
  rownames(Spatial.Deconv.Data$location.real) <- colnames(Spatial.Deconv.Data$comData$T_test)
  Spatial.Deconv.Data$P_train <- Spatial.Deconv.Data$comData$P_train
  Spatial.Deconv.Data$P_test <- Spatial.Deconv.Data$comData$P_test
  return(Spatial.Deconv.Data)
  
}
##window is test simulation win
mnnCorrection <- function(Spatial.Deconv.Data,norm_plot = "None", cos.norm.in = TRUE, cos.norm.out = FALSE,window = 10, plot = FALSE, onlySim = F, norm_corr = "zero",useCombat=F) {
  ##correct expression batch effect for bulk data using MNN
  library(batchelor)
  ##use single cell
  B1 <- Spatial.Deconv.Data$comData$Train$data
  B2 <- Spatial.Deconv.Data$comData$Test$data
  
  # B1 <- Spatial.Deconv.Data$spatial.pseudo
  # B2 <- Spatial.Deconv.Data$spatial.real
  ##input data with cpm and log2-transformation
  B1_cpm <- calculateCPM(B1)
  B1_log <- log2(B1_cpm + 1)
  
  interGenes <- intersect(rownames(B1),rownames(B2))
  B1_marker <- B1_log[interGenes,]
  B2_cpm <- calculateCPM(B2)
  B2_log <- log2(B2_cpm + 1)
  B2_marker <- B2_log[interGenes,]
  
  B1 <- as.matrix(B1_marker)
  B2 <- as.matrix(B2_marker)
  data_combine <- cbind(B1,B2)
  if(plot){
    mnnPlot(data_combine,B1,B2,norm = norm_plot)
  }
  
  

  m.out <- mnnCorrect(B1, B2, cos.norm.in = cos.norm.in,cos.norm.out = cos.norm.out)
  exprCorr <- assay(m.out, "corrected")
  if(norm_corr == "minmax"){
    exprCorr <- (exprCorr - min(exprCorr)) / (max(exprCorr) - min(exprCorr))

  }else if(norm_corr == "zero"){
    exprCorr[exprCorr < 0] = 0
  }else if(norm_corr == "add") {
    exprCorr <- exprCorr + abs(min(exprCorr))
  }
  
  if(plot){
    mnnPlot(exprCorr,B1,B2,norm = norm_plot)
  }
  
  comData_ <- Spatial.Deconv.Data$comData
  comData_$Train$data <- exprCorr[,1:ncol(B1)]
  comData_$Test$data <- exprCorr[,(ncol(B1) + 1):ncol(exprCorr)]
  Spatial.Deconv.Data$comData <- simulateSpots(Spatial.Deconv.Data$comData, win = window)
  comData_ <- simulateSpots(comData_, win = window)
  comData_ <- simulateSpotsSampling(comData_)
  
  Spatial.Deconv.Data$comData <- simulateSpotsSampling(Spatial.Deconv.Data$comData)
  Spatial.Deconv.Data$P_train <- comData_$P_train
  Spatial.Deconv.Data$P_test <- comData_$P_test
  ###batch correct(cpm + log)
  Spatial.Deconv.Data$spatial.pseudo <- comData_$T_train[Spatial.Deconv.Data$markerGenes,]
  Spatial.Deconv.Data$spatial.real <- comData_$T_test[Spatial.Deconv.Data$markerGenes,]
  
  ###calculate c
  getcTrain <- comData_$Train$data
  getcTest <- comData_$Test$data
  colnames(getcTrain) <- comData_$Train$pData$cellType
  colnames(getcTest) <- comData_$Test$pData$cellType
  if(!onlySim){
    Spatial.Deconv.Data$C_ref <- as.matrix(Get_C(getcTrain,norm = "None")[['C_ref']])[Spatial.Deconv.Data$markerGenes,]
    Spatial.Deconv.Data$C <- as.matrix(Get_C(getcTest, norm = "None")[['C_ref']])[Spatial.Deconv.Data$markerGenes,]
    
  }

  Spatial.Deconv.Data$location.real <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(comData_$T_test),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(comData_$T_test),split="_"),"[",3)))
  rownames(Spatial.Deconv.Data$location.real) <- colnames(comData_$T_test)
  Spatial.Deconv.Data$location.pseudo <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(comData_$T_train),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(comData_$T_train),split="_"),"[",3)))
  
  if(!all(is.na(Spatial.Deconv.Data$location.pseudo$x))) {
    rownames(Spatial.Deconv.Data$location.pseudo) <- colnames(comData_$T_train)
  } else {
    print("There are no locations in T_train! continue..........................")
  }
  
  
  
  ###combat
  if(useCombat){
    message("Use combat.............")
    library(sva)
    library(scater)
    set.seed(44)
    B1 <- Spatial.Deconv.Data$spatial.pseudo
    B2 <- Spatial.Deconv.Data$spatial.real
    ##input data with cpm 
    B1_marker <- B1[Spatial.Deconv.Data$markerGenes,]
    B2_marker <- B2[Spatial.Deconv.Data$markerGenes,]
    
    B1 <- B1_marker
    B2 <- B2_marker
    data_combine <- cbind(B1,B2)
    
    mnnPlot(data_combine,B1,B2,norm = norm_plot)
    
    batch <- c(rep(1, ncol(B1)), rep(2, ncol(B2))) 
    exprCorr <- ComBat(data_combine, batch=batch)
    Spatial.Deconv.Data$BatchCorr <- exprCorr
    exprCorr[exprCorr < 0] <- 0
    mnnPlot(exprCorr,B1,B2,norm = norm_plot)
    
    
    
    
    
    Spatial.Deconv.Data$T_combine <- exprCorr

  }else{
    Spatial.Deconv.Data$BatchCorr <- exprCorr[Spatial.Deconv.Data$markerGenes,]
    Spatial.Deconv.Data$T_combine <- cbind(Spatial.Deconv.Data$spatial.pseudo,Spatial.Deconv.Data$spatial.real)
    Spatial.Deconv.Data$spatial.real <- Spatial.Deconv.Data$comData$T_test[Spatial.Deconv.Data$markerGenes,]
    Spatial.Deconv.Data$BatchCorr <- Spatial.Deconv.Data$T_combine
    Spatial.Deconv.Data$spatial.pseudo.cpm.marker <- Spatial.Deconv.Data$spatial.real[Spatial.Deconv.Data$markerGenes,]

  }
  
  
  # # Spatial.Deconv.Data$BatchCorr <- exprCorr[Spatial.Deconv.Data$markerGenes,]
  # Spatial.Deconv.Data$T_combine <- cbind(Spatial.Deconv.Data$spatial.pseudo,Spatial.Deconv.Data$spatial.real)
  Spatial.Deconv.Data$spatial.real <- Spatial.Deconv.Data$comData$T_test[Spatial.Deconv.Data$markerGenes,]
  # Spatial.Deconv.Data$BatchCorr <- Spatial.Deconv.Data$T_combine
  # Spatial.Deconv.Data$spatial.pseudo.cpm.marker <- Spatial.Deconv.Data$spatial.real[Spatial.Deconv.Data$markerGenes,]
  # 
  if(onlySim) {
    Spatial.Deconv.Data$T_combine <- NULL
    Spatial.Deconv.Data$spatial.pseudo <- Spatial.Deconv.Data$comData$T_train[Spatial.Deconv.Data$markerGenes,]

    
     }
  return(Spatial.Deconv.Data)
  
}


combatSeqCorrection <- function(Spatial.Deconv.Data,norm_plot = "cpm") {
  ##correct expression batch effect for bulk data using MNN
  library(sva)
  library(scater)
  set.seed(44)
  B1 <- Spatial.Deconv.Data$spatial.pseudo
  B2 <- Spatial.Deconv.Data$spatial.real
  ##input data with cpm 
  B1_marker <- B1[,]
  B2_marker <- B2[,]
  
  B1 <- B1_marker
  B2 <- B2_marker
  data_combine <- cbind(B1,B2)
  
  mnnPlot(data_combine,B1,B2,norm = norm_plot)
  
  batch <- c(rep(1, ncol(B1)), rep(2, ncol(B2))) 
  exprCorr <- ComBat_seq(data_combine, batch=batch, group=NULL)
  
  mnnPlot(exprCorr,B1,B2,norm = norm_plot)
  
  
  
  
  Spatial.Deconv.Data$BatchCorr <- calculateCPM(exprCorr)[Spatial.Deconv.Data$markerGenes,]
  Spatial.Deconv.Data$T_combine <- Spatial.Deconv.Data$BatchCorr
  return(Spatial.Deconv.Data)
  
}

ComBatSeq_ref <- function (counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE, 
                          shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL, 
                          ref.batch = NULL, BPPARAM = BiocParallel::bpparam("SerialParam")) 
{
  library(edgeR)
  
  batch <- as.factor(batch)
  if (any(table(batch) <= 1)) {
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }
  keep_lst <- lapply(levels(batch), function(b) {
    which(apply(counts[, batch == b], 1, function(x) {
      !all(x == 0)
    }))
  })
  keep <- Reduce(intersect, keep_lst)
  rm <- setdiff(1:nrow(counts), keep)
  countsOri <- counts
  counts <- counts[keep, ]
  dge_obj <- DGEList(counts = counts)
  

  n_batch <- nlevels(batch)
  batches_ind <- lapply(1:n_batch, function(i) {
    which(batch == levels(batch)[i])
  })
  n_batches <- sapply(batches_ind, length)
  n_sample <- sum(n_batches)
  cat("Found", n_batch, "batches\n")
  batchmod <- model.matrix(~-1 + batch)
  if (!is.null(ref.batch)) {
    if (!(ref.batch %in% levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    message("Using batch = ", ref.batch, " as a reference batch (this batch won't change)")
    ref <- which(levels(batch) == ref.batch)
    # batchmod[, ref] <- 1
  } else {
    ref <- NULL
  }
  group <- as.factor(group)
  if (full_mod & nlevels(group) > 1) {
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }
  else {
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data = as.data.frame(t(counts)))
  }
  if (!is.null(covar_mod)) {
    if (is.data.frame(covar_mod)) {
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), 
                                         function(i) {
                                           model.matrix(~covar_mod[, i])
                                         }))
    }
    covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x) {
      all(x == 1)
    })]
  }
  mod <- cbind(mod, covar_mod)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n_batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")
    }
    if (ncol(design) > (n_batch + 1)) {
      if ((qr(design[, -c(1:n_batch)])$rank < ncol(design[, 
                                                          -c(1:n_batch)]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      }
      else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")
      }
    }
  }
  NAs = any(is.na(counts))
  if (NAs) {
    cat(c("Found", sum(is.na(counts)), "Missing Data Values\n"), 
        sep = " ")
  }
  cat("Estimating dispersions\n")
  disp_common <- sapply(1:n_batch, function(i) {
    if ((n_batches[i] <= ncol(design) - ncol(batchmod) + 
         1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)) {
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                                   design = NULL, subset = nrow(counts)))
    }
    else {
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                                   design = mod[batches_ind[[i]], ], subset = nrow(counts)))
    }
  })
  genewise_disp_lst <- lapply(1:n_batch, function(j) {
    if ((n_batches[j] <= ncol(design) - ncol(batchmod) + 
         1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)) {
      return(rep(disp_common[j], nrow(counts)))
    }
    else {
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], 
                                    design = mod[batches_ind[[j]], ], dispersion = disp_common[j], 
                                    prior.df = 0))
    }
  })
  names(genewise_disp_lst) <- paste0("batch", levels(batch))
  phi_matrix <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (k in 1:n_batch) {
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], 
                                              n_batches[k])
  }
  cat("Fitting the GLM model\n")
  glm_f <- glmFit(dge_obj, design = design, dispersion = phi_matrix, 
                  prior.count = 1e-04)
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample)
  new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) + 
    vec2mat(alpha_g, ncol(counts))
  glm_f2 <- glmFit.default(dge_obj$counts, design = design, 
                           dispersion = phi_matrix, offset = new_offset, prior.count = 1e-04)
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  

  
  
  if (shrink) {
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii) {
      if (ii == 1) {
        mcres <- mcint_fun(dat = counts[, batches_ind[[ii]]], 
                           mu = mu_hat[, batches_ind[[ii]]], gamma = gamma_hat[, 
                                                                               ii], phi = phi_hat[, ii], gene.subset.n = gene.subset.n)
      }
      else {
        invisible(capture.output(mcres <- mcint_fun(dat = counts[, 
                                                                 batches_ind[[ii]]], mu = mu_hat[, batches_ind[[ii]]], 
                                                    gamma = gamma_hat[, ii], phi = phi_hat[, ii], 
                                                    gene.subset.n = gene.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0("batch", levels(batch))
    gamma_star_mat <- lapply(monte_carlo_res, function(res) {
      res$gamma_star
    })
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res) {
      res$phi_star
    })
    phi_star_mat <- do.call(cbind, phi_star_mat)
    if (!shrink.disp) {
      cat("Apply shrinkage to mean only\n")
      phi_star_mat <- phi_hat
    }
  }
  else {
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }
  # 参考批次处理
  if (!is.null(ref)) {
    gamma_star_mat[, -ref] <- gamma_star_mat[, -ref] - gamma_star_mat[, ref]
    gamma_star_mat[, ref] <- 0
    
    phi_star_mat[, -ref] <- phi_star_mat[, ref]
  }
  
  mu_star <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (jj in 1:n_batch) {
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]]) - 
                                          vec2mat(gamma_star_mat[, jj], n_batches[jj]))
  }
  phi_star <- rowMeans(phi_star_mat)
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (kk in 1:n_batch) {
    counts_sub <- counts[, batches_ind[[kk]]]
    old_mu <- mu_hat[, batches_ind[[kk]]]
    old_phi <- phi_hat[, kk]
    new_mu <- mu_star[, batches_ind[[kk]]]
    new_phi <- phi_star
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub = counts_sub, 
                                                          old_mu = old_mu, old_phi = old_phi, new_mu = new_mu, 
                                                          new_phi = new_phi)
  }
  if (!is.null(ref.batch)) {
    adjust_counts[, batches_ind[[ref]]] <- counts[, batches_ind[[ref]]]
  }
  adjust_counts_whole <- matrix(NA, nrow = nrow(countsOri), 
                                ncol = ncol(countsOri))
  dimnames(adjust_counts_whole) <- dimnames(countsOri)
  adjust_counts_whole[keep, ] <- adjust_counts
  adjust_counts_whole[rm, ] <- countsOri[rm, ]
  return(list(counts = adjust_counts_whole,gamma=gamma_star_mat))
}
ComBatSeq_My <- function (counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE, 
          shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL) 
{
  library(edgeR)
  
  batch <- as.factor(batch)
  if (any(table(batch) <= 1)) {
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }
  keep_lst <- lapply(levels(batch), function(b) {
    which(apply(counts[, batch == b], 1, function(x) {
      !all(x == 0)
    }))
  })
  keep <- Reduce(intersect, keep_lst)
  rm <- setdiff(1:nrow(counts), keep)
  countsOri <- counts
  counts <- counts[keep, ]
  dge_obj <- DGEList(counts = counts)
  n_batch <- nlevels(batch)
  batches_ind <- lapply(1:n_batch, function(i) {
    which(batch == levels(batch)[i])
  })
  n_batches <- sapply(batches_ind, length)
  n_sample <- sum(n_batches)
  cat("Found", n_batch, "batches\n")
  batchmod <- model.matrix(~-1 + batch)
  group <- as.factor(group)
  if (full_mod & nlevels(group) > 1) {
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }
  else {
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data = as.data.frame(t(counts)))
  }
  if (!is.null(covar_mod)) {
    if (is.data.frame(covar_mod)) {
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), 
                                         function(i) {
                                           model.matrix(~covar_mod[, i])
                                         }))
    }
    covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x) {
      all(x == 1)
    })]
  }
  mod <- cbind(mod, covar_mod)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n_batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")
    }
    if (ncol(design) > (n_batch + 1)) {
      if ((qr(design[, -c(1:n_batch)])$rank < ncol(design[, 
                                                          -c(1:n_batch)]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      }
      else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")
      }
    }
  }
  NAs = any(is.na(counts))
  if (NAs) {
    cat(c("Found", sum(is.na(counts)), "Missing Data Values\n"), 
        sep = " ")
  }
  cat("Estimating dispersions\n")
  disp_common <- sapply(1:n_batch, function(i) {
    if ((n_batches[i] <= ncol(design) - ncol(batchmod) + 
         1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)) {
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                                   design = NULL, subset = nrow(counts)))
    }
    else {
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                                   design = mod[batches_ind[[i]], ], subset = nrow(counts)))
    }
  })
  genewise_disp_lst <- lapply(1:n_batch, function(j) {
    if ((n_batches[j] <= ncol(design) - ncol(batchmod) + 
         1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)) {
      return(rep(disp_common[j], nrow(counts)))
    }
    else {
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], 
                                    design = mod[batches_ind[[j]], ], dispersion = disp_common[j], 
                                    prior.df = 0))
    }
  })
  names(genewise_disp_lst) <- paste0("batch", levels(batch))
  phi_matrix <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (k in 1:n_batch) {
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], 
                                              n_batches[k])
  }
  cat("Fitting the GLM model\n")
  glm_f <- glmFit(dge_obj, design = design, dispersion = phi_matrix, 
                  prior.count = 1e-04)
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample)
  new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) + 
    vec2mat(alpha_g, ncol(counts))
  glm_f2 <- glmFit.default(dge_obj$counts, design = design, 
                           dispersion = phi_matrix, offset = new_offset, prior.count = 1e-04)
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  if (shrink) {
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii) {
      if (ii == 1) {
        mcres <- mcint_fun(dat = counts[, batches_ind[[ii]]], 
                           mu = mu_hat[, batches_ind[[ii]]], gamma = gamma_hat[, 
                                                                               ii], phi = phi_hat[, ii], gene.subset.n = gene.subset.n)
      }
      else {
        invisible(capture.output(mcres <- mcint_fun(dat = counts[, 
                                                                 batches_ind[[ii]]], mu = mu_hat[, batches_ind[[ii]]], 
                                                    gamma = gamma_hat[, ii], phi = phi_hat[, ii], 
                                                    gene.subset.n = gene.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0("batch", levels(batch))
    gamma_star_mat <- lapply(monte_carlo_res, function(res) {
      res$gamma_star
    })
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res) {
      res$phi_star
    })
    phi_star_mat <- do.call(cbind, phi_star_mat)
    if (!shrink.disp) {
      cat("Apply shrinkage to mean only\n")
      phi_star_mat <- phi_hat
    }
  }
  else {
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }
  mu_star <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (jj in 1:n_batch) {
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]]) - 
                                          vec2mat(gamma_star_mat[, jj], n_batches[jj]))
  }
  phi_star <- rowMeans(phi_star_mat)
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (kk in 1:n_batch) {
    counts_sub <- counts[, batches_ind[[kk]]]
    old_mu <- mu_hat[, batches_ind[[kk]]]
    old_phi <- phi_hat[, kk]
    new_mu <- mu_star[, batches_ind[[kk]]]
    new_phi <- phi_star
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub = counts_sub, 
                                                          old_mu = old_mu, old_phi = old_phi, new_mu = new_mu, 
                                                          new_phi = new_phi)
  }
  adjust_counts_whole <- matrix(NA, nrow = nrow(countsOri), 
                                ncol = ncol(countsOri))
  dimnames(adjust_counts_whole) <- dimnames(countsOri)
  adjust_counts_whole[keep, ] <- adjust_counts
  adjust_counts_whole[rm, ] <- countsOri[rm, ]
  missing_rows <- setdiff(rownames(adjust_counts_whole), rownames(gamma_star_mat))
  if (length(missing_rows) > 0) {
    zero_rows <- matrix(0, 
                        nrow = length(missing_rows),
                        ncol = ncol(gamma_star_mat),
                        dimnames = list(missing_rows))
    colnames(zero_rows) <- colnames(gamma_star_mat)
    gamma_star_mat <- rbind(gamma_star_mat, zero_rows)
  }
  return(list(counts = adjust_counts_whole,gamma=gamma_star_mat))
}


####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  weights <- pos_res <- list()
  for(i in 1:nrow(dat)){
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    #LH <- sapply(1:G_sub, function(j){sum(log2(dnbinom(x, mu=m[j,], size=1/phi_sub[j])+1))})  
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i]))
    }else{
      pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
    }
    
    weights[[i]] <- as.matrix(LH/sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  weights <- do.call(cbind, weights)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"], weights=weights)	
  return(res)
} 


####  Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(counts_sub[a, b] <= 1){
        new_counts_sub[a,b] <- counts_sub[a, b]
      }else{
        tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
        if(abs(tmp_p-1)<1e-4){
          new_counts_sub[a,b] <- counts_sub[a, b]  
          # for outlier count, if p==1, will return Inf values -> use original count instead
        }else{
          new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}



mapDisp <- function(old_mu, new_mu, old_phi, divider){
  new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
  for(a in 1:nrow(old_mu)){
    for(b in 1:ncol(old_mu)){
      old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
      new_var <- old_var / (divider[a, b]^2)
      new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
    }
  }
  return(new_phi)
}





batchEffectCorrectionUsingCombatSeq <- function(Spatial.Deconv.Data, seed, plotmnn=T) {
  set.seed(seed)
  B1 <- Spatial.Deconv.Data$comData$T_train
  B2 <- Spatial.Deconv.Data$comData$T_test
  gene <- intersect(rownames(B1),rownames(B2))
  B1 <- B1[gene,]
  B2 <- B2[gene,]
  data_combine <- cbind(B1,B2)
  if(plotmnn){
    mnnPlot(data_combine,B1,B2,norm = "None")
  } 
  batch <- c(rep(1, ncol(B1)), rep(2, ncol(B2))) 
  exprCorr <- ComBatSeq_My(data_combine, batch=batch)
  if(plotmnn){
    mnnPlot(exprCorr$counts,B1,B2,norm = "None")
  } 
  
  Spatial.Deconv.Data$BatchCorr.original <- exprCorr$counts
  Spatial.Deconv.Data$BatchCorr <- calculateCPM(exprCorr)[Spatial.Deconv.Data$markerGenes,]
  if(plotmnn){
    mnnPlot(Spatial.Deconv.Data$BatchCorr,B1,B2,norm = "None")
  } 
  Spatial.Deconv.Data$T_combine <- Spatial.Deconv.Data$BatchCorr
  
  return(exprCorr)
  
}

