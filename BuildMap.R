




buildMap <- function(Spatial.Deconv.Data, Methods = c("RCTD","STdeconvolve","spatialDWLS","SPOTlight","Redeconve","TOAST/-P","CARD","Deconf","Linseed","RNASieve","Cell2location","DSA","ssKL","ssFrobenius"),UMI_min_sigma = 300,k.score = 30,useAllGenes= F,k.weight=50) {
 ###run deconvolution methods as benchmarking
  library(peakRAM)
  
  library(scater)
  if(!all((c("spatial.real.cpm.marker","C_ref","markerslist","spatial.real.count.marker","sc","location.real","C") %in% names(Spatial.Deconv.Data)))){
    stop("Spatial.Deconv.Data is not complete!")
  }
  ST <- as.matrix(Spatial.Deconv.Data$spatial.real.cpm.marker)
  cref <- (Spatial.Deconv.Data$C_ref)
  C_ref_count <- Spatial.Deconv.Data$C_ref.count.marker
  c_count <- Spatial.Deconv.Data$C.count.marker
  mgs <- Spatial.Deconv.Data$markerslist
  ST_count <- Spatial.Deconv.Data$spatial.real.count.marker
  index = which(colSums(Spatial.Deconv.Data$sc$data) == 0)
  if(length(index)!=0){
    Spatial.Deconv.Data$sc$data <- Spatial.Deconv.Data$sc$data[,-index]
    Spatial.Deconv.Data$sc$meta <- Spatial.Deconv.Data$sc$meta[-index,]
  }
  if(useAllGenes){
    ref_rna.count.marker <- Spatial.Deconv.Data$sc$data
    # 
    ref_rna <- calculateCPM(Spatial.Deconv.Data$sc$data)###intercept matrix with marker genes after cpm standadization
  }else{
    ref_rna <- calculateCPM(Spatial.Deconv.Data$sc$data)[unique(Spatial.Deconv.Data$markerslist$gene),]###intercept matrix with marker genes after cpm standadization
    ref_rna.count.marker <- Spatial.Deconv.Data$sc$data[unique(Spatial.Deconv.Data$markerslist$gene),] ###count datas,with marker genes
  }
  # 

  if(is.null(Spatial.Deconv.Data$batchCorrMethod)){
    stop("is.null(Spatial.Deconv.Data$batchCorrMethod)")
  }
  if(is.matrix(Spatial.Deconv.Data$BatchCorr) & Spatial.Deconv.Data$batchCorrMethod == "MNN") {
    print("use cpm normalized & log2-transformed sc reference!")
    ref_rna <- log2(calculateCPM(Spatial.Deconv.Data$sc$data) + 1)[unique(Spatial.Deconv.Data$markerslist$gene),]###intercept matrix with marker genes after cpm standadization and log2-transformation
    
  }
  
  ref_rna_meta <- Spatial.Deconv.Data$sc$meta
  location <- Spatial.Deconv.Data$location.real 
  rownames(ref_rna_meta) <- colnames(ref_rna)
  markers <- Spatial.Deconv.Data$markerslist$gene
  
  real_P <- Spatial.Deconv.Data$P_test
  real_C <- Spatial.Deconv.Data$C
  
  
  RESULTS <- Spatial.Deconv.Data$RESULTSMap$RESULTS
  pcc_rmse_p <- Spatial.Deconv.Data$RESULTSMap$pcc_rmse_p
  pcc_rmse_c <- Spatial.Deconv.Data$RESULTSMap$pcc_rmse_c
  time <- Spatial.Deconv.Data$RESULTSMap$time
  mem <- Spatial.Deconv.Data$RESULTSMap$mem
  if(is.null(Methods)) {
    stop("Methods name is NULL,Please set deconvolution method name!")
  }
  if("all" %in% Methods) {
    Methods = c("RCTD","STdeconvolve","spatialDWLS","SPOTlight","Redeconve","TOAST/-P","CARD","Deconf","Linseed","RNASieve","Cell2location","DSA","ssKL","ssFrobenius")
  }
  if("cellMix" %in% Methods) {
    Methods = c("DSA","ssKL","ssFrobenius","Deconf")
  }
  if("R4.4.1" %in% Methods) {
    Methods = c("RCTD","CARD","Tangram","Seurat","TOAST/-P","Linseed")
  }
  tryCatch({
    print(paste("---------------------Run deconvolution methods:",Methods,"------------------------------",sep = " "))
    if("RCTD" %in% Methods) {
       # start_time <- Sys.time()
       # mem_before <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
       RES <- run_rctd_my(ref_rna.count.marker, ref_rna_meta, ST_count, location, mode = "doublet", UMI_min_sigma = UMI_min_sigma)
       RESULTS$RCTD <- RES[[1]]
       # end_time <- Sys.time()
       # mem_after <- 
       if(!is.null(real_P)){
        pcc_rmse_p$RCTD <- getPearsonRMSE(NORM(RESULTS$RCTD),real_P)
       }
       time$RCTD <- RES[[3]]
       mem$RCTD <- RES[[2]]
      print("####----------------------------RCTD--------------------------------------------")
      print(pcc_rmse_p$RCTD)    
    }
   
    if("Seurat" %in% Methods){
      
      RES <- my_Seurat(ref_rna.count.marker,ref_rna_meta,ST_count,k.score = k.score,k.weight = k.weight)
      RESULTS$Seurat <- RES[[1]]
      # end_time <- Sys.time()
      # mem_after <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
      time$Seurat <- RES[[3]]
      mem$Seurat <- RES[[2]]
      if(!is.null(real_P)){
        pcc_rmse_p$Seurat <- getPearsonRMSE(NORM(RESULTS$Seurat),real_P)
      }
      print("####----------------------------Seurat--------------------------------------------")
      print(pcc_rmse_p$Seurat)   
      
    }
    
  
    if("spatialDWLS" %in% Methods) {
      # start_time <- SysC.time()
      # mem_before <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
      RES <- DWLS_my(count = ST_count, ref_rna = ref_rna.count.marker, ref_rna_meta = ref_rna_meta, Sig = NULL, n_cells = 20,python_pth = NULL)
      RESULTS$spatialDWLS <- RES[[1]]
      # end_time <- Sys.time()
      # mem_after <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
      time$spatialDWLS <- RES[[3]]
      mem$spatialDWLS <- RES[[2]]
      colnames(RESULTS$spatialDWLS) <- colnames(ST_count)
      if(!is.null(real_P)){
        pcc_rmse_p$spatialDWLS <- getPearsonRMSE(NORM(RESULTS$spatialDWLS),real_P)
        
      }
      print("####----------------------------spatialDWLS------------------------------------")
      print(pcc_rmse_p$spatialDWLS)    
    }
    
  
    
    
    if("TOAST/-P" %in% Methods) {
      start_time <- Sys.time()
      RES <- toast_my_para(st = ST,ref = cref,marker = mgs, cellTypeName = "CT", geneName = "gene")
      RESULTS[["TOAST/-P"]] <- RES[[1]]
 
      time[["TOAST/-P"]] <- RES[[3]]
      mem[["TOAST/-P"]] <- RES[[2]]
      if(!is.null(real_P)){
        pcc_rmse_p[["TOAST/-P"]] <- getPearsonRMSE(NORM(RESULTS[["TOAST/-P"]]$result_p),real_P)
        
      }
      pcc_rmse_c[["TOAST/-P"]] <- getPearsonRMSE(RESULTS[["TOAST/-P"]]$result_c,real_C)
      print("####----------------------------TOAST/-P-------------------------------------")
      print(pcc_rmse_p[["TOAST/-P"]])
      print(pcc_rmse_c[["TOAST/-P"]])    
    }
    
    if("CARD" %in% Methods) {
      # start_time <- Sys.time()
      # mem_before <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
      RES <- card_my(ST_count,ref_rna.count.marker,ref_rna_meta,location)
      RESULTS$CARD = RES[[1]]
      
      # end_time <- Sys.time()
      # mem_after <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
      time$CARD <- RES[[3]]
      mem$CARD <- RES[[2]]
            if(!is.null(real_P)){
        pcc_rmse_p$CARD <- getPearsonRMSE(NORM(RESULTS$CARD),real_P)
        
      }
      print("####----------------------------CARD--------------------------------------")
      print(pcc_rmse_p$CARD)    
    }
    
    if("Deconf" %in% Methods) {
      # start_time <- Sys.time()
      # mem_before <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
      RES <- deconf_my(ST,mgs)
      RESULTS$Deconf <- RES[[1]]
      # end_time <- Sys.time()
      # mem_after <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
      time$Deconf <- RES[[3]]
      mem$Deconf <- RES[[2]]
      
      pcc_rmse_c$Deconf <- getPearsonRMSE(RESULTS$Deconf$c,real_C)
      if(!is.null(real_P)){
        pcc_rmse_p$Deconf <- getPearsonRMSE(NORM(RESULTS$Deconf$p),real_P)
        
      }
      print("####----------------------------Deconf--------------------------------------")
      print(pcc_rmse_p$Deconf)
      print(pcc_rmse_c$Deconf)    
    }
    
    
    if("Linseed" %in% Methods) {
      # start_time <- Sys.time()
      # mem_before <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
      RES <- linseed_my(ST,cref)
      RESULTS$Linseed <- RES[[1]]
      # end_time <- Sys.time()
      # mem_after <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
      time$Linseed <- RES[[3]]
      mem$Linseed <- RES[[2]]
      
      if(!is.null(real_P)){
        pcc_rmse_p$Linseed <- getPearsonRMSE(NORM(RESULTS$Linseed$p),real_P)
        
      }
      pcc_rmse_c$Linseed <- getPearsonRMSE((RESULTS$Linseed$c),real_C)
      print("####----------------------------Linseed--------------------------------------")
      print(pcc_rmse_p$Linseed)
      print(pcc_rmse_c$Linseed)    
    }
    
    if("RNASieve" %in% Methods) {
      # start_time <- Sys.time()
      # mem_before <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
      
      RES <- RNASievePy(Spatial.Deconv.Data)
      RESULTS$RNASieve <- list(result_p = RES[[1]],result_c = RES[[2]])
      time$RNASieve <- RES[[4]]
      mem$RNASieve <- RES[[3]]
      
      if(!is.null(real_P)){
        pcc_rmse_p$RNASieve <- getPearsonRMSE(NORM(RES[[1]]),real_P)
        
      }
      pcc_rmse_c$RNASieve <- getPearsonRMSE((RES[[2]]),real_C)
      print("####----------------------------RNASieve--------------------------------------")
      print(pcc_rmse_p$RNASieve)
      print(pcc_rmse_c$RNASieve)    
    }
    if("Tangram" %in% Methods){

      RES <- TangramPy(Spatial.Deconv.Data)
      RESULTS$Tangram <- RES[[1]]
      time$Tangram <- RES[[3]]
      mem$Tangram <- RES[[2]]
      
      RESULTS$Tangram <- as.data.frame(t(RESULTS$Tangram))
      RESULTS$Tangram <- as.matrix(RESULTS$Tangram)
      if(!is.null(real_P)){
        pcc_rmse_p$Tangram <- getPearsonRMSE(NORM(RESULTS$Tangram),real_P)
        
      }
      print("####----------------------------Tangram--------------------------------------")
      print(pcc_rmse_p$Tangram)
    }
  
    
    if("DSA" %in% Methods) {
      # start_time <- Scys.time()
      # mem_before <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
      ####"DSA" "ssKL" and "ssFrobenius" need Package CellMix which need R version <= 4.1.0
      RES <- dsa_my(ST,mgs)
      RESULTS$DSA <- RES[[1]]
      # end_time <- Sys.time()
      # mem_after <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
      time$DSA <- RES[[3]]
      mem$DSA <- RES[[2]]
      
      pcc_rmse_c$DSA <- getPearsonRMSE(RESULTS$DSA$c,real_C)
      if(!is.null(real_P)){
        pcc_rmse_p$DSA <- getPearsonRMSE(NORM(RESULTS$DSA$p),real_P)
        
      }
      print("####----------------------------DSA--------------------------------------")
      print(pcc_rmse_p$DSA)
      print(pcc_rmse_c$DSA)
    }
    
    if("ssKL" %in% Methods) {
 
      RES <- sskl_my(ST,mgs)
      RESULTS$ssKL <- RES[[1]]
      time$ssKL <- RES[[3]]
      mem$ssKL <- RES[[2]]
      
      pcc_rmse_c$ssKL <- getPearsonRMSE(RESULTS$ssKL$c,real_C)
      if(!is.null(real_P)){
        pcc_rmse_p$ssKL <- getPearsonRMSE(NORM(RESULTS$ssKL$p),real_P)
        
      }
      print("####----------------------------ssKL--------------------------------------")
      print(pcc_rmse_p$ssKL)
      print(pcc_rmse_c$ssKL)
      
    }
    
    if("ssFrobenius" %in% Methods) {
      # start_ticme <- Sys.time()
      # mem_before <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
      RES <- ssFrobenius_my(ST,mgs)
      RESULTS$ssFrobenius <- RES[[1]]
      time$ssFrobenius <- RES[[3]]
      mem$ssFrobenius <- RES[[2]]
      
      pcc_rmse_c$ssFrobenius <- getPearsonRMSE(RESULTS$ssFrobenius$c,real_C)
      if(!is.null(real_P)){
        pcc_rmse_p$ssFrobenius <- getPearsonRMSE(NORM(RESULTS$ssFrobenius$p),real_P)
        
      }
      print("####----------------------------ssFrobenius--------------------------------------")
      print(pcc_rmse_p$ssFrobenius)
      print(pcc_rmse_c$ssFrobenius)
      
    }
    
    
    if("MyMethod" %in% Methods) {
      
    }
    
    
    
  }, error = function(e){
    print('--------------------------There is something wrong!----------------------')
    message(e$message)
    Spatial.Deconv.Data$RESULTSMap <- list(RESULTS = RESULTS, pcc_rmse_c = pcc_rmse_c, pcc_rmse_p = pcc_rmse_p)
    if(!is.null(Spatial.Deconv.Data$RESULTSMap$pcc_rmse_p) && !is.null(Spatial.Deconv.Data$RESULTSMap$pcc_rmse_c)){
      Spatial.Deconv.Data <- combineResults(Spatial.Deconv.Data)
    }
    
    return(Spatial.Deconv.Data)
  },finally = {
    
    Spatial.Deconv.Data$RESULTSMap <- list(RESULTS = RESULTS, pcc_rmse_c = pcc_rmse_c, pcc_rmse_p = pcc_rmse_p,time=time,mem = mem)
    if(!is.null(Spatial.Deconv.Data$RESULTSMap$pcc_rmse_p) || !is.null(Spatial.Deconv.Data$RESULTSMap$pcc_rmse_c)){
      Spatial.Deconv.Data <- combineResults(Spatial.Deconv.Data)
    }
    print(paste("Tasks of", Methods, "are done!",sep = " "))
    return(Spatial.Deconv.Data)
  })
 
  
  
  

  
}



combineResults <- function(Spatial.Deconv.Data) {
  library(dplyr)
  pcc_rmse_p = Spatial.Deconv.Data$RESULTSMap$pcc_rmse_p
  pcc_rmse_c = Spatial.Deconv.Data$RESULTSMap$pcc_rmse_c
  time = Spatial.Deconv.Data$RESULTSMap$time
  mem = Spatial.Deconv.Data$RESULTSMap$mem
  
  map_p <- NULL
  for (i in pcc_rmse_p) {
    map_p <- rbind(map_p,i)
  }
  rownames(map_p) <- names(pcc_rmse_p)
  if(!is.null(map_p)){
    colnames(map_p) <- c("RMSE of P","Pearson of P")
  }
  
  map_c <- NULL
  
  for (i in pcc_rmse_c) {
  
    map_c <- rbind(map_c,i)
  }

  rownames(map_c) <- names(pcc_rmse_c)
  if(!is.null(map_c)){
    colnames(map_c) <- c("RMSE of C","Pearson of C")
    
  }
  
  map_time <- NULL
  for (i in time) {
    
    map_time <- rbind(map_time,i)
  }
  rownames(map_time) <- names(time)
  colnames(map_time) <- c("Time(s)")

  map_mem <- NULL
  for (i in mem) {
    map_mem <- rbind(map_mem,i)
  }
  rownames(map_mem) <- names(mem)
  colnames(map_mem) <- c("Memory(MB)")
  
  map_p$key <- rownames(map_p)
  map_c$key <- rownames(map_c)
  if(!is.null(map_c) & !is.null(map_p)){
    map <- left_join(map_p, map_c, by = "key")  
    
  }else{
    map <- map_p
  }
  rownames(map) <- map$key
  map <- map[,!(colnames(map) %in% c("key"))]
  
  map_time <- map_time[rownames(map),]
  map_mem <- map_mem[rownames(map),]
  
  Spatial.Deconv.Data$RESULTSMap$MAP = cbind(map,map_time,map_mem)
  

  return(Spatial.Deconv.Data)
}



combineResults.aripurity <- function(Spatial.Deconv.Data) {
  library(dplyr)
  aripurity.kmeans = Spatial.Deconv.Data$RESULTSMap$aripurity.kmeans
  aripurity.dominant = Spatial.Deconv.Data$RESULTSMap$aripurity.dominant
  
  map.kmeans <- NULL
  if(!(is.null(aripurity.kmeans))){
      for (i in aripurity.kmeans) {
    map.kmeans <- rbind(map.kmeans,i)
  }
  rownames(map.kmeans) <- names(aripurity.kmeans)
  colnames(map.kmeans) <- c("ARI.kmeans","Purity.kmeans")
  map.kmeans <- as.data.frame(map.kmeans)
  map.kmeans$key <- rownames(map.kmeans)
  
  
  }

  map.dominant <- NULL
  if(!(is.null(aripurity.dominant))){
    for (i in aripurity.dominant) {
      map.dominant <- rbind(map.dominant,i)
    }
    rownames(map.dominant) <- names(aripurity.dominant)
    colnames(map.dominant) <- c("ARI.dominant","Purity.dominant")
    map.dominant <- as.data.frame(map.dominant)
    map.dominant$key <- rownames(map.dominant)
  }

  
  
  
  if(!(is.null(map.kmeans)) && !(is.null(map.dominant))) {
    map <- left_join(map.kmeans, map.dominant, by = "key")  
  }else if(!(is.null(map.kmeans))){
    map <- map.kmeans
  }else {
    map <- map.dominant
  }
  
  rownames(map) <- map$key
  map <- map[,!(colnames(map) %in% c("key"))]
  if(!is.null(Spatial.Deconv.Data$RESULTSMap$MAP.airpurity) && length(colnames(Spatial.Deconv.Data$RESULTSMap$MAP.airpurity)) == length(colnames(map))){
    Spatial.Deconv.Data$RESULTSMap$MAP.airpurity = rbind(Spatial.Deconv.Data$RESULTSMap$MAP.airpurity,map)
  } else {
    Spatial.Deconv.Data$RESULTSMap$MAP.airpurity = map
  }
  return(Spatial.Deconv.Data)
}

combineMyMethod <- function(Spatial.Deconv.Data, rowNum = -1) {
  ##combine map with parafine result
  if(is.null(Spatial.Deconv.Data$PARAS)){
    stop("Spatial.Deconv.Data$PARAS is NULL,run DSSC_spatial_with_only_lambda first!")
  }
  if(rowNum == -1) {
    sorted_df <- Spatial.Deconv.Data$PARAS[order(Spatial.Deconv.Data$PARAS$rank_average, decreasing = TRUE), ]  
    first_row <- sorted_df[1, ] 
    Spatial.Deconv.Data$RESULTSMap$MyMethod <- first_row
    first_row <- cbind(first_row$RMSE_of_P,first_row$Pearson_of_P, first_row$RMSE_of_C, first_row$Pearson_of_C)
    colnames(first_row) <- c("RMSE of P","Pearson of P","RMSE of C","Pearson of C")
    rownames(first_row) <- "My Method"
    Spatial.Deconv.Data$RESULTSMap$MAPMy <- rbind(first_row, Spatial.Deconv.Data$RESULTSMap$MAP)    
  } else {
    select_row <- Spatial.Deconv.Data$PARAS[rowNum,]
    Spatial.Deconv.Data$RESULTSMap$MyMethod <- select_row
    select_row <- cbind(select_row$RMSE_of_P,select_row$Pearson_of_P, select_row$RMSE_of_C, select_row$Pearson_of_C)
    colnames(select_row) <- c("RMSE of P","Pearson of P","RMSE of C","Pearson of C")
    rownames(select_row) <- "My Method"
    Spatial.Deconv.Data$RESULTSMap$MAPMy <- rbind(select_row, Spatial.Deconv.Data$RESULTSMap$MAP)    
  }

  return(Spatial.Deconv.Data)
}

Cell2locationPy <- function(Spatial.Deconv.Data) {
  ##use python script to run Cell2location deconvolution
  
  library(reticulate)
  ST_count = Spatial.Deconv.Data$spatial.real.count.marker
  ST_count = t(ST_count)
  ST_count <- as.data.frame(ST_count)
  ref_rna.count.marker <- Spatial.Deconv.Data$sc$data[unique(Spatial.Deconv.Data$markerslist$gene),] ###count datas,with marker genes
  ref_rna.count.marker = t(ref_rna.count.marker)
  ref_rna.count.marker <- as.data.frame(ref_rna.count.marker)
  ref_rna_meta <- Spatial.Deconv.Data$sc$meta
  location <- Spatial.Deconv.Data$location.real
  use_condaenv("D:/softwares/anaconda/envs/c2l_new/python.exe")
  source_python("H:/LinYifan/spatialDSSC-BE/R/Cell2locationMy.py")
  result <- py$Cell2locationDeconv(ST_count, location, ref_rna.count.marker, ref_rna_meta)
  return(result)
  
  
}

RNASievePy <- function(Spatial.Deconv.Data) {
  ##use python script to run RNASieve deconvolution
  library(scater)
  ST <- Spatial.Deconv.Data$spatial.real.cpm.marker
  ref_rna <- calculateCPM(Spatial.Deconv.Data$sc$data)[unique(Spatial.Deconv.Data$markerslist$gene),]###intercept matrix with marker genes after cpm standadization
  ref_rna_meta <- Spatial.Deconv.Data$sc$meta
  gene_name <- rownames(ST)
  gene_name <- as.data.frame(gene_name)
  colnames(gene_name) <- c('x')
  colnames(ref_rna) <- ref_rna_meta$cellType
  sc_groupby_class <- lapply(unique(ref_rna_meta$cellType),function(cellType) {
    as.matrix(ref_rna[ ,colnames(ref_rna) %in% cellType])
  })
  names(sc_groupby_class) <- unique(ref_rna_meta$cellType) 
  
  library(reticulate)
  use_condaenv("D:/softwares/anaconda/envs/rnasieve/python.exe")
  source_python("H:/LinYifan/spatialDSSC-BE/R/RNASieveMy.py")
  
  # start_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  # start_time = Sys.time()
  results <- py$RNASieveDeconve(ST, gene_name, sc_groupby_class)
  # end_time = Sys.time()
  # end = as.numeric(bench_process_memory()[1][1]) / (1024)^2 
  results[[1]] <- t(results[[1]])
  results[[1]] <- as.matrix(results[[1]])
  results[[2]] <- as.matrix(results[[2]])
  ###rna-sieve results like : Bulk1 Bulk2....
  colnames(results[[1]]) <- colnames(ST)
  return(results)
  
  
  
  
}



run_rctd_my <- function(scrna_ref,scrna_ref_meta,count_matrix,coords,mode = "multi",require_int = TRUE,
                        gene_cutoff = 0, fc_cutoff = 0, gene_cutoff_reg = 0, fc_cutoff_reg = 0,UMI_min_sigma = 300){
  
  
  set.seed(44)
  library(spacexr)
  library(Matrix)
  library(doParallel)
  
  counts <- scrna_ref
  meta_data <- scrna_ref_meta
  meta_data$cellType <- sub("/","_",meta_data$cellType)
  meta_data$cellType <- sub("/","_",meta_data$cellType)
  cell_counts = table(meta_data$cellType)
  meta_work <- as.data.frame(meta_data[,c("cellID","cellType")])
  colnames(meta_work) <- c("barcode", "clusters")
  meta_data <- meta_work
  cell_types <- meta_data$clusters
  names(cell_types) <- meta_data$barcode # create cell_types named list
  cell_types <- as.factor(cell_types)
  #cell_types <- cell_types[!(libsize <= 0)]
  meta_data$nUMI <- colSums(as.matrix(counts))
  nUMI <- meta_data$nUMI
  names(nUMI) <- meta_data$barcode

  nUMI_sp <- colSums(count_matrix)

  reference <- Reference(as.matrix(counts), cell_types, nUMI,require_int = require_int,min_UMI = 0)
  ### Create SpatialRNA object------------------
  puck <- SpatialRNA(coords, count_matrix, nUMI_sp,require_int = require_int)

  myRCTD <- create.RCTD(puck, reference, max_cores = 8, gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg,UMI_min = 1, UMI_max = 2e+07, counts_MIN = 1,UMI_min_sigma = UMI_min_sigma)
  # start_mem <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
  # start_time <- Sys.time()
  # mem = as.numeric(bench_memory(myRCTD <- run.RCTD(myRCTD, doublet_mode = mode))[1]) / (1024)^2
  peak <- peakRAM(myRCTD <- run.RCTD(myRCTD, doublet_mode = mode))
  
  
  # end_time <- Sys.time()
  # end_mem <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
  results <- myRCTD@results
  norm_weights = normalize_weights(results$weights) 
  result_p <- as.matrix(norm_weights)
  result_p <- t(result_p)
  
  return(list(result_p = result_p,mem = peak$Peak_RAM_Used_MiB * 1.048576,time = peak$Elapsed_Time_sec))
}
#count_matrix = ST_count
#####这个方法只接受int数据
restrictCorpus_My <- function(counts, removeAbove = 1, removeBelow = 0.05, alpha = 0.05, 
                              nTopOD = 1000, plot = FALSE, verbose = TRUE) {

  
    vi <- rowSums(as.matrix(counts) > 0) >= ncol(counts) * removeAbove
    if (verbose) {
      message(paste0("Removing ", sum(vi), " genes present in ", 
                     removeAbove * 100, "% or more of pixels..."))
    }
    counts <- counts[!vi, ]
    if (verbose) {
      message(paste0(nrow(counts), " genes remaining..."))
    }
    vi <- rowSums(as.matrix(counts) > 0) <= ncol(counts) * removeBelow
    if (verbose) {
      message(paste0("Removing ", sum(vi), " genes present in ", 
                     removeBelow * 100, "% or less of pixels..."))
    }
    counts <- counts[!vi, ]
    if (verbose) {
      message(paste0(nrow(counts), " genes remaining..."))
    }
    # if (verbose) {
    #   message(paste0("Restricting to overdispersed genes with alpha = ", 
    #                  alpha, "..."))
    # }
    # OD <- getOverdispersedGenes(counts, alpha = alpha, plot = plot, 
    #                             details = TRUE, verbose = verbose)
    # if (!is.na(nTopOD)) {
    #   if (verbose) {
    #     message(" Using top ", nTopOD, " overdispersed genes.", 
    #             "\n")
    #   }
    #   OD_filt <- OD$df[OD$ods, ]
    #   if (dim(OD_filt)[1] < nTopOD) {
    #     if (verbose) {
    #       message(" number of top overdispersed genes available: ", 
    #               dim(OD_filt)[1], "\n")
    #     }
    #     od_genes <- rownames(OD_filt)
    #   }
    #   else {
    #     od_genes <- rownames(OD_filt[order(OD_filt$lpa), 
    #     ][1:nTopOD, ])
    #   }
    # }
    # else {
    #   od_genes <- OD$ods
    # }
    countsFiltRestricted <- counts
    if (plot) {
      par(mfrow = c(1, 2), mar = rep(5, 4))
      hist(log10(Matrix::colSums(countsFiltRestricted) + 1), 
           breaks = 20, main = "Genes Per Pixel")
      hist(log10(Matrix::rowSums(countsFiltRestricted) + 1), 
           breaks = 20, main = "Pixels Per Gene")
    }
    if (dim(countsFiltRestricted)[1] > 1000) {
      message("Genes in corpus > 1000 (", dim(countsFiltRestricted)[1], 
              "). This may cause model fitting to take a while. Consider reducing the number of genes.", 
              "\n")
    }
    return(countsFiltRestricted)
  
  
}


STd_my <- function(count_matrix, cref, betaScale = 1){
  set.seed(44)
  library(STdeconvolve)
  K = length(colnames(cref))
  ####count数据，而且使用marker基因
  #genes <- unique(inter_hard_3_new$markers_$genes)
  cd <- count_matrix
  ## remove pixels with too few genes
  counts <- cleanCounts(cd, min.lib.size = 1)
  ## feature select for genes
  corpus <- restrictCorpus_My(counts, removeAbove=10, removeBelow = 0)
  #counts <- cleanCounts(as.matrix(cd), min.lib.size = 1)
  ldas <- fitLDA(t(as.matrix(cd)), Ks = K)#seq(2, 9, by = 1)
  optLDA <- optimalModel(models = ldas, opt = "min")
  results <- getBetaTheta(optLDA,perc.filt = 0.05,betaScale = betaScale)
  deconProp <- t(results$theta)
  deconGexp <- t(results$beta)
  result_p <- deconProp
  result_c <- deconGexp
  # getwd()
  # setwd("H:/LinYifan/")
  # source("CDSC.R")
  # source("function_help.R")
  # # cref <- cref
  # result_c <- result_c[intersect(rownames(cref),rownames(result_c)),]
  # cref <- cref[intersect(rownames(cref),rownames(result_c)),]
  # ct <- Row_label(result_c,cref)
  # colnames(result_c) <- ct
  # rownames(result_p) <- ct
  # RESULTS = NULL
  # result_C_PCCRMSE = getPearsonRMSE(result_c, cref)
  return(list(result_c = result_c,result_p = result_p))
}




DWLS_my <- function(count, ref_rna, ref_rna_meta, n_cells = 5, Sig = NULL,python_pth = NULL){
  set.seed(44)
  library(Giotto)
  if(is.null(python_pth)){
    my_python_path= "C:/Users/Dell/AppData/Local/r-miniconda/envs/giotto_env/python.exe"
    # my_python_path= "D:/softwares/anaconda/envs/gst/python.exe"
    
  }else{
    my_python_path= python_pth
    
  }
  instrs = createGiottoInstructions(python_path = my_python_path)
  st_data <- createGiottoObject(
    raw_exprs = count,
    instructions = instrs
  )
  
  st_data <- normalizeGiotto(gobject = st_data)
  st_data <- calculateHVG(gobject = st_data)
  gene_metadata = fDataDT(st_data)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  # if (length(featgenes) <= 15) {
  #   featgenes = gene_metadata$gene_ID
  # }
  featgenes = gene_metadata$gene_ID
  st_data <- runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(st_data, genes_to_use = featgenes, scale_unit = F)
  st_data <- runUMAP(st_data, dimensions_to_use = 1:5, n_neighbors = 5)
  st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
  st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)
  sc_data <- createGiottoObject(
    raw_exprs = ref_rna,
    instructions = instrs
  )
  sc_data <- normalizeGiotto(gobject = sc_data)
  sc_data <- calculateHVG(gobject = sc_data)
  gene_metadata = fDataDT(sc_data)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  # if (length(featgenes) <= 15) {
  #   featgenes = gene_metadata$gene_ID
  # }
  featgenes = gene_metadata$gene_ID
  sc_data <- runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
  sc_data@cell_metadata$leiden_clus <- as.character(ref_rna_meta[,c("cellType")])
  # scran_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
  #                                                    method = 'scran',
  #                                                    expression_values = 'normalized',
  #                                                    cluster_column = 'leiden_clus')
  # Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 100)])
  norm_exp<-2^(sc_data@norm_expr)-1
  id<-sc_data@cell_metadata$leiden_clus
  ExprSubset<-norm_exp
  Sig_exp<-NULL
  for (i in unique(id)){
    Sig_exp<-cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig_exp)<-unique(id)

  # mem = as.numeric(bench_memory(st_data <- runDWLSDeconv(st_data,sign_matrix = Sig_exp, n_cell = n_cells))[1]) / (1024)^2
  peak <- peakRAM(st_data <- runDWLSDeconv(st_data,sign_matrix = Sig_exp, n_cell = n_cells))
  # end_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  spatial_DWLS_result <- st_data@spatial_enrichment$DWLS
  cell_id <- spatial_DWLS_result$cell_ID
  spatial_DWLS_result <- as.data.frame(spatial_DWLS_result)
  
  spatial_DWLS_result <- t(spatial_DWLS_result)
  colnames(spatial_DWLS_result) <- cell_id
  ct <- colnames(Sig_exp)
  DWLS <- as.matrix(spatial_DWLS_result[ct,])
  num_matrix <- as.numeric(DWLS)
  num_matrix <- matrix(num_matrix, nrow = nrow(DWLS))
  colnames(num_matrix) <- colnames(DWLS)
  # num_matrix <- NORM(num_matrix)
  rownames(num_matrix) <- ct
  return(list(result = num_matrix,mem = peak$Peak_RAM_Used_MiB * 1.048576  , time = peak$Elapsed_Time_sec))
}





DWLS_my_test <- function(count, ref_rna, ref_rna_meta, n_cells = 5, Sig = NULL,python_pth = NULL){
  set.seed(44)
  library(Giotto)
  if(is.null(python_pth)){
    my_python_path= "D:/softwares/anaconda/envs/rnasieve/python.exe"
    
  }else{
    my_python_path= python_pth
    
  }
  instrs = createGiottoInstructions(python_path = my_python_path)
  st_data <- createGiottoObject(
    raw_exprs = count,
    instructions = instrs
  )
  st_data <- normalizeGiotto(gobject = st_data,library_size_norm = FALSE)
  st_data <- calculateHVG(gobject = st_data)
  gene_metadata = fDataDT(st_data)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  st_data <- runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(st_data, genes_to_use = featgenes, scale_unit = F)
  st_data <- runUMAP(st_data, dimensions_to_use = 1:10)
  st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
  st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)
  sc_data <- createGiottoObject(
    raw_exprs = ref_rna,
    norm_expr = ref_rna,
    instructions = instrs
  )
  # sc_data <- normalizeGiotto(gobject = sc_data)
  sc_data <- calculateHVG(gobject = sc_data)
  gene_metadata = fDataDT(sc_data)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  sc_data <- runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
  sc_data@cell_metadata$leiden_clus <- as.character(ref_rna_meta[,"cellType"])
  # scran_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
  #                                                    method = 'scran',
  #                                                    expression_values = 'normalized',
  #                                                    cluster_column = 'leiden_clus')
  # Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 100)])
  norm_exp<-2^(sc_data@norm_expr)-1
  id<-sc_data@cell_metadata$leiden_clus
  ExprSubset<-norm_exp
  Sig_exp<-NULL
  for (i in unique(id)){
    Sig_exp<-cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig_exp)<-unique(id)
  # Sig_exp = Sig
  st_data <- runDWLSDeconv(st_data,sign_matrix = Sig_exp, n_cell = n_cells)
  
  spatial_DWLS_result <- st_data@spatial_enrichment$DWLS
  cell_id <- spatial_DWLS_result$cell_ID
  spatial_DWLS_result <- as.data.frame(spatial_DWLS_result)
  
  spatial_DWLS_result <- t(spatial_DWLS_result)
  colnames(spatial_DWLS_result) <- cell_id
  ct <- colnames(Sig)
  DWLS <- as.matrix(spatial_DWLS_result[ct,])
  num_matrix <- as.numeric(DWLS)
  num_matrix <- matrix(num_matrix, nrow = nrow(DWLS))
  colnames(num_matrix) <- colnames(DWLS)
  # num_matrix <- NORM(num_matrix)
  rownames(num_matrix) <- ct
  return(num_matrix)
}

my_music <- function(sccount,scmeta,count){
  scData = NULL
  library(SingleCellExperiment)
  scData$sc$C <- as.matrix(sccount)
  ncells <- length(colnames(scData$sc$C))
  if(!("sampleInfo" %in% colnames(scmeta))){
    scmeta$sampleInfo <- scmeta$cellID
  }
  
  scData$sc$phenoDataC <- as.data.frame(scmeta)
  
  #name <- sub("^X(\\d+)_(.\\d+)$", "\\1", colnames(scrnaseq))
  #scData$sc$phenoDataC$SubjectName <- name
  
  rownames(scData$sc$phenoDataC) = scData$sc$phenoDataC$cellID
  u <- matrix(as.matrix(scData$sc$C),ncol = ncells)
  rownames(u) <- rownames(scData$sc$C)
  scData$sc$C.eset <- SingleCellExperiment(assays = list(counts = u),colData = scData$sc$phenoDataC)
  #---------MuSiC-----------
  library(MuSiC)
  select.ct = unique(scData$sc$phenoDataC$cellType)
  #select.ct = intersect(unique(as.character(scData$sc$phenoDataC$cellType)),unique(ComDate1$scData_2$full_phenoData$cellType))
  result_music = t(music_prop(bulk.mtx = count, 
                              sc.sce = scData$sc$C.eset,
                              clusters = "cellType",
                              samples = 'cellID',
                              select.ct = select.ct,
                              
  )$Est.prop.weighted)
  #result_music <- result_music[intersect(rownames(result_music),rownames(inter_hard_3_new$indata$P_test)),]
  #result$MuSiC$result = getPearsonRMSE(NORM(result_music),inter_hard_3_new$indata$P_test)
  #result$MuSiC$result
  return(NORM(result_music))
}




spotlight_my <- function(sccount,scemeta,count){
  set.seed(44)
  library(SPOTlight)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(scater)
  library(scran)
  meta_data <- scemeta
  meta_data$free_annotation = meta_data$cellType
  sce <- SingleCellExperiment(assays = list(counts = sccount,logcounts = sccount), colData = meta_data)
  # sce <- sce #reference sc
  # names(assays(sce))
  # # Get vector indicating which genes are neither ribosomal or mitochondrial
  genes <- rownames(sce)
  dec <- modelGeneVar(sce, subset.row = genes)
  # Get the top 3000 genes.
  hvg <- getTopHVGs(dec, n = 3000)
  Brain_ST <-  as.matrix(count)
  colLabels(sce) <- colData(sce)$free_annotation
  # Compute marker genes
  mgs <- scoreMarkers(sce, subset.row = genes)
  mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)
  # split cell indices by identity
  idx <- split(seq(ncol(sce)), sce$cellType)
  ### Deconvolution
  res <- SPOTlight(
    x = sce, 
    y = Brain_ST,
    groups = as.character(meta_data$cellType), # 也可以是cluster，
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")
  #Extract data from `SPOTlight`:
  str(res) #查看结果类型 
  decon_mtrx <- res$mat
  result_p <- decon_mtrx
  result_p <- t(result_p)
  return(result_p)
}




rede_my <- function(st,sc,annotations){
  set.seed(44)
  library(Redeconve)
  st <- st
  #ref <- ref
  res <- deconvoluting(sc, st, genemode = "def", hpmode = "def", dopar = T, ncores = 4)
  # ref <- get.ref(sc, annotations, dopar = T)
  # result <- deconvoluting(
  #   ref = ref,
  #   st = st,
  #   genemode = "default",
  #   hpmode = "autoselection",
  #   normalize = T,
  #   thre = 1e-10,
  #   dopar = T,
  #   ncores = 8,
  #   realtime = F
  # )
  ####result is spots with different cells
  annotations <- droplevels(annotations)
  res.ctmerge <- sc2type(res, annotations)
  res.prop <- to.proportion(res.ctmerge)
  return(as.matrix(res.prop))##res.pop
  
}

toast_my <- function(st,marker,ref,cellTypeName = "cluster", geneName = "genes"){
  
  ###TOAST/P-
  ####PRF DECONVOLUTION METHOD
  # 
  library(TOAST)
  #   
  K=ncol(ref)
  set.seed(1234)
  nmarker = nrow(st)
  start_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  start_time = Sys.time()
  result <- csDeconv(as.matrix(st), K = K, TotalIter = 30, bound_negative = TRUE,nMarker = nmarker) 
  end_time = Sys.time()
  end_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  
  
  rownames(result$W) <- unlist(selmarker)
  result$W <- result$W[unique(unlist(selmarker)),]
  ct <- Row_label(result$W,ref)
  colnames(result$W) <- ct
  rownames(result$H) <- ct
  ct <- intersect(as.character(ct),as.character(colnames(ref)))
  gene <- intersect(rownames(ref),rownames(result$W))
  
  result_c <- (result$W)[gene,ct]
  real_c <- ref[gene,ct]
  #result_pcc_rmse = getPearsonRMSE(result_c,real_c)
  
  colnames(result$H) <- colnames(st)
  result_p <- result$H[ct,]
  return(list(result = list(result_c = result_c,result_p = result_p),mem = (end_mem - start_mem) ,time = end_time - start_time))
}



toast_my_para <- function(st,marker,ref,cellTypeName = "cluster", geneName = "genes"){
  
  ###TOAST/P-
  ####PRF DECONVOLUTION METHOD
  # 
  Ymat <- st
  library(TOAST)
  #   
  
  # 
  SelMarker <- marker
  #   # selmarker <- NULL
  cell_name <- unique(SelMarker[,cellTypeName])
  # 
  selmarker <- split(SelMarker[,geneName],SelMarker[,cellTypeName])
  selmarker <- selmarker[sapply(selmarker, length) > 1]
  # start_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  start_time = Sys.time()
  # mem = as.numeric(bench_memory( result <- MDeconv(Ymat, selmarker,
  #                                                  alpha = NULL, sigma = NULL,
  #                                                  epsilon = 0.001,
  #                                                  verbose = TRUE))[1]) / (1024)^2
  
  peak <- peakRAM(result <- MDeconv(Ymat, selmarker,
                                    alpha = NULL, sigma = NULL,
                                    epsilon = 0.001,
                                    verbose = TRUE))
  ###TOAST/-P
  end_time = Sys.time()
  # end_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  
  rownames(result$W) <- unlist(selmarker)
  result$W <- result$W[unique(unlist(selmarker)),]
  colnames(ref) <- gsub("[^a-zA-Z0-9]",".",colnames(ref))
  names(selmarker) <- gsub("[^a-zA-Z0-9]",".",names(selmarker))
  ct <- Row_label(result$W,ref[rownames(result$W),names(selmarker)])
  colnames(result$W) <- ct
  rownames(result$H) <- ct
  ct <- intersect(as.character(ct),as.character(colnames(ref)))
  gene <- intersect(rownames(ref),rownames(result$W))
  
  result_c <- (result$W)[gene,ct]
  real_c <- ref[gene,ct]
  #result_pcc_rmse = getPearsonRMSE(result_c,real_c)
  
  colnames(result$H) <- colnames(st)
  result_p <- result$H[ct,]
  return(list(results = list(result_c = result_c,result_p = result_p),mem = peak$Peak_RAM_Used_MiB * 
                                                                               1.048576, time = peak$Elapsed_Time_sec))
}




dssc_my <- function(ST,cref,P_real,C_real,T_train,P_train,useARI = TRUE,
                    lambda1 = c(10^-3,10^-1,10^-2),
                    lambda2 = c(10^-3,10^-1,10^-2),
                    lambdaC = c(10^3,3000,10000),
                    m_values = c(0),
                    num_iter_max = 3000,
                    l = c(1000),
                    seedd = 44,
                    TerCondition = 10^-8
                    
){
  

  source("CDSC.R")
  source("function_help.R")
  
  
  #----------find Paramater ///----------
  result <- list()
  lambda1 <- c(10^-3,10^-1,10^-2)
  lambda2 <- c(10^-3,10^-1,10^-2)
  lambdaC <- c(10^3,3000,10000)
  
  m_values <- c(0)
  num_iter_max <- 3000
  l <-c(0)
  result$CDSC4$seedd = seedd
  
  intra_sc_Data <- NULL
  intra_sc_Data$Indata$TerCondition = TerCondition
  
  
  
  ###使用参考数据集产生marker基因和C_ref矩阵
  
  
  new_bulk_T <- ST
  #genes <- unique(data$markers$genes)
  
  
  
  #ss_cal_output(spatial.pseudo$pseudo.data,spatial.pseudo2$pseudo.data, usefov = FALSE)
  
  
  
  #------------cyclic to do deconvolution----------
  pb <- txtProgressBar(style = 3)
  star_time <- Sys.time()
  result$CDSC4$Ss <- SM(t(new_bulk_T))
  result$CDSC4$Sg <- SM(new_bulk_T)
  result_para_c = NULL
  result_para_p = NULL
  pearson_para_c = NULL
  number_iter = NULL
  num = 1
  g_values <- c(0)
  g_i = 1
  para_find = NULL 
  library(dplyr)
  
  for (l_i in 1:length(l)){
    #LAMDA <- caclulate_lamda(ComDate1$indata$T_test, ComDate1$markers, as.matrix(ComDate1$indata$C_ref), l = l[l_i])
    for (dir_i in 1:length(lambda1)){
      for (dir_j in 1:length(lambda2)){
        for (dir_k in 1:length(lambdaC)){
          result_CDSC = CDSC(new_bulk_T, as.matrix(cref), dim(cref)[2], 
                             lambda1 = lambda1[dir_i], lambda2 = lambda2[dir_j], lambdaC = lambdaC[dir_k],
                             10^-8,result$CDSC4$seedd,
                             result$CDSC4$Ss,result$CDSC4$Sg,all_number = num_iter_max)
          # result_CDSC = CDSC(ComDate1$indata$T_test, as.matrix(ComDate1$indata$C_ref), dim(ComDate1$indata$C_ref)[2], 
          #                    lambda1 = lambda1, lambda2 = lambda2, lambdaC = lambdaC,
          #                    10^-8,result$CDSC4$seedd,
          #                    result$CDSC4$Ss,SM(ComDate1$indata$T_test),all_number = num_iter_max)
          
          #result_CDSC_original <- result_CDSC
          for(m_i in 1:length(m_values)){
            result_para_c <- result_CDSC[[1]]
            
            
            result1 = NULL
            if(!all(is.na(result_para_c) == FALSE) ){
              break
            }
            
            
            result_para_p_all <- result_CDSC[[2]]
            #result_para_p_all <- filter_p(result_para_p_all, LAMDA, m = m_values[m_i])
            
            
            #real_P_mixname <- (ncol(intra_sc_Data$Indata$P_train) + 1):(ncol(intra_sc_Data$Indata$P_train) + ncol(intra_sc_Data$Indata$P_test))
            #pre_P_mixname <- 1:ncol(intra_sc_Data$Indata$P_train)
            result_para_p_real <- result_para_p_all
            
            
            
            result_para_p_pre <- NULL
            number_iter <- result_CDSC[[3]]
            # result1 <- calculate_result_1(result_c = result_para_c,C_ref = data$indata$C_ref[genes,],
            #                                                      lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
            #                                                      number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = 10^-8,
            #                                                      leastnum = 3)
            # 
            result1 <- calculate_result_1(result_c = result_para_c,result_p = result_para_p_real,
                                          T = ST,C = C_real,C_ref = cref,P = P_real, 
                                          lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
                                          number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = intra_sc_Data$Indata$TerCondition,
                                          leastnum = 3,ctlabels=NULL,result_p_1 = result_para_p_pre, P_1 = P_train, T_1 = T_train, result_p_combine = result_para_p_all,T_combine = new_bulk_T, l = l[l_i], m = m_values[m_i], g = g_values[g_i])
            if(useARI == TRUE){
              mymethod <- draw_kmeans_result(result_para_p_real)
              result1 <- cbind(result1,ari(mymethod,ground_truth = read.table("D:/locations_real.txt")),purityyy(mymethod,ground_truth = read.table("D:/locations_real.txt")))
              para_find =  rbind(para_find,result1)
            }
            para_find =  rbind(para_find,result1)
            
            num = num + 1
            setTxtProgressBar(pb, num/(length(lambda1)*length(lambda2)*length(lambdaC)*length(l)*length(m_values)))
          }
        }
      }
    }
  }
  end_time <- Sys.time()
  close(pb)
  print(star_time)
  print(end_time)
  end_time-star_time
  return((para = para_find))
}



dssc_my_dir <- function(ST,cref,P_real,C_real,T_train,P_train,useARI = TRUE,l1 = 0.001,l2 = 0.001,lc = 100){
  
  getwd()
  setwd("D:/大四上/毕设/CDSC-master/code/CDSC-master/CDSC-master/code")
  
  source("CDSC.R")
  source("function_help.R")
  
  
  
  #----------find Paramater ///----------
  result <- list()
  
  
  m_values <- c(0)
  num_iter_max <- 3000
  l <-c(0)
  result$CDSC4$seedd = 44
  intra_sc_Data <- NULL
  intra_sc_Data$Indata$TerCondition = 10^-8
  
  
  
  ###使用参考数据集产生marker基因和C_ref矩阵
  
  
  new_bulk_T <- ST
  
  #------------cyclic to do deconvolution----------
  pb <- txtProgressBar(style = 3)
  star_time <- Sys.time()
  result$CDSC4$Ss <- SM(t(new_bulk_T))
  result$CDSC4$Sg <- SM(new_bulk_T)
  result_para_c = NULL
  result_para_p = NULL
  pearson_para_c = NULL
  number_iter = NULL
  num = 1
  g_values <- c(0)
  g_i = 1
  para_find = NULL 
  library(dplyr)
  
  
  result_CDSC = CDSC(new_bulk_T, as.matrix(cref), dim(cref)[2], 
                     lambda1 = l1, lambda2 = l2, lambdaC = lc,
                     10^-8,result$CDSC4$seedd,
                     result$CDSC4$Ss,result$CDSC4$Sg,all_number = num_iter_max)
  # result_CDSC = CDSC(ComDate1$indata$T_test, as.matrix(ComDate1$indata$C_ref), dim(ComDate1$indata$C_ref)[2], 
  #                    lambda1 = lambda1, lambda2 = lambda2, lambdaC = lambdaC,
  #                    10^-8,result$CDSC4$seedd,
  #                    result$CDSC4$Ss,SM(ComDate1$indata$T_test),all_number = num_iter_max)
  
  #result_CDSC_original <- result_CDSC
  for(m_i in 1:length(m_values)){
    result_para_c <- result_CDSC[[1]]
    
    
    result1 = NULL
    if(!all(is.na(result_para_c) == FALSE) ){
      break
    }
    
    
    result_para_p_all <- result_CDSC[[2]]
    #result_para_p_all <- filter_p(result_para_p_all, LAMDA, m = m_values[m_i])
    
    
    #real_P_mixname <- (ncol(intra_sc_Data$Indata$P_train) + 1):(ncol(intra_sc_Data$Indata$P_train) + ncol(intra_sc_Data$Indata$P_test))
    #pre_P_mixname <- 1:ncol(intra_sc_Data$Indata$P_train)
    result_para_p_real <- result_para_p_all
    
    
    
    result_para_p_pre <- NULL
    number_iter <- result_CDSC[[3]]
    result1 <- calculate_result_1(result_c = result_para_c,result_p = result_para_p_real,
                                  T = ST,C = C_real,C_ref = cref,P = P_real, 
                                  lambda1=l1, lambda2=l2, lambdaC=lc,
                                  number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = intra_sc_Data$Indata$TerCondition,
                                  leastnum = 3,ctlabels=NULL,result_p_1 = result_para_p_pre, P_1 = P_train, T_1 = T_train, result_p_combine = result_para_p_all,T_combine = new_bulk_T, l = l[1], m = m_values[1], g = g_values[1])
    
    
    
    num = num + 1
  }
  
  
  return(list(result_p = result_para_p_real,result_c_pcc_rmse = result1,result_c = result_para_c))
}



my_method <- function(ST,cref,P_real,C_real,T_train,P_train,useARI = TRUE,Ss,mgs,gsva_mat){
  
  
  getwd()
  setwd("D:/大四上/毕设/CDSC-master/code/CDSC-master/CDSC-master/code")
  
  source("CDSC.R")
  source("function_help.R")
  
  
  #----------find Paramater ///----------
  result <- list()
  lambda1 <- c(1)
  lambda2 <- c(0.1,0.01,0.001)
  lambdaC <- c(0.1,0.01,0.001)
  
  
  
  
  m_values <- seq(0,0.1,0.02)
  num_iter_max <- 3000
  l <-c(0)
  g_values <- c(0)
  result$CDSC4$seedd = 44
  intra_sc_Data <- NULL
  intra_sc_Data$Indata$TerCondition = 10^-8
  
  
  
  ###使用参考数据集产生marker基因和C_ref矩阵
  
  #genes <- unique(inter_hard_3_new$markers_$genes)
  new_bulk_T <- cbind(T_train,ST)
  
  
  
  
  #ss_cal_output(spatial.pseudo$pseudo.data,spatial.pseudo2$pseudo.data, usefov = FALSE)
  
  
  
  #------------cyclic to do deconvolution----------
  pb <- txtProgressBar(style = 3)
  star_time <- Sys.time()
  result$CDSC4$Ss <- Ss
  result$CDSC4$Sg <- SM(new_bulk_T)
  result_para_c = NULL
  result_para_p = NULL
  pearson_para_c = NULL
  number_iter = NULL
  num = 1
  library(scater)
  
  para_find = NULL 
  library(dplyr)
  
  #gsva_result <- gsva_apply(new_bulk_T,inter_hard_3_new$markers_,"Poisson")
  for (l_i in 1:length(l)){
    for (g_i in 1:length(g_values)){
      LAMBDA <- caclulate_lamda(new_bulk_T, mgs, as.matrix(cref), l = l[l_i], g = g_values[g_i],gsva_result = gsva_mat)
      for (dir_i in 1:length(lambda1)){
        for (dir_j in 1:length(lambda2)){
          for (dir_k in 1:length(lambdaC)){
            result_CDSC = CDSC_1(new_bulk_T, as.matrix(cref), dim(cref)[2], 
                                 lambda1 = lambda1[dir_i], lambda2 = lambda2[dir_j], lambdaC = lambdaC[dir_k],LAMBDA = LAMBDA,
                                 10^-8,result$CDSC4$seedd,
                                 result$CDSC4$Ss,result$CDSC4$Sg,all_number = num_iter_max)
            # result_CDSC = CDSC(ComDate1$indata$T_test, as.matrix(ComDate1$indata$C_ref), dim(ComDate1$indata$C_ref)[2], 
            #                    lambda1 = lambda1, lambda2 = lambda2, lambdaC = lambdaC,
            #                    10^-8,result$CDSC4$seedd,
            #                    result$CDSC4$Ss,SM(ComDate1$indata$T_test),all_number = num_iter_max)
            
            #result_CDSC_original <- result_CDSC
            for(m_i in 1:length(m_values)){
              result_para_c <- result_CDSC[[1]]
              
              
              result1 = NULL
              if(!all(is.na(result_para_c) == FALSE) ){
                break
              }
              
              
              result_para_p_all <- result_CDSC[[2]]
              result_para_p_all <- filter_p(result_para_p_all, LAMBDA, m = m_values[m_i])
              
              
              real_P_mixname <- (ncol(P_train) + 1):ncol(result_para_p_all)
              pre_P_mixname <- 1:ncol(P_train)
              result_para_p_real <- result_para_p_all[,real_P_mixname]
              
              
              
              result_para_p_pre <- result_para_p_all[,pre_P_mixname]
              number_iter <- result_CDSC[[3]]
              # result1 <- calculate_result_1(result_c = result_para_c,C_ref = data$indata$C_ref[genes,],
              #                                                      lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
              #                                                      number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = 10^-8,
              #                                                      leastnum = 3)
              # 
              
              result1 <- calculate_result_1(result_c = result_para_c,result_p = result_para_p_real,
                                            T = ST,C = C_real,C_ref = cref,P = P_real,
                                            lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
                                            number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = intra_sc_Data$Indata$TerCondition,
                                            leastnum = 3,ctlabels=NULL,result_p_1 = result_para_p_pre, P_1 = P_train, T_1 = T_train, result_p_combine = result_para_p_all,T_combine = new_bulk_T, l = l[l_i], m = m_values[m_i], g = g_values[g_i])
              
              if(useARI == TRUE){
                mymethod <- draw_kmeans_result(result_para_p_real)
                result1 <- cbind(result1,ari(mymethod,ground_truth = read.table("D:/locations_real.txt")),purityyy(mymethod,ground_truth = read.table("D:/locations_real.txt")))
                para_find =  rbind(para_find,result1)
              }
              
              
              
              para_find =  rbind(para_find,result1)
              
              
              num = num + 1
              setTxtProgressBar(pb, num/(length(lambda1)*length(lambda2)*length(lambdaC)*length(l)*length(m_values)*length(g_values)))
            }
          }
        }
      }
    }
  }
  end_time <- Sys.time()
  close(pb)
  print(star_time)
  print(end_time)
  end_time-star_time
  
  
  return(para_find)
}


getARIPurity <- function(clu,ground_truth,anno = c("Layer")) {
  return(list(ARI = ari(clu,ground_truth,anno = c("Layer")),Purity = purityyy(clu,ground_truth,anno = c("Layer"))))
}

getARIPurity.kmeans <- function(result_p,location,anno = c("Layer"), methodTitle,ground_truth, draw.plot = TRUE, seed = 44) {
  clu <- draw.pie.cluster(result_p,location,methodTitle = methodTitle, draw.plot = draw.plot, seed = seed)
  
  return(list(ARI.kmeans = ari(clu,ground_truth,anno = c("Layer")),Purity.kmeans = purityyy(clu,ground_truth,anno = c("Layer"))))
}

getARIPurity.dominant <- function(result_p,location,anno = c("Layer"), methodTitle,ground_truth, draw.plot = TRUE, seed = 44) {
  clu <- draw.pie.main(result_p,location,methodTitle = methodTitle, draw.plot = draw.plot , seed = seed)
  return(list(ARI.dominant = ari(clu,ground_truth,anno = c("Layer")),Purity.dominant = purityyy(clu,ground_truth,anno = c("Layer"))))
}

getARIPurity.kmeans.oap <- function(result_p,location,anno = c("Layer"), methodTitle,ground_truth, draw.plot = TRUE, seed = 44, fileName = "H:/LinYifan/OAP/new.pptx", colors = NULL, font_size = 15) {
  library(eoffice)
  result <- draw.pie.cluster.oap(result_p,location,methodTitle = methodTitle, draw.plot = draw.plot, seed = seed, colors = colors,font_size = font_size)
  topptx(result$plot, filename = fileName)
  return(list(ARI = ari(result$clu,ground_truth,anno = c("Layer")),Purity = purityyy(result$clu,ground_truth,anno = c("Layer"))))
}

getARIPurity.dominant.oap <- function(result_p,location,anno = c("Layer"), methodTitle,ground_truth, draw.plot = TRUE, seed = 44, fileName = "H:/LinYifan/OAP/new.pptx", colors = NULL, font_size = 15) {
  library(eoffice)
  result <- draw.pie.main.oap(result_p,location,methodTitle = methodTitle, draw.plot = draw.plot , seed = seed, colors = colors, font_size = font_size)
  topptx(result$plot, filename = fileName)
  return(list(ARI = ari(result$clu,ground_truth,anno = c("Layer")),Purity = purityyy(result$clu,ground_truth,anno = c("Layer"))))
}




ari <- function(clu,ground_truth,anno = c("Layer")){
  #ground_truth <- layer_manual_PDAC
  library(aricode)
  library(gtools)
  g1 <- as.matrix(clu)
  g2 <- as.data.frame(ground_truth[[anno]])
  rownames(g2) <- rownames(ground_truth)
  g2 <- as.matrix(g2)
  g1 <- g1[mixedsort(rownames(g1)),]
  g2 <- g2[mixedsort(rownames(g2)),]
  rownames(g1) <- rownames(g2)
  # # 
  # # spot_order <- rownames(g1)
  # # rownames(g2) <- rownames(ground_truth)
  # # g2 <- g2[spot_order,]
  # g1 <- as.numeric(g1)
  # #g2 <- as.matrix(g2)
  # # g2 <- as.numeric(g2)
  # 
  # library(flexclust)
  # tab <- table(g2,g1)
  # #print(tab)
  # similarity <- randIndex(table(g2,g1))
  # #similarity <- compare(g2, g1, method="adjusted.rand")
  return(ARI(g1,g2))
}


ClusterPurity <- function(true_labels, pred_labels) {
  sum(apply(table(pred_labels, true_labels), 1, max)) / length(true_labels)
}

# 自定义F-score函数（宏平均）
calculate_fscore <- function(true_labels, pred_labels) {
  # 生成混淆矩阵
  confusion_matrix <- table(true_labels, pred_labels)
  
  # 计算每个类别的Precision和Recall
  precision <- diag(confusion_matrix) / colSums(confusion_matrix)
  recall <- diag(confusion_matrix) / rowSums(confusion_matrix)
  
  # 处理NA（除零情况）
  precision[is.na(precision)] <- 0
  recall[is.na(recall)] <- 0
  
  # 计算每个类别的F1
  f1_scores <- 2 * (precision * recall) / (precision + recall)
  f1_scores[is.na(f1_scores)] <- 0  # 处理NaN
  
  # 宏平均
  mean(f1_scores)
}


purityyy <- function(clu,ground_truth,anno = c("Layer")){
  library(funtimes)
  #ground_truth <- layer_manual_PDAC
  g1 <- as.matrix(clu)
  g2 <- as.data.frame(ground_truth[[anno]])
  rownames(g2) <- rownames(ground_truth)
  g2 <- as.matrix(g2)
  g1 <- g1[mixedsort(rownames(g1)),]
  g2 <- g2[mixedsort(rownames(g2)),] 
  # spot_order <- rownames(g1)
  rownames(g1) <- rownames(g2)
  # rownames(g2) <- rownames(ground_truth)
  # g2 <- g2[spot_order,]
  # g2 <- as.numeric(g2)
  # 
  # library(flexclust)
  # tab <- table(g2,g1)
  # #print(tab)
  similarity <- purity(g2,g1)
  # similarity <- sum(apply(tab, 2, max)) / sum(tab)
  #similarity <- compare(g2, g1, method="adjusted.rand")
  return(similarity$pur)
}










my_method_directly <- function(ST,cref,P_real,C_real,T_train,P_train,useARI = TRUE,Ss,mgs,gsva_mat,l1,l2,lc,l,m=0.02,g=-0.1,num = 5000){
  
  
  getwd()
  setwd("D:/大四上/毕设/CDSC-master/code/CDSC-master/CDSC-master/code")
  
  source("CDSC.R")
  source("function_help.R")
  
  #----------find Paramater ///----------
  result <- list()
  
  
  
  lambda1 <- c(l1)
  lambda2 <- c(l2)
  lambdaC <- c(lc)
  
  
  m_values <- c(m)
  num_iter_max <- num
  l <-c(l)
  g_values <- c(g)
  result$CDSC4$seedd = 44
  intra_sc_Data <- NULL
  intra_sc_Data$Indata$TerCondition = 10^-10
  
  
  
  ###使用参考数据集产生marker基因和C_ref矩阵
  
  new_bulk_T <- cbind(T_train,ST)
  
  
  
  
  #------------cyclic to do deconvolution----------
  pb <- txtProgressBar(style = 3)
  star_time <- Sys.time()
  result$CDSC4$Ss <- Ss
  result$CDSC4$Sg <- SM(new_bulk_T)
  result_para_c = NULL
  result_para_p = NULL
  pearson_para_c = NULL
  number_iter = NULL
  num = 1
  library(scater)
  
  para_find = NULL 
  library(dplyr)
  
  #gsva_result <- gsva_apply(new_bulk_T,inter_hard_3_new$markers_,"Poisson")
  for (l_i in 1:length(l)){
    for (g_i in 1:length(g_values)){
      LAMBDA <- caclulate_lamda(new_bulk_T, mgs, as.matrix(cref), l = l[l_i], g = g_values[g_i],gsva_result = gsva_mat)
      for (dir_i in 1:length(lambda1)){
        for (dir_j in 1:length(lambda2)){
          for (dir_k in 1:length(lambdaC)){
            result_CDSC = CDSC_1(new_bulk_T, as.matrix(cref), dim(cref)[2], 
                                 lambda1 = lambda1[dir_i], lambda2 = lambda2[dir_j], lambdaC = lambdaC[dir_k],LAMBDA = LAMBDA,
                                 10^-8,result$CDSC4$seedd,
                                 result$CDSC4$Ss,result$CDSC4$Sg,all_number = num_iter_max)
            
            for(m_i in 1:length(m_values)){
              result_para_c <- result_CDSC[[1]]
              
              
              result1 = NULL
              if(!all(is.na(result_para_c) == FALSE) ){
                break
              }
              
              
              result_para_p_all <- result_CDSC[[2]]
              result_para_p_all <- filter_p(result_para_p_all, LAMBDA, m = m_values[m_i])
              
              
              real_P_mixname <- (ncol(P_train) + 1):ncol(result_para_p_all)
              pre_P_mixname <- 1:ncol(P_train)
              result_para_p_real <- result_para_p_all[,real_P_mixname]
              
              
              
              result_para_p_pre <- result_para_p_all[,pre_P_mixname]
              number_iter <- result_CDSC[[3]]
              # result1 <- calculate_result_1(result_c = result_para_c,C_ref = data$indata$C_ref[genes,],
              #                                                      lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
              #                                                      number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = 10^-8,
              #                                                      leastnum = 3)
              # 
              
              # result1 <- calculate_result_1(result_c = result_para_c,result_p = result_para_p_real,
              #                               T = ST,C = C_real,C_ref = cref,P = NULL, 
              #                               lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
              #                               number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = intra_sc_Data$Indata$TerCondition,
              #                               leastnum = 3,ctlabels=NULL,result_p_1 = result_para_p_pre, P_1 = P_train, T_1 = T_train, result_p_combine = result_para_p_all,T_combine = new_bulk_T, l = l[l_i], m = m_values[m_i], g = g_values[g_i])
              result1 <- calculate_result_1(result_c = result_para_c,result_p = result_para_p_real,
                                            T = ST,C = C_real,C_ref = cref,P = P_real,
                                            lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
                                            number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = intra_sc_Data$Indata$TerCondition,
                                            leastnum = 3,ctlabels=NULL,result_p_1 = result_para_p_pre, P_1 = P_train, T_1 = T_train, result_p_combine = result_para_p_all,T_combine = new_bulk_T, l = l[l_i], m = m_values[m_i], g = g_values[g_i])
              
              para_find =  rbind(para_find,result1)
              
              num = num + 1
              setTxtProgressBar(pb, num/(length(lambda1)*length(lambda2)*length(lambdaC)*length(l)*length(m_values)))
            }
          }
        }
      }
    }
  }
  end_time <- Sys.time()
  close(pb)
  print(star_time)
  print(end_time)
  end_time-star_time
  
  ct <- Row_label(result_para_c,cref)
  colnames(result_para_c) <- ct
  rownames(result_para_p_real) <- ct
  
  return(list(para_find = para_find,result_c = result_para_c,result_p = result_para_p_real))
}

#a <- generate_binary_matrix(g1, result_dominant$card)

###画kmeans聚类图 注意这里的location值不可以太大，不然画不出来
draw_kmeans <- function(mat,colors = c("#FCD4C7","#A8AF7F","#F47C50"),lables = c("cluster2","cluster3","cluster1"),k=2,scale_size = 1){
  library(factoextra)
  library(cluster)
  getwd()
  setwd("D:/大四上/毕设/CDSC-master/code/CDSC-master/CDSC-master/code")
  
  source("CDSC.R")
  source("function_help.R")
  mat <- NORM(mat)
  
  location <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(mat),split="x"),"[",1)),y=as.numeric(sapply(strsplit(colnames(mat),split="x"),"[",2)))
  rownames(location) = colnames(mat)
  set.seed(44)
  mat <- t(mat)
  # 调用kmeans聚类算法 k = 4
  n <- nrow(mat)  
  
  km <- kmeans(mat, centers = k)
  clu <- as.matrix(km$cluster)
  
  # labels 是标签向量，spots 是 spot 数据
  labels <- as.character(clu)
  spots <- mat
  
  
  # 使用函数生成 0/1 矩阵
  binary_matrix <- generate_binary_matrix(labels, spots)
  binary_matrix <- t(binary_matrix)
  rownames(binary_matrix) <- lables
  binary_matrix <- binary_matrix[,rownames(location)]
  #("#CFA6B2","#FCD4C7","#6D889F","#F47C50","#A8AF7F"
  location <- location / scale_size
  library(CARD)
  p_cdsc <- CARD.visualize.pie(
    proportion = t((binary_matrix)),
    spatial_location = (location),
    colors = colors,
    radius = 0.37)
  print(p_cdsc)
  
  
  #library(ggplot2)
  # 
  # 
  # data <- cbind(location,clu)
  # # 创建散点图
  # ggplot(data, aes(x = x, y = y, color = clu)) +
  #   geom_point(size = 4) +
  #   scale_color_manual(values = c("purple","yellow","orange")) +
  #   labs(title = "Scatter Plot with Colored Labels", x = "X Axis Label", y = "Y Axis Label")
  # 
  # 
  
  return(clu)
}

my_method_inter <- function(data_all){
  
  
  getwd()
  setwd("D:/大四上/毕设/CDSC-master/code/CDSC-master/CDSC-master/code")
  
  source("CDSC.R")
  source("function_help.R")
  
  inter_hard_3_new <- data_all
  
  colnames(inter_hard_3_new$indata$P_train) <- sub("X","",colnames(inter_hard_3_new$indata$P_train))
  colnames(inter_hard_3_new$indata$P_test) <- sub("X","",colnames(inter_hard_3_new$indata$P_test))
  
  #----------find Paramater ///----------
  result <- list()
  lambda1 <- c(0.01,0.1,0.001)
  lambda2 <- c(0.001,0.1,0.01)
  lambdaC <- c(0.1,1,10)
  
  
  
  
  m_values <- seq(0,0.2,0.02)
  num_iter_max <- 3000
  l <-c(1,10,100)
  g_values <- c(-0.1,-0.2)
  result$CDSC4$seedd = 44
  intra_sc_Data <- NULL
  intra_sc_Data$Indata$TerCondition = 10^-8
  
  
  
  ###使用参考数据集产生marker基因和C_ref矩阵
  
  #genes <- unique(inter_hard_3_new$markers_$genes)
  new_bulk_T <- cbind(inter_hard_3_new$indata$T_train,inter_hard_3_new$indata$T_test)
  
  
  
  
  #ss_cal_output(spatial.pseudo$pseudo.data,spatial.pseudo2$pseudo.data, usefov = FALSE)
  
  
  
  #------------cyclic to do deconvolution----------
  pb <- txtProgressBar(style = 3)
  star_time <- Sys.time()
  result$CDSC4$Ss <- inter_hard_3_new$indata$ss
  result$CDSC4$Sg <- SM(new_bulk_T)
  result_para_c = NULL
  result_para_p = NULL
  pearson_para_c = NULL
  number_iter = NULL
  num = 1
  library(scater)
  
  para_find = NULL 
  library(dplyr)
  
  #gsva_result <- gsva_apply(new_bulk_T,inter_hard_3_new$markers_,"Poisson")
  for (l_i in 1:length(l)){
    for (g_i in 1:length(g_values)){
      LAMBDA <- caclulate_lamda(new_bulk_T, inter_hard_3_new$markers, as.matrix(inter_hard_3_new$indata$C_ref), l = l[l_i], g = g_values[g_i],gsva_result = inter_hard_3_new$indata$gsva_mat_combine)
      for (dir_i in 1:length(lambda1)){
        for (dir_j in 1:length(lambda2)){
          for (dir_k in 1:length(lambdaC)){
            result_CDSC = CDSC_1(new_bulk_T, as.matrix(inter_hard_3_new$indata$C_ref), dim(inter_hard_3_new$indata$C_ref)[2], 
                                 lambda1 = lambda1[dir_i], lambda2 = lambda2[dir_j], lambdaC = lambdaC[dir_k],LAMBDA = LAMBDA,
                                 10^-8,result$CDSC4$seedd,
                                 result$CDSC4$Ss,result$CDSC4$Sg,all_number = num_iter_max)
            
            for(m_i in 1:length(m_values)){
              result_para_c <- result_CDSC[[1]]
              
              
              result1 = NULL
              if(!all(is.na(result_para_c) == FALSE) ){
                break
              }
              
              
              result_para_p_all <- result_CDSC[[2]]
              result_para_p_all <- filter_p(result_para_p_all, LAMBDA, m = m_values[m_i])
              
              
              real_P_mixname <- (ncol(inter_hard_3_new$indata$P_train) + 1):ncol(result_para_p_all)
              pre_P_mixname <- 1:ncol(inter_hard_3_new$indata$P_train)
              result_para_p_real <- result_para_p_all[,real_P_mixname]
              
              
              
              result_para_p_pre <- result_para_p_all[,pre_P_mixname]
              number_iter <- result_CDSC[[3]]
              # result1 <- calculate_result_1(result_c = result_para_c,C_ref = data$indata$C_ref[genes,],
              #                                                      lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
              #                                                      number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = 10^-8,
              #                                                      leastnum = 3)
              # 
              
              result1 <- calculate_result_1(result_c = result_para_c,result_p = result_para_p_real,
                                            T = inter_hard_3_new$indata$T_test,C = inter_hard_3_new$indata$C,C_ref = inter_hard_3_new$indata$C_ref,P = inter_hard_3_new$indata$P_test,
                                            lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
                                            number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = intra_sc_Data$Indata$TerCondition,
                                            leastnum = 3,ctlabels=NULL,result_p_1 = result_para_p_pre, P_1 = inter_hard_3_new$indata$P_train, T_1 = inter_hard_3_new$indata$T_train, result_p_combine = result_para_p_all,T_combine = new_bulk_T, l = l[l_i], m = m_values[m_i], g = g_values[g_i])
              
              if(!is.null(result_para_p_real)){
                #mymethod <- draw_kmeans_result(result_para_p_real)
                #result1 <- cbind(result1,ari(mymethod))
                para_find =  rbind(para_find,result1)
              }
              
              
              num = num + 1
              setTxtProgressBar(pb, num/(length(lambda1)*length(lambda2)*length(lambdaC)*length(l)*length(m_values)))
            }
          }
        }
      }
    }
  }
  end_time <- Sys.time()
  close(pb)
  print(star_time)
  print(end_time)
  end_time-star_time
  
  
  return(para_find)
}




my_method_inter_directly <- function(data_all){
  
  
  getwd()
  setwd("D:/大四上/毕设/CDSC-master/code/CDSC-master/CDSC-master/code")
  
  source("CDSC.R")
  source("function_help.R")
  
  inter_hard_3_new <- data_all
  
  colnames(inter_hard_3_new$indata$P_train) <- sub("X","",colnames(inter_hard_3_new$indata$P_train))
  colnames(inter_hard_3_new$indata$P_test) <- sub("X","",colnames(inter_hard_3_new$indata$P_test))
  
  #----------find Paramater ///----------
  result <- list()
  lambda1 <- c(0.1)
  lambda2 <- c(0.001)
  lambdaC <- c(1)
  
  
  
  
  m_values <- seq(0)
  num_iter_max <- 3000
  l <-c(1)
  g_values <- c(-0.2)
  result$CDSC4$seedd = 44
  intra_sc_Data <- NULL
  intra_sc_Data$Indata$TerCondition = 10^-8
  
  
  
  ###使用参考数据集产生marker基因和C_ref矩阵
  
  #genes <- unique(inter_hard_3_new$markers_$genes)
  new_bulk_T <- cbind(inter_hard_3_new$indata$T_train[unique(inter_hard_3_new$markers$genes),],inter_hard_3_new$indata$T_test[unique(inter_hard_3_new$markers$genes),])
  
  
  
  
  #ss_cal_output(spatial.pseudo$pseudo.data,spatial.pseudo2$pseudo.data, usefov = FALSE)
  
  
  
  #------------cyclic to do deconvolution----------
  pb <- txtProgressBar(style = 3)
  star_time <- Sys.time()
  result$CDSC4$Ss <- inter_hard_3_new$indata$ss
  result$CDSC4$Sg <- SM(new_bulk_T)
  result_para_c = NULL
  result_para_p = NULL
  pearson_para_c = NULL
  number_iter = NULL
  num = 1
  library(scater)
  
  para_find = NULL 
  library(dplyr)
  
  #gsva_result <- gsva_apply(new_bulk_T,inter_hard_3_new$markers_,"Poisson")
  for (l_i in 1:length(l)){
    for (g_i in 1:length(g_values)){
      LAMBDA <- caclulate_lamda(new_bulk_T, inter_hard_3_new$markers, as.matrix(inter_hard_3_new$indata$C_ref), l = l[l_i], g = g_values[g_i],gsva_result = inter_hard_3_new$indata$gsva_mat_combine)
      for (dir_i in 1:length(lambda1)){
        for (dir_j in 1:length(lambda2)){
          for (dir_k in 1:length(lambdaC)){
            result_CDSC = CDSC_1(new_bulk_T, as.matrix(inter_hard_3_new$indata$C_ref), dim(inter_hard_3_new$indata$C_ref)[2], 
                                 lambda1 = lambda1[dir_i], lambda2 = lambda2[dir_j], lambdaC = lambdaC[dir_k],LAMBDA = LAMBDA,
                                 10^-8,result$CDSC4$seedd,
                                 result$CDSC4$Ss,result$CDSC4$Sg,all_number = num_iter_max)
            
            for(m_i in 1:length(m_values)){
              result_para_c <- result_CDSC[[1]]
              
              
              result1 = NULL
              if(!all(is.na(result_para_c) == FALSE) ){
                break
              }
              
              
              result_para_p_all <- result_CDSC[[2]]
              result_para_p_all <- filter_p(result_para_p_all, LAMBDA, m = m_values[m_i])
              
              
              real_P_mixname <- (ncol(inter_hard_3_new$indata$P_train) + 1):ncol(result_para_p_all)
              pre_P_mixname <- 1:ncol(inter_hard_3_new$indata$P_train)
              result_para_p_real <- result_para_p_all[,real_P_mixname]
              
              
              
              result_para_p_pre <- result_para_p_all[,pre_P_mixname]
              number_iter <- result_CDSC[[3]]
              # result1 <- calculate_result_1(result_c = result_para_c,C_ref = data$indata$C_ref[genes,],
              #                                                      lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
              #                                                      number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = 10^-8,
              #                                                      leastnum = 3)
              # 
              
              result1 <- calculate_result_1(result_c = result_para_c,result_p = result_para_p_real,
                                            T = inter_hard_3_new$indata$T_test,C = inter_hard_3_new$indata$C,C_ref = inter_hard_3_new$indata$C_ref,P = inter_hard_3_new$indata$P_test,
                                            lambda1=lambda1[dir_i], lambda2=lambda2[dir_j], lambdaC=lambdaC[dir_k],
                                            number_iter=number_iter,seedd = result$CDSC4$seedd,TerCondition = intra_sc_Data$Indata$TerCondition,
                                            leastnum = 3,ctlabels=NULL,result_p_1 = result_para_p_pre, P_1 = inter_hard_3_new$indata$P_train, T_1 = inter_hard_3_new$indata$T_train, result_p_combine = result_para_p_all,T_combine = new_bulk_T, l = l[l_i], m = m_values[m_i], g = g_values[g_i])
              
              if(!is.null(result_para_p_real)){
                #mymethod <- draw_kmeans_result(result_para_p_real)
                #result1 <- cbind(result1,ari(mymethod))
                para_find =  rbind(para_find,result1)
              }
              
              
              num = num + 1
              setTxtProgressBar(pb, num/(length(lambda1)*length(lambda2)*length(lambdaC)*length(l)*length(m_values)))
            }
          }
        }
      }
    }
  }
  end_time <- Sys.time()
  close(pb)
  print(star_time)
  print(end_time)
  end_time-star_time
  
  
  return(para_find)
}


# 
# ###画dominant图
# draw_domi_my <- function(mat){
#   mat <- NORM(mat)
#   
#   for (i in 1:ncol(mat)) {
#     max_val <- max(mat[,i])  # 获取每行的最大值
#     mat[,i] <- ifelse(((mat[,i] == max_val)), 1, 0)  #  将最大值置为1，其他置为0
#   }
#   
#   colors = c("#CFA6B2","#FCD4C7","#6D889F","#F47C50","#A8AF7F")#("#CFA6B2","#FCD4C7","#6D889F","#F47C50","#A8AF7F"
#   location <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(mat),split="x"),"[",1)),y=as.numeric(sapply(strsplit(colnames(mat),split="x"),"[",2)))
#   rownames(location) = colnames(mat)
#   p_cdsc <- CARD.visualize.pie(
#     proportion = t((mat)),
#     spatial_location = (location),
#     colors = colors,
#     radius = 0.4)
#   print(p_cdsc)
#   return(mat)
# }

# 
# ####画比例pie图
# #mat形式celltype*spot
# draw_pie_my <- function(mat,colors = c("#CFA6B2","#FCD4C7","#6D889F","#F47C50","#A8AF4F","pink","red","yellow","gray"),scale_size = 100){
#   mat <- NORM(mat)
#   #colors = c("#CFA6B2","#FCD4C7","#6D889F","#F47C50","#A8AF4F","#A8BF7F","#A8AF7F","#B8AF7F","#C8AF7F")#("#CFA6B2","#FCD4C7","#6D889F","#F47C50","#A8AF7F"
#   location <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(mat),split="x"),"[",1)),y=as.numeric(sapply(strsplit(colnames(mat),split="x"),"[",2)))
#   rownames(location) = colnames(mat)
#   location <- (location / scale_size)
#   p_cdsc <- CARD.visualize.pie(
#     proportion = t((mat)),
#     spatial_location = (location),
#     colors = colors,
#     radius = 0.37)
#   print(p_cdsc)
# }



card_my <- function(MOB_raw, sce_counts, sce_meta, location){
  
  library(SingleCellExperiment)
  library(pbmcapply)
  library(CARD)
  ct.varname = "cellType"
  sample.varname = "sampleInfo"
  if(!("sampleInfo" %in% colnames(sce_meta))){
    sce_meta$sampleInfo <- rep("sample0",nrow(sce_meta))
  }
  ct.select = as.character(unique(sce_meta$cellType))
  rownames(location) = colnames(MOB_raw)
  ##### create CARD object
  colnames(location) <- c("x","y")

  CARD_obj = createCARDObject(
    sc_count = sce_counts,
    sc_meta = sce_meta,
    spatial_count = MOB_raw,
    spatial_location = location,
    ct.varname = ct.varname,
    ct.select = ct.select,
    sample.varname = sample.varname,
    minCountGene = 0,
    minCountSpot = 0) 
  # start_time <- Sys.time()
  # start_mem <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
  # mem = as.numeric(bench_memory(CARD_obj <-  CARD_deconv_My(CARD_obj))[1]) / (1024)^2
  peak <- peakRAM(CARD_obj <-  CARD_deconv_My(CARD_obj))
  
  # end_time <- Sys.time()
  # end_mem <- as.numeric(bench_process_memory()[1][1]) / (1024)^2
  return(list(result =t(CARD_obj@Proportion_CARD),mem = peak$Peak_RAM_Used_MiB * 
                1.048576, time = peak$Elapsed_Time_sec))
}


CARD_deconv_My <- function(CARD_object) {
  library(fields)
  library(gtools)
  ct.select = CARD_object@info_parameters$ct.select
  ct.varname = CARD_object@info_parameters$ct.varname
  sample.varname = CARD_object@info_parameters$sample.varname
  cat(paste0("## create reference matrix from scRNASeq...\n"))
  sc_eset = CARD_object@sc_eset
  Basis_ref = createscRef(sc_eset, ct.select, ct.varname, 
                          sample.varname)
  Basis = Basis_ref$basis
  Basis = Basis[, colnames(Basis) %in% ct.select]
  Basis = Basis[, match(ct.select, colnames(Basis))]
  spatial_count = CARD_object@spatial_countMat
  commonGene = intersect(rownames(spatial_count), rownames(Basis))
  commonGene = commonGene[!(commonGene %in% commonGene[grep("mt-", 
                                                            commonGene)])]
  cat(paste0("## Select Informative Genes! ...\n"))
  common = commonGene
  Xinput = spatial_count
  rm(spatial_count)
  B = Basis
  rm(Basis)
  Xinput = Xinput[order(rownames(Xinput)), ]
  B = B[order(rownames(B)), ]
  B = B[rownames(B) %in% common, ]
  Xinput = Xinput[rownames(Xinput) %in% common, ]
  Xinput = Xinput[rowSums(Xinput) > 0, ]
  Xinput = Xinput[, colSums(Xinput) > 0]
  colsumvec = colSums(Xinput)
  Xinput_norm = sweep(Xinput, 2, colsumvec, "/")
  B = B[rownames(B) %in% rownames(Xinput_norm), ]
  B = B[match(rownames(Xinput_norm), rownames(B)), ]
  spatial_location = CARD_object@spatial_location
  spatial_location = spatial_location[rownames(spatial_location) %in% 
                                        colnames(Xinput_norm), ]
  spatial_location = spatial_location[match(colnames(Xinput_norm), 
                                            rownames(spatial_location)), ]
  norm_cords = spatial_location[, c("x", "y")]
  norm_cords$x = norm_cords$x - min(norm_cords$x)
  norm_cords$y = norm_cords$y - min(norm_cords$y)
  scaleFactor = max(norm_cords$x, norm_cords$y)
  norm_cords$x = norm_cords$x/scaleFactor
  norm_cords$y = norm_cords$y/scaleFactor
  ED <- rdist(as.matrix(norm_cords))
  cat(paste0("## Deconvolution Starts! ...\n"))
  set.seed(20200107)
  Vint1 = as.matrix(rdirichlet(ncol(Xinput_norm), rep(10, 
                                                      ncol(B))))
  colnames(Vint1) = colnames(B)
  rownames(Vint1) = colnames(Xinput_norm)
  b = rep(0, length(ct.select))
  isigma = 0.1
  epsilon = 1e-04
  phi = c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  kernel_mat <- exp(-ED^2/(2 * isigma^2))
  diag(kernel_mat) <- 0
  rm(ED)
  rm(Xinput)
  rm(norm_cords)
  gc()
  mean_X = mean(Xinput_norm)
  mean_B = mean(B)
  Xinput_norm = Xinput_norm * 0.1/mean_X
  B = B * 0.1/mean_B
  gc()
  ResList = list()
  Obj = c()
  for (iphi in 1:length(phi)) {
    res = CARDref(XinputIn = as.matrix(Xinput_norm), UIn = as.matrix(B), 
                  WIn = kernel_mat, phiIn = phi[iphi], max_iterIn = 1000, 
                  epsilonIn = epsilon, initV = Vint1, initb = rep(0, 
                                                                  ncol(B)), initSigma_e2 = 0.1, initLambda = rep(10, 
                                                                                                                 length(ct.select)))
    rownames(res$V) = colnames(Xinput_norm)
    colnames(res$V) = colnames(B)
    ResList[[iphi]] = res
    Obj = c(Obj, res$Obj)
  }
  Optimal = which(Obj == max(Obj))
  Optimal = Optimal[length(Optimal)]
  OptimalPhi = phi[Optimal]
  OptimalRes = ResList[[Optimal]]
  cat(paste0("## Deconvolution Finish! ...\n"))
  CARD_object@info_parameters$phi = OptimalPhi
  CARD_object@Proportion_CARD = sweep(OptimalRes$V, 1, rowSums(OptimalRes$V), 
                                      "/")
  CARD_object@algorithm_matrix = list(B = B * mean_B/0.1, 
                                      Xinput_norm = Xinput_norm * mean_X/0.1, Res = OptimalRes)
  CARD_object@spatial_location = spatial_location
  return(CARD_object)
  
}


dsa_my <- function(st,marker){
  require(CellMix)
  result <- NULL
  if(is.null(marker)){
    ML <- getMarkers(st)
  }else{
    md = marker
    ML = CellMix::MarkerList()
    ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
  }
  # start_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  # start_time = Sys.time()
  # mem = as.numeric(bench_memory(  result <- CellMix::ged(as.matrix(st), ML, method = "DSA",verbose=TRUE,
  #                       sscale = TRUE, exact = FALSE, maxIter=500, log=FALSE))[1]) / (1024)^2
  
  peak <- peakRAM(result <- CellMix::ged(as.matrix(st), ML, method = "DSA",verbose=TRUE,
                                 sscale = TRUE, exact = FALSE, maxIter=500, log=FALSE))
  # end_time = Sys.time()
  # end_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  result$c = result@fit@W
  result$p = result@fit@H
  return(list(result = result, mem = peak$Peak_RAM_Used_MiB * 
                1.048576, time = peak$Elapsed_Time_sec))
}

sskl_my <- function(st,marker){
  require(CellMix)
  ML = CellMix::MarkerList()
  md = marker
  ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
  # start_mem <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
  # start_time <- Sys.time()
  result <- CellMix::ged(as.matrix(st), ML, method = "ssKL",verbose=TRUE,
                         sscale = FALSE, maxIter=500, log = FALSE)
  # mem <- as.numeric(bench_memory(result <- CellMix::ged(as.matrix(st), ML, method = "ssKL",verbose=TRUE,
  #                                                       sscale = FALSE, maxIter=500, log = FALSE))[1]) / (1024)^2  # 重置并获取初始内存
  peak <- peakRAM(result <- CellMix::ged(as.matrix(st), ML, method = "ssKL",verbose=TRUE,
                                         sscale = FALSE, maxIter=500, log = FALSE))
  # end_time <- Sys.time()


  result$c = result@fit@W
  result$p = result@fit@H
  
  return(list(result = result,mem = peak$Peak_RAM_Used_MiB * 
                1.048576, time = peak$Elapsed_Time_sec))
}

ssFrobenius_my <- function(st,marker){
  require(CellMix)
  md = marker
  ML = CellMix::MarkerList()
  ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
  # start_mem <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
  # start_time <- Sys.time() 
  # mem = as.numeric(bench_memory(  result <- CellMix::ged(as.matrix(st), ML, method = "ssFrobenius",verbose=TRUE,
  #                        sscale = FALSE, maxIter = 500, log = FALSE))[1]) / (1024)^2
  
  peak <- peakRAM(result <- CellMix::ged(as.matrix(st), ML, method = "ssFrobenius",verbose=TRUE,
                                         sscale = FALSE, maxIter = 500, log = FALSE))
  # end_time <- Sys.time()
  # end_mem <- as.numeric(bench_process_memory()[1][1]) / (1024)^2  # 重置并获取初始内存
  
  result$c = result@fit@W
  result$p = result@fit@H
  
  return(list(result = result, mem = peak$Peak_RAM_Used_MiB * 
                1.048576, time = peak$Elapsed_Time_sec))
}


deconf_my <- function(st,marker){
  require(CellMix)
  md = marker
  ML = CellMix::MarkerList()
  ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
  # start_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  # start_time = Sys.time()
  # mem = as.numeric(bench_memory(  result <- CellMix::ged(as.matrix(st), ML, method = "deconf"))[1]) / (1024)^2
  peak <- peakRAM(result <- CellMix::ged(as.matrix(st), ML, method = "deconf"))
  # end_time = Sys.time()
  # end_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  result$c = result@fit@W
  result$p = result@fit@H
  return(list(results = result,mem = peak$Peak_RAM_Used_MiB * 
                1.048576, time = peak$Elapsed_Time_sec))
  
}

linCor <- function(yref,
                   iters = 100,
                   pval = .01,
                   n.types = NULL,
                   scree = 'cumvar',
                   logTransform = F){
  library(linseed)
  if (all(class(yref) != c('matrix'))){
    stop("matrix not supplied in yref")
  }
  
  if (logTransform){
    yref = log2(yref+1)
  }
  
  row.means = rowSums(yref)
  yref = yref[row.means != 0,]
  
  lo <- linseed::LinseedObject$new(yref)
  lo$calculatePairwiseLinearity()
  lo$calculateSpearmanCorrelation()
  lo$calculateSignificanceLevel(iters)
  lo$filterDatasetByPval(pval)
  
  if (all(lo$genes$pvals > pval)) {
    return(list(prop = matrix(rep(100,n.types * nrow(yref)),
                              ncol = n.types),
                sigs = NA))
  }
  
  if (is.null(n.types)){
    n.types = findNumberCells(yref,scree = scree)
  }
  
  lo$setCellTypeNumber(n.types)
  lo$project("full")
  lo$project("filtered")
  lo$smartSearchCorners(dataset="filtered", error="norm")
  lo$deconvolveByEndpoints()
  
  
  return(list(prop = lo$proportions,
              sig = lo$signatures))
  
}
linseed_my <- function(st,cref) {
  # start_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  start_time = Sys.time()
  # mem = as.numeric(bench_memory(  linseed.rs <-  linCor(yref = st,
  #                                                     iters = 100,
  #                                                     pval = 1000,
  #                                                     n.types = length(colnames(cref)),
  #                                                     scree = 'drop',
  #                                                     logTransform = F))[1]) / (1024)^2
  peak <- peakRAM(linseed.rs <-  linCor(yref = st,iters = 100, pval = 1000,n.types = length(colnames(cref)), scree = 'drop',logTransform = F))
  end_time = Sys.time()
  # end_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  
  result <- NULL
  result$p <- linseed.rs[[1]]
  result$c <- linseed.rs[[2]]
  cref <- cref[intersect(rownames(cref),rownames(result$c)),]
  ct <- Row_label(result$c,cref)
  colnames(result$c) <- ct
  rownames(result$p) <- ct
  return(list(result = result,mem = peak$Peak_RAM_Used_MiB * 
                1.048576, time = peak$Elapsed_Time_sec))
  
}


cellid_my <- function(st,ref_c) {
  require(CellDistinguisher)
  result <- NULL
  result$CellDistinguisher <- CellDistinguisher::gecd_CellDistinguisher(
    st,
    genesymb=rownames(st),
    numCellClasses = dim(ref_c)[2],
    minDistinguisherAlternatives=1,
    maxDistinguisherAlternatives=100,
    minAlternativesLengthsNormalized=0.5,
    expressionQuantileForScale = 0.75,
    expressionQuantileForFilter=0.999,
    expressionConcentrationRatio=0.333,
    probesWithGenesOnly = FALSE,
    verbose=0)
  result$CellDist.deconv <-
    tryCatch(CellDistinguisher::gecd_DeconvolutionByDistinguishers(
      as.matrix(bulkData$indata$T),
      result$CellDistinguisher$bestDistinguishers,
      nonNegativeOnly = FALSE,
      convexSolution = FALSE,
      verbose = 0),
      error = function(e) return(list(sampleCompositions =
                                        matrix(rep(1/dim(bulkData$indata$C)[2],
                                                   dim(bulkData$indata$C)[2]*ncol(bulkData$indata$T)),
                                               ncol=dim(bulkData$indata$C)[2]))))
  result$CellDistinguisher$p <- result$CellDist.deconv$sampleComposition
  result$CellDistinguisher$c <- result$CellDist.deconv$cellSubclassSignatures
  # labels <- Row_label((result$CellDistinguisher$c),ref_c)
  # rownames(result$CellDistinguisher$p) <- labels
  # colnames(result$CellDistinguisher$c) <- labels
  return(result$CellDistinguisher)
}

rnasieve_my <- function(st,sc,sc_meta,STRING_name,path){
  library(jsonlite)
  colnames(sc) <- sc_meta$cellType
  sc_groupby_class <- lapply(unique(sc_meta$cellType),function(cellType) {
    as.matrix(sc[ ,colnames(sc) %in% cellType])
  })
  names(sc_groupby_class) <- unique(sc_meta$cellType) 
  output_path_sc <- paste(path, "/sc_",STRING_name,".json",sep = "")
  output_path_bulk <- paste(path, "/bulk_",STRING_name,".json",sep = "")
  write_json(sc_groupby_class,output_path_sc)
  write_json(st,output_path_bulk)
  gene_name <- rownames(sc)
  write.csv(gene_name,paste(path, "/genename_",STRING_name,".txt",sep = ""))
  
}

rna_sieve_read <- function(STRING_name,real_p)
{
  result_path = "H:/LinYifan/rnasieve/" 
  result = NULL
  result$p <- read.csv(paste(result_path,"result_p_",STRING_name,".csv",sep = ""),header = T,row.names = 1)
  result$c <- read.csv(paste(result_path,"result_c_",STRING_name,".csv",sep = ""),header = T,row.names = 1)
  result$p <- t(result$p)
  result$p <- as.matrix(result$p)
  result$c <- as.matrix(result$c)
  ###rna-sieve results like : Bulk1 Bulk2....
  colnames(result$p) <- colnames(real_p)
  return(result)
}


cell2location_output_my <- function(ST, sc, sc_meta, spatial, STRINGNAME, path) {
  
  write.csv(ST, paste(path, STRINGNAME, "_ST.csv", sep = ""))
  write.csv(sc, paste(path, STRINGNAME, "_sc.csv", sep = ""))
  write.csv(sc_meta, paste(path, STRINGNAME, "_sc_meta.csv", sep = ""))
  write.csv(spatial, paste(path, STRINGNAME, "_sp.csv", sep = ""))
  
  
}




getClusterIndex.dominant.map <- function(Spatial.Deconv.Data,Methods = NULL ,anno = c("Layer")){
  if(is.null(Spatial.Deconv.Data$RESULTSMap)){
    if(is.null(Spatial.Deconv.Data$runSpatialDSSCBEResult)){
      stop("error")
    }
    
  }else{
    if(is.null(Methods)){
      Methods = names(Spatial.Deconv.Data$RESULTSMap$RESULTS)
    }
  }
  
  if(!is.null(Spatial.Deconv.Data$runSpatialDSSCBEResult)){
    Methods = c(Methods,"My Method")
  }
  aripurity.dominant <- NULL
  for (method in Methods) {
    if(is.matrix(Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]])){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]
      
    }else if(method == "Deconf" || method == "DSA" || method == "ssFrobenius" || method == "ssKL"){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]$p
      
    }else if(method == "STdeconvolve" || method == "TOAST/-P" || method == "TOAST"|| method == "TOAST/P-"){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]$result_p
      
    }else if(method == "My Method"){
      result_p = Spatial.Deconv.Data$runSpatialDSSCBEResult$result_para_p_real
    } else {
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]][[1]]
    }
    ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData
    aripurity.dominant[[method]] <- getClusterIndex.dominant(result_p,ground_truth,anno)
    
  }
  if(is.null(Spatial.Deconv.Data$RESULTSMap$ClusterIndex.dominant.MAP)){
    Spatial.Deconv.Data$RESULTSMap$ClusterIndex.dominant.MAP <- round(as.matrix(do.call(rbind, lapply(aripurity.dominant,unlist))),4)
  }else{
    Spatial.Deconv.Data$RESULTSMap$ClusterIndex.dominant.MAP <- unique(rbind(Spatial.Deconv.Data$RESULTSMap$ClusterIndex.dominant.MAP,round(as.matrix(do.call(rbind, lapply(aripurity.dominant,unlist))),4)))
  }
  return(Spatial.Deconv.Data)
}

getClusterIndex.dominant <- function(result_p,ground_truth,anno = c("Layer")){
  ###单独计算优势细胞类型聚类指标的函数
  ##1.得到每个spot的优势细胞类型
  ##input : groundtruth 为整个anno
  library(aricode)
  library(gtools)
  library(dendextend)
  clu <- getDominant(result_p)
  if(nrow(clu) != nrow(ground_truth)){
    stop("Nrows don't match")
  }
  # clu <- clu[mixedsort(rownames(clu)),]
  ground_truth <- ground_truth[rownames(clu),anno]
  clu <- c(clu$main)
  
  
  conf_mat <- confusion_matrix_clustering(ground_truth,clu)
  
  return(list(ARI = ARI(ground_truth,clu),NMI = NMI(ground_truth,clu),Purity = ClusterPurity(ground_truth,clu),FMI = FM_index(ground_truth,clu)[1],F1 = f1_score(conf_mat)))
  
  
  
}


getDominant <- function(proportion){
  library(gtools)
  
  P = as.data.frame(proportion)
  P = P[, mixedsort(colnames(P))]
  P = P[mixedsort(rownames(P)),]##mix the celltypes and the samples
  P <- t(P)
  P <- as.data.frame(P)
  max <- apply(P, 1, max, na.rm = TRUE)
  dominant = as.data.frame(max)
  colnames(dominant) <- c("main")
  
  for (i in 1:nrow(P)) {
    max.index = which.max(P[i,])
    dominant[i,] <- colnames(P)[max.index]
  }
  rownames(dominant) <- rownames(P)
  return(dominant)
}


# 输入：真实标签 (truelabels) 和预测标签 (predlabels)
confusion_matrix_clustering <- function(truelabels, predlabels) {
  n <- length(truelabels)
  
  # 初始化计数
  TP <- 0  # 真实同类，预测同类
  FP <- 0  # 真实不同类，预测同类
  FN <- 0  # 真实同类，预测不同类
  TN <- 0  # 真实不同类，预测不同类
  
  # 遍历所有样本对
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # 真实是否同类
      true_same <- (truelabels[i] == truelabels[j])
      # 预测是否同类
      pred_same <- (predlabels[i] == predlabels[j])
      
      # 更新混淆矩阵
      if (true_same && pred_same) {
        TP <- TP + 1
      } else if (!true_same && pred_same) {
        FP <- FP + 1
      } else if (true_same && !pred_same) {
        FN <- FN + 1
      } else {
        TN <- TN + 1
      }
    }
  }
  
  # 返回混淆矩阵
  matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE,
         dimnames = list(c("True Same", "True Different"),
                         c("Pred Same", "Pred Different")))
}

precision <- function(conf_mat) {
  TP <- conf_mat[1, 1]
  FP <- conf_mat[2, 1]
  return(TP / (TP + FP))
}

recall <- function(conf_mat) {
  TP <- conf_mat[1, 1]
  FN <- conf_mat[1, 2]
  return(TP / (TP + FN))
}

f1_score <- function(conf_mat) {
  p <- precision(conf_mat)
  r <- recall(conf_mat)
  return(2 * (p * r) / (p + r))
}

my_Seurat <- function(sccount,scmeta,count,k.score = 30,k.weight = 50) {
  library(Seurat)
  library(dplyr)
  library(SeuratDisk)
  
  sc_rna <- CreateSeuratObject(
    counts = sccount,          # 输入计数矩阵
    project = "scRNA_project",  # 项目名称（自定义）
    meta.data = scmeta,         # 添加细胞注释
    assay = "RNA"
  )
  spatial <- CreateSeuratObject(
    counts = count,
    project = "Spatial_project",
    assay = "RNA"  # 指定为空间数据
  )
  celltype_final = "cellType"
  
  sc_rna <- SCTransform(sc_rna)
  spatial <- SCTransform(spatial)
 
  # mem = as.numeric(bench_memory(anchors <- FindTransferAnchors(reference=sc_rna, query = spatial, dims = 1:30, normalization.method = 'SCT',k.score = k.score))[1]) / (1024)^2
  # mem = as.numeric(bench_memory(predictions <- TransferData(anchorset = anchors, refdata = sc_rna@meta.data[,celltype_final], dims = 1:30,k.weight = k.weight))[1]) / (1024)^2
  peak <- peakRAM({anchors <- FindTransferAnchors(reference=sc_rna, query = spatial, dims = 1:30, normalization.method = 'SCT',k.score = k.score)
  predictions <- TransferData(anchorset = anchors, refdata = sc_rna@meta.data[,celltype_final], dims = 1:30,k.weight = k.weight)})
  end_time = Sys.time()
  # end_mem = as.numeric(bench_process_memory()[1][1]) / (1024)^2
  predictions <- t(predictions)
  predictions <- predictions[2:(nrow(predictions) - 1),]
  pre <- apply(predictions,2,as.numeric)
  rownames(pre) <- rownames(predictions)
  predictions <- pre
  rownames(predictions) <- gsub("prediction.score.","",rownames(predictions)) 
  return(list(result = predictions,mem = peak$Peak_RAM_Used_MiB * 1.048576, time = peak$Elapsed_Time_sec))
  }


TangramPy <- function(Spatial.Deconv.Data) {

  library(reticulate)
  ST_count = Spatial.Deconv.Data$spatial.real.cpm.marker
  ST_count = t(ST_count)
  ST_count <- as.data.frame(ST_count)
  ref_rna.count.marker <- calculateCPM(Spatial.Deconv.Data$sc$data)[unique(Spatial.Deconv.Data$markerslist$gene),] ###count datas,with marker genes
  ref_rna.count.marker = t(ref_rna.count.marker)
  ref_rna.count.marker <- as.data.frame(ref_rna.count.marker)
  ref_rna_meta <- Spatial.Deconv.Data$sc$meta
  location <- Spatial.Deconv.Data$location.real
  use_condaenv("D:/softwares/anaconda/envs/gst/python.exe")
  source_python("H:/LinYifan/spatialDSSC-BE/R/TangramMy.py")

  result <- py$TangramDeconv(ST_count, location, ref_rna.count.marker, ref_rna_meta)
  
  return(result)
  
  
}






get.resultPs <- function(Spatial.Deconv.Data,Methods = NULL){
  if(is.null(Spatial.Deconv.Data$RESULTSMap)){
    if(is.null(Spatial.Deconv.Data$runSpatialDSSCBEResult)){
      stop("error")
    }
    
  }else{
    if(is.null(Methods)){
      Methods = names(Spatial.Deconv.Data$RESULTSMap$RESULTS)
    }
  }
  
  if(!is.null(Spatial.Deconv.Data$runSpatialDSSCBEResult)){
    Methods = c(Methods,"My Method")
  }
  for (method in Methods) {
    if(is.matrix(Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]])){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]
      
    }else if(method == "Deconf" || method == "DSA" || method == "ssFrobenius" || method == "ssKL"){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]$p
      
    }else if(method == "STdeconvolve" || method == "TOAST/-P" || method == "TOAST"|| method == "TOAST/P-"){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]$result_p
      
    }else if(method == "My Method"){
      result_p = Spatial.Deconv.Data$runSpatialDSSCBEResult$result_para_p_real
    } else {
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]][[1]]
    }
    # result_p = as.data.frame(result_p)
    result_p = result_p[, mixedsort(colnames(result_p))]
    result_p = result_p[mixedsort(rownames(result_p)),]##mix the celltypes and the samples
    result_p <- NORM(result_p)
    Spatial.Deconv.Data$RESULTSMap$resultPs[[method]] <- result_p
    
  }

  return(Spatial.Deconv.Data)
}



getresultPs.dominant <- function(Spatial.Deconv.Data,Methods = NULL) {
  if(is.null(Spatial.Deconv.Data$RESULTSMap)){
    if(is.null(Spatial.Deconv.Data$runSpatialDSSCBEResult)){
      stop("error")
    }
    
  }else{
    if(is.null(Methods)){
      Methods = names(Spatial.Deconv.Data$RESULTSMap$RESULTS)
    }
  }
  
  if(!is.null(Spatial.Deconv.Data$runSpatialDSSCBEResult)){
    Methods = c(Methods,"My Method")
  }
  for (method in Methods) {
    if(is.matrix(Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]])){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]
      
    }else if(method == "Deconf" || method == "DSA" || method == "ssFrobenius" || method == "ssKL"){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]$p
      
    }else if(method == "STdeconvolve" || method == "TOAST/-P" || method == "TOAST"|| method == "TOAST/P-"){
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]]$result_p
      
    }else if(method == "My Method"){
      result_p = Spatial.Deconv.Data$runSpatialDSSCBEResult$result_para_p_real
    } else {
      result_p = Spatial.Deconv.Data$RESULTSMap$RESULTS[[method]][[1]]
    }
    library(gtools)
    
    # result_p = as.data.frame(result_p)
    result_p = result_p[, mixedsort(colnames(result_p))]
    result_p = result_p[mixedsort(rownames(result_p)),]##mix the celltypes and the samples
    binary_mat <- matrix(1e-5, 
                         nrow = nrow(result_p), 
                         ncol = ncol(result_p),
                         dimnames = dimnames(result_p))
    
    # 找到每列最大值的索引
    max_indices <- apply(result_p, 2, which.max)
    
    # 通过矩阵索引赋值
    binary_mat[cbind(max_indices, 1:ncol(result_p))] <- 1
    Spatial.Deconv.Data$RESULTSMap$resultPs.dominant[[method]] <- binary_mat
    
  }
  
  return(Spatial.Deconv.Data)
}
