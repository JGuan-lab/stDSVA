


simulateSpotsSpatialSc <- function(sp.exp, sp.meta, sp.loc, window = 500, CoordinateX = 'X', CoordinateY = 'Y', CoordinateZ = NULL) {
###simulate spots use python
###input:spatial data with single-cell resolution
  library(reticulate)
  source_python("simulate.py")  
  pd <- import("pandas")
  sp.exp <- pd$DataFrame(sp.exp,index = sp.meta$cellID)
  sp.meta <- pd$DataFrame(sp.meta,index = sp.meta$cellID)
  sp.loc <- pd$DataFrame(sp.loc,index = sp.meta$cellID)
  resultsCombine = py$Simulated(sp.exp, sp.meta, sp.loc, CoordinateX, CoordinateY, window, CoordinateZ) 
  return(resultsCombine)
}




simulateSpots <- function(comData, CoordinateX = 'X', CoordinateY = 'Y', window = 500) {
  if("FOV" %in% colnames(comData$Train$pData)){
    comData$Test$simulateData <- simulateSpotsSpatialSc(t(comData$Test$data),comData$Test$pData,comData$Test$pData[,c(CoordinateX, CoordinateY, "FOV")],window = window,CoordinateZ= "FOV")
  }else{
    comData$Test$simulateData <- simulateSpotsSpatialSc(t(comData$Test$data),comData$Test$pData,comData$Test$pData[,c(CoordinateX, CoordinateY)],window = window)
    
    }
  ###simulate spots using DSSC object
  
  comData <- getSimulateSpotsReal(comData)
  index_nonezero = which(colSums(comData$T_test) != 0)
  comData$T_test <- comData$T_test[, colSums(comData$T_test) != 0]
  
  comData$P_test <- comData$P_test[, index_nonezero]
  return(comData)
  
}


simulateSpotsGridPseudo <- function(comData, CoordinateX = 'X', CoordinateY = 'Y', window = 500) {
  ###simulate data by GRID using object
  if("FOV" %in% colnames(comData$Train$pData)){
    comData$Train$simulateData <- simulateSpotsSpatialSc(t(comData$Train$data),comData$Train$pData,comData$Train$pData[,c(CoordinateX, CoordinateY, "FOV")],window = window,CoordinateZ= "FOV")
    comData <- getSimulateSpotsPseudo(comData)
    index_nonezero = which(colSums(comData$T_train) != 0)
    comData$T_train <- comData$T_train[, colSums(comData$T_train) != 0]
    
    comData$P_train <- comData$P_train[, index_nonezero]
  }else{

    comData$Train$simulateData <- simulateSpotsSpatialSc(t(comData$Train$data),comData$Train$pData,comData$Train$pData[,c(CoordinateX, CoordinateY)],window = window)
    comData <- getSimulateSpotsPseudo(comData)
    }

  return(comData)
}


getSimulateSpotsReal <- function(comData) {
  
  ###get T and P from simulation results
  combined_cell_counts <-  comData$Test$simulateData[[1]]
  combined_spot_loc <- comData$Test$simulateData[[2]]
  combined_spot_exp <- comData$Test$simulateData[[3]]
  combined_spot_clusters <- comData$Test$simulateData[[4]]
  
  
  comData$T_test <- t(combined_spot_exp)
  dim(comData$T_test)
  dim(combined_spot_loc)
  colnames(comData$T_test) <-  paste("real",combined_spot_loc$x,combined_spot_loc$y,sep = "_")
  rownames(comData$T_test) <- rownames(comData$Test$data)
  comData$P_test <- apply(combined_spot_clusters,1,function(x) x / sum(x))
  colnames(comData$P_test) <- colnames(comData$T_test)
  dim(comData$P_test)
  return(comData)
  
}


getSimulateSpotsPseudo <- function(comData) {
  
  ###get T and P from simulation results
  combined_cell_counts <-  comData$Train$simulateData[[1]]
  combined_spot_loc <- comData$Train$simulateData[[2]]
  combined_spot_exp <- comData$Train$simulateData[[3]]
  combined_spot_clusters <- comData$Train$simulateData[[4]]
  
  
  comData$T_train <- t(combined_spot_exp)
  dim(comData$T_train)
  dim(combined_spot_loc)
  if("FOV" %in% colnames(comData$Train$pData)){
    colnames(comData$T_train) <-  paste("pseudo",combined_spot_loc$fov,combined_spot_loc$x,combined_spot_loc$y,sep = "_")
    
  }else{
    colnames(comData$T_train) <-  paste("pseudo",combined_spot_loc$x,combined_spot_loc$y,sep = "_")
    
  }
  rownames(comData$T_train) <- rownames(comData$Train$data)
  comData$P_train <- apply(combined_spot_clusters,1,function(x) x / sum(x))
  colnames(comData$P_train) <- colnames(comData$T_train)
  dim(comData$P_train)
  return(comData)
  
}


simulateSpotsMy <- function(scRNAseq,scRNA.meta,meta.var = "cellType",cellNum = 1:20,
                            n = 1000,
                            verbose = TRUE,UMI=2000) {
  # simulate spots using single-cell RNAseq data
  # simulate by sampling
  # Check variables
  set.seed(44)
  if (!is.character(meta.var)) stop("ERROR: meta.var must be a character string!")
  if (!is.numeric(n)) stop("ERROR: n must be an integer!")
  if (!is.logical(verbose)) stop("ERROR: verbose must be a logical object!")
  
  suppressMessages(require(DropletUtils)) # For the downsampling
  suppressMessages(require(purrr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tidyr))
  
  scRNA.meta[, meta.var] <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                 x = scRNA.meta[, meta.var],
                                 perl = TRUE)
  has_coords <- all(c("X", "Y") %in% colnames(scRNA.meta))
  
  print("Generating synthetic test spots...")
  start_gen <- Sys.time()
  # create progress bar
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  # Save count matrix
  #count_mtrx <- as.matrix(se_obj@assays$RNA@counts)
  count_mtrx <- scRNAseq
  ds_spots <- lapply(seq_len(n), function(i) {
    
    # Select between 2 and 10 cells randomly from the count matrix
    cell_pool <- sample(colnames(count_mtrx), sample(x = cellNum, size = 1))
    
    # Determine the weight each cell will have on the synthetic spot
    # weigh <-runif(length(cell_pool))
    # weigh <- weigh/sum(weigh)
    
    # We're not going to sum the reads as is bc spots are **enriched**
    # so we'll add up the counts and downsample to the ~depth of a typical spot.
    
    # Create a name for the spot with the necessary info to deconvolute it
    pos <- which(colnames(count_mtrx) %in% cell_pool)
    tmp_ds <- scRNA.meta[pos, ] %>% mutate(weight = 1)
    # tmp_ds[, "weight"] <- weigh
    name_simp <- paste("spot_", i, sep = "")
    
    spot_ds <- tmp_ds %>%
      dplyr::select(all_of(meta.var), weight) %>%
      # dplyr::mutate(clust_vr = paste("clust_",
      #                                tmp_ds[, clust_vr], sep = "")) %>%
      dplyr::group_by(!! sym(meta.var)) %>%
      dplyr::summarise(sum_weights = sum(weight)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = all_of(meta.var),
                         values_from = sum_weights) %>%
      dplyr::mutate(name = name_simp)
    
    
    # If coordinates are available, calculate average position
    spot_coords <- NULL
    if (has_coords) {
      spot_coords <- tmp_ds %>%
        summarise(X = mean(X, na.rm = TRUE),
                  Y = mean(Y, na.rm = TRUE))
    }
    
    # Generate synthetic spot
    
    ## Here we multiply each vector by its weight
    # weighted <- lapply(1:length(cell_pool), function(ii) expr <- as.integer(round(weigh[[ii]]*count_mtrx[hvg,cell_pool[[ii]]],0)))
    
    ## Next step we add all the vectors by position
    # syn_spot <- Reduce(`+`,weighted)
    # ret_ds <- data.frame(gene=hvg, tmp=syn_spot)
    
    ## Here we add up the counts of each cell
    syn_spot <- rowSums(as.matrix(count_mtrx[, cell_pool])); sum(syn_spot)
    names_genes <- names(syn_spot)
    
    ## Downsample
    ### 25k is a bit above average 20k UMIs observed in spatial transcriptomics data then downsample to 20k
    if (sum(syn_spot) > UMI * 1.25) {
      syn_spot_sparse <- DropletUtils::downsampleMatrix(Matrix::Matrix(syn_spot, sparse = T),
                                                        prop = UMI / sum(syn_spot))
      # syn_spot_sparse <- Matrix::Matrix(syn_spot, sparse = T)
      } else {
      syn_spot_sparse <- Matrix::Matrix(syn_spot, sparse = T)
    }
    
    rownames(syn_spot_sparse) <- names_genes
    colnames(syn_spot_sparse) <- name_simp
    
    # update progress bar
    setTxtProgressBar(pb, i)
    
    return(list(syn_spot_sparse, spot_ds, spot_coords))
  })
  
  ds_syn_spots <- purrr::map(ds_spots, 1) %>%
    base::Reduce(function(m1, m2) cbind(unlist(m1), unlist(m2)), .)
  
  # Generate dataframe of spot characteristic
  ds_spots_metadata <- purrr::map(ds_spots, 2) %>%
    dplyr::bind_rows() %>%
    data.frame()
  
  ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
  
  # change column order so that its progressive
  lev_mod <- gsub("[\\+|\\ |\\/]", ".", unique(scRNA.meta[, meta.var]))
  colnames(ds_spots_metadata) <- gsub("[\\+|\\ |\\/]", ".", colnames(ds_spots_metadata))
  # all_cn <- c(paste("clust_", lev_mod, sep = ""), "name") # This was set to deal when cluster names are numeric
  
  # Check if there are missing columns (Cell types not selected) and add them as all 0s
  if (sum(lev_mod %in% colnames(ds_spots_metadata)) == (length(unique(scRNA.meta[, meta.var])) + 1)) {
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  } else {
    
    missing_cols <- lev_mod[which(!lev_mod %in% colnames(ds_spots_metadata))]
    ds_spots_metadata[missing_cols] <- 0
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  }
  # Generate dataframe of spot coordinates if available
  ds_spots_coords <- NULL
  if (has_coords) {
    ds_spots_coords <- purrr::map(ds_spots, 3) %>%
      dplyr::bind_rows() %>%
      data.frame()
  }
  # Close progress bar
  close(pb)
  
  print(sprintf("Generation of %s test spots took %s mins", n,
                round(difftime(Sys.time(), start_gen, units = "mins"), 2)))
  # Prepare output list
  output_list <- list(
    topic_profiles = ds_syn_spots,
    cell_composition = ds_spots_metadata
  )
  
  # Add coordinates to output if available
  if (has_coords) {
    output_list$spot_coordinates <- ds_spots_coords
    print("output consists of a list with three elements: 1) weighted count matrix, 2) metadata for each spot, 3) coordinates for each spot")
  } else {
    print("output consists of a list with two elements: 1) weighted count matrix, 2) metadata for each spot")
  }
  
  return(output_list)
  # print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  # 
  # return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}





simulateSpotsSampling <- function(comData,n=NULL,cellNum=NULL,UMI=2000) {
  ###simulate data by sampling using object
  if(is.null(n)){
    trainSimulate <- simulateSpotsMy(comData$Train$data,comData$Train$pData,n = length(colnames(comData$T_test)),cellNum=cellNum,UMI=UMI)
    
  }else{
    trainSimulate <- simulateSpotsMy(comData$Train$data,comData$Train$pData,n = n,cellNum=cellNum,UMI=UMI)
    
  }
  combined_spot_exp <- as.matrix(trainSimulate$topic_profiles)
  combined_spot_clusters <- trainSimulate$cell_composition
  comData$T_train <- combined_spot_exp
  colnames(combined_spot_exp) <-  paste("pseudo",colnames(combined_spot_exp),sep = "_")
  colnames(comData$T_train) <- colnames(combined_spot_exp)
  comData$P_train <- apply(combined_spot_clusters,1,function(x) x / sum(x))
  colnames(comData$P_train) <- colnames(comData$T_train)
  return(comData)
  
  
}

min_max_scale <- function(x, new_min, new_max) {
  (x - min(x)) / (max(x) - min(x)) * (new_max - new_min) + new_min
}

simulateSpotsSampling.spatial <- function(comData,n=NULL,cellNum=NULL,UMI=2000) {
  ###simulate data by sampling using DSSC object
  if(is.null(n)){
    trainSimulate <- simulateSpotsMy(comData$Train$data,comData$Train$pData,n = length(colnames(comData$T_test)),cellNum=cellNum,UMI=UMI)
    
  }else{
    trainSimulate <- simulateSpotsMy(comData$Train$data,comData$Train$pData,n = n,cellNum=cellNum,UMI=UMI)
    
  }
  combined_spot_exp <- as.matrix(trainSimulate$topic_profiles)
  combined_spot_clusters <- trainSimulate$cell_composition
  positions <- trainSimulate$spot_coordinates
  
  comData$T_train <- combined_spot_exp
  colnames(combined_spot_exp) <-  paste("pseudo",positions$X,positions$Y,sep = "_")
  colnames(comData$T_train) <- colnames(combined_spot_exp)
  comData$P_train <- apply(combined_spot_clusters,1,function(x) x / sum(x))
  colnames(comData$P_train) <- colnames(comData$T_train)
  return(comData)
  
  
}
