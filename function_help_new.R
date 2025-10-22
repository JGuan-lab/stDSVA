


#' This functions takes in a seurat object and creates different mixtures resembling spots at different proportions
#'
#' @param se_obj seurat object.
#' @param clust_vr Name of the variable containing the cell clustering
#' @param n number of spots to generate
#' @param verbose Name of the variable containing the cell clustering
#' @return This function returnsa dataframe where each column is the proportion of each spot and a sparse matrix with the synthetic spots with the genenames as rownames and colnames being the spot names
#' @export
#' @examples
#'

test_spot_fun <- function(se_obj,
                          clust_vr,
                          n = 1000,
                          verbose = TRUE) {
  
  # Check variables
  if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.numeric(n)) stop("ERROR: n must be an integer!")
  if (!is.logical(verbose)) stop("ERROR: verbose must be a logical object!")
  
  suppressMessages(require(DropletUtils)) # For the downsampling
  suppressMessages(require(purrr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tidyr))
  
  se_obj@meta.data[, clust_vr] <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                       x = se_obj@meta.data[, clust_vr],
                                       perl = TRUE)
  print("Generating synthetic test spots...")
  start_gen <- Sys.time()
  # create progress bar
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  # Save count matrix
  #count_mtrx <- as.matrix(se_obj@assays$RNA@counts)
  count_mtrx <- GetAssayData(object = se_obj, slot = "counts")
  ds_spots <- lapply(seq_len(n), function(i) {
    
    # Select between 2 and 10 cells randomly from the count matrix
    cell_pool <- sample(colnames(count_mtrx), sample(x = 1:20, size = 1))
    
    # Determine the weight each cell will have on the synthetic spot
    # weigh <-runif(length(cell_pool))
    # weigh <- weigh/sum(weigh)
    
    # We're not going to sum the reads as is bc spots are **enriched**
    # so we'll add up the counts and downsample to the ~depth of a typical spot.
    
    # Create a name for the spot with the necessary info to deconvolute it
    pos <- which(colnames(count_mtrx) %in% cell_pool)
    tmp_ds <- se_obj@meta.data[pos, ] %>% mutate(weight = 1)
    # tmp_ds[, "weight"] <- weigh
    name_simp <- paste("spot_", i, sep = "")
    
    spot_ds <- tmp_ds %>%
      dplyr::select(all_of(clust_vr), weight) %>%
      # dplyr::mutate(clust_vr = paste("clust_",
      #                                tmp_ds[, clust_vr], sep = "")) %>%
      dplyr::group_by(!! sym(clust_vr)) %>%
      dplyr::summarise(sum_weights = sum(weight)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = all_of(clust_vr),
                         values_from = sum_weights) %>%
      dplyr::mutate(name = name_simp)
    
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
    if (sum(syn_spot) > 25000) {
      syn_spot_sparse <- DropletUtils::downsampleMatrix(Matrix::Matrix(syn_spot, sparse = T),
                                                        prop = 20000 / sum(syn_spot))
    } else {
      syn_spot_sparse <- Matrix::Matrix(syn_spot, sparse = T)
    }
    
    rownames(syn_spot_sparse) <- names_genes
    colnames(syn_spot_sparse) <- name_simp
    
    # update progress bar
    setTxtProgressBar(pb, i)
    
    return(list(syn_spot_sparse, spot_ds))
  })
  
  ds_syn_spots <- purrr::map(ds_spots, 1) %>%
    base::Reduce(function(m1, m2) cbind(unlist(m1), unlist(m2)), .)
  
  # Generate dataframe of spot characteristic
  ds_spots_metadata <- purrr::map(ds_spots, 2) %>%
    dplyr::bind_rows() %>%
    data.frame()
  
  ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
  
  # change column order so that its progressive
  lev_mod <- gsub("[\\+|\\ |\\/]", ".", unique(se_obj@meta.data[, clust_vr]))
  colnames(ds_spots_metadata) <- gsub("[\\+|\\ |\\/]", ".", colnames(ds_spots_metadata))
  # all_cn <- c(paste("clust_", lev_mod, sep = ""), "name") # This was set to deal when cluster names are numeric
  
  # Check if there are missing columns (Cell types not selected) and add them as all 0s
  if (sum(lev_mod %in% colnames(ds_spots_metadata)) == (length(unique(se_obj@meta.data[, clust_vr])) + 1)) {
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  } else {
    
    missing_cols <- lev_mod[which(!lev_mod %in% colnames(ds_spots_metadata))]
    ds_spots_metadata[missing_cols] <- 0
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  }
  
  # Close progress bar
  close(pb)
  
  print(sprintf("Generation of %s test spots took %s mins", n,
                round(difftime(Sys.time(), start_gen, units = "mins"), 2)))
  print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}


##---------------------simulate spots through scRNAseq---------------------------
##

simulateSpotsMy <- function(scRNAseq,scRNA.meta,meta.var = "cellType",cellNum = 1:20,
                          n = 1000,
                          verbose = TRUE) {
  
  # Check variables
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
    if (sum(syn_spot) > 25000) {
      syn_spot_sparse <- DropletUtils::downsampleMatrix(Matrix::Matrix(syn_spot, sparse = T),
                                                        prop = 20000 / sum(syn_spot))
    } else {
      syn_spot_sparse <- Matrix::Matrix(syn_spot, sparse = T)
    }
    
    rownames(syn_spot_sparse) <- names_genes
    colnames(syn_spot_sparse) <- name_simp
    
    # update progress bar
    setTxtProgressBar(pb, i)
    
    return(list(syn_spot_sparse, spot_ds))
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
  
  # Close progress bar
  close(pb)
  
  print(sprintf("Generation of %s test spots took %s mins", n,
                round(difftime(Sys.time(), start_gen, units = "mins"), 2)))
  print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}

compute_kl_matrix <- function(X){
  library(entropy)
  num_cols <- ncol(X)
  kl_matrix <- matrix(0, nrow = num_cols, ncol = num_cols)
  for (i in 1:num_cols){
    for (j in 1:num_cols){
    P <- X[, i] + 0.01
    Q <- X[, j] + 0.01
    kl_div <- KL.plugin(P, Q)
    kl_matrix[i, j] <- kl_div
    }
  }
  return(kl_matrix)
}
compute_dis_one_vs_all <- function(loc,i){
  x_a <- loc$x[i]
  y_a <- loc$y[i]
  dis <- matrix(0, nrow = 1, ncol = nrow(loc))
  
  for (j in 1:nrow(loc)) {
    x_b <- loc$x[j]
    y_b <- loc$y[j]
    dis_ab = sqrt((x_a - x_b)^2 + (y_a - y_b)^2)
    dis[,j] <- dis_ab
    
  }
  return(dis)
}
compute_loc_matrix <- function(loc){
  ##compute relative distance between sample1 and sample2
  n_samples <- nrow(loc)
  dis_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  for (i in 1:n_samples) {
    for (j in 1:n_samples) {
      ##
      sample_a_vs_all = compute_dis_one_vs_all(loc,i)
      sample_b_vs_all = compute_dis_one_vs_all(loc,j)
      dis_matrix[i,j] = sum((sample_a_vs_all - sample_b_vs_all)^2)
    }
    
  }
  
  return(dis_matrix)
}
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
compute_similarity_PASTE <- function(kl,dis,alpha = 0.1){
  kl <- min_max_norm(kl)
  dis <- min_max_norm(dis)
  kl <- 1 / (1 + kl)
  dis <- 1 / (1 + dis)
  result_matrix <- (1 - alpha) * kl + alpha * dis

  return(result_matrix)
  
  
}


get_ss <- function(x = NULL, coor = NULL,alpha = 0.1){
  if(is.null(x)){
    kl <- matrix(0, nrow = nrows(coor), ncol = nrows(coor))
  }else if(is.null(coor)){
    dis <- matrix(0, nrow = ncols(x), ncol = ncols(x))
  }else{
    kl <- compute_kl_matrix(x)
    dis <- compute_loc_matrix(coor)
  }
  ss <- compute_similarity_PASTE(kl, dis, alpha)
  colnames(ss) <- colnames(x)
  rownames(ss) <- colnames(x)
  return(ss)
}


getCs <- function(comData, Spatial.Deconv.Data){
  traindata = comData$Train$data
  colnames(traindata) <- comData$Train$pData$cellType
  c_ref <- Get_C(traindata,norm = "cpm")[["C_ref"]]
  testdata = comData$Test$data
  if("cellType" %in% colnames(comData$Test$pData)){
    ###if the ground-truth cell type annotations of test data is known
    colnames(testdata) <- comData$Test$pData$cellType
    c <- Get_C(testdata,norm = "cpm")[["C_ref"]]
    Spatial.Deconv.Data$C <- as.matrix(c)[Spatial.Deconv.Data$markerGenes,]
    
  }else{
    Spatial.Deconv.Data$C <- NULL
  }
  Spatial.Deconv.Data$C_ref <- as.matrix(c_ref)[Spatial.Deconv.Data$markerGenes,]
  Spatial.Deconv.Data$C_ref.cpm.marker <- Spatial.Deconv.Data$C_ref
  Spatial.Deconv.Data$C.cpm.marker <- Spatial.Deconv.Data$C
  Spatial.Deconv.Data$sc <- list(data = comData$Train$data[Spatial.Deconv.Data$markerGenes,], meta = comData$Train$pData)
  Spatial.Deconv.Data$C_ref.count.marker <- as.matrix(comData$C_ref[unique(Spatial.Deconv.Data$markerslist$gene),])
  Spatial.Deconv.Data$C.count.marker <- as.matrix(comData$C[unique(Spatial.Deconv.Data$markerslist$gene),])
  return(Spatial.Deconv.Data)
}


getLocations <- function(Spatial.Deconv.Data){
  Spatial.Deconv.Data$location.real <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_test),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_test),split="_"),"[",3)))
  if("X" %in% Spatial.Deconv.Data$comData$pData){
    Spatial.Deconv.Data$location.pseudo <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_train),split="_"),"[",2)), y = as.numeric(sapply(strsplit(colnames(Spatial.Deconv.Data$comData$T_train),split="_"),"[",3)))
    rownames(Spatial.Deconv.Data$location.pseudo) <- colnames(Spatial.Deconv.Data$comData$T_train)
  }else{
    Spatial.Deconv.Data$location.pseudo <- NULL
  }
  rownames(Spatial.Deconv.Data$location.real) <- colnames(Spatial.Deconv.Data$comData$T_test)
  Spatial.Deconv.Data$location.real.original <- Spatial.Deconv.Data$location.real
  Spatial.Deconv.Data$location.pseudo.original <- Spatial.Deconv.Data$location.pseudo
  return(Spatial.Deconv.Data)
}

  getCountCPMMarkerData <- function(Spatial.Deconv.Data) {
  Spatial.Deconv.Data$T_combine <- Spatial.Deconv.Data$BatchCorr
  Spatial.Deconv.Data$spatial.real.count.marker <- Spatial.Deconv.Data$comData$T_test[Spatial.Deconv.Data$markerGenes,]###矫正之前的counts数据
  Spatial.Deconv.Data$spatial.pseudo.count.marker <- Spatial.Deconv.Data$comData$T_train[Spatial.Deconv.Data$markerGenes,]###矫正之前的counts数据
  Spatial.Deconv.Data$spatial.real.cpm.marker <- calculateCPM(Spatial.Deconv.Data$comData$T_test)[Spatial.Deconv.Data$markerGenes,]###先标准化再截取
  Spatial.Deconv.Data$spatial.pseudo.cpm.marker <- calculateCPM(Spatial.Deconv.Data$comData$T_train)[Spatial.Deconv.Data$markerGenes,]###先标准化再截取
  return(Spatial.Deconv.Data)
  
}
