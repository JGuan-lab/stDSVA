


draw.pie <- function(proportion, spatial_location, colors = NULL, radius = NULL, seed = NULL, methodTitle = NULL, font_size= 15,ratio_rate=1) {
  ##proportion : celltypes * spots
  ##spatial_location: spots * (x,y)
  library(gtools)
  library(ggplot2)
  library(scatterpie)
  P = as.data.frame(proportion)
  P = P[, mixedsort(colnames(P))]
  P = P[mixedsort(rownames(P)),]##mix the celltypes and the samples
  P <- t(P)
  P <- as.data.frame(P)
  location = as.data.frame(spatial_location)
  if(all(c("x","y") %in% colnames(location))) {
    colnames(location) <- c("X","Y")
  } else{
    if(!all(c("X","Y") %in% colnames(location))) {
      stop("Location's form not true!")
    }
  }
  if (length(rownames(P)) != length(rownames(location))) {
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorDict <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
              "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", 
              "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", 
              "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
              "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#ffffb3", 
              "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", 
              "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
              "#8c6bb1", "#ff9da7", "#66c2a5", "#fc8d62", "#8da0cb", 
              "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
              "#b2df8a", "#cab2d6", "#6a3d9a", "#ffff99", "#33a02c",
              "#fb9a99")
  
  if (is.null(colors)) {
    if (ncol(P) > length(colorDict)) {
      colors = colorRampPalette(colorDict)(ncol(P))
    }
    else {
      if (is.null(seed)) {
        iseed = 44
      }
      else {
        iseed = seed
      }
      set.seed(iseed)
      colors = colorDict[sample(1:length(colorDict), 
                                     ncol(P))]
    }
  }
  else {
    colors = colors
  }
  
  
  location <- location[rownames(P),]
  data = cbind(P, location)
  ct.select = colnames(P)
  if (is.null(radius)) {
    radius = (max(data$X) - min(data$X)) * (max(data$Y) - 
                                              min(data$Y))
    radius = radius/nrow(data)
    radius = radius/pi
    radius = sqrt(radius) * 0.85
  }
  else {
    radius = radius
  }
  p = suppressMessages(ggplot() + geom_scatterpie(aes(x = X, 
                                                      y = Y, r = radius), data = data, cols = ct.select, color = NA) + 
                         coord_fixed(ratio = ratio_rate * max(data$X)/max(data$Y)) + scale_fill_manual(values = colors) + 
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), 
                               panel.background = element_blank(), plot.background = element_blank(), 
                               panel.border = element_rect(colour = "grey89", fill = NA, 
                                                           size = 0.5), axis.text = element_blank(), axis.ticks = element_blank(), 
                               axis.title = element_blank(), legend.title = element_text(size = 16, 
                                                                                         face = "bold"), legend.text = element_text(size = font_size), 
                               legend.key = element_rect(colour = "transparent", 
                                                         fill = "white"), legend.key.size = unit(0.45, 
                                                                                                 "cm"), strip.text = element_text(size = 16, 
                                                                                                                                  face = "bold"),
                               plot.title = element_text(hjust = 0.5,  
                                                         vjust = 1,     
                                                         size = 15,   
                                                         face = "bold", 
                                                         color = "black"),legend.position = "bottom") + 
                         guides(fill = guide_legend(title = "Cell Type")) + ggtitle(methodTitle))
  return(p)
}




draw.pie.oap <- function(proportion, spatial_location, colors = NULL, radius = NULL, seed = NULL, methodTitle = NULL, fileName = "H:/LinYifan/OAP/new.pptx", font_size = 16) {
  ##proportion : celltypes * spots
  ##spatial_location: spots * (x,y)
  library(gtools)
  library(ggplot2)
  library(scatterpie)
  P = as.data.frame(proportion)
  P = P[, mixedsort(colnames(P))]
  P = P[mixedsort(rownames(P)),]##mix the celltypes and the samples
  P <- t(P)
  P <- as.data.frame(P)
  location = as.data.frame(spatial_location)
  if(all(c("x","y") %in% colnames(location))) {
    colnames(location) <- c("X","Y")
  } else{
    if(!all(c("X","Y") %in% colnames(location))) {
      stop("Location's form not true!")
    }
  }
  if (length(rownames(P)) != length(rownames(location))) {
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorDict <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                 "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                 "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", 
                 "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", 
                 "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#ffffb3", 
                 "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", 
                 "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
                 "#8c6bb1", "#ff9da7", "#66c2a5", "#fc8d62", "#8da0cb", 
                 "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
                 "#b2df8a", "#cab2d6", "#6a3d9a", "#ffff99", "#33a02c",
                 "#fb9a99")
  
  if (is.null(colors)) {
    if (ncol(P) > length(colorDict)) {
      colors = colorRampPalette(colorDict)(ncol(P))
    }
    else {
      if (is.null(seed)) {
        iseed = 44
      }
      else {
        iseed = seed
      }
      set.seed(iseed)
      colors = colorDict[sample(1:length(colorDict), 
                                ncol(P))]
    }
  }
  else {
    colors = colors
  }
  
  
  location <- location[rownames(P),]
  data = cbind(P, location)
  ct.select = colnames(P)
  if (is.null(radius)) {
    radius = (max(data$X) - min(data$X)) * (max(data$Y) - 
                                              min(data$Y))
    radius = radius/nrow(data)
    radius = radius/pi
    radius = sqrt(radius) * 0.85
  }
  else {
    radius = radius
  }
  p = suppressMessages(ggplot() + geom_scatterpie(aes(x = X, 
                                                      y = Y, r = radius), data = data, cols = ct.select, color = NA) + 
                         coord_fixed(ratio = 1 * max(data$X)/max(data$Y)) + scale_fill_manual(values = colors) + 
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), 
                               panel.background = element_blank(), plot.background = element_blank(), 
                               panel.border = element_rect(colour = "grey89", fill = NA, 
                                                           size = 0.5), axis.text = element_blank(), axis.ticks = element_blank(), 
                               axis.title = element_blank(), legend.title = element_text(size = 16, 
                                                                                         face = "bold"), legend.text = element_text(size = font_size), 
                               legend.key = element_rect(colour = "transparent", 
                                                         fill = "white"), legend.key.size = unit(0.45, 
                                                                                                 "cm"), strip.text = element_text(size = 16, 
                                                                                                                                  face = "bold"),
                               plot.title = element_text(hjust = 0.5,  
                                                         vjust = 1,     
                                                         size = 15,   
                                                         face = "bold", 
                                                         color = "black"),legend.position = "bottom") + 
                         guides(fill = guide_legend(title = "Cell Type")) + ggtitle(methodTitle))
  library(eoffice)
  topptx(p,fileName)
  return(p)
  
}

draw.pie.main.oap <- function(proportion, spatial_location, colors = NULL, radius = NULL, seed = NULL, methodTitle = NULL, draw.plot = TRUE, font_size = 16) {
  
  ##draw the max celltype plot
  ##proportion : celltypes * spots
  ##spatial_location: spots * (x,y)
  library(gtools)
  library(ggplot2)
  library(scatterpie)
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
    P[i,][P[i,] < max[i]] = 10e-8
    P[i,][P[i,] == max[i]] = 1
    
  }
  rownames(dominant) <- rownames(P)
  location = as.data.frame(spatial_location)
  if(all(c("x","y") %in% colnames(location))) {
    colnames(location) <- c("X","Y")
  } else{
    if(!all(c("X","Y") %in% colnames(location))) {
      stop("Location's form not true!")
    }
  }
  if (length(rownames(P)) != length(rownames(location))) {
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorDict <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                 "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                 "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", 
                 "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", 
                 "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#ffffb3", 
                 "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", 
                 "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
                 "#8c6bb1", "#ff9da7", "#66c2a5", "#fc8d62", "#8da0cb", 
                 "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
                 "#b2df8a", "#cab2d6", "#6a3d9a", "#ffff99", "#33a02c",
                 "#fb9a99")
  
  if (is.null(colors)) {
    if (ncol(P) > length(colorDict)) {
      colors = colorRampPalette(colorDict)(ncol(P))
    }
    else {
      if (is.null(seed)) {
        iseed = 44
      }
      else {
        iseed = seed
      }
      set.seed(iseed)
      colors = colorDict[sample(1:length(colorDict), 
                                ncol(P))]
    }
  }
  else {
    colors = colors
  }
  
  
  location <- location[rownames(P),]
  data = cbind(P, location)
  ct.select = colnames(P)
  if (is.null(radius)) {
    radius = (max(data$X) - min(data$X)) * (max(data$Y) - 
                                              min(data$Y))
    radius = radius/nrow(data)
    radius = radius/pi
    radius = sqrt(radius) * 0.85
  }
  else {
    radius = radius
  }
  
  
  if(draw.plot) {
    p = suppressMessages(ggplot() + geom_scatterpie(aes(x = X, 
                                                        y = Y, r = radius), data = data, cols = ct.select, color = NA) + 
                           coord_fixed(ratio = 1 * max(data$X)/max(data$Y)) + scale_fill_manual(values = colors) + 
                           theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), 
                                 panel.background = element_blank(), plot.background = element_blank(), 
                                 panel.border = element_rect(colour = "grey89", fill = NA, 
                                                             size = 0.5), axis.text = element_blank(), axis.ticks = element_blank(), 
                                 axis.title = element_blank(), legend.title = element_text(size = 16, 
                                                                                           face = "bold"), legend.text = element_text(size = font_size), 
                                 legend.key = element_rect(colour = "transparent", 
                                                           fill = "white"), legend.key.size = unit(0.45, 
                                                                                                   "cm"), strip.text = element_text(size = 16, 
                                                                                                                                    face = "bold"),
                                 plot.title = element_text(hjust = 0.5,  
                                                           vjust = 1,     
                                                           size = 15,   
                                                           face = "bold", 
                                                           color = "black"), legend.position = "bottom") + 
                           guides(fill = guide_legend(title = "Cell Type")) + ggtitle(methodTitle))
    print(p)
  }
  return(list(plot = p,clu = dominant))
  
}


draw.pie.main <- function(proportion, spatial_location, colors = NULL, radius = NULL, seed = NULL, methodTitle = NULL, draw.plot = TRUE,font_size = 16) {
  
  ##draw the max celltype plot
  ##proportion : celltypes * spots
  ##spatial_location: spots * (x,y)
  library(gtools)
  library(ggplot2)
  library(scatterpie)
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
    P[i,][P[i,] < max[i]] = 10e-8
    P[i,][P[i,] == max[i]] = 1
    
  }
  rownames(dominant) <- rownames(P)
  location = as.data.frame(spatial_location)
  if(all(c("x","y") %in% colnames(location))) {
    colnames(location) <- c("X","Y")
  } else{
    if(!all(c("X","Y") %in% colnames(location))) {
      stop("Location's form not true!")
    }
  }
  if (length(rownames(P)) != length(rownames(location))) {
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorDict <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                 "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                 "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", 
                 "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", 
                 "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#ffffb3", 
                 "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", 
                 "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
                 "#8c6bb1", "#ff9da7", "#66c2a5", "#fc8d62", "#8da0cb", 
                 "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
                 "#b2df8a", "#cab2d6", "#6a3d9a", "#ffff99", "#33a02c",
                 "#fb9a99")
  
  if (is.null(colors)) {
    if (ncol(P) > length(colorDict)) {
      colors = colorRampPalette(colorDict)(ncol(P))
    }
    else {
      if (is.null(seed)) {
        iseed = 44
      }
      else {
        iseed = seed
      }
      set.seed(iseed)
      colors = colorDict[sample(1:length(colorDict), 
                                ncol(P))]
    }
  }
  else {
    colors = colors
  }
  
  
  location <- location[rownames(P),]
  data = cbind(P, location)
  ct.select = colnames(P)
  if (is.null(radius)) {
    radius = (max(data$X) - min(data$X)) * (max(data$Y) - 
                                              min(data$Y))
    radius = radius/nrow(data)
    radius = radius/pi
    radius = sqrt(radius) * 0.85
  }
  else {
    radius = radius
  }
  
  
  if(draw.plot) {
      p = suppressMessages(ggplot() + geom_scatterpie(aes(x = X, 
                                                      y = Y, r = radius), data = data, cols = ct.select, color = NA) + 
                         coord_fixed(ratio = 1 * max(data$X)/max(data$Y)) + scale_fill_manual(values = colors) + 
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), 
                               panel.background = element_blank(), plot.background = element_blank(), 
                               panel.border = element_rect(colour = "grey89", fill = NA, 
                                                           size = 0.5), axis.text = element_blank(), axis.ticks = element_blank(), 
                               axis.title = element_blank(), legend.title = element_text(size = 16, 
                                                                                         face = "bold"), legend.text = element_text(size = font_size), 
                               legend.key = element_rect(colour = "transparent", 
                                                         fill = "white"), legend.key.size = unit(0.45, 
                                                                                                 "cm"), strip.text = element_text(size = 16, 
                                                                                                                                  face = "bold"),
                               plot.title = element_text(hjust = 0.5,  
                                                         vjust = 1,     
                                                         size = 15,   
                                                         face = "bold", 
                                                         color = "black"), legend.position = "bottom") + 
                         guides(fill = guide_legend(title = "Cell Type")) + ggtitle(methodTitle))
  print(p)
  }

  return(dominant)
}


draw.pie.cluster <- function(proportion, location, colors = NULL,clusterName = NULL, k = 4, methodTitle = NULL, draw.plot = TRUE, seed = 44){
  library(factoextra)
  library(cluster)
  proportion <- NORM(proportion)
  if(length(rownames(location)) != length(colnames(proportion))) {
    stop("Error!location and proportion donot match!")
  }
  set.seed(44)
  P <- t(proportion)
  # kmeans
  n <- nrow(P)  
  location <- location[rownames(P),]
  
  km <- kmeans(P, centers = k)
  clu <- as.matrix(km$cluster)
  labels <- as.character(clu)
  spots <- P
  
  # if(is.null(colors)) {
  #   
  #   colors = c("#8dd3c7", "#ffffb3", 
  #              "#bebada", "#fb8072")
  # }
  if(is.null(clusterName)) {
    clusterName <- paste("cluster",unique(labels),sep = "")
  }
  binary_matrix <- generate_binary_matrix(labels, spots)
  binary_matrix <- t(binary_matrix)
  rownames(binary_matrix) <- clusterName
  binary_matrix <- binary_matrix[,rownames(location)]
  if(draw.plot) {
    plot <- draw.pie(
      proportion = binary_matrix,
      spatial_location = (location),
      colors = colors,
      methodTitle = methodTitle, seed = seed)
    print(plot)
  }


  return(clu)
}

draw.pie.cluster.oap <- function(proportion, location, colors = NULL,clusterName = NULL, k = 4, methodTitle = NULL, draw.plot = TRUE, seed = 44, font_size = 15){
  library(factoextra)
  library(cluster)
  proportion <- NORM(proportion)
  if(length(rownames(location)) != length(colnames(proportion))) {
    stop("Error!location and proportion donot match!")
  }
  set.seed(44)
  P <- t(proportion)
  # kmeans
  n <- nrow(P)  
  location <- location[rownames(P),]
  
  km <- kmeans(P, centers = k)
  clu <- as.matrix(km$cluster)
  labels <- as.character(clu)
  spots <- P
  
  # if(is.null(colors)) {
  #   
  #   colors = c("#8dd3c7", "#ffffb3", 
  #              "#bebada", "#fb8072")
  # }
  if(is.null(clusterName)) {
    clusterName <- paste("cluster",unique(labels),sep = "")
  }
  binary_matrix <- generate_binary_matrix(labels, spots)
  binary_matrix <- t(binary_matrix)
  rownames(binary_matrix) <- clusterName
  binary_matrix <- binary_matrix[,rownames(location)]
  if(draw.plot) {
    plot <- draw.pie(
      proportion = binary_matrix,
      spatial_location = (location),
      colors = colors,
      methodTitle = methodTitle, seed = seed, font_size = font_size)
    print(plot)
  }
  
  
  return(list(plot = plot,clu = clu))
}
generate_binary_matrix <- function(labels, spots) {
  unique_labels <- unique(labels)
  binary_matrix <- sapply(unique_labels, function(l) {
    ifelse(labels == l, 1, 10e-8)
  })
  colnames(binary_matrix) <- unique_labels
  rownames(binary_matrix) <- rownames(spots)
  return(binary_matrix)
}


draw_plot <- function(data){
  library(ggplot2)
  if(all(c("x","y") %in% colnames(data))){
    data$X <- data$x
    data$Y <- data$Y
    data <- data[,!colnames(data) %in% c("x","y")]
    
  } else if(all(c("X","Y") %in% colnames(data))){
  ggplot(data, aes(x = X, y = Y, color = Layer)) +
    geom_point(size = 3.5) +
    theme_minimal()
  } else {
    stop("missing locations!")
  }
}



draw.pie.all <- function(Spatial.Deconv.Data,Methods = c("CARD","RCTD","STdeconvolve")) {
  if(is.null(Spatial.Deconv.Data$RESULTSMap) || is.null(Methods)){
    stop("error")
  }
  aripurity.dominant <- NULL
  aripurity.kmeans <- NULL
  
  if("CARD" %in% Methods) {
    
    print(draw.pie(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD"))
    aripurity.dominant$CARD <- getARIPurity.dominant(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    aripurity.kmeans$CARD <- getARIPurity.kmeans(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    
  }
  if("RCTD" %in% Methods) {
    print(draw.pie(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD"))
    aripurity.dominant$RCTD <- getARIPurity.dominant(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    aripurity.kmeans$RCTD <- getARIPurity.kmeans(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    
  }
  
  if("STdeconvolve" %in% Methods) {
    print(draw.pie(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve"))
    aripurity.dominant$STdeconvolve <- getARIPurity.dominant(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    aripurity.kmeans$STdeconvolve <- getARIPurity.kmeans(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    
  }
  if("RNASieve" %in% Methods) {
    print(draw.pie(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve"))
    aripurity.dominant$RNASieve <- getARIPurity.dominant(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    aripurity.kmeans$RNASieve <- getARIPurity.kmeans(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    
  }
  
  if("Cell2location" %in% Methods) {
    print(draw.pie(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location"))
    aripurity.dominant$Cell2location <- getARIPurity.dominant(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    aripurity.kmeans$Cell2location <- getARIPurity.kmeans(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData)
    
  }  
  Spatial.Deconv.Data$RESULTSMap$aripurity.dominant <- aripurity.dominant
  Spatial.Deconv.Data$RESULTSMap$aripurity.kmeans <- aripurity.kmeans
  Spatial.Deconv.Data <- combineResults.aripurity(Spatial.Deconv.Data)
  return(Spatial.Deconv.Data)

}




draw.pie.all.oap <- function(Spatial.Deconv.Data,Methods = c("CARD","RCTD","STdeconvolve"),dataName = "mob_loc",seed1 = 44, seed2 = 44, seed3 = 44, type = "all", colors = NULL,font_size = 15) {
  if(is.null(Spatial.Deconv.Data$RESULTSMap) || is.null(Methods)){
    if(is.null(realPDAC_noloc_2$runSpatialDSSCBEResult)){
      stop("error")
    }
    
  }
  aripurity.dominant <- NULL
  aripurity.kmeans <- NULL
  if(type == "all"){
    if("My Method" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data[["runSpatialDSSCBEResult"]][["result_para_p_real"]],Spatial.Deconv.Data$location.real,methodTitle = "Method",fileName = paste("H:/LinYifan/OAP/", dataName ,"_MyMethod.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$MyMethod <- getARIPurity.dominant.oap(Spatial.Deconv.Data[["runSpatialDSSCBEResult"]][["result_para_p_real"]],Spatial.Deconv.Data$location.real,methodTitle = "My Method",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_MyMethod_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$MyMethod <- getARIPurity.kmeans.oap(Spatial.Deconv.Data[["runSpatialDSSCBEResult"]][["result_para_p_real"]],Spatial.Deconv.Data$location.real,methodTitle = "My Method",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_MyMethod_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("CARD" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD",fileName = paste("H:/LinYifan/OAP/", dataName ,"_CARD.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$CARD <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_CARD_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$CARD <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_CARD_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    
    if("TOAST/P-" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS[['TOAST/P-']]$result_p,Spatial.Deconv.Data$location.real,methodTitle = "TOAST/P-",fileName = paste("H:/LinYifan/OAP/", dataName ,"_TOASTP-.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant[['TOAST/P-']] <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS[['TOAST/P-']]$result_p,Spatial.Deconv.Data$location.real,methodTitle = "TOAST/P-",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_TOASTP-_dominant.pptx",sep = ""),seed = seed2,font_size = font_size)
      aripurity.kmeans[['TOAST/P-']] <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS[['TOAST/P-']]$result_p,Spatial.Deconv.Data$location.real,methodTitle = "TOAST/P-",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_TOASTP-_kmeans.pptx",sep = ""),seed = seed3,font_size = font_size)
      
    }
    
    if("Redeconve" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Redeconve,Spatial.Deconv.Data$location.real,methodTitle = "Redeconve",fileName = paste("H:/LinYifan/OAP/", dataName ,"_Redeconve.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$Redeconve <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Redeconve,Spatial.Deconv.Data$location.real,methodTitle = "Redeconve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Redeconve_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$Redeconve <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Redeconve,Spatial.Deconv.Data$location.real,methodTitle = "Redeconve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Redeconve_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    if("SPOTlight" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$SPOTlight,Spatial.Deconv.Data$location.real,methodTitle = "SPOTlight",fileName = paste("H:/LinYifan/OAP/", dataName ,"_SPOTlight.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$SPOTlight <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$SPOTlight,Spatial.Deconv.Data$location.real,methodTitle = "SPOTlight",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_SPOTlight_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$SPOTlight <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$SPOTlight,Spatial.Deconv.Data$location.real,methodTitle = "SPOTlight",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_SPOTlight_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("Linseed" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Linseed$p,Spatial.Deconv.Data$location.real,methodTitle = "Linseed",fileName = paste("H:/LinYifan/OAP/", dataName ,"_Linseed.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$Linseed <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Linseed$p,Spatial.Deconv.Data$location.real,methodTitle = "Linseed",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Linseed_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$Linseed <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Linseed$p,Spatial.Deconv.Data$location.real,methodTitle = "Linseed",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Linseed_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("DSA" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$DSA$p,Spatial.Deconv.Data$location.real,methodTitle = "DSA",fileName = paste("H:/LinYifan/OAP/", dataName ,"_DSA.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$DSA <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$DSA$p,Spatial.Deconv.Data$location.real,methodTitle = "DSA",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_DSA_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$DSA <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$DSA$p,Spatial.Deconv.Data$location.real,methodTitle = "DSA",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_DSA_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("ssKL" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssKL$p,Spatial.Deconv.Data$location.real,methodTitle = "ssKL",fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssKL.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$ssKL <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssKL$p,Spatial.Deconv.Data$location.real,methodTitle = "ssKL",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssKL_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$ssKL <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssKL$p,Spatial.Deconv.Data$location.real,methodTitle = "ssKL",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssKL_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    
    if("ssFrobenius" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssFrobenius$p,Spatial.Deconv.Data$location.real,methodTitle = "ssFrobenius",fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssFrobenius.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$ssFrobenius <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssFrobenius$p,Spatial.Deconv.Data$location.real,methodTitle = "ssFrobenius",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssFrobenius_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$ssFrobenius <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssFrobenius$p,Spatial.Deconv.Data$location.real,methodTitle = "ssFrobenius",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssFrobenius_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    } 
    
    if("Deconf" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Deconf$p,Spatial.Deconv.Data$location.real,methodTitle = "Deconf",fileName = paste("H:/LinYifan/OAP/", dataName ,"_Deconf.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$Deconf <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Deconf$p,Spatial.Deconv.Data$location.real,methodTitle = "Deconf",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Deconf_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$Deconf <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Deconf$p,Spatial.Deconv.Data$location.real,methodTitle = "Deconf",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Deconf_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    } 
    
    if("RCTD" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD",fileName = paste("H:/LinYifan/OAP/", dataName ,"_RCTD.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$RCTD <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_RCTD_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$RCTD <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_RCTD_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("STdeconvolve" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve",fileName = paste("H:/LinYifan/OAP/", dataName ,"_STdeconvolve.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$STdeconvolve <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_STdeconvolve_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$STdeconvolve <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_STdeconvolve_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    if("RNASieve" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve",fileName = paste("H:/LinYifan/OAP/", dataName ,"_RNASieve.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$RNASieve <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_RNASieve_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$RNASieve <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_RNASieve_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("Cell2location" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location",fileName = paste("H:/LinYifan/OAP/", dataName ,"_Cell2location.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$Cell2location <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Cell2location_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$Cell2location <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Cell2location_kmeans.pptx",sep = ""),,seed = seed3, colors = colors,font_size = font_size)
      
    }  
    if("spatialDWLS" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$spatialDWLS,Spatial.Deconv.Data$location.real,methodTitle = "spatialDWLS",fileName = paste("H:/LinYifan/OAP/", dataName ,"_spatialDWLS.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
      aripurity.dominant$spatialDWLS <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$spatialDWLS,Spatial.Deconv.Data$location.real,methodTitle = "spatialDWLS",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_spatialDWLS_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)
      aripurity.kmeans$spatialDWLS <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$spatialDWLS,Spatial.Deconv.Data$location.real,methodTitle = "spatialDWLS",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_spatialDWLS_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
  }
  
  if(type == "kmeans"){
    if("My Method" %in% Methods) {
      
      aripurity.kmeans$MyMethod <- getARIPurity.kmeans.oap(Spatial.Deconv.Data[["runSpatialDSSCBEResult"]][["result_para_p_real"]],Spatial.Deconv.Data$location.real,methodTitle = "My Method",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_MyMethod_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("CARD" %in% Methods) {
      
     
      aripurity.kmeans$CARD <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_CARD_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    
    if("TOAST/P-" %in% Methods) {
      
      aripurity.kmeans[['TOAST/P-']] <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS[['TOAST/P-']]$result_p,Spatial.Deconv.Data$location.real,methodTitle = "TOAST/P-",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_TOASTP-_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("Redeconve" %in% Methods) {
     
      aripurity.kmeans$Redeconve <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Redeconve,Spatial.Deconv.Data$location.real,methodTitle = "Redeconve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Redeconve_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    if("SPOTlight" %in% Methods) {
      
      aripurity.kmeans$SPOTlight <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$SPOTlight,Spatial.Deconv.Data$location.real,methodTitle = "SPOTlight",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_SPOTlight_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("Linseed" %in% Methods) {
      
      aripurity.kmeans$Linseed <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Linseed$p,Spatial.Deconv.Data$location.real,methodTitle = "Linseed",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Linseed_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("DSA" %in% Methods) {
      
      aripurity.kmeans$DSA <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$DSA$p,Spatial.Deconv.Data$location.real,methodTitle = "DSA",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_DSA_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("ssKL" %in% Methods) {
      
      aripurity.kmeans$ssKL <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssKL$p,Spatial.Deconv.Data$location.real,methodTitle = "ssKL",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssKL_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    
    if("ssFrobenius" %in% Methods) {
      
      aripurity.kmeans$ssFrobenius <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssFrobenius$p,Spatial.Deconv.Data$location.real,methodTitle = "ssFrobenius",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssFrobenius_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    } 
    
    if("Deconf" %in% Methods) {

      aripurity.kmeans$Deconf <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Deconf$p,Spatial.Deconv.Data$location.real,methodTitle = "Deconf",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Deconf_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    } 
    
    if("RCTD" %in% Methods) {
      
      aripurity.kmeans$RCTD <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_RCTD_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("STdeconvolve" %in% Methods) {
      
      aripurity.kmeans$STdeconvolve <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_STdeconvolve_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    if("RNASieve" %in% Methods) {
      
      aripurity.kmeans$RNASieve <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_RNASieve_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }
    
    if("Cell2location" %in% Methods) {
      
      aripurity.kmeans$Cell2location <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Cell2location_kmeans.pptx",sep = ""),seed = seed3, colors = colors,font_size = font_size)
      
    }  
    if("spatialDWLS" %in% Methods) {
      
      aripurity.kmeans$spatialDWLS <- getARIPurity.kmeans.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$spatialDWLS,Spatial.Deconv.Data$location.real,methodTitle = "spatialDWLS",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_spatialDWLS_kmeans.pptx",sep = ""),,seed = seed3, colors = colors,font_size = font_size)
      
    }
  }
  
  if(type == "dominant"){
    if("My Method" %in% Methods) {
      aripurity.dominant$MyMethod <- getARIPurity.dominant.oap(Spatial.Deconv.Data[["runSpatialDSSCBEResult"]][["result_para_p_real"]],Spatial.Deconv.Data$location.real,methodTitle = "My Method",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_MyMethod_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    if("CARD" %in% Methods) {
      
      aripurity.dominant$CARD <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_CARD_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    
    if("TOAST/P-" %in% Methods) {
      
      aripurity.dominant[['TOAST/P-']] <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS[['TOAST/P-']]$result_p,Spatial.Deconv.Data$location.real,methodTitle = "TOAST/P-",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_TOASTP-_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    if("Redeconve" %in% Methods) {
      
      aripurity.dominant$Redeconve <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Redeconve,Spatial.Deconv.Data$location.real,methodTitle = "Redeconve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Redeconve_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    if("SPOTlight" %in% Methods) {
      
      aripurity.dominant$SPOTlight <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$SPOTlight,Spatial.Deconv.Data$location.real,methodTitle = "SPOTlight",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_SPOTlight_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    if("Linseed" %in% Methods) {
      
      aripurity.dominant$Linseed <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Linseed$p,Spatial.Deconv.Data$location.real,methodTitle = "Linseed",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Linseed_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    if("DSA" %in% Methods) {
      
      aripurity.dominant$DSA <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$DSA$p,Spatial.Deconv.Data$location.real,methodTitle = "DSA",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_DSA_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    if("ssKL" %in% Methods) {
      
      aripurity.dominant$ssKL <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssKL$p,Spatial.Deconv.Data$location.real,methodTitle = "ssKL",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssKL_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    
    if("ssFrobenius" %in% Methods) {
      
      aripurity.dominant$ssFrobenius <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssFrobenius$p,Spatial.Deconv.Data$location.real,methodTitle = "ssFrobenius",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssFrobenius_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    } 
    
    if("Deconf" %in% Methods) {
      
      aripurity.dominant$Deconf <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Deconf$p,Spatial.Deconv.Data$location.real,methodTitle = "Deconf",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Deconf_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    } 
    
    if("RCTD" %in% Methods) {
      aripurity.dominant$RCTD <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_RCTD_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    if("STdeconvolve" %in% Methods) {
      aripurity.dominant$STdeconvolve <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_STdeconvolve_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    if("RNASieve" %in% Methods) {
      aripurity.dominant$RNASieve <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_RNASieve_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
    
    if("Cell2location" %in% Methods) {
      aripurity.dominant$Cell2location <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_Cell2location_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }  
    if("spatialDWLS" %in% Methods) {
      aripurity.dominant$spatialDWLS <- getARIPurity.dominant.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$spatialDWLS,Spatial.Deconv.Data$location.real,methodTitle = "spatialDWLS",ground_truth = Spatial.Deconv.Data$comData$Test$full_phenoData,fileName = paste("H:/LinYifan/OAP/", dataName ,"_spatialDWLS_dominant.pptx",sep = ""),seed = seed2, colors = colors,font_size = font_size)

    }
  }
  if(type == "proportions"){
    if("My Method" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data[["runSpatialDSSCBEResult"]][["result_para_p_real"]],Spatial.Deconv.Data$location.real,methodTitle = "Method",fileName = paste("H:/LinYifan/OAP/", dataName ,"_MyMethod.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
   
    }
    
    if("CARD" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$CARD,Spatial.Deconv.Data$location.real,methodTitle = "CARD",fileName = paste("H:/LinYifan/OAP/", dataName ,"_CARD.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
 
    }
    
    
    if("TOAST/P-" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS[['TOAST/P-']]$result_p,Spatial.Deconv.Data$location.real,methodTitle = "TOAST/P-",fileName = paste("H:/LinYifan/OAP/", dataName ,"_TOASTP-.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
 
    }
    
    if("Redeconve" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Redeconve,Spatial.Deconv.Data$location.real,methodTitle = "Redeconve",fileName = paste("H:/LinYifan/OAP/", dataName ,"_Redeconve.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))

    }
    if("SPOTlight" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$SPOTlight,Spatial.Deconv.Data$location.real,methodTitle = "SPOTlight",fileName = paste("H:/LinYifan/OAP/", dataName ,"_SPOTlight.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))

    }
    
    if("Linseed" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Linseed$p,Spatial.Deconv.Data$location.real,methodTitle = "Linseed",fileName = paste("H:/LinYifan/OAP/", dataName ,"_Linseed.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
  
    }
    
    if("DSA" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$DSA$p,Spatial.Deconv.Data$location.real,methodTitle = "DSA",fileName = paste("H:/LinYifan/OAP/", dataName ,"_DSA.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
  
    }
    
    if("ssKL" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssKL$p,Spatial.Deconv.Data$location.real,methodTitle = "ssKL",fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssKL.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
 
    }
    
    
    if("ssFrobenius" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$ssFrobenius$p,Spatial.Deconv.Data$location.real,methodTitle = "ssFrobenius",fileName = paste("H:/LinYifan/OAP/", dataName ,"_ssFrobenius.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
    
    } 
    
    if("Deconf" %in% Methods) {
      
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Deconf$p,Spatial.Deconv.Data$location.real,methodTitle = "Deconf",fileName = paste("H:/LinYifan/OAP/", dataName ,"_Deconf.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
   
    } 
    
    if("RCTD" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RCTD,Spatial.Deconv.Data$location.real,methodTitle = "RCTD",fileName = paste("H:/LinYifan/OAP/", dataName ,"_RCTD.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
  
    }
    
    if("STdeconvolve" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$STdeconvolve$result_p,Spatial.Deconv.Data$location.real,methodTitle = "STdeconvolve",fileName = paste("H:/LinYifan/OAP/", dataName ,"_STdeconvolve.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
   
    }
    if("RNASieve" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$RNASieve[[1]],Spatial.Deconv.Data$location.real,methodTitle = "RNASieve",fileName = paste("H:/LinYifan/OAP/", dataName ,"_RNASieve.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
    
    }
    
    if("Cell2location" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$Cell2location,Spatial.Deconv.Data$location.real,methodTitle = "Cell2location",fileName = paste("H:/LinYifan/OAP/", dataName ,"_Cell2location.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
 
    }  
    if("spatialDWLS" %in% Methods) {
      print(draw.pie.oap(Spatial.Deconv.Data$RESULTSMap$RESULTS$spatialDWLS,Spatial.Deconv.Data$location.real,methodTitle = "spatialDWLS",fileName = paste("H:/LinYifan/OAP/", dataName ,"_spatialDWLS.pptx",sep = ""),seed = seed1, colors = colors,font_size = font_size))
  
    }
  }

  if(is.null(Spatial.Deconv.Data$RESULTSMap$aripurity.dominant) && is.null(Spatial.Deconv.Data$RESULTSMap$aripurity.kmeans)){
    Spatial.Deconv.Data$RESULTSMap$aripurity.dominant <- aripurity.dominant
    Spatial.Deconv.Data$RESULTSMap$aripurity.kmeans <- aripurity.kmeans
  } else {
    Spatial.Deconv.Data$RESULTSMap$aripurity.dominant <- rbind(Spatial.Deconv.Data$RESULTSMap$aripurity.dominant,aripurity.dominant)
    Spatial.Deconv.Data$RESULTSMap$aripurity.kmeans <- rbind(Spatial.Deconv.Data$RESULTSMap$aripurity.kmeans,aripurity.kmeans)
    
  }
  Spatial.Deconv.Data <- combineResults.aripurity(Spatial.Deconv.Data)
  return(Spatial.Deconv.Data)
  
}



#' 绘制多方法反卷积结果的空间饼图拼图
#'
#' @param RESULTS 包含多种方法反卷积结果的列表，每个元素应为细胞类型比例矩阵（celltypes × spots）
#' @param location 空间坐标数据框，需包含x/y列（spots × (x,y)）
#' @param colors 可选的颜色向量，长度需≥细胞类型数（默认使用draw.pie内置调色板）
#' @param radius 饼图半径（默认NULL自动计算）
#' @param ncol 拼图的列数（默认根据方法数自动计算）
#' @param titles 各子图标题（默认使用RESULTS的命名）
#' @param font_size 图例文字大小（默认15）
#' @param legend_position 图例位置（"bottom", "right", "none"等，默认"bottom"）
#' 
#' @return 返回拼合后的ggplot对象，可通过print()显示或ggsave()保存
#'
#' @examples
plot_multiple_decon_pies <- function(RESULTS, location,  layer = NULL,
                                     colors = NULL, radius = NULL, colors_layer = NULL,
                                     ncol = NULL, titles = NULL,
                                     font_size = 15, 
                                     legend_font_size = 12,
                                     legend_position = "bottom",
                                     show_tags = TRUE,
                                     reorder_methods = NULL,
                                     family = "Times New Roman",seed = NULL) {
  location <- location[mixedsort(rownames(location)),]
  # 检查依赖包
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("请先安装patchwork包：install.packages('patchwork')")
  }
  
  # 参数校验
  if (!is.list(RESULTS) || length(RESULTS) == 0) {
    stop("RESULTS必须是非空列表")
  }
  
  # 方法重排序（将"My Method"置顶）
  # 方法重命名为 stDSVA
  if ("My Method" %in% names(RESULTS)) {
    names(RESULTS)[names(RESULTS) == "My Method"] <- "stDSVA"
  }
  if ("Deconf" %in% names(RESULTS)) {
    names(RESULTS)[names(RESULTS) == "Deconf"] <- "deconf"
  }

      if('ssKL' %in% names(RESULTS)){
        ordered_names <- c('stDSVA',
                           'deconf',
                           'DSA',
                           'Linseed',
                           'RNASieve',
                           'ssFrobenius',
                           'ssKL',
                           'TOAST/-P',
                           'CARD',
                           'RCTD',
                           'spatialDWLS',
                           'Seurat',
                           'Tangram')
      }else {
        ordered_names <- c('stDSVA',
                           'deconf',
                           'DSA',
                           'Linseed',
                           'RNASieve',
                           'ssFrobenius',
                       
                           'TOAST/-P',
                           'CARD',
                           'RCTD',
                           'spatialDWLS',
                           'Seurat',
                           'Tangram')
      }

      RESULTS <- RESULTS[ordered_names]

  if (is.null(names(RESULTS))) {
    names(RESULTS) <- paste0("Method", seq_along(RESULTS))
  }
  if (is.null(titles)) {
    titles <- names(RESULTS)
  } else if (length(titles) != length(RESULTS)) {
    stop("titles长度必须与RESULTS相同")
  }
  
  # 自定义主题函数（统一字体和图例）
  custom_theme <- function(base_size = font_size, legend_size = legend_font_size) {
    ggplot2::theme(
      text = element_text(family = family),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = legend_size),
      plot.title = element_text(
        family = family,
        size = base_size + 1,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      strip.text = element_text(family = family, size = base_size, face = "bold")
    )
  }
  # 生成layer标签图（如果提供了layer）
  if (!is.null(layer)) {
    # 检查layer和location的匹配性
    common_spots <- intersect(rownames(location), names(layer))
    if (length(common_spots) == 0) stop("layer和location无匹配的spot名称")
    
    # 创建layer的"伪比例矩阵"（每个spot属于一个layer）
    layer_matrix <- model.matrix(~ 0 + factor(layer[common_spots]))
    colnames(layer_matrix) <- levels(factor(layer[common_spots]))
    rownames(layer_matrix) <- common_spots
    
    # 绘制layer饼图
    layer_pie <- draw.pie(
      proportion = t(layer_matrix),  # 转置以符合draw.pie的输入格式
      spatial_location = location[common_spots, , drop = FALSE],
      colors = colors_layer,
      radius = radius,
      methodTitle = "Layer Annotations",
      font_size = 6
    ) + 
      custom_theme() +
      guides(fill = guide_legend(
        title = NULL,
        title.position = "top",
        ncol = length(unique(layer)),
        keywidth = unit(0.4, "lines"),
        keyheight = unit(0.4, "lines")
        )) +
          theme(legend.position = "bottom", 
                lengend.text = element_text(size = 8),
                legend.spacing = unit(0.1, "cm"),
                legend.margin = margin(0, 0, 0, 0))
  }
  # 生成各方法的饼图
  ct_names <- rownames(RESULTS[[1]])
  ct_names_clean <- gsub("\\.", " ", ct_names)
  plot_list <- lapply(seq_along(RESULTS), function(i) {
    print(names(RESULTS))
    method_name <- names(RESULTS)[i]
    cat("Method:", method_name, "\n")  # 检查方法名
    prop_matrix <- RESULTS[[i]]
    rownames(prop_matrix) <- gsub("\\.", " ", rownames(prop_matrix))
    cat("prop_matrix colnames:", head(colnames(prop_matrix)), "\n")  # 检查列名
    
    cat("location rownames:", head(rownames(location)), "\n")  # 检查行名
    
    common_spots <- intersect(rownames(location), colnames(prop_matrix))
    cat("Common spots count:", length(common_spots), "\n")  # 检查匹配数
    
    if (length(common_spots) == 0) 
      stop(paste(method_name, "无匹配的spot名称"))
    
    # 数据对齐
    common_spots <- intersect(rownames(location), colnames(prop_matrix))
    if (length(common_spots) == 0) stop(paste(method_name, "无匹配的spot名称"))
    prop_matrix <- prop_matrix[, common_spots, drop = FALSE]
    location_sub <- location[common_spots, , drop = FALSE]
    
    # 绘制单方法饼图
    p <- draw.pie(
      proportion = prop_matrix[ct_names_clean,],
      spatial_location = location_sub,
      colors = colors,
      radius = radius,
      methodTitle = titles[i],
      font_size = font_size,seed=seed
    ) + 
      custom_theme() +
      guides(fill = guide_legend(
        title = NULL,
        title.position = "top",
        ncol = ceiling(sqrt(ncol(prop_matrix))),  # 图例自动分列
        keywidth = unit(0.5,"lines"),
        keyheight = unit(0.5,"lines")
      ))
      
      # 控制图例显示
      if (i < length(RESULTS)) {
        p <- p + theme(legend.position = "none")
      } else {
        p <- p + 
          theme(legend.position = "bottom") +
          guides(fill = guide_legend(
            title = NULL,
            title.position = "bottom",
            ncol = ncol,          # 设置图例列数
            byrow = TRUE       # 按行填充图例项
          ))
      }
      
      return(p)
  })
  
  # 组合所有图形
  if (!is.null(layer)) {
    # 将layer图作为第一个图
    all_plots <- c(list(layer_pie), plot_list)
    titles <- c("Layer Annotations", titles)
  } else {
    all_plots <- plot_list
  }
    # 拼合参数
    # if (is.null(ncol)) ncol <- ceiling(sqrt(length(RESULTS)))
  # 计算最优列数以形成近似矩形布局
  if (is.null(ncol)) {
    n_plots <- length(all_plots)
    # 尝试找到一个接近平方根的整数列数
    ncol <- ceiling(sqrt(n_plots))
    # 检查是否能形成更紧凑的矩形
    possible_cols <- ceiling(n_plots / floor(sqrt(n_plots)))
    if (abs(n_plots/ncol - floor(n_plots/ncol)) > 
        abs(n_plots/possible_cols - floor(n_plots/possible_cols))) {
      ncol <- possible_cols
    }
  }
    combined_plot <- patchwork::wrap_plots(all_plots, ncol = ncol)
    
    # 添加子图标签
    if (show_tags) {
      combined_plot <- combined_plot + 
        patchwork::plot_annotation(
          tag_levels = 'A',
          theme = theme(
            plot.tag = element_text(
              size = font_size + 2, 
              face = "bold",
              family = family
            )
          )
        )
    }
    
    return(combined_plot)
}



plot_dim_reduction <- function(expr_matrix, meta, 
                               method = c("pca", "tsne", "umap"),
                               n_top_genes = 2000,
                               point_size = 1.5,
                               label_clusters = TRUE,
                               palette = NULL,
                               plot_title = NULL,
                               return_model = FALSE) {
  
  # 检查必要包
  required_pkgs <- c("ggplot2", "Rtsne", "umap", "dplyr", "ggrepel")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("请先安装", pkg, "包：install.packages('", pkg, "')", sep = ""))
    }
  }
  
  # 参数校验
  method <- match.arg(method)
  if (!"cellType" %in% colnames(meta)) {
    stop("meta数据框中必须包含cellType列")
  }
  if (!all(colnames(expr_matrix) %in% rownames(meta))) {
    stop("expr_matrix的列名必须与meta的行名匹配")
  }
  
  # 数据预处理
  meta <- meta[colnames(expr_matrix), , drop = FALSE]  # 对齐样本顺序
  expr_matrix <- as.matrix(expr_matrix)
  
  # 选择高变基因（HVG）
  if (nrow(expr_matrix) > n_top_genes) {
    gene_vars <- matrixStats::rowVars(expr_matrix)
    top_genes <- order(gene_vars, decreasing = TRUE)[1:n_top_genes]
    expr_matrix <- expr_matrix[top_genes, ]
  }
  
  # 转置矩阵（细胞×基因）
  data_for_dimred <- t(expr_matrix)
  
  # 执行降维
  set.seed(42)  # 保证可重复性
  if (method == "pca") {
    dimred <- prcomp(data_for_dimred, scale. = TRUE)$x[, 1:2]
    colnames(dimred) <- c("Dim1", "Dim2")
  } else if (method == "tsne") {
    dimred <- Rtsne::Rtsne(data_for_dimred, perplexity = 30, 
                           check_duplicates = FALSE)$Y
    colnames(dimred) <- c("Dim1", "Dim2")
  } else if (method == "umap") {
    dimred <- umap::umap(data_for_dimred)$layout
    colnames(dimred) <- c("Dim1", "Dim2")
  }
  
  # 合并降维结果与元数据
  plot_data <- data.frame(
    Sample = rownames(dimred),
    Dim1 = dimred[, 1],
    Dim2 = dimred[, 2],
    CellType = meta$cellType
  )
  
  # 计算各类别中心位置（用于标签）
  label_pos <- plot_data %>%
    group_by(CellType) %>%
    summarise(
      Dim1 = median(Dim1),
      Dim2 = median(Dim2)
    )
  
  # 颜色设置
  n_types <- length(unique(meta$cellType))
  if (is.null(palette)) {
    if (n_types <= 8) {
      palette <- RColorBrewer::brewer.pal(n_types, "Set2")
    } else {
      palette <- scales::hue_pal()(n_types)
    }
  }
  
  # 绘制图形
  p <- ggplot(plot_data, aes(x = Dim1, y = Dim2, color = CellType)) +
    geom_point(size = point_size, alpha = 0.7) +
    scale_color_manual(values = palette) +
    labs(
      x = paste0(toupper(method), " 1"),
      y = paste0(toupper(method), " 2"),
      color = "Cell Type",
      title = plot_title
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  # 添加类别标签
  if (label_clusters) {
    p <- p + ggrepel::geom_label_repel(
      data = label_pos,
      aes(label = CellType),
      color = "black",
      size = 3,
      box.padding = 0.5,
      show.legend = FALSE
    )
  }
  
  # 返回结果
  if (return_model) {
    return(list(plot = p, model = dimred))
  } else {
    return(p)
  }
}

# 
# 
# pcaPlot <- function(data_combine, B1, B2, norm = "None", output_pptx = "PCA_plot.pptx") {
#   library(ggplot2)
#   library(scater)
#   library(officer)
#   
#   # 数据标准化（可选）
#   if(norm %in% c("cpm", "CPM")) {
#     data_combine <- calculateCPM(data_combine)
#   }
#   
#   # PCA计算
#   pca_result <- prcomp(t(data_combine), scale. = TRUE)
#   pca_df <- as.data.frame(pca_result$x[, 1:2])
#   pca_df$label <- c(rep("Pseudo-spots", ncol(B1)), rep("Real-spots", ncol(B2)))
#   colnames(pca_df) <- c("PC1", "PC2", "label")
#   
#   # 计算方差解释率
#   var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
#   
#   # 绘图
#   p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = label)) +
#     geom_point(size = 3, alpha = 0.8) +
#     scale_color_manual(values = c("#B291B5", "#3BA997")) + 
#     # xlab(paste0("PC1 (", var_explained[1], "%)")) +
#     # ylab(paste0("PC2 (", var_explained[2], "%)")) +
#     xlab(paste0("PC1")) +
#     ylab(paste0("PC2")) +
#     theme_classic() +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "right",
#       legend.title = element_blank(),
#       axis.line = element_line(size = 0.5, color = "black"),
#       axis.text = element_text(color = "black", size = 10),
#       axis.title = element_text(color = "black", size = 12)
#     ) +
#     guides(color = guide_legend(override.aes = list(size = 4)))
#   
#   # 输出到PPTX
#   doc <- read_pptx()
#   doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
#   doc <- ph_with(doc, value = p, location = ph_location_fullsize())
#   print(doc, target = output_pptx)
#   
#   return(p)  # 返回ggplot对象以便进一步调整
# }


pcaPlot <- function(data_combine, B1, B2, norm = "None", 
                    point_size = 3, point_alpha = 0.8,
                    output_pptx = "PCA_plot.pptx") {
  # Load required packages
  library(ggplot2)
  library(scater)
  library(officer)
  
  # Data normalization (optional)
  if(norm %in% c("cpm", "CPM")) {
    data_combine <- calculateCPM(data_combine)
  }
  
  # PCA calculation
  pca_result <- prcomp(t(data_combine), scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$label <- c(rep("Pseudo-spots", ncol(B1)), rep("Real-spots", ncol(B2)))
  colnames(pca_df) <- c("PC1", "PC2", "label")
  
  # Calculate variance explained
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  
  # Plotting with four borders
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = label)) +
    # Add points with customizable size and alpha
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = c("#B291B5", "#3BA997")) + 
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    # Custom theme with four borders
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      # Add all four borders
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Output to PPTX if filename is provided
  if(!is.null(output_pptx)) {
    # doc <- read_pptx()
    # doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
    # doc <- ph_with(doc, value = p, location = ph_location_fullsize())
    # print(doc, target = output_pptx)
    library(eoffice)
    topptx(p,output_pptx,width = 10,height = 8)
  }
  
  return(p)  # Return ggplot object for further adjustments
}

pcaPlotWithCellType <- function(data_combine, B1, B2, celltype_props = NULL,
                                norm = "None", point_size = 3, point_alpha = 0.8,
                                celltype_colors = NULL, output_pptx = NULL) {
  # Load required packages
  library(ggplot2)
  library(scater)
  library(scales)
  
  # Data normalization (optional)
  if(norm %in% c("cpm", "CPM")) {
    data_combine <- calculateCPM(data_combine)
  }
  
  # PCA calculation
  pca_result <- prcomp(t(data_combine), scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  
  # Add spot type labels
  pca_df$label <- c(rep("Pseudo-spots", ncol(B1)), rep("Real-spots", ncol(B2)))
  
  # Add dominant celltype if provided
  if(!is.null(celltype_props)) {
    # Check if celltype_props matches our data
    if(nrow(celltype_props) != nrow(pca_df)) {
      stop("Celltype proportions matrix doesn't match the number of spots")
    }
    
    # Get dominant celltype for each spot
    dominant_ct <- apply(celltype_props, 1, function(x) {
      colnames(celltype_props)[which.max(x)]
    })
    pca_df$dominant_celltype <- factor(dominant_ct)
    
    # Set default colors if not provided
    if(is.null(celltype_colors)) {
      n_ct <- length(levels(pca_df$dominant_celltype))
      celltype_colors <- hue_pal()(n_ct)
      names(celltype_colors) <- levels(pca_df$dominant_celltype)
    } else {
      # Check if all celltypes have colors
      missing_ct <- setdiff(levels(pca_df$dominant_celltype), names(celltype_colors))
      if(length(missing_ct) > 0) {
        warning("Some cell types missing colors: ", paste(missing_ct, collapse = ", "))
        # Add default colors for missing types
        n_missing <- length(missing_ct)
        default_colors <- hue_pal()(n_missing)
        names(default_colors) <- missing_ct
        celltype_colors <- c(celltype_colors, default_colors)
      }
    }
  }
  
  # Calculate variance explained
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  
  # Create base plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.2, "cm"),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  
  # Add points with different aesthetics based on available data
  if(!is.null(celltype_props)) {
    p <- p + 
      geom_point(aes(shape = label, fill = dominant_celltype),
                 size = point_size, alpha = point_alpha) +
      scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24),
                         name = "Spot type") +
      scale_fill_manual(values = celltype_colors,
                        name = "Dominant cell type") +
      guides(
        shape = guide_legend(order = 1),
        fill = guide_legend(override.aes = list(shape = 21), order = 2)
      )
  } else {
    p <- p +
      geom_point(aes(shape = label, color = label),
                 size = point_size, alpha = point_alpha) +
      scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24),
                         name = "Spot type") +
      scale_color_manual(values = c("Pseudo-spots" = "#B291B5", "Real-spots" = "#3BA997"),
                         name = "Spot type") +
      guides(shape = guide_legend(order = 1),
             color = guide_legend(order = 1))
  }
  
  # Output to PPTX if requested
  if(!is.null(output_pptx)) {
    if(!requireNamespace("eoffice", quietly = TRUE)) {
      install.packages("eoffice")
    }
    eoffice::topptx(p, output_pptx, width = 10, height = 8)
  }
  
  return(p)
}


compare_celltype_pcc <- function(true_matrix, pred_batch_matrix, pred_nobatch_matrix,colorset = c("#4DBBD5","#E64B35")) {
  # 输入验证
  if (!all(dim(true_matrix) == dim(pred_batch_matrix)) || 
      !all(dim(true_matrix) == dim(pred_nobatch_matrix))) {
    stop("All input matrices must have the same dimensions")
  }
  
  # if (!identical(rownames(true_matrix), rownames(pred_batch_matrix)) || 
  #     !identical(rownames(true_matrix), rownames(pred_nobatch_matrix))) {
  #   stop("All input matrices must have identical row names (cell types)")
  # }
  ct <- rownames(true_matrix)
  true_matrix <- true_matrix[ct,]
  pred_batch_matrix <- pred_batch_matrix[ct,]
  pred_nobatch_matrix <- pred_nobatch_matrix[ct,]
  # 计算每个样本的PCC
  calculate_sample_pcc <- function(mat1, mat2) {
    sapply(1:ncol(mat1), function(i) {
      cor(mat1[,i], mat2[,i], method = "pearson")
    })
  }
  
  # 为每个样本计算PCC
  pcc_batch <- calculate_sample_pcc(true_matrix, pred_batch_matrix)
  pcc_nobatch <- calculate_sample_pcc(true_matrix, pred_nobatch_matrix)
  
  # 正确创建数据框
  plot_data <- data.frame(
    PCC = c(pcc_batch, pcc_nobatch),
    Method = factor(rep(c("With Batch Correction", "Without Batch Correction"), 
                        each = length(pcc_batch))),
    SampleID = rep(colnames(true_matrix), 2)
  )
  
  # 执行t检验
  t_test_result <- t.test(pcc_batch, pcc_nobatch, paired = TRUE)
  
  # 定义颜色
  colors <- c("With Batch Correction" = colorset[1], 
              "Without Batch Correction" = colorset[2])
  
  # 绘制箱线图
  p <- ggplot(plot_data, aes(x = Method, y = PCC)) +
    # 透明箱体，彩色边框
    geom_boxplot(aes(color = Method), 
                 fill = "transparent", 
                 width = 0.6,
                 outlier.shape = NA,
                 size = 0.8) +
    # 彩色散点
    geom_jitter(aes(color = Method), 
                width = 0.15, 
                size = 2.5, 
                alpha = 0.7) +
    # 设置颜色
    scale_color_manual(values = colors) +
    # 标题和标签
    labs(title = "Comparison of Prediction Accuracy",
         subtitle = paste0("Paired t-test p-value = ", format.pval(t_test_result$p.value, digits = 3)),
         x = "",
         y = "Pearson Correlation Coefficient (PCC)") +
    # 主题设置：添加四边形外框
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 0.8),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(color = "black", size = 11),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      legend.position = "none",
      aspect.ratio = 0.8
    ) +
    # 添加显著性标记
    stat_compare_means(
      comparisons = list(c("With Batch Correction", "Without Batch Correction")),
      method = "t.test", 
      paired = TRUE,
      label = "p.signif", 
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")),
      size = 5,
      vjust = 0.5
    ) +
    # 添加y轴范围扩展
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  # 返回结果
  result <- list(
    plot = p,
    pcc_values = plot_data,
    t_test = t_test_result,
    summary_stats = aggregate(PCC ~ Method, data = plot_data, 
                              FUN = function(x) c(mean = mean(x), sd = sd(x)))
  )
  
  return(result)
}






plot_proportion_fit <- function(P_real, P_infer, title = "Real vs. Inferred Proportions") {
  
  library(ggplot2)
  library(ggpubr)
  # 确保输入矩阵的行名（细胞类型）和列名（spots）匹配
  if (!identical(rownames(P_real), rownames(P_infer))) {
    stop("Row names (cell types) of P_real and P_infer must match!")
  }
  
  # 将矩阵转换为长格式数据框
  melt_matrix <- function(mat, value_name) {
    as.data.frame.table(mat, responseName = value_name) %>%
      setNames(c("CellType", "Spot", value_name))
  }
  
  df_real <- melt_matrix(P_real, "Real")
  df_infer <- melt_matrix(P_infer, "Inferred")
  
  # 合并数据
  df_combined <- merge(df_real, df_infer, by = c("CellType", "Spot"))
  
  # 计算拟合优度指标（R-squared）
  r_squared <- round(cor(df_combined$Real, df_combined$Inferred)^2, 3)
  
  # 绘图
  p <- ggplot(df_combined, aes(x = Real, y = Inferred, color = CellType)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # 添加y=x参考线
    labs(
      x = "Real Proportion",
      y = "Inferred Proportion",
      title = title,
      subtitle = paste0("R-squared = ", r_squared)
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_color_brewer(palette = "Set1") +  # 使用Set1调色板区分细胞类型
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")  # 添加相关系数
  
  return(p)
}



#' Visualize cell type proportions as a spatial heatmap
#'
#' @param prop_matrix A matrix/dataframe of cell type proportions (spots x cell types)
#' @param cell_type Character string specifying the cell type to visualize
#' @param locations A dataframe with spot coordinates (must contain 'x' and 'y' columns)
#' @param title Plot title (optional)
#' @param point_size Size of spots in the plot (default = 3)
#' @param palette Color palette (default = "viridis")
#' @param legend_title Title for the color legend (optional)
#'
#' @return A ggplot object displaying the spatial heatmap
#'
#' @examples
#' # prop_matrix <- read.csv("cell_proportions.csv", row.names = 1)
#' # locations <- data.frame(x = runif(100), y = runif(100))
#' # visualize_spatial_proportion(prop_matrix, "T_cells", locations)
visualize_spatial_proportion <- function(prop_matrix, cell_type, locations, 
                                         title = NULL, point_size = 3,
                                         palette = "viridis", 
                                         legend_title = NULL) {
  prop_matrix <- t(prop_matrix)
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("viridis", quietly = TRUE) && palette == "viridis") {
    message("Using default palette instead of viridis (package not installed)")
    palette <- "RdYlBu"
  }
  
  # Validate inputs
  if (!cell_type %in% colnames(prop_matrix)) {
    stop(paste("Cell type", cell_type, "not found in proportion matrix"))
  }
  if (nrow(prop_matrix) != nrow(locations)) {
    stop("Number of spots in proportion matrix and locations don't match")
  }
  if (!all(c("x", "y") %in% colnames(locations))) {
    stop("Locations dataframe must contain 'x' and 'y' columns")
  }
  
  # Prepare data
  plot_data <- data.frame(
    x = locations$x,
    y = locations$y,
    proportion = prop_matrix[, cell_type]
  )
  
  # Set default title if none provided
  if (is.null(title)) {
    title <- cell_type
  }
  
  # Set default legend title if none provided
  if (is.null(legend_title)) {
    legend_title <- "Proportion"
  }
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = proportion)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_gradient(low = "white", high = "red", 
                                  limits = c(0, 1),  # 确保比例在0-1范围内
                                  na.value = "grey50" ) +  # 处理NA值
    ggplot2::labs(title = title, x = "X coordinate", y = "Y coordinate", color = legend_title) +
    ggplot2::theme_minimal() +
    ggplot2::coord_fixed() 
  print(p)
  return(p)
}

locationPlot <- function(loc_B1, loc_B2, norm = "None", 
                         point_size = 3, point_alpha = 0.8, axis_ratio = "auto", 
                         output_pptx = NULL) {
  # Load required packages
  library(ggplot2)
  library(scater)
  library(officer)
  
  # Prepare location data
  if("x" %in% colnames(loc_B1)){
    loc_df <- rbind(
      data.frame(X = loc_B1$x, Y = loc_B1$y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$x, Y = loc_B2$y, label = "Real-spots")
    )
  } else {
    loc_df <- rbind(
      data.frame(X = loc_B1$X, Y = loc_B1$Y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$X, Y = loc_B2$Y, label = "Real-spots")
    )
  }
  
  # 计算坐标范围比例
  x_range <- diff(range(loc_df$X))
  y_range <- diff(range(loc_df$Y))
  range_ratio <- y_range / x_range
  
  # 自动调整图形宽高比
  if(axis_ratio == "auto") {
    plot_ratio <- ifelse(range_ratio > 2, range_ratio * 0.6, 
                         ifelse(range_ratio < 0.5, 1/(range_ratio * 0.6), 1))
  } else {
    plot_ratio <- axis_ratio
  }
  
  # 定义箭头样式
  arrow_style <- arrow(
    angle = 15,      # 箭头角度
    length = unit(0.2, "cm"),  # 箭头长度
    type = "closed"  # 封闭箭头
  )
  
  # 绘图（添加箭头）
  p <- ggplot(loc_df, aes(x = X, y = Y, color = label)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = c("#ffbc98", "#95B0B0")) + 
    xlab("X Coordinate") +
    ylab("Y Coordinate") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      panel.border = element_blank(),
      # 添加箭头（关键修改部分）
      axis.line.x = element_line(
        color = "black", 
        size = 0.5, 
        arrow = arrow_style  # x轴箭头
      ),
      axis.line.y = element_line(
        color = "black", 
        size = 0.5, 
        arrow = arrow_style  # y轴箭头
      ),
      aspect.ratio = 1/plot_ratio
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    coord_fixed()
  
  # 输出到PPTX（可选）
  if(!is.null(output_pptx)) {
    library(eoffice)
    topptx(p, output_pptx, width = 6, height = 7)
  }
  
  return(p)
}


locationPlot_2 <- function(loc_B1, loc_B2, celltype_props = NULL, 
                           point_size = 3, point_alpha = 0.8, axis_ratio = "auto", 
                           celltype_colors = NULL, output_pptx = NULL,scale_y = c(320,550),scale_x=c(-150,-120)) {
  
  # Load required packages
  library(ggplot2)
  library(scales)  # 用于 hue_pal 颜色生成
  library(ggbreak)
  
  # Prepare location data
  if("x" %in% colnames(loc_B1)) {
    loc_df <- rbind(
      data.frame(X = loc_B1$x, Y = loc_B1$y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$x, Y = loc_B2$y, label = "Real-spots")
    )
  } else {
    loc_df <- rbind(
      data.frame(X = loc_B1$X, Y = loc_B1$Y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$X, Y = loc_B2$Y, label = "Real-spots")
    )
  }
  
  # 添加形状信息
  loc_df$shape <- ifelse(loc_df$label == "Pseudo-spots", "circle", "triangle")
  
  # 如果有细胞类型比例矩阵
  if(!is.null(celltype_props)) {
    # 验证行名匹配
    if(!all(rownames(loc_B1) %in% rownames(celltype_props)) && 
       !all(rownames(loc_B2) %in% rownames(celltype_props))) {
      stop("Cell type proportion matrix rownames don't match location data rownames")
    }
    
    # 获取优势细胞类型
    dominant_ct <- apply(celltype_props, 1, function(x) {
      colnames(celltype_props)[which.max(x)]
    })
    
    loc_df$dominant_celltype <- factor(c(
      dominant_ct[rownames(loc_B1)],
      dominant_ct[rownames(loc_B2)]
    ))
    
    # 设置默认颜色（如果未提供）
    if(is.null(celltype_colors)) {
      n_ct <- length(levels(loc_df$dominant_celltype))
      celltype_colors <- hue_pal()(n_ct)
      names(celltype_colors) <- levels(loc_df$dominant_celltype)
    } else {
      # 检查提供的颜色是否覆盖所有细胞类型
      missing_ct <- setdiff(levels(loc_df$dominant_celltype), names(celltype_colors))
      if(length(missing_ct) > 0) {
        warning("Some cell types missing colors: ", paste(missing_ct, collapse = ", "))
        # 为缺失的细胞类型添加默认颜色
        n_missing <- length(missing_ct)
        default_colors <- hue_pal()(n_missing)
        names(default_colors) <- missing_ct
        celltype_colors <- c(celltype_colors, default_colors)
      }
    }
  }
  
  # 计算坐标范围比例
  x_range <- diff(range(loc_df$X))
  y_range <- diff(range(loc_df$Y))
  range_ratio <- y_range / x_range
  
  # 自动调整图形宽高比
  if(axis_ratio == "auto") {
    plot_ratio <- ifelse(range_ratio > 2, range_ratio * 0.6, 
                         ifelse(range_ratio < 0.5, 1/(range_ratio * 0.6), 1))
  } else {
    plot_ratio <- axis_ratio
  }
  
  # 定义箭头样式
  arrow_style <- arrow(angle = 15, length = unit(0.2, "cm"), type = "closed")
  
  # 创建基础绘图
  p <- ggplot(loc_df, aes(x = X, y = Y)) +
    xlab("X Coordinate") + 
    ylab("Y Coordinate") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.2, "cm"),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.5, arrow = arrow_style),
      axis.line.y = element_line(color = "black", size = 0.5, arrow = arrow_style),)+scale_y_break(scale_y)
      # aspect.ratio = 1/plot_ratio)
    # ) +
    # coord_fixed()
  
  
  
  ###添加截断
  # p <- p + scale_y_break(scale_y)
  
  # p <- p + scale_x_break(c(-150,-120))
  
  # 添加几何对象
  if (!is.null(celltype_props)) {
    p <- p + 
      geom_point(aes(shape = label, fill = dominant_celltype),
                 size = point_size, alpha = point_alpha) +
      scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24),
                         name = "Spot type") +
      scale_fill_manual(values = celltype_colors,
                        name = "Dominant cell type") +
      guides(
        shape = guide_legend(order = 1),
        fill = guide_legend(override.aes = list(shape = 21), order = 2)
      )
  } else {
    p <- p +
      geom_point(aes(shape = label, color = label),
                 size = point_size, alpha = point_alpha) +
      scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24),
                         name = "Spot Type") +
      scale_color_manual(values = c("Pseudo-spots" = "#ffbc98", "Real-spots" = "#95B0B0"),
                         name = "Spot Type") +
      guides(shape = guide_legend(order = 1),
             color = guide_legend(order = 1))
  }
  
  
  # 输出到PPTX（可选）
  if(!is.null(output_pptx)) {
    if(!requireNamespace("eoffice", quietly = TRUE)) {
      install.packages("eoffice")
    }
    eoffice::topptx(p, output_pptx, width = 6, height = 7)
  }
  
  return(p)
}


locationPlot_break <- function(loc_B1, loc_B2, celltype_props = NULL, 
                           point_size = 3, point_alpha = 0.8, axis_ratio = "auto", 
                           celltype_colors = NULL, output_pptx = NULL) {
  
  # 加载必要的包
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(gridExtra)
  
  # 准备位置数据
  if("x" %in% colnames(loc_B1)) {
    loc_df <- rbind(
      data.frame(X = loc_B1$x, Y = loc_B1$y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$x, Y = loc_B2$y, label = "Real-spots")
    )
  } else {
    loc_df <- rbind(
      data.frame(X = loc_B1$X, Y = loc_B1$Y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$X, Y = loc_B2$Y, label = "Real-spots")
    )
  }
  
  # 添加形状信息
  loc_df$shape <- ifelse(loc_df$label == "Pseudo-spots", "circle", "triangle")
  
  # 如果有细胞类型比例矩阵
  if(!is.null(celltype_props)) {
    # 验证行名匹配
    if(!all(rownames(loc_B1) %in% rownames(celltype_props)) && 
       !all(rownames(loc_B2) %in% rownames(celltype_props))) {
      stop("Cell type proportion matrix rownames don't match location data rownames")
    }
    
    # 获取优势细胞类型
    dominant_ct <- apply(celltype_props, 1, function(x) {
      colnames(celltype_props)[which.max(x)]
    })
    
    loc_df$dominant_celltype <- factor(c(
      dominant_ct[rownames(loc_B1)],
      dominant_ct[rownames(loc_B2)]
    ))
    
    # 设置默认颜色（如果未提供）
    if(is.null(celltype_colors)) {
      n_ct <- length(levels(loc_df$dominant_celltype))
      celltype_colors <- hue_pal()(n_ct)
      names(celltype_colors) <- levels(loc_df$dominant_celltype)
    } else {
      # 检查提供的颜色是否覆盖所有细胞类型
      missing_ct <- setdiff(levels(loc_df$dominant_celltype), names(celltype_colors))
      if(length(missing_ct) > 0) {
        warning("Some cell types missing colors: ", paste(missing_ct, collapse = ", "))
        # 为缺失的细胞类型添加默认颜色
        n_missing <- length(missing_ct)
        default_colors <- hue_pal()(n_missing)
        names(default_colors) <- missing_ct
        celltype_colors <- c(celltype_colors, default_colors)
      }
    }
  }
  
  # 定义箭头样式
  arrow_style <- arrow(angle = 15, length = unit(0.2, "cm"), type = "closed")
  
  # 创建基础绘图（无图例）
  base_plot <- function(show_x = TRUE, show_y = TRUE) {
    p <- ggplot(loc_df, aes(x = X, y = Y)) +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12),
        # panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
        # panel.border =  element_blank()
      )
    
    # 控制坐标轴显示
    if (!show_x) {
      p <- p + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    } else {
      p <- p + xlab("") + 
        theme(axis.line.x = element_line(color = "black", size = 0.5, arrow = arrow_style))
    }
    
    if (!show_y) {
      p <- p + theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    } else {
      p <- p + ylab("") + 
        theme(axis.line.y = element_line(color = "black", size = 0.5, arrow = arrow_style))
    }
    
    # 添加几何对象
    if (!is.null(celltype_props)) {
      p <- p + 
        geom_point(aes(shape = label, fill = dominant_celltype),
                   size = point_size, alpha = point_alpha) +
        scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24))
    } else {
      p <- p +
        geom_point(aes(shape = label, color = label),
                   size = point_size, alpha = point_alpha) +
        scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24)) +
        scale_color_manual(values = c("Pseudo-spots" = "#ffbc98", "Real-spots" = "#95B0B0"))
    }
    return(p)
  }
  
  # 创建四个子图
  plot_top_left <- base_plot(show_x = FALSE, show_y = TRUE) +
    coord_cartesian(xlim = c(-180, -150), ylim = c(550, 580)) +
    theme(plot.margin = unit(c(0.1, 0, 0, 0.1), "cm"))
  
  plot_bottom_left <- base_plot(show_x = TRUE, show_y = TRUE) +
    coord_cartesian(xlim = c(-180, -150), ylim = c(300, 330)) +
    theme(plot.margin = unit(c(0, 0, 0.1, 0.1), "cm"))
  
# 右上空白图 - 完全空白无边框
  plot_top_right <- ggplot() + 
    theme_void() +
    theme(plot.margin = unit(c(0.1, 0.1, 0, 0), "cm"),
          panel.background = element_blank())
  
  plot_bottom_right <- base_plot(show_x = TRUE, show_y = FALSE) +
    coord_cartesian(xlim = c(-110, -80), ylim = c(300, 330)) +
    theme(plot.margin = unit(c(0, 0.1, 0.1, 0), "cm"))
  
  # 修复图例提取问题
  if (!is.null(celltype_props)) {
    # 创建专门的图例绘图
    legend_data <- unique(loc_df[, c("label", "dominant_celltype")])
    
    legend_plot <- ggplot(legend_data, aes(x = 1, y = 1)) +
      geom_point(aes(shape = label, fill = dominant_celltype), size = 0) +
      scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24)) +
      scale_fill_manual(values = celltype_colors) +
      guides(
        shape = guide_legend(
          title = NULL, 
          override.aes = list(size = 4, fill = "gray50")
        ),
        fill = guide_legend(
          title = NULL, 
          override.aes = list(shape = 21, size = 4)
        )
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.5, "cm"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
      )
    
    # 安全提取图例
    legend <- tryCatch({
      get_legend(legend_plot)
    }, error = function(e) {
      message("图例提取错误: ", e$message)
      NULL
    })
  } else {
    # 无细胞类型情况的图例
    legend_plot <- ggplot(loc_df, aes(x = X, y = Y)) +
      geom_point(aes(shape = label, color = label), size = point_size) +
      scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24)) +
      scale_color_manual(values = c("Pseudo-spots" = "#ffbc98", "Real-spots" = "#95B0B0")) +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm")
      )
    
    legend <- get_legend(legend_plot)
  }
  
  # 组合四个子图
  plots <- arrangeGrob(
    plot_top_left, plot_top_right,
    plot_bottom_left, plot_bottom_right,
    nrow = 2, ncol = 2,
    widths = c(1, 1), heights = c(1, 1)
  )
  
  # 创建包含图例的完整图形
  if (!is.null(legend)) {
    combined_plot <- arrangeGrob(
      plots,
      legend,
      nrow = 2,
      heights = c(10, 1)
    )
  } else {
    combined_plot <- plots
    warning("图例创建失败，继续无图例输出")
  }
  
  # 添加全局坐标轴标题
  final_plot <- grid.arrange(
    combined_plot,
    left = textGrob("Y Coordinate", rot = 90, 
                    gp = gpar(fontsize = 14, fontface = "bold", col = "black")),
    bottom = textGrob("X Coordinate", 
                      gp = gpar(fontsize = 14, fontface = "bold", col = "black")),
    padding = unit(1, "line")
  )
  
  # 修复PPTX输出问题
  if(!is.null(output_pptx)) {
    if(!requireNamespace("officer", quietly = TRUE)) {
      install.packages("officer")
    }
    if(!requireNamespace("rvg", quietly = TRUE)) {
      install.packages("rvg")
    }
    
    # 创建PPTX文档
    doc <- officer::read_pptx()
    doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
    
    # 创建临时绘图文件
    temp_plot <- tempfile(fileext = ".png")
    png(temp_plot, width = 2000, height = 1500, res = 300)
    grid.draw(final_plot)
    dev.off()
    
    # 将绘图添加到PPTX
    doc <- officer::ph_with(
      doc, 
      value = officer::external_img(temp_plot),
      location = officer::ph_location_type(type = "body")
    )
    
    # 保存文档
    print(doc, target = output_pptx)
    message("成功保存PPTX到: ", output_pptx)
    
    # 清理临时文件
    unlink(temp_plot)
  }
  
  return(final_plot)
}
locationPlot_break1 <- function(loc_B1, loc_B2, celltype_props = NULL, 
                               point_size = 3, point_alpha = 0.8, axis_ratio = "auto", 
                               celltype_colors = NULL, output_pptx = NULL) {
  
  # 加载必要的包
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(gridExtra)
  library(grid)
  
  # 准备位置数据
  if("x" %in% colnames(loc_B1)) {
    loc_df <- rbind(
      data.frame(X = loc_B1$x, Y = loc_B1$y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$x, Y = loc_B2$y, label = "Real-spots")
    )
  } else {
    loc_df <- rbind(
      data.frame(X = loc_B1$X, Y = loc_B1$Y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$X, Y = loc_B2$Y, label = "Real-spots")
    )
  }
  
  # 添加形状信息
  loc_df$shape <- ifelse(loc_df$label == "Pseudo-spots", "circle", "triangle")
  
  # 如果有细胞类型比例矩阵
  if(!is.null(celltype_props)) {
    # 验证行名匹配
    if(!all(rownames(loc_B1) %in% rownames(celltype_props)) || 
       !all(rownames(loc_B2) %in% rownames(celltype_props))) {
      stop("Cell type proportion matrix rownames don't match location data rownames")
    }
    
    # 获取优势细胞类型
    dominant_ct <- apply(celltype_props, 1, function(x) {
      colnames(celltype_props)[which.max(x)]
    })
    
    loc_df$dominant_celltype <- factor(c(
      dominant_ct[rownames(loc_B1)],
      dominant_ct[rownames(loc_B2)]
    ))
    
    # 设置默认颜色（如果未提供）
    if(is.null(celltype_colors)) {
      n_ct <- length(levels(loc_df$dominant_celltype))
      celltype_colors <- hue_pal()(n_ct)
      names(celltype_colors) <- levels(loc_df$dominant_celltype)
    } else {
      # 检查提供的颜色是否覆盖所有细胞类型
      missing_ct <- setdiff(levels(loc_df$dominant_celltype), names(celltype_colors))
      if(length(missing_ct) > 0) {
        warning("Some cell types missing colors: ", paste(missing_ct, collapse = ", "))
        # 为缺失的细胞类型添加默认颜色
        n_missing <- length(missing_ct)
        default_colors <- hue_pal()(n_missing)
        names(default_colors) <- missing_ct
        celltype_colors <- c(celltype_colors, default_colors)
      }
    }
  }
  
  # 定义箭头样式
  arrow_style <- arrow(angle = 15, length = unit(0.2, "cm"), type = "closed")
  
  # 创建基础绘图（无图例）
  base_plot <- function(show_x = TRUE, show_y = TRUE) {
    p <- ggplot(loc_df, aes(x = X, y = Y)) +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_blank(),  # 移除所有轴标题
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
    
    # 控制坐标轴显示
    if (!show_x) {
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    } else {
      p <- p + 
        theme(axis.line.x = element_line(color = "black", size = 0.5, arrow = arrow_style))
    }
    
    if (!show_y) {
      p <- p + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    } else {
      p <- p + 
        theme(axis.line.y = element_line(color = "black", size = 0.5, arrow = arrow_style))
    }
    
    # 添加几何对象
    if (!is.null(celltype_props)) {
      p <- p + 
        geom_point(aes(shape = label, fill = dominant_celltype),
                   size = point_size, alpha = point_alpha) +
        scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24)) +
        scale_fill_manual(values = celltype_colors)
    } else {
      p <- p +
        geom_point(aes(shape = label, color = label),
                   size = point_size, alpha = point_alpha) +
        scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24)) +
        scale_color_manual(values = c("Pseudo-spots" = "#ffbc98", "Real-spots" = "#95B0B0"))
    }
    return(p)
  }
  
  # 创建四个子图
  plot_top_left <- base_plot(show_x = FALSE, show_y = TRUE) +
    coord_cartesian(xlim = c(-180, -150), ylim = c(550, 580)) +
    theme(plot.margin = unit(c(0.1, 0, 0, 0.1), "cm"))
  
  plot_bottom_left <- base_plot(show_x = TRUE, show_y = TRUE) +
    coord_cartesian(xlim = c(-180, -150), ylim = c(300, 330)) +
    theme(plot.margin = unit(c(0, 0, 0.1, 0.1), "cm"))
  
  # 右上空白图 - 完全空白无边框
  plot_top_right <- ggplot() + 
    theme_void() +
    theme(plot.margin = unit(c(0.1, 0.1, 0, 0), "cm"),
          panel.background = element_blank())
  
  plot_bottom_right <- base_plot(show_x = TRUE, show_y = FALSE) +
    coord_cartesian(xlim = c(-110, -80), ylim = c(300, 330)) +
    theme(plot.margin = unit(c(0, 0.1, 0.1, 0), "cm"))
  
  # 创建图例
  if (!is.null(celltype_props)) {
    # 创建图例专用绘图
    legend_data <- unique(loc_df[, c("label", "dominant_celltype")])
    
    legend_plot <- ggplot(legend_data, aes(x = 1, y = 1)) +
      geom_point(aes(shape = label, fill = dominant_celltype), size = 0) +
      scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24)) +
      scale_fill_manual(values = celltype_colors) +
      guides(
        shape = guide_legend(
          title = NULL, 
          override.aes = list(size = 4, fill = "gray50")
        ),
        fill = guide_legend(
          title = NULL, 
          override.aes = list(shape = 21, size = 4)
        )
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.5, "cm"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
      )
    
    # 提取图例
    legend <- get_legend(legend_plot)
  } else {
    # 无细胞类型情况的图例
    legend_plot <- ggplot(loc_df, aes(x = X, y = Y)) +
      geom_point(aes(shape = label, color = label), size = point_size) +
      scale_shape_manual(values = c("Pseudo-spots" = 21, "Real-spots" = 24)) +
      scale_color_manual(values = c("Pseudo-spots" = "#ffbc98", "Real-spots" = "#95B0B0")) +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm")
      )
    
    legend <- get_legend(legend_plot)
  }
  
  # 使用cowplot组合图形
  # 组合四个子图
  plot_matrix <- plot_grid(
    plot_top_left, plot_top_right,
    plot_bottom_left, plot_bottom_right,
    nrow = 2, ncol = 2,
    align = "hv"
  )
  
  # 添加图例
  plot_with_legend <- plot_grid(
    plot_matrix,
    legend,
    nrow = 2,
    rel_heights = c(10, 1)
  )
  
  # 添加全局坐标轴标题
  final_plot <- ggdraw() +
    cowplot::draw_plot(plot_with_legend, x = 0.05, y = 0.05, width = 0.9, height = 0.9) +
    draw_label("Y Coordinate", x = 0.02, y = 0.5, angle = 90, 
               size = 14, fontface = "bold") +
    draw_label("X Coordinate", x = 0.5, y = 0.02, 
               size = 14, fontface = "bold")
  
  # 输出到PPTX
  if(!is.null(output_pptx)) {
    if(!requireNamespace("eoffice", quietly = TRUE)) {
      install.packages("eoffice")
    }
    # 直接使用topptx输出
    eoffice::topptx(final_plot, output_pptx, width = 8, height = 6)
    message("成功保存PPTX到: ", output_pptx)
  }
  
  return(final_plot)
}

locationPlot_break2 <- function(loc_B1, loc_B2, celltype_props = NULL, 
                                point_size = 3, point_alpha = 0.8, axis_ratio = "auto", 
                                celltype_colors = NULL, output_pptx = NULL) {
  
  # 加载必要的包
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(gridExtra)
  library(grid)
  
  # 准备位置数据
  if("x" %in% colnames(loc_B1)) {
    loc_df <- rbind(
      data.frame(X = loc_B1$x, Y = loc_B1$y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$x, Y = loc_B2$y, label = "Real-spots")
    )
  } else {
    loc_df <- rbind(
      data.frame(X = loc_B1$X, Y = loc_B1$Y, label = "Pseudo-spots"),
      data.frame(X = loc_B2$X, Y = loc_B2$Y, label = "Real-spots")
    )
  }
  
  # 如果有细胞类型比例矩阵
  if(!is.null(celltype_props)) {
    # 验证行名匹配
    if(!all(rownames(loc_B1) %in% rownames(celltype_props)) || 
       !all(rownames(loc_B2) %in% rownames(celltype_props))) {
      stop("Cell type proportion matrix rownames don't match location data rownames")
    }
    
    # 获取优势细胞类型
    dominant_ct <- apply(celltype_props, 1, function(x) {
      colnames(celltype_props)[which.max(x)]
    })
    
    loc_df$dominant_celltype <- factor(c(
      dominant_ct[rownames(loc_B1)],
      dominant_ct[rownames(loc_B2)]
    ))
    
    # 设置默认颜色（如果未提供）
    if(is.null(celltype_colors)) {
      n_ct <- length(levels(loc_df$dominant_celltype))
      celltype_colors <- hue_pal()(n_ct)
      names(celltype_colors) <- levels(loc_df$dominant_celltype)
    } else {
      # 检查提供的颜色是否覆盖所有细胞类型
      missing_ct <- setdiff(levels(loc_df$dominant_celltype), names(celltype_colors))
      if(length(missing_ct) > 0) {
        warning("Some cell types missing colors: ", paste(missing_ct, collapse = ", "))
        # 为缺失的细胞类型添加默认颜色
        n_missing <- length(missing_ct)
        default_colors <- hue_pal()(n_missing)
        names(default_colors) <- missing_ct
        celltype_colors <- c(celltype_colors, default_colors)
      }
    }
  }
  
  # 定义箭头样式
  arrow_style <- arrow(angle = 15, length = unit(0.2, "cm"), type = "closed")
  
  # 创建基础绘图（无图例）
  base_plot <- function(show_x = TRUE, show_y = TRUE) {
    p <- ggplot(loc_df, aes(x = X, y = Y)) +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_blank(),  # 移除所有轴标题
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
    
    # 控制坐标轴显示
    if (!show_x) {
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    } else {
      p <- p + 
        theme(axis.line.x = element_line(color = "black", size = 0.5, arrow = arrow_style))
    }
    
    if (!show_y) {
      p <- p + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    } else {
      p <- p + 
        theme(axis.line.y = element_line(color = "black", size = 0.5, arrow = arrow_style))
    }
    
    # 添加几何对象
    if (!is.null(celltype_props)) {
      p <- p + 
        geom_point(aes(fill = dominant_celltype),
                   shape = 21, size = point_size, alpha = point_alpha) +
        scale_fill_manual(values = celltype_colors)
    } else {
      p <- p +
        geom_point(aes(shape = label, color = label),
                   size = point_size, alpha = point_alpha) +
        scale_shape_manual(values = c("Pseudo-spots" = 16, "Real-spots" = 16)) +
        scale_color_manual(values = c("Pseudo-spots" = "#ffbc98", "Real-spots" = "#95B0B0"))
    }
    return(p)
  }
  
  # 创建四个子图
  plot_top_left <- base_plot(show_x = FALSE, show_y = TRUE) +
    coord_cartesian(xlim = c(-180, -150), ylim = c(550, 580)) +
    theme(plot.margin = unit(c(0.1, 0, 0, 0.1), "cm"))
  
  plot_bottom_left <- base_plot(show_x = TRUE, show_y = TRUE) +
    coord_cartesian(xlim = c(-180, -150), ylim = c(300, 330)) +
    theme(plot.margin = unit(c(0, 0, 0.1, 0.1), "cm"))
  
  # 右上空白图 - 完全空白无边框
  plot_top_right <- ggplot() + 
    theme_void() +
    theme(plot.margin = unit(c(0.1, 0.1, 0, 0), "cm"),
          panel.background = element_blank())
  
  plot_bottom_right <- base_plot(show_x = TRUE, show_y = FALSE) +
    coord_cartesian(xlim = c(-110, -80), ylim = c(300, 330)) +
    theme(plot.margin = unit(c(0, 0.1, 0.1, 0), "cm"))
  
  # # 创建图例
  # if (!is.null(celltype_props)) {
  #   # 创建图例专用绘图
  #   legend_data <- unique(loc_df[, c("dominant_celltype")])
  #   
  #   legend_plot <- ggplot(legend_data, aes(x = 1, y = 1)) +
  #     geom_point(aes(fill = dominant_celltype), shape = 21, size = 0) +
  #     scale_fill_manual(values = celltype_colors) +
  #     guides(
  #       fill = guide_legend(
  #         title = NULL, 
  #         override.aes = list(shape = 21, size = 4)
  #       )
  #     ) +
  #     theme_void() +
  #     theme(
  #       legend.position = "bottom",
  #       legend.box = "horizontal",
  #       legend.direction = "horizontal",
  #       legend.spacing.x = unit(0.5, "cm"),
  #       legend.key = element_rect(fill = "white", color = "white"),
  #       legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  #     )
  #   
  #   # 提取图例
  #   legend <- get_legend(legend_plot)
  # } else {
  #   # 无细胞类型情况的图例
  #   legend_plot <- ggplot(loc_df, aes(x = X, y = Y)) +
  #     geom_point(aes(fill = label), shape = 21, size = point_size) +
  #     scale_fill_manual(values = c("Pseudo-spots" = "#ffbc98", "Real-spots" = "#95B0B0")) +
  #     theme(
  #       legend.position = "bottom",
  #       legend.box = "horizontal",
  #       legend.title = element_blank(),
  #       legend.spacing.x = unit(0.5, "cm")
  #     )
  #   
  #   legend <- get_legend(legend_plot)
  # }
  
  # 使用cowplot组合图形
  # 组合四个子图
  plot_matrix <- plot_grid(
    plot_top_left, plot_top_right,
    plot_bottom_left, plot_bottom_right,
    nrow = 2, ncol = 2,
    align = "hv"
  )
  
  # # 添加图例
  # plot_with_legend <- plot_grid(
  #   plot_matrix,
  #   legend,
  #   nrow = 2,
  #   rel_heights = c(10, 1)
  # )
  
  # 添加全局坐标轴标题
  final_plot <- ggdraw() +
    cowplot::draw_plot(plot_matrix, x = 0.05, y = 0.05, width = 0.9, height = 0.9) +
    draw_label("Y Coordinate", x = 0.02, y = 0.5, angle = 90, 
               size = 14, fontface = "bold") +
    draw_label("X Coordinate", x = 0.5, y = 0.02, 
               size = 14, fontface = "bold")
  
  # 输出到PPTX
  if(!is.null(output_pptx)) {
    if(!requireNamespace("eoffice", quietly = TRUE)) {
      install.packages("eoffice")
    }
    # 直接使用topptx输出
    eoffice::topptx(final_plot, output_pptx, width = 8, height = 6)
    message("成功保存PPTX到: ", output_pptx)
  }
  
  return(final_plot)
}


plotRealLayer <- function(labels, coordinates, 
                          title = "Layer annotations",
                          point.size = 3.5, 
                          layer_colors = c("GCL" = "#FCDABA", 
                                           "GL" = "#A7D2BA", 
                                           "ONL" = "#C39398", 
                                           "MCL" = "#ABC6E4"),
                          padding = 0.1,
                          legend.cex = 0.8,
                          legend.spacing = 1.5) {  # 新增legend.spacing控制行间距
  
  # 确保输入数据长度匹配
  if (length(labels) != nrow(coordinates)) {
    stop("Length of labels must match number of rows in coordinates")
  }
  
  # 创建颜色向量
  point_colors <- layer_colors[labels]
  
  # 设置绘图参数
  par(bty = "o", lwd = 2, col = "lightgray", xpd = TRUE, 
      mar = c(3, 3, 3, 10))  # 右侧边距增大到10以容纳图例
  
  # 计算坐标范围（增加padding比例的留白）
  x_range <- range(coordinates[, 1])
  y_range <- range(coordinates[, 2])
  x_padding <- diff(x_range) * padding
  y_padding <- diff(y_range) * padding
  
  # 绘制点图
  plot(coordinates, 
       pch = 16,
       cex = point.size,
       col = point_colors,
       main = "",
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty = "o",
       col.axis = "lightgray",
       lwd = 2,
       xlim = x_range + c(-x_padding, x_padding),
       ylim = y_range + c(-y_padding, y_padding))
  
  # 添加居中标题（更靠近图形）
  title(main = title, line = 0.5, cex.main = 1.2,cex.main = 1.5)
  
  # 添加图例（增大行间距）
  legend(x = par("usr")[2] * 1.05,  # 稍微右移
         y = mean(par("usr")[3:4]),
         legend = names(layer_colors), 
         col = layer_colors, 
         pch = 16,
         pt.cex = legend.cex,
         cex = 0.9,
         bty = "n",
         title = "Layer",
         title.adj = 0,               # 标题左对齐
         y.intersp = legend.spacing,  # 使用参数控制行间距
         text.width = max(strwidth(names(layer_colors))), # 统一文本宽度
         x.intersp = 0.8,
         text.col = "black",title.cex = 1.3,# 图例标题字号
         title.font = 2)  # 符号和文字的间距
  
  grDevices::recordPlot()
}




plotDominantCellType <- function(proportion_matrix, coordinates, 
                                 title = "Dominant Cell Type Annotations",
                                 point.size = 3.5, 
                                 cell_colors = c("EPL-IN" = "#c4d8e9", 
                                                 "GC" = "#ffe2ce", 
                                                 "M-TC" = "#bebebe", 
                                                 "OSNs" = "#ffadad",
                                                  "PGC" = "#d7bde2"),
                                 padding = 0.1,
                                 legend.cex = 1.2,
                                 legend.spacing = 1.5,
                                 legend.title.cex = 1.3,
                                 title.cex = 1.5,
                                 legend_add = F) {
  
  # 检查输入
  if (ncol(proportion_matrix) != nrow(coordinates)) {
    stop("Number of columns in proportion_matrix must match number of rows in coordinates")
  }
  
  # 确定每个spot的优势细胞类型
  dominant_types <- apply(proportion_matrix, 2, function(x) {
    rownames(proportion_matrix)[which.max(x)]
  })
  
  # 检查颜色是否覆盖所有细胞类型
  missing_types <- setdiff(unique(dominant_types), names(cell_colors))
  if (length(missing_types) > 0) {
    stop(paste("Missing color definitions for cell types:", paste(missing_types, collapse = ", ")))
  }
  
  # 创建颜色向量
  point_colors <- cell_colors[dominant_types]
  
  # 设置绘图参数
  par(bty = "o", lwd = 2, col = "lightgray", xpd = TRUE, 
      mar = c(3, 3, 4, 10))  # 上边距增大到4以适应大标题
  
  # 计算坐标范围
  x_range <- range(coordinates[, 1])
  y_range <- range(coordinates[, 2])
  x_padding <- diff(x_range) * padding
  y_padding <- diff(y_range) * padding
  
  # 绘制点图
  plot(coordinates, 
       pch = 16,
       cex = point.size,
       col = point_colors,
       main = "",
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty = "o",
       col.axis = "lightgray",
       lwd = 2,
       xlim = x_range + c(-x_padding, x_padding),
       ylim = y_range + c(-y_padding, y_padding))
  
  # 添加主标题（加大字号）
  title(main = title, line = 1, cex.main = title.cex, font.main = 2)
  names(cell_colors) <- gsub("\\."," ",names(cell_colors))
  
  
  if(legend_add){
    legend(
      x = par("usr")[2] * 1.02,      # 放在图形右侧（稍微超出右边界）
      y = par("usr")[4],             # 顶部对齐
      legend = names(cell_colors), 
      col = cell_colors, 
      pch = 16,
      horiz = T,                 # 垂直排列（默认）
      # ncol = 3,
      # nrow =2,                      # 分两列显示
      pt.cex = 3,                  # 加大符号大小（原1.3改为2.0）
      cex = 1.2,                     # 加大字体大小（原0.7改为1.2）
      bty = "n",                     # 无边框
      title = "Dominant cell types",
      title.adj = 0,                 # 标题左对齐
      title.cex = 1.3,               # 加大标题字号（原0.8改为1.3）
      title.font = 2,                # 标题加粗
      xpd = NA,                      # 允许在图外绘制
      xjust = 0,                     # 从x坐标向右绘制
      yjust = 1,                     # 从y坐标向下绘制（顶部对齐）
      text.width = NULL,             # 自动计算宽度
      x.intersp = 1,                 # 调整水平间距
      y.intersp = 1.5,               # 调整垂直间距（加大行距）
      text.col = "black"
    )
  }

  # # 添加图例（增大行间距）
  # legend(x = par("usr")[2] * 1.05,  # 稍微右移
  #        y = mean(par("usr")[3:4]),
  #        legend = names(layer_colors), 
  #        col = layer_colors, 
  #        pch = 16,
  #        pt.cex = legend.cex,
  #        cex = 0.9,
  #        bty = "n",
  #        title = "Layer",
  #        title.adj = 0,               # 标题左对齐
  #        y.intersp = legend.spacing,  # 使用参数控制行间距
  #        text.width = max(strwidth(names(layer_colors))), # 统一文本宽度
  #        x.intersp = 0.8,
  #        text.col = "black",title.cex = 1.3,# 图例标题字号
  #        title.font = 2)  # 符号和文字的间距
  
  grDevices::recordPlot()
}