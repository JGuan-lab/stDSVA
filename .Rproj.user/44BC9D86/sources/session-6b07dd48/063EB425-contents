


simulateSpotsSpatialSc <- function(sp.exp, sp.meta, sp.loc, window = 500, CoordinateX = 'X', CoordinateY = 'Y', CoordinateZ = NULL) {
  library(reticulate)
  source_python("H:/LinYifan/spatialDSSC-BE/R/simulate.py")  
  pd <- import("pandas")
  sp.exp <- pd$DataFrame(sp.exp,index = sp.meta$cellID)
  sp.meta <- pd$DataFrame(sp.meta,index = sp.meta$cellID)
  sp.loc <- pd$DataFrame(sp.loc,index = sp.meta$cellID)


  resultsCombine = py$Simulated(sp.exp, sp.meta, sp.loc, CoordinateX, CoordinateY, window, CoordinateZ) 
  return(resultsCombine)
  
  
  
}



