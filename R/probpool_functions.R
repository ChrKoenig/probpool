disp.pool <- function(occurrence.surfaces, disp.ability, method=c("negexp","fattail"),cond.surfaces=NULL, longlat=TRUE) {
  # TODO: check extent of rasters, type of data etc.
  occurrences <- rasterToPoints(occurrence.surfaces)
  
  if (is.null(cond.surfaces)){
    
    distances <- spDists(occurrences[,c(1:2)], occurrences[,c(1:2)], longlat = longlat)
    
    dispersal <- lapply(1:dim(occurrence.surfaces)[3], function(x) {
      
      dispersal.x <- sapply(1:nrow(occurrences), function(y) calcD(occurrences[,x+2],distances[,y], ifelse(length(disp.ability)==1,disp.ability,disp.ability[x]), method = method[1]))
      dispersal.x.rst <- occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] <- dispersal.x
      return(dispersal.x.rst)
    })
    
  } else {
    
    # TODO add extra margin around conductance layers?
    
    # replace NAs by 0
    cond.surfaces[is.na(cond.surfaces)] <- 1/100
    # truncation: replace <1/100 by 0                                   
    cond.surfaces[cond.surfaces<1/100] <- 1/100
    
    if (class(cond.surfaces)=="RasterLayer") {
      # create transition raster
      spec.trans <- transition(cond.surfaces, mean, 8)
      # geocorrection
      spec.trans <- geoCorrection(spec.trans, type="c")
      
      # calculate commute distances
      distances <- commuteDistance(spec.trans, occurrences[,c(1:2)])/2
      distances <- as.matrix(distances)
    }
    
    
    dispersal <- lapply(1:dim(occurrence.surfaces)[3], function(x) {
      
      if (class(cond.surfaces) %in% c("RasterStack","RasterBrick")) {
        # create transition raster
        spec.trans <- transition(cond.surfaces[[x]], mean, 8)
        # geocorrection
        spec.trans <- geoCorrection(spec.trans, type="c")
        
        # calculate commute distances
        distances <- commuteDistance(spec.trans, occurrences[,c(1:2)])/2
        distances <- as.matrix(distances)
      }
      
      dispersal.x <- sapply(1:nrow(occurrences), function(y) calcD(occurrences[,x+2],distances[,y], ifelse(length(disp.ability)==1,disp.ability,disp.ability[x]), method = method[1]))  # check 5000
      dispersal.x.rst <- occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] <- dispersal.x
      return(dispersal.x.rst)
    })
  }
  
  dispersal <- stack(dispersal)      
  
  return(dispersal)
}