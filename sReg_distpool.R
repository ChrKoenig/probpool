load("dispersal_ability.RData")
load("occ_rst_stack.RData")
load("suit_rst_stack.RData")

library(raster); library(gdistance)


##############################################
# calcD is a function that calculates the 
# dispersal probabilities for a species to 
# all locations on the landscale from a focal
# cell. It takes a vector containing the species 
# presence/absence information, a vector containing
# the distances from a focal cell to other others
# and a constant describing the species' dispersal
# ability. It returns a single probability value
# describing the chance of dispersing into that cell

# Also, allow two methods of calculating the 
# dispersal kernel (negative exponential 
# and fat-tailed)

calcD = function(occupancy, distance, k, method = "negexp")
{
  index = which(occupancy > 0)
  if(method == "negexp") {distFunction = function(d,k){1 - prod(1-exp(-1*d/k))}}
  if(method == "fattail") {distFunction = function(d,k){1 - prod(1/d^k)}}
  return(distFunction(distance[index],k))
}




disp_pool <- function(occurrence.surfaces, disp.ability=NULL, cost.surfaces=NULL) {
  
  occurrences <- rasterToPoints(occurrence.surfaces)
  
  if (is.null(cost.surfaces)){
    
    distances <- spDists(occurrences[,c(1:2)], occurrences[,c(1:2)], longlat = TRUE)
    
    dispersal <- lapply(1:dim(occurrence.surfaces)[3], function(x) {
      
      dispersal.x <- sapply(1:nrow(occurrences), function(y) calcD(occurrences[,x+2],distances[,y]/50, ifelse(is.null(disp.ability),1,disp.ability[x]), method = "negexp"))  # check 50
      dispersal.x.rst <- occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] <- dispersal.x
      return(dispersal.x.rst)
    })
    
  } else {
    
    # replace NAs by 0
    cost.surfaces[is.na(cost.surfaces)] <- 1/100
    # truncation: replace NAs by 0
    cost.surfaces[cost.surfaces<1/100] <- 1/100
    
    dispersal <- lapply(1:dim(occurrence.surfaces)[3], function(x) {
      
      # create transition raster
      spec.trans <- transition(cost.surfaces[[x]], mean, 8)
      # geocorrection
      spec.trans <- geoCorrection(spec.trans, type="c")
      
      # calculate commute distances
      distances <- commuteDistance(spec.trans, occurrences[,c(1:2)])
      distances <- as.matrix(distances)
      
      dispersal.x <- sapply(1:nrow(occurrences), function(y) calcD(occurrences[,x+2],distances[,y]/5000, ifelse(is.null(disp.ability),1,disp.ability[x]), method = "negexp"))  # check 5000
      dispersal.x.rst <- occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] <- dispersal.x
      return(dispersal.x.rst)
    })
  }
  
  
  dispersal <- stack(dispersal,layers=names(occurrence.surfaces))      
  
  return(dispersal)
}


par(mfrow=c(2,2))

disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack[[1:5]])
plot(disp.rst.stack[[5]])

disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack[[1:5]], cost.surfaces=suit.rst.stack[[1:5]])
plot(disp.rst.stack[[5]])

disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack[[1:5]], disp.ability=dispersal.ability[1:5])
plot(disp.rst.stack[[5]])

disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack[[1:5]], disp.ability=dispersal.ability[1:5], cost.surfaces=suit.rst.stack[[1:5]])
plot(disp.rst.stack[[5]])






      