# Validity check function
prob_pool_check = function(object){
  if(all(sapply(object@pools, is.null))){ 
    return("Please provide at least one of the following arguments: env_pool, disp_pool, occurrences")
  }
  if(is.null(object@interaction_matrix) && is.null(object@pools$env_pool) && is.null(object@pools$disp_pool)){
    return("interaction_matrix is missing.")
  }
  if(is.null(object@interaction_method)){
    return("Unknown interaction_method. Choose '1' for modification or '2' for multiplication")
  }
  if(!all(sapply(object@pools, function(pool) {extends(class(pool), "Raster") || is.null(pool)}))){ # check types
    return("Invalid argument. Please provide a raster object.")
  }
  pool.dims = matrix(sapply(object@pools, dim)) # Check dimensions
  if(any(apply(pool.dims, 1, function(x) length(unique(x)) != 1))){
    return("All raster objects need to have the same dimensions.")
  }
  if(!all(sapply(object@pools, function(x){is.null(x) || all.equal(names(object@pools$prob_pool), names(x))}))){ # check species names
    return("Species names do not match.")
  }
  if(object@interaction_method == 2 && !is.null(object@interaction_matrix) && is.null(object@pools$occurrences)){
    return("Multiplication approach (interaction_method = 2) requires species occurrences.")
  }
  # TODO: Check for abundance data and rescale if prob > 1
  return(TRUE)
}


#####################################################################
########################## HELPER FUNCTIONS #########################
multiply.pools = function(env.pool, disp.pool){
  if(is.null(env.pool)){
    return(disp.pool)
  } else if(is.null(disp.pool)){
    return(env.pool)
  } else {
    return(env.pool * disp.pool)
  }
}

calc.prob = function(probabilities, interaction.matrix, interaction.method, occurrences = NULL){
  interaction.matrix = as.matrix(interaction.matrix) # in case of dist object being provided
  if(!is.null(occurrences)){
    tmp.probs = occurrences # Only exists when no env/disp pool present --> calc interactions from occurrences
  } else {
    tmp.probs = probabilities # calc interactions from probabilities
  }
  
  # Create "interaction pool"
  interaction.pool = calc(tmp.probs, function(prob.cell){
    sapply(1:length(prob.cell), function(species.index){
      prob.cell[species.index] * mean((prob.cell * interaction.matrix[,species.index]))
    })
  })
  
  # Calculate prob.pool
  if(interaction.method == 1){
    int.pos = calc(interaction.pool, fun = function(x){x[x <= 0] = NA; x}) # only positive local interactions
    int.neg = calc(interaction.pool, fun = function(x){x[x > 0] = NA; x})  # only negative local interactions
    for(i in 1:nlayers(probabilities)){
      tmp = probabilities[[i]]
      tmp.pos = tmp + (1-tmp) * int.pos[[i]]
      tmp.neg = tmp + (tmp * int.neg[[i]])
      tmp[!is.na(tmp.pos)] = tmp.pos[!is.na(tmp.pos)] 
      tmp[!is.na(tmp.neg)] = tmp.neg[!is.na(tmp.neg)]
      probabilities[[i]] = tmp
    }
  } else if(interaction.method == 2){
    warning("This interaction method is outdated. Results are not interpretable as probabilies anymore")
    probabilities = probabilities * ((interaction.pool + 1) / 2)
  }
  return(probabilities)
}

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


# occurence.surfaces needs to be a raster stack including rasters of species occurances or abundances with values of 0 for absences and values > 0 for occurances. Abundances are not considered though.

# disp.ability needs to be a vector of dispersal abilities. Scaling range? Larger values = higher dispersal abilities

# cond.surfaces needs to be a raster stack including rasters of conductance values between 0 and one
# NAs will be replaced by 0 and 0s will be replaced by small value; set variable for this?


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


# allD is a wrapper function for calcD that applies the calculation
# of calcD to all cells and all species. It takes the full occupancy
# matrix, the full pairwise distance matrix and a vector of each
# species' dispersal ability. As with calcD, you specify the dispersal
# function by specifying "method"
allD = function(occupancy,distance,k,method = "negexp")
{
  Pd = sapply(1:ncol(occupancy),function(i){sapply(1:nrow(occupancy), function(j){calcD(occupancy[,i],distance[j,],k[i],method)})})
  return(Pd)
}

# allE is a function that calculates the environmental suitability
# for all species across all sites. There are two methods available:
# 1. consider the distance between a focal site and the closest (in 
# environmental space) occupied site, and rank it against all distances
# 2. use beals smoothing
calcE = function(occupancy, envdist,site)
{
  index = which(occupancy > 0)
  es = 1-length(which(envdist[,site] < min(envdist[index,site])))/256
  return(es)
}

allE = function(occupancy,environment=NULL, method)
{
  if(method == "mindist")
  {
    dists = as.matrix(dist(environment))
    Pe = sapply(1:ncol(occupancy),function(i){sapply(1:nrow(occupancy), function(j){calcE(d[,i],dists,j)})})
  }
  if(method == "beals")
  {
    require(vegan)
    Pe = beals(occupancy)	
  }
  return(Pe)
}

#TODO: Function to convert species by sites to raster stack

# occurence.surfaces needs to be a raster stack including rasters of species occurences or abundances with values of 0 for absences and values > 0 for occurances. Values wil be scaled to range from 0 to 1
# # or the disp.pool, the env.pool or the product of both
# 
# # interaction.matrix is a species by species matrix (may be asymmetric) with interactions assumed to be directed from the species in the row to the species in the colums
# 
# bio_pool <- function(occurrences, interaction.matrix, abundance=TRUE) {
#   occurrences <- values(occurences)
#   occurrences <- occurrences[complete.cases(occurrences),]
#   
#   if(abundance){ 
#     occurrences <- occurrences/max(occurrences)
#   } else {
#     occurrences[occurrences > 0] <- 1 
#   }
#   
#   # multiply the incoming interactions of each species x (columns in int.matrix)
#   # with the occurrence/probability of all other species for the given site y
#   interactions <- lapply(1:dim(occurences)[3], function(x) {
#     interactions.x <- t(sapply(1:nrow(occurrences), function(y) occurrences[y,] * interaction.matrix[,x]))
#     interactions.x <- rowMeans(interactions.x)
#     interactions.x.rst <- occurences[[x]]
#     interactions.x.rst[!is.na(interactions.x.rst)] <- interactions.x
#     return(interactions.x.rst)
#   })
#   
#   interactions <- stack(interactions)     
#   return(interactions)
# }


