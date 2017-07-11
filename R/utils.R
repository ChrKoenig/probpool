########################## HELPER FUNCTIONS #########################
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

mult_pools = function(pool1, pool2){
  if(is.null(pool1)){
    return(pool2)
  } else if(is.null(pool2)){
    return(pool1)
  } else {
    return(pool1 * pool2)
  }
}

calc_prob = function(probabilities, interaction_matrix, interaction_method, occurrences = NULL){
  interaction_matrix = as.matrix(interaction_matrix) # in case of dist object being provided
  if(!is.null(occurrences)){
    tmp_probs = occurrences # Only exists when no env/disp pool present --> calc interactions based on occurrences
  } else {
    tmp_probs = probabilities # calculate interactions based on probabilities
  }
  
  # Create "interaction pool"
  interaction_pool = raster::calc(tmp_probs, fun = function(prob_cell){
    sapply(1:length(prob_cell), function(species_index){
      prob_cell[species_index] = mean((prob_cell * interaction_matrix[,species_index]))
    })
  })
  
  # Calculate prob.pool
  if(interaction_method == 1){
    int_pos = raster::calc(interaction_pool, fun = function(x){x[x <= 0] = NA; x}) # only positive local interactions
    int_neg = raster::calc(interaction_pool, fun = function(x){x[x > 0] = NA; x})  # only negative local interactions
    for(i in 1:raster::nlayers(probabilities)){
      tmp = probabilities[[i]]
      tmp_pos = tmp + (1-tmp) * int_pos[[i]]
      tmp_neg = tmp + (tmp * int_neg[[i]])
      tmp[!is.na(tmp_pos)] = tmp_pos[!is.na(tmp_pos)] 
      tmp[!is.na(tmp_neg)] = tmp_neg[!is.na(tmp_neg)]
      probabilities[[i]] = tmp
    }
  } else if(interaction_method == 2){
    warning("This interaction method is deprecated. Results are not interpretable as probabilies anymore")
    probabilities = probabilities * ((interaction_pool + 1) / 2)
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

##############################################
# allD is a wrapper function for calcD that applies the calculation
# of calcD to all cells and all species. It takes the full occupancy
# matrix, the full pairwise distance matrix and a vector of each
# species' dispersal ability. As with calcD, you specify the dispersal
# function by specifying "method"
allD = function(occupancy,distance,k,method = "negexp"){
  Pd = sapply(1:ncol(occupancy),function(i){sapply(1:nrow(occupancy), function(j){calcD(occupancy[,i],distance[j,],k[i],method)})})
  return(Pd)
}
##############################################
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

allE = function(occupancy,environment=NULL, method){
  if(method == "mindist"){
    dists = as.matrix(dist(environment))
    Pe = sapply(1:ncol(occupancy),function(i){sapply(1:nrow(occupancy), function(j){calcE(d[,i],dists,j)})})
  }
  if(method == "beals"){
    Pe = beals(occupancy)	
  }
  return(Pe)
}