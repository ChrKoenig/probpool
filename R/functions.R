#' @import raster
NULL

#' @include probpool_class.R
NULL


#########################################
# Constructor function for Probpool
#' Create a Probpool object
#' 
#' Constructor function for class \code{Probpool}
#' 
#' [DETAILS GO HERE]
#' 
#' @param env_pool A Raster* object containing environmental probabilities
#' @param disp_pool A Raster* object containing dispersal probabilies
#' @param occurences A Raster* object containing occurrence data for individual species. Each species is represented by a binary RasterLayer, where cell values of 0 indicate absence and 1 indicates presence. 
#' @param interaction_matrix An interaction matrix
#' @param interaction_method A number specifying the method used for incorporating species interacitons in the calculation of the probabilistic species pool. (see details)
#' @return An object of class \code{Probpool}
#' @seealso \code{\link{disppool}}, \code{\link{Probpool-class}}
#' @references Karger, D. N. et al. (2016), Delineating probabilistic species pools in ecology and biogeography. Global Ecology and Biogeography, 25: 489-501. doi:10.1111/geb.12422
#' @examples example_probpool = probpool(env_pool = example1_env, disp_pool = example1_disp)
#' @export
probpool = function(env_pool = NULL, disp_pool = NULL, occurrences = NULL,
                    interaction_matrix = NULL, interaction_method = 1){
  prob_pool = NULL
  if(is.null(interaction_matrix) && (!is.null(env_pool) || !is.null(disp_pool))){ # No interactions, multiply env * disp
    prob_pool = mult_pools(env_pool, disp_pool)
    interaction_method = 3
  } else { # Interactions present
    if(!is.null(env_pool) || !is.null(disp_pool)){ # Base probabilities from env/disp layer
      prob_pool_raw = mult_pools(env_pool, disp_pool)
      prob_pool = calc_prob(prob_pool_raw, interaction_matrix, interaction_method)
    } else { # Estimate base probabilities from occurence layer
      # occurrences = raster::calc(occurrences, function(x){x[x>1] = 1; x}) # Convert abundance to occurrence
      n.total = raster::nlayers(occurrences)
      n.mean = mean(raster::values(sum(occurrences)), na.rm = T)
      base.prob = n.mean/n.total
      prob_pool_raw = raster::calc(occurrences, function(x){ # Assign uniform distribution to raster
        x[!is.na(x)] = base.prob
        return(x)
      }) 
      prob_pool = calc_prob(prob_pool_raw, interaction_matrix, interaction_method, occurrences)
    }
  } 
  
  # Create Probpool object
  pools = list(occurrences = occurrences,
               env_pool = env_pool, 
               disp_pool = disp_pool, 
               prob_pool = prob_pool)
  pools = pools[!sapply(pools, is.null)] # remove empty pools
  species_richness = lapply(pools, FUN = function(x) sum(x))
  result = new("Probpool", 
               pools = pools,
               interaction_matrix = interaction_matrix,
               interaction_method = c("Modification", "Multiplication", "None")[interaction_method],
               species_names = names(prob_pool),
               species_total = length(names(prob_pool)),
               species_mean = round(mean(raster::values(sum(prob_pool)), 1, na.rm = T)),
               species_richness = species_richness,
               slots = c("pools","interaction_matrix","interaction_method", "species_names", "species_total", "species_mean", "species_richness", "slots"))
  class(result) = "Probpool"
  return(result)
}

#' Create a probabilistic dispersal pool.
#' 
#' Creates a RasterBrick of dispersal probabilities that can be used as argument in \code{\link{probpool}}.
#' 
#' [DETAILS GO HERE]
#' 
#' @param disp_ability A numeric vector of dispersal abilities. [WHAT IS THE UNIT HERE?]
#' @param occurence.surfaces A Raster* object containing occurrence data for individual species. Each species is represented by a binary RasterLayer, where cell values of 0 indicate absence and 1 indicates presence. 
#' @param cost.surfaces Optional. A Raster* object representing landscape conductance for dispersal. Takes values between 0 and 1, where 0 indicates complete resistance (\eqn{P(dispersal to raster cell) = 0}) and 1 indicates complete conductance (\eqn{P(dispersal to raster cell) = p}). If a single RasterLayer is provided, the conductance matrix will be recycled for each species. 
#' @param method A number specifying the kernel used for estimating dispersal distance. (see details)
#' @param longlat Logical. Should calculations be based on long/lat?
#' @param transitionFunction Function. transitionFunction in gdistance::transition
#' @return A raster stack containing species-specific probabilities of dispersal from areas of occurrence to a certain raster cell within a given time frame.
#' @seealso \code{\link{probpool}}
#' @references Tamme, R., Götzenberger, L., Zobel, M., Bullock, J. M., Hooftman, D. A. P., Kaasik, A. and Pärtel, M. (2014), Predicting species' maximum dispersal distances from simple plant traits. Ecology, 95: 505–513. doi:10.1890/13-1000.1
#' @examples [EXAMPLE GOES HERE]
#' @export
disppool = function(disp_ability, occurrence.surfaces, cost.surfaces=NULL, method=c("negexp","fattail"), longlat=TRUE, transitionFunction=function(x) 1/mean(x)) { # method only one attribute?
  # TODO?: check extent of rasters, type of data etc.
  
  # Dispersal ability k is the distance a species can disperse in the given time frame
  # If longlat = TRUE the unit for k is km, if longlat= FALSE it is the unit of the raster coordinate system (for example no of cells).
  # If cost.surfaces or a single cost.surface are/is supplied, longlat is ignored and the output is scaled to no of cells, so k should be in number of cells.
  
  
  # calcD is a function that calculates the 
  # dispersal probabilities for a species to 
  # all locations on the landscale from a focal
  # cell. It takes a vector containing the species 
  # presence/absence information, a vector containing
  # the distances from a focal cell to others
  # and a constant describing the species' dispersal
  # ability. It returns a single probability value
  # describing the chance of dispersing into that cell
  # It allows two methods of calculating the 
  # dispersal kernel (negative exponential 
  # and fat-tailed)
  
  
  calcD = function(occupancy, distance, k, method = "negexp")
  {
    index = which(occupancy > 0)
    if(method == "negexp") {distFunction = function(d,k){1 - prod(1-exp(-1*d/k))}}
    if(method == "fattail") {distFunction = function(d,k){1 - prod(1/d^k)}}
    return(distFunction(distance[index],k))
  }
  
  occurrences = raster::rasterToPoints(occurrence.surfaces)
  if (is.null(cost.surfaces)){
    distances = sp::spDists(occurrences[,c(1:2)], occurrences[,c(1:2)], longlat = longlat)
    dispersal = lapply(1:dim(occurrence.surfaces)[3], function(x){
      dispersal.x = sapply(1:nrow(occurrences), function(y){
        calcD(occurrences[,x+2],distances[,y], ifelse(length(disp_ability)==1,disp_ability,disp_ability[x]), method = method[1])
      })
      dispersal.x.rst = occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] = dispersal.x
      return(dispersal.x.rst)
    })
  } else { 
    if (class(cost.surfaces)=="RasterLayer") {
      spec.trans = gdistance::transition(cost.surfaces, transitionFunction, 8) # create transition raster
      spec.trans = gdistance::geoCorrection(spec.trans, type="c", scl=TRUE)
      distances = gdistance::costDistance(spec.trans, occurrences[,c(1:2)], occurrences[,c(1:2)]) # calculate cost distance
      distances = as.matrix(distances)
    }
    dispersal = lapply(1:dim(occurrence.surfaces)[3], function(x) {
      if (class(cost.surfaces) %in% c("RasterStack","RasterBrick")){
        spec.trans = gdistance::transition(cost.surfaces[[x]], transitionFunction, 8) # create transition raster
        spec.trans = gdistance::geoCorrection(spec.trans, type="c", scl=TRUE)
        distances = gdistance::costDistance(spec.trans, occurrences[,c(1:2)], occurrences[,c(1:2)]) # calculate cost distance
        distances = as.matrix(distances)
      }
      dispersal.x = sapply(1:nrow(occurrences), function(y){
        calcD(occurrences[,x+2],distances[,y], ifelse(length(disp_ability)==1,disp_ability,disp_ability[x]), method = method[1])
      })
      dispersal.x.rst = occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] = dispersal.x
      return(dispersal.x.rst)
    })
  }
  dispersal = raster::stack(dispersal)      
  return(dispersal)
}