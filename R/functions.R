#' @import raster
NULL

#' @include probpool_class.R
NULL


#########################################
# Constructor function for class Probpool
#' Create a probabilistic species pool
#' 
#' Description goes here
#' 
#' @param env_pool A Raster* object containing environmental probabilities
#' @param disp_pool A Raster* object containing dispersal probabilies
#' @param occurences A Raster* object containing occurrence data for individual species. Each species is represented by a binary RasterLayer, where cell values of 0 indicate absence and 1 indicates presence. 
#' @param interaction_matrix An interaction matrix
#' @param interaction_method A number specifying the method used for incorporating species interacitons in the calculation of the probabilistic species pool. (see details)
#' @return A raster stack containing species-specific probabilities of occurrence
#' @seealso \code{\link{Probpool}}
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

#' Create a probabilistic dispersal pool
#' 
#' @param disp_ability A numeric vector of dispersal abilities. [WHAT IS THE UNIT HERE?]
#' @param occurences A Raster* object containing occurrence data for individual species. Each species is represented by a binary RasterLayer, where cell values of 0 indicate absence and 1 indicates presence. 
#' @param conductance Optional. A Raster* object representing landscape conductance for dispersal. Takes values between 0 and 1, where 0 indicates complete resistance (\eqn{P(dispersal to raster cell) = 0}) and 1 indicates complete conductance (\eqn{P(dispersal to raster cell) = p}). If a single RasterLayer is provided, the conductance matrix will be recycled for each species. 
#' @param method A number specifying the kernel used for estimating dispersal distance. (see details)
#' @param longlat Logical. Should calculations be based on long/lat?
#' @return A raster stack containing species-specific probabilities of dispersal from areas of occurrence to a certain raster cell within a given time frame.
#' @seealso \code{\link{probpool}}
#' #' @export
disppool = function(disp_ability, occurrences, conductance, method=c("negexp","fattail"), longlat=TRUE) {
  # TODO: check extent of rasters, type of data etc.
  occurrences = raster::rasterToPoints(occurrence.surfaces)
  if (is.null(conductance)){
    distances = spDists(occurrences[,c(1:2)], occurrences[,c(1:2)], longlat = longlat)
    dispersal = lapply(1:dim(occurrence.surfaces)[3], function(x){
      dispersal.x = sapply(1:nrow(occurrences), function(y){
        calcD(occurrences[,x+2],distances[,y], ifelse(length(disp_ability)==1,disp_ability,disp_ability[x]), method = method[1])
      })
      dispersal.x.rst = occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] = dispersal.x
      return(dispersal.x.rst)
    })
  } else { # TODO add extra margin around conductance layers?
    conductance[is.na(conductance) | conductance < 0.01] = 0.01 # Replace NAs and very small values (incl. 0) by 0.01
    if (class(conductance)=="RasterLayer") {
      spec.trans = gdistance::transition(conductance, mean, 8) # create transition raster
      spec.trans = gdistance::geoCorrection(spec.trans, type="c")
      distances = gdistance::commuteDistance(spec.trans, occurrences[,c(1:2)])/2 # calculate commute distances
      distances = as.matrix(distances)
    }
    dispersal = lapply(1:dim(occurrence.surfaces)[3], function(x) {
      if (class(conductance) %in% c("RasterStack","RasterBrick")){
        spec.trans = transition(conductance[[x]], mean, 8) # create transition raster
        spec.trans = geoCorrection(spec.trans, type="c")
        distances = commuteDistance(spec.trans, occurrences[,c(1:2)])/2 # calculate commute distances
        distances = as.matrix(distances)
      }
      dispersal.x = sapply(1:nrow(occurrences), function(y){
        calcD(occurrences[,x+2],distances[,y], ifelse(length(disp_ability)==1,disp_ability,disp_ability[x]), method = method[1])
      })  # check 5000
      dispersal.x.rst = occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] = dispersal.x
      return(dispersal.x.rst)
    })
  }
  dispersal = stack(dispersal)      
  return(dispersal)
}