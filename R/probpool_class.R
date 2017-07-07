######################### CLASS DEFINITION ########################
#' @import methods
#' @import raster
NULL

#' Class "Prob_pool"
#' 
#' Core class of the \code{probpool} package. An object of class Prob_pool contains 
#' 
#' @slot "pools" Rasterbrick objects, containing 
#' @slot "interaction_matrix" A species by species matrix (may be asymmetric) with interactions directed from the species in the row to the species in the colums,
#' @slot "interaction_method" Method to be used for incorporating biotic interaction in the calculation of the probabilistic species pool
#' @slot "species_names" Names of species included in the Prob_pool
#' @slot "species_total" Total number of species included in the Prob_pool
#' @slot "species_mean" Mean number of species per raster cell
#' @slot "species_richness" Sum of Prob_pool@pools, representing the species richness of individual pools
#' @seealso \code{\link{probpool}}
#' @examples test_probpool = prob_pool(env_pool = env, 
#'  disp_pool = disp, occurrences = occ, interaction_matrix = int, interaction_method = 1)
setClass("Probpool",
         slots = c(pools = "list", 
                      interaction_matrix = "ANY", # can be empty
                      interaction_method = "character",
                      species_names = "character",
                      species_total = "numeric",
                      species_mean = "numeric",
                      species_richness = "list",
                      slots = "character"),
                  validity = prob_pool_check)

# Constructor function
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
      occurrences = raster::calc(occurrences, function(x){x[x>1] = 1; x}) # Convert abundance to occurrence
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
}

#####################################################################
######################### METHOD DEFINITIONS ########################
# print
setMethod(f = "print", signature = "Probpool", definition =  function(x){
  cat("Probabilistic species pool \n\n")
  cat(paste("Pools              : ", paste(names(x@pools), collapse = ", "), sep = ""), "\n")
  cat(paste("Species (total)    : ", x@species_total, "\n", sep = ""))
  cat(paste("Species (mean)     : ", x@species_mean, "\n", sep = ""))
  cat(paste("Interaction method : ", x@interaction_method, "\n", sep = ""))
  cat(paste("Resolution         : ", paste(round(raster::res(x@pools$prob_pool), 3), collapse = " x "), " (x,y)\n", 
                                     sep = ""))
  cat(paste("Extent             : ", paste(round(x@pools$prob_pool@extent[1:4], 3), 
                                           collapse = ", "), " (xmin, xmax, ymin, ymax)", sep = ""))
  return(NULL)
})
#' @rdname summary
#' @export
print = function(x) UseMethod("print", x)

#-------------------------------------------------------------------------------------------
# summary
setMethod("summary", "Probpool",  function(object){
  print("!AFAPSF")
  print(object)
  smry = list(pools = paste(names(object@pools), collapse = ", "),
       species_total = object@species_total,
       species_mean =  object@species_mean,
       interaction_method = object@interaction_method,
       resolution = round(res(object@pools$prob_pool), 3),
       extent = round(object@pools$prob_pool@extent[1:4], 3)
       )
  return(smry)
})
#' Summarizing a Probpool object
#' 
#' Prints and returns a \code{list} of properties for an object of class \code{Probpool}
#' @param x An object of class \code{Probpool}
#' @return A \code{list} of properties of \code{object}
#' @export
summary <- function(x) UseMethod("summary", x)

#-------------------------------------------------------------------------------------------
# plot
setMethod("plot", "Probpool", function(x, focal_species = NULL, focal_unit = NULL, ...){
  par_old = par()
  on.exit(par(par_old))
  moreargs = eval(substitute(list(...)))
  if(is.null(focal.species) && is.null(focal.unit)){
    # Plot species richness maps
    color.theme <- rasterVis::rasterTheme(region = rev(terrain.colors(100)))
    richness.maps = raster::stack(x@species_richness[!is.na(x@species_richness)])
    do.call(rasterVis::levelplot, c(x = quote(richness.maps), par.settings = quote(color.theme), moreargs))
  } else if(!is.null(focal.species) && is.null(focal.unit)){
    # Plot probability pools for focal species
    if(!focal.species %in% x@species_names){stop("Species not found.")}
    color.theme <- rasterVis::rasterTheme(region = rev(terrain.colors(100)))
    species.maps = raster::stack(lapply(x@pools, FUN = "[[", i = focal.species))
    do.call(rasterVis::levelplot, c(x = quote(species.maps), par.settings = quote(color.theme), moreargs))
  } else if(is.null(focal.species) && !is.null(focal.unit)){
    # Plot probability distributions for focal unit
    focal.probs = sapply(x@pools[names(x@pools) != "occurrences"], raster::extract, y = focal.unit)
    focal.pdf = apply(focal.probs, 2, poinbin::ppoibin, kk = 1:x@species_total)
    focal.pdf = reshape2::melt(focal.pdf)
    focal.pdf$type = "Probability density"
    focal.cdf = apply(focal.probs, 2, poinbin::dpoibin, kk = 1:x@species_total)
    focal.cdf = reshape2::melt(focal.cdf)
    focal.cdf$type = "Cumulative probability density"
    
    focal.all = rbind(focal.pdf, focal.cdf)
    lattice::xyplot(value ~ Var1 | factor(Var2) + factor(type), data = focal.all, ylab = "probability", xlab= "species", type = "l")
  } else {
    stop("Please provide only one argument (species/focal.unit)")
  }
})

#' Plot a Probpool object. 
#' 
#' Method to plot and visualize a \code{Probpool} object. Depending on the arguments provided, different plotting behaviour (see Details). Additional parameters are passed to the  the '...' argument
#'
#' @param x An object of class \code{Probpool}
#' @param focal.species The name of a particular species in x. If provided, probability pools for that ...
#' @param focal.unit An object comapatible with the \code{\link[raster]{raster::extract}} function, i.e 
#' \itemize{
#'   \item points represented by a two-column matrix or data.frame
#'   \item a numeric vector representing cell numbers
#'   \item a SpatialPoints*; SpatialPolygons*; SpatialLines object
#'   \item an \code{\link[raster]{raster::Extent}} object
#' }
#' @param ... Additional arguments for \code{\link[lattice]{lattice::levelplot}} or \code{\link[lattice]{lattice::xyplot}}
#' @return None
#' @examples
#' my_prob_pool = prob_pool(env_pool = env, disp_pool = disp)
#' plot(my_prob_pool)
#' plot(my_prob_pool, focal.species = "Olea europaea")
#' plot(my_prob_pool, focal.unit = 132)) # Cell number
#' plot(my_prob_pool, focal.unit = c(24,26))) # Cell index
#' plot(my_prob_pool, focal.unit = extent(c(7,8,49,53))) # Object of class raster::Extent
#' @export
plot = function(x, focal_species = NULL, focal_unit = NULL, ...) UseMethod("plot")
