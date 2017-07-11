######################### CLASS DEFINITION ########################
#' @import methods
#' @import raster
NULL

#' @include utils.R
NULL

#' "Prob_pool" Class definition
#' 
#' Core class of the \code{probpool} package.
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
                   slots = "character"))

setValidity("Probpool", prob_pool_check)

#####################################################################
######################### METHOD DEFINITIONS ########################
print.probpool = function(x){
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
}
#' @rdname summary
#' @export
setMethod("print",
          signature(x = "Probpool"),
          function(x){print.probpool(x)})

show.probpool = function(object){
  print.probpool(object)
}
#' @rdname summary
#' @export
setMethod("show",
          signature(object = "Probpool"),
          function(object){show.probpool(object)})

summary.probpool = function(object,...){
  smry = list(pools = names(object@pools),
              species_total = object@species_total,
              species_mean =  object@species_mean,
              interaction_method = object@interaction_method,
              resolution = round(raster::res(object@pools$prob_pool), 3),
              extent = round(object@pools$prob_pool@extent[1:4], 3)
  )
  return(smry)
}

#' Summarize a Probpool object.
#' 
#' Summary functions for an object of class Probpool. Prints and returns a \code{list} of properties for an object of class \code{Probpool}
#' 
#' @param x An object of class \code{Probpool}
#' @return A \code{list} of properties of \code{object}
#' @export
setMethod("summary",
          signature(object = "Probpool"),
          function(object,...){summary.probpool(object,...)})

#-------------------------------------------------------------------------------------------
#' Plot a Probpool object. 
#' 
#' Method to plot and visualize a \code{Probpool} object. Depending on the arguments provided, different plotting behaviour (see Details). Additional parameters are passed to the  the '...' argument
#'
#' @param x An object of class \code{Probpool}
#' @param focal_species The name of a particular species in x. If provided, probability pools for that ...
#' @param focal_unit An object comapatible with the \code{\link[raster]{raster::extract}} function, i.e 
#' \itemize{
#'   \item points represented by a two-column matrix or data.frame
#'   \item a numeric vector representing cell numbers
#'   \item a SpatialPoints*; SpatialPolygons*; SpatialLines object
#'   \item an \code{\link[raster]{raster::Extent}} object
#' }
#' @param ... Additional arguments for \code{\link[lattice]{lattice::levelplot}} or \code{\link[lattice]{lattice::xyplot}}
#' @return None
#' @examples
#' my_prob_pool = probpool(env_pool = env, disp_pool = disp)
#' plot(my_prob_pool)
#' plot(my_prob_pool, focal_species = "Olea europaea")
#' plot(my_prob_pool, focal_unit = 132)) # Cell number
#' plot(my_prob_pool, focal_unit = c(24,26))) # Cell index
#' plot(my_prob_pool, focal_unit = extent(c(7,8,49,53))) # Object of class raster::Extent
#' @export
plot.probpool = function(x, focal_species = NULL, focal_unit = NULL, ...){
  par_old = par()
  on.exit(suppressWarnings(par(par_old)))
  moreargs = eval(substitute(list(...)))
  if(is.null(focal_species) && is.null(focal_unit)){
    # Plot species richness maps
    color_theme <- rasterVis::rasterTheme(region = rev(terrain.colors(100)))
    richness_maps = raster::stack(x@species_richness[!is.na(x@species_richness)])
    do.call(rasterVis::levelplot, c(x = quote(richness_maps), par.settings = quote(color_theme), moreargs))
  } else if(!is.null(focal_species) && is.null(focal_unit)){
    # Plot probability pools for focal species
    if(!focal_species %in% x@species_names){stop("Species not found.")}
    color_theme <- rasterVis::rasterTheme(region = rev(terrain.colors(100)))
    species_maps = raster::stack(lapply(x@pools, FUN = "[[", i = focal_species))
    do.call(rasterVis::levelplot, c(x = quote(species_maps), par.settings = quote(color_theme), moreargs))
  } else if(is.null(focal_species) && !is.null(focal_unit)){
    # Plot probability distributions for focal unit
    focal_probs = sapply(x@pools[names(x@pools) != "occurrences"], raster::extract, y = focal_unit)
    focal_pdf = apply(focal_probs, 2, poibin::ppoibin, kk = 1:x@species_total)
    focal_pdf = reshape2::melt(focal_pdf)
    focal_pdf$type = "Probability density"
    focal_cdf = apply(focal_probs, 2, poibin::dpoibin, kk = 1:x@species_total)
    focal_cdf = reshape2::melt(focal_cdf)
    focal_cdf$type = "Cumulative probability density"
    
    focal_all = rbind(focal_pdf, focal_cdf)
    do.call(lattice::xyplot, c(x = formula(value ~ Var1 | factor(Var2) + factor(type)), 
                               data = quote(focal_all), type = quote("l"), moreargs))
  } else {
    stop("Please provide only one argument (species/focal_unit)")
  }
}
setMethod("plot",
          signature(x = "Probpool"),
          function(x){plot.probpool(x)})