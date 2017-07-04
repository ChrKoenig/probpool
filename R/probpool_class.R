library(methods)
library(raster)
library(rasterVis)
library(latticeExtra)
library(poibin)
library(reshape2)

###################################################################
# Probabilistic species pool class definition
# author: Christian König
# date:   June 2017
#
# current issues:
# - no abundance support in occurrences data
# - not thoroughly tested
# - disp.pool creation should be included in the function call instead of a separate function
# 
###################################################################
######################### CLASS DEFINITION ########################
# Validity check function
is.valid.prob.pool = function(object){
  if(all(sapply(object@pools, is.null))){ 
    return("Please provide at least one of the following arguments: env.pool, disp.pool, occurrences")
  }
  if(is.null(object@interaction.matrix) && is.null(object@pools$env.pool) && is.null(object@pools$disp.pool)){
    return("Interaction.matrix is missing.")
  }
  if(is.na(object@interaction.method)){
    return("Unknown interaction.method. Choose '1' for modification or '2' for multiplication")
  }
  if(!all(sapply(object@pools, function(pool) {extends(class(pool), "Raster") || is.null(pool)}))){ # check types
    return("Invalid argument. Please provide a raster object.")
  }
  pool.dims = matrix(sapply(object@pools, dim)) # Check dimensions
  if(any(apply(pool.dims, 1, function(x) length(unique(x)) != 1))){
    return("All raster objects need to have the same dimensions.")
  }
  if(!all(sapply(object@pools, function(x){is.null(x) || all.equal(names(object@pools$prob.pool), names(x))}))){ # check species names
    return("Species names do not match.")
  }
  if(object@interaction.method == 2 && !is.null(object@interaction.matrix) && is.null(object@pools$occurrences)){
    return("Multiplication approach (interaction.method = 2) requires species occurrences.")
  }
  # TODO: Check for abundance data and rescale if prob > 1
  return(TRUE)
}

# Class skeleton
setClass("prob.pool",
         slots = c("pools", 
                   "interaction.matrix",
                   "interaction.method",
                   "species.names",
                   "species.total",
                   "species.mean",
                   "species.richness",
                   "slots"),
         validity = is.valid.prob.pool)


# Constructor function
prob.pool = function(env.pool = NULL, disp.pool = NULL, occurrences = NULL,
                     interaction.matrix = NULL, interaction.method = 1){
  
  prob.pool = NULL
  if(is.null(interaction.matrix) && (!is.null(env.pool) || !is.null(disp.pool))){ # No interactions, multiply env * disp
    prob.pool = multiply.pools(env.pool, disp.pool)
    interaction.method = NA
  } else { # Interactions present
    if(!is.null(env.pool) || !is.null(disp.pool)){ # Base probabilities from env/disp layer
      prob.pool.raw = multiply.pools(env.pool, disp.pool)
      prob.pool = calc.prob(prob.pool.raw, interaction.matrix, interaction.method)
    } else { # Estimate base probabilities from occurence layer
      occurrences = calc(occurrences, function(x){x[x>1] = 1; x}) # Convert abundance to occurrence
      n.total = dim(occurrences)[3]
      n.mean = mean(raster::values(sum(occurrences)), na.rm = T)
      base.prob = n.mean/n.total
      prob.pool.raw = calc(occurrences, function(x){ # Assign uniform distribution to raster
        x[!is.na(x)] = base.prob
      }) 
      prob.pool = calc.prob(prob.pool.raw, interaction.matrix, interaction.method, occurrences)
    }
  } 
  
  # Create prob.pool object
  pools = list(occurrences = occurrences,
               env.pool = env.pool, 
               disp.pool = disp.pool, 
               prob.pool = prob.pool)
  pools = pools[!sapply(pools, is.null)] # remove empty pools
  species.richness = lapply(pools, FUN = function(x) sum(x))
  result = new("prob.pool", 
               pools = pools,
               interaction.matrix = interaction.matrix,
               interaction.method = c("Modification", "Multiplication", "None")[interaction.method],
               species.names = names(prob.pool),
               species.total = length(names(prob.pool)),
               species.mean = round(mean(raster::values(sum(prob.pool)), 1, na.rm = T)),
               species.richness = species.richness,
               slots = c("pools","interaction.matrix","interaction.method", "species.names", "species.total", "species.mean", "species.richness"))
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

calc.prob = function(probabilities, interaction.matrix, interaction.method, occurrences){
  interaction.matrix = as.matrix(interaction.matrix) # in case of dist object being provided
  if(exists("occurrences")){
    tmp.probs = occurrences
  } else {
    tmp.probs = probabilities
  }
  
  # Create "interaction pool"
  interaction.pool = calc(tmp.probs, function(prob.cell){
    sapply(1:length(prob.cell), function(species.index){
      prob.cell[species.index] * mean((prob.cell * interaction.matrix[,species.index]))
    })
  })
  
  # Calculate prob.pool
  if(interaction.method == 1){
    int.pos = calc(interaction.pool, fun = function(x){x[x <= 0] = 0; x})
    int.neg = calc(interaction.pool, fun = function(x){x[x > 0] = 0; x})
    probabilities = probabilities + (1-probabilities) * int.pos # positive interactions
    probabilities = probabilities + (probabilities * int.neg) # negative interactions
  } else if(interaction.method == 2){
    stop("not implemented")
    interaction = (interaction + 1) / 2
    prob.cell.new = species.cell * interaction
  }
  
  return(probabilities)
}

#####################################################################
######################### METHOD DEFINITIONS ########################
setMethod("summary", "prob.pool",  function(object){
  cat("Probabilistic species pool \n\n")
  cat(paste("Pools              : ", paste(names(object@pools), collapse = ", "), sep = ""), "\n")
  cat(paste("Species (total)    : ", object@species.total, "\n", sep = ""))
  cat(paste("Species (mean)     : ", object@species.mean, "\n", sep = ""))
  cat(paste("Interaction method : ", object@interaction.method, "\n", sep = ""))
  cat(paste("Resolution         : ", paste(round(res(object@pools$prob.pool), 3), collapse = " x "), " (x,y)\n", sep = ""))
  cat(paste("Extent             : ", paste(round(object@pools$prob.pool@extent[1:4], 3), 
                                           collapse = ", "), " (xmin, xmax, ymin, ymax)", sep = ""))
})

setMethod("print", "prob.pool", function(x){
  summary(x)
})


setMethod("plot", "prob.pool", function(x, focal.species = NULL, focal.unit = NULL, ...){
  moreargs = eval(substitute(list(...)))
  if(is.null(focal.species) && is.null(focal.unit)){
    # Plot species richness maps
    color.theme <- rasterTheme(region = rev(terrain.colors(30)))
    richness.maps = stack(x@species.richness[!is.na(x@species.richness)])
    do.call(levelplot, c(x = quote(richness.maps), par.settings = quote(color.theme), moreargs))
  } else if(!is.null(focal.species) & is.null(focal.unit)){
    # Plot probability pools for focal species
    if(!focal.species %in% x@species.names){stop("Species not found.")}
    color.theme <- rasterTheme(region = rev(terrain.colors(30)))
    species.maps = stack(lapply(x@pools, FUN = "[[", i = focal.species))
    do.call(levelplot, c(x = quote(species.maps), par.settings = quote(color.theme), moreargs))
  } else if(is.null(focal.species) & !is.null(focal.unit)){
    # Plot probability distributions for focal unit
    focal.probs = sapply(x@pools[names(x@pools) != "occurrences"], extract, y = focal.unit)
    focal.pdf = apply(focal.probs, 2, ppoibin, kk = 1:x@species.total)
    focal.pdf = melt(focal.pdf)
    focal.pdf$type = "Probability density"
    focal.cdf = apply(focal.probs, 2, dpoibin, kk = 1:x@species.total)
    focal.cdf = melt(focal.cdf)
    focal.cdf$type = "Cumulative probability density"
    
    focal.all = rbind(focal.pdf, focal.cdf)
    xyplot(value ~ Var1 | factor(Var2) + factor(type), data = focal.all, ylab = "probability", xlab= "species", col = "gray", type = "l")
  } else {
    stop("Please provide only one argument (species/focal.unit)")
  }
})