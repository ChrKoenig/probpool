library(methods)
library(raster)
library(rasterVis)#
library(poibin)

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
    interaction.method = 3
  } else { # Interactions present
    if(!is.null(env.pool) || !is.null(disp.pool)){ # Base probabilities from env/disp layer
      prob.pool.raw = multiply.pools(env.pool, disp.pool)
      prob.pool = calc.prob(prob.pool.raw, interaction.matrix, interaction.method)
    } else { # Estimate base probabilities from occurence layer
      occurrences = calc(occurrences, function(x){x[x>1] = 1; x}) # Convert abundance to occurrence
      n.total = dim(occurrences)[3]
      n.mean = mean(raster::values(sum(occurrences)), na.rm = T)
      prob.pool.raw = calc(occurrences, function(x){ # Assign uniform distribution to raster
        x[!is.na(x)] = n.mean/n.total
        return(x)
      }) 
      prob.pool = calc.prob(prob.pool.raw, interaction.matrix, interaction.method)
    }
  } 
  
  # Create prob.pool object
  result = new("prob.pool", 
               pools = list(occurrences = occurrences,
                            env.pool = env.pool, 
                            disp.pool = disp.pool, 
                            prob.pool = prob.pool),
               interaction.matrix = interaction.matrix,
               interaction.method = c("Modification", "Multiplication", "None")[interaction.method],
               species.names = names(prob.pool),
               species.total = length(names(prob.pool)),
               species.mean = round(mean(raster::values(sum(prob.pool)), 1, na.rm = T)),
               species.richness = list(occurrences = if(is.null(occurrences)){NA} else {sum(occurrences)},
                                       env.pool = if(is.null(env.pool)){NA} else {sum(env.pool)}, 
                                       disp.pool = if(is.null(disp.pool)){NA} else {sum(disp.pool)},
                                       prob.pool = sum(prob.pool)),
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

calc.prob = function(probabilities, interaction.matrix, interaction.method){
  interaction.matrix = as.matrix(interaction.matrix) # in case of dist object being provided
  
  if(interaction.method == 2){
    warning("Caution: Results are not interpretable as probabilities using this interaction.method.")
  }
  
  modify = function(probability, interaction){
    if(is.na(probability) | is.na(interaction)){return(NA)}
    if(interaction > 0){return(probability + (1-probability) * interaction)}
    if(interaction <= 0){return(probability + (probability * interaction))}
  }
  
  prob.pool = calc(probabilities, function(prob.cell){
    interaction = sapply(1:length(prob.cell), function(species.index){ 
      prob.cell * interaction.matrix[,species.index]
    })
    interaction = rowMeans(interaction)
    if(interaction.method == 1){
      prob.cell.new = mapply(modify, prob.cell, interaction)
    } else if(interaction.method == 2){
      interaction = (interaction + 1) / 2
      prob.cell.new = species.cell * interaction
    }
    return(prob.cell.new)
  })
  return(prob.pool)
}

#####################################################################
######################### METHOD DEFINITIONS ########################
setMethod("summary", "prob.pool",  function(object){
  cat("Probabilistic species pool \n\n")
  cat(paste("Pools              : ", paste(names(which(!sapply(object@pools, is.null))), collapse = ", "), sep = ""), "\n")
  cat(paste("Species (total)    : ", object@species.total, "\n", sep = ""))
  cat(paste("Species (avg/cell) : ", object@species.mean, "\n", sep = ""))
  cat(paste("Interaction method : ", object@interaction.method, "\n", sep = ""))
  cat(paste("Resolution         : ", paste(round(res(object@pools$prob.pool), 3), collapse = " x "), " (x,y)\n", sep = ""))
  cat(paste("Extent             : ", paste(round(object@pools$prob.pool@extent[1:4], 3), 
                                           collapse = ", "), " (xmin, xmax, ymin, ymax)", sep = ""))
})

setMethod("print", "prob.pool", function(x){
  summary(x)
})


setMethod("plot", "prob.pool", function(x, species = NULL, focal.unit = NULL){
  if(is.null(species) && is.null(focal.unit)){
    color_theme <- rasterTheme(region = rev(terrain.colors(30)))
    richness_maps = stack(x@species.richness[!is.na(x@species.richness)])
    levelplot(richness_maps, par.settings = color_theme)
  } else if(!is.null(species) && is.null(focal.unit)){
    color_theme <- rasterTheme(region = rev(terrain.colors(30)))
    richness_maps = stack(x@species.richness[!is.na(x@species.richness)])
    levelplot(richness_maps, par.settings = color_theme)
  } else if(is.null(species) && !is.null(focal.unit)){
    focal.probs = extract(prob.pool@pools[[pool]], focal.unit)
    barplot(dpoibin(1:prob.pool@species.total, focal.probs), names.arg = 1:prob.pool@species.total, ylab = "Probabilities", xlab= "Species", main = pool)
    barplot(ppoibin(1:prob.pool@species.total, focal.probs), names.arg = 1:prob.pool@species.total, ylab = "Probabilities", xlab= "Species", main = pool)
  } else {
    stop("Please provide only one argument (species/focal.unit)")
  }
})

# Function to plot pool probability density function
prob.pool.PDF = function(prob.pool, pool, focal.unit){
  focal.probs = extract(prob.pool@pools[[pool]], focal.unit)
  barplot(ppoibin(1:prob.pool@species.total, focal.probs), names.arg = 1:prob.pool@species.total, ylab = "Probabilities", xlab= "Species", main = pool)
}

# Function to plot pool cumulative density function
plotPoolCDF<-function(prob.pool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(prob.pool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(ppoibin(seq(1,length(prob.pool@pools[[pool]][loc[1]][1,]),1), prob.pool@pools[[pool]][loc[1]][1,]), ylim=c(0,1),
          names.arg = prob.pool@species.names, ylab = "Probabilities", xlab= "Species", main = pool)
}

plotFacComp<-function(prob.pool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(prob.pool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(rev(sort(prob.pool@pools[[pool]][loc[1]][1,])), ylim=c(-1,1), ylab = "Probabilities", xlab= "Species", main = pool, col="red")
}
