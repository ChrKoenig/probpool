library(methods)
library(raster)

###################################################################
# Probabilistic species pool class definition
# author: Christian König
# date:   March 2017
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
  if(is.null(object@interaction.matrix) & is.null(object@pools$env.pool) & is.null(object@pools$disp.pool)){
    return("Interaction.matrix is missing.")
  }
  if(is.na(object@interaction.method)){
    return("Unknown interaction.method. Choose '1' for modification or '2' for multiplication")
  }
  if(!all(sapply(object@pools, function(pool) {extends(class(pool), "Raster") | is.null(pool)}))){ # check types
    return("Invalid argument. Please provide a raster object.")
  }
  pool.dims = matrix(sapply(object@pools, dim)) # Check dimensions
  if(any(apply(pool.dims, 1, function(x) length(unique(x)) != 1))){
    return("All raster objects need to have the same dimensions.")
  }
  if(!all(sapply(object@pools, function(x){is.null(x) || all.equal(names(object@pools$prob.pool), names(x))}))){ # check species names
    return("Species names do not match.")
  }
  if(object@interaction.method == 2 & !is.null(object@interaction.matrix) & is.null(object@pools$occurrences)){
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
                   "species.tot",
                   "species.avg",
                   "species.names",
                   "PSI",
                   "slots"),
         validity = is.valid.prob.pool)


# Constructor function
prob.pool = function(env.pool = NULL, disp.pool = NULL, occurrences = NULL,
                     interaction.matrix = NULL, interaction.method = 1){
  
  prob.pool = NULL
  if(is.null(interaction.matrix) & !is.null(env.pool) | !is.null(disp.pool)){ # No interactions, multiply env * disp
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
      prob.pool.raw = occurrences 
      values(prob.pool.raw) = n.mean/n.total # Assign uniform distribution to raster
      prob.pool = calc.prob(prob.pool.raw, interaction.matrix, interaction.method)
    }
  } 
  
  # Create prob.pool object
  result = new("prob.pool", 
               pools = list(env.pool = env.pool, 
                            disp.pool = disp.pool, 
                            occurrences = occurrences, 
                            prob.pool = prob.pool),
               interaction.matrix = interaction.matrix,
               interaction.method = c("Modification", "Multiplication", "None")[interaction.method],
               species.tot = length(names(prob.pool)),
               species.avg = round(mean(raster::values(sum(prob.pool)), 1, na.rm = T)),
               species.names = names(prob.pool),
               PSI = list(env.pool = if(is.null(env.pool)){NA} else {sum(env.pool)},
                          disp.pool = if(is.null(disp.pool)){NA} else {sum(disp.pool)},
                          occurrences = if(is.null(occurrences)){NA} else {sum(occurrences)},
                          prob.pool = sum(prob.pool)),
               slots = c("pools","interaction.matrix","interaction.method", "species.tot", "species.avg", "species","PSI"))
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
  if(interaction.method == 2){
    warning("Caution: Results are not interpretable as probabilities using this interaction.method.")
  }
  
  modify = function(probability, interaction){
    if(is.na(probability) | is.na(interaction)){return(NA)}
    if(interaction > 0){return(probability + (1-probability)*interaction)}
    if(interaction <= 0){return(probability + (probability*interaction))}
  }
  
  prob.pool = calc(probabilities, function(species.cell){
    interaction = sapply(1:length(species.cell), function(species.index){ 
      species.cell * interaction.matrix[,species.index]
    })
    interaction = rowMeans(interaction)
    if(interaction.method == 1){
      warning("Test")
      prob.cell = mapply(modify, species.cell, interaction)
    } else if(interaction.method == 2){
      interaction = (interaction + 1) / 2
      prob.cell = species.cell * interaction
    }
    return(prob.cell)
  })
}

#####################################################################
######################### METHOD DEFINITIONS ########################
setMethod("summary", "prob.pool",  function(object, ...){
  cat("Probabilistic species pool \n\n")
  cat(paste("Pools              : ", paste(names(which(!sapply(object@pools, is.null))), collapse = ", "), sep = ""), "\n")
  cat(paste("Species (total)    : ", object@species.tot, "\n", sep = ""))
  cat(paste("Species (avg/cell) : ", object@species.avg, "\n", sep = ""))
  cat(paste("Interaction method : ", object@interaction.method, "\n", sep = ""))
  cat(paste("Resolution         : ", paste(round(res(object@pools$prob.pool), 3), collapse = " x "), " (x,y)\n", sep = ""))
  cat(paste("Extent             : ", paste(round(object@pools$prob.pool@extent[1:4], 3), collapse = ", "), " (xmin, xmax, ymin, ymax)", sep = ""))
})

setMethod("print", "prob.pool", function(x){
  summary(x)
})




setMethod("plot", c("prob.pool"),
          function(x, focalunit=FALSE,...)
          {
            if(length(focalunit)==2)
            {
              par(mfrow=c(x@pool.count,3))
              if(is.null(x@pools$disp.pool)==FALSE)
              {
                plotPoolProbs(x,"disp.pool",focalunit)
                plotPoolPDF(x,"disp.pool",focalunit)
                plotPoolCDF(x,"disp.pool",focalunit)
              }
              if(is.null(x@pools$env.pool)==FALSE)
              {
                plotPoolProbs(x,"env.pool",focalunit)
                plotPoolPDF(x,"env.pool",focalunit)
                plotPoolCDF(x,"env.pool",focalunit)
              }
              if(is.null(x@pools$bio.pool)==FALSE)
              {
                if(x@method=="simple_multiplication")
                {
                  plotPoolProbs(x,"bio.pool",focalunit)
                  plotPoolPDF(x,"bio.pool",focalunit)
                  plotPoolCDF(x,"bio.pool",focalunit)  
                }else{
                  plotFacComp(x,"bio.pool",focalunit)
                  plotPoolPDF(x,"bio.pool",focalunit)
                  plotPoolCDF(x,"bio.pool",focalunit)
                }
              }
              if(is.null(x@pools$comb.pool)==FALSE)
              {
                plotPoolProbs(x,"comb.pool",focalunit)
                plotPoolPDF(x,"comb.pool",focalunit)
                plotPoolCDF(x,"comb.pool",focalunit)
              }
            }else{
              par(mfrow=c(1,x@pool.count))    
              if(is.null(x@PSI$env.pool)==FALSE)
              {
                plotRasterPool(x,"disp.pool")
              }
              if(is.null(x@PSI$env.pool)==FALSE)
              {
                plotRasterPool(x,"env.pool")
              }
              if(is.null(x@PSI$env.pool)==FALSE)
              {
                plotRasterPool(x,"bio.pool")
              }
              if(is.null(x@PSI$env.pool)==FALSE)
              {
                plotRasterPool(x,"comb.pool")
              }
            }
          })
