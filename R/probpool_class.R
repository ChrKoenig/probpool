library(methods)
library(raster)


###################################################################
######################### CLASS DEFINITION ########################
# Constructor function
probpool = function(env.pool = NULL, disp.pool = NULL, bio.pool = NULL){
  pools = list(env.pool = env.pool, disp.pool = disp.pool, bio.pool = bio.pool)
  
  # Calculate comb.pool
  comb.pool = disp.pool * env.pool # multiplied probabilitie  # modify probabilities
  
  pool.count = length(which(!sapply(pools, is.null)))
  new("probpool", 
      pools = pools,
      pool.count = pool.count,
      species = names(pools[[min(which(!sapply(pools, is.null)))]]),
      PSI = list(env.pool = sum(env.pool), disp.pool = sum(disp.pool), bio.pool = sum(bio.pool)),
      slots = c("pools", "pool.count", "species", "PSI")
  )
}

# Validity function
is.valid.probpool = function(object){
  errors = character()
  if(all(sapply(object@pools, is.null))){ # Check arguments
    errors = c(errors,"no probabilities provided")
  }
  for(pool in object@pools){ # Check types
    if(!(extends(class(pool), "Raster") | is.matrix(pool) | is.null(pool))){
      errors = c(errors, "Invalid pool type. Please provide a raster object or a matrix")
    }
  }
  pool.dims = matrix(sapply(object@pools, dim)) # Check dimensions
  if(any(apply(pool.dims, 1, function(x) length(unique(x)) != 1))){
    errors = c(errors, "All pools need to have the same dimensions")
  }
  # check species names
  
  return(ifelse(length(errors) == 0, TRUE, errors))
}

# Class
setClass("probpool",
         slots = c(pools = "list",
                   pool.count = "numeric",
                   species = "character",
                   PSI = "list",
                   slots = "character"),
         validity = is.valid.probpool)

###################################################################
######################### METHOD DEFINITIONS ########################
setMethod("summary", "probpool",  function(object, ...){
  cat("Probabilistic species pool \n\n")
  cat(paste("Pools:\t", paste(names(which(!sapply(object@pools, is.null))), collapse = ", "), sep = ""), "\n")
  cat(paste("Species present:\t", length(object@species), "\n", sep = ""))
})

setMethod("print", "probpool", function(x){
  summary(x)
})

setMethod("plot", c("probpool"),
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
                plotPoolProbs(x,"bio.pool",focalunit)
                plotPoolPDF(x,"bio.pool",focalunit)
                plotPoolCDF(x,"bio.pool",focalunit)
              }
              if(is.null(x@pools$comb.pool)==FALSE)
              {
                plotPoolProbs(x,"comb.pool",focalunit)
                plotPoolPDF(x,"comb.pool",focalunit)
                plotPoolCDF(x,"comb.pool",focalunit)
              }
            }else{
              par(mfrow=c(1,x@pool.count))    
              if(is.null(x@PSI$disp.pool)==FALSE)  
              {
                plotRasterPool(x,"disp.pool")
              }
              if(is.null(x@PSI$env.pool)==FALSE)  
              {
                plotRasterPool(x,"env.pool")
              }
              if(is.null(x@PSI$bio.pool)==FALSE)  
              {
                plotRasterPool(x,"bio.pool")
              }
              if(is.null(x@PSI$comb.pool)==FALSE)  
              {
                plotRasterPool(x,"comb.pool")
              }
            }
          })
