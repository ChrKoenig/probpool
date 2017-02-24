library(methods)
library(raster)

###################################################################
######################### CLASS DEFINITION ########################
# Class
setClass("prob.pool",
         slots = c(pools = "list",
                   pool.count = "numeric",
                   species = "character",
                   PSI = "list",
                   slots = "character"),
         validity = is.valid.prob.pool)

# Constructor function
prob.pool = function(env.pool = NULL, disp.pool = NULL, 
                     interaction.matrix = NULL, interaction.method = 1, occurences = NULL){
  pools = list(env.pool = env.pool, disp.pool = disp.pool, bio.pool = bio.pool, comb.pool = NULL)
  
  # Calculate comb.pool
  if(is.null(env.pool)){env.pool = 1} 
  if(is.null(disp.pool)){disp.pool = 1}
  disp.env.pool = disp.pool * env.pool # multiply probabilities, equals 1 if none is provided

  # 2 facilitation_competition:  
  # modify probabilities
  # raster > 1 works
  # values(raster) extracts all the stuff
  pool.count = length(which(!sapply(pools, is.null)))
  new("prob.pool", 
      pools = pools,
      pool.count = pool.count,
      species = names(pools[[min(which(!sapply(pools, is.null)))]]),
      PSI = list(env.pool = sum(pools$env.pool), disp.pool = sum(pools$disp.pool), 
                 bio.pool = sum(pools$bio.pool), comb.pool = sum(pools$comb.pool)),
      slots = c("pools", "pool.count", "species", "PSI")
  )
  
  
  if(is.numeric(comb.pool.raw)){ # disp.pool and env.pool not provided, bio.pool is guaranteed to be present due to validity checks
    if(bio.pool.method == "modify"){
      pools$comb.pool = (bio.pool+1)/2
    } else if(!is.null(bio.pool) & bio.pool.method == "multiply"){
      warning("Multiplication approach is not reliable when positive interactions are present.")
      pools$comb.pool = comb.pool.raw * (bio.pool+1)/2
    } else if(!is.null(bio.pool) & bio.pool.method == "modify"){
      pools$comb.pool = "XXXXX NOT IMPLEMENTED XXXXXX"
    } else {
      warning("No bio.pool provided. Ignoring bio.pool.method")
      pools$comb.pool = comb.pool.raw
    }
  }
}

# Function for validity check
is.valid.prob.pool = function(object){
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

###################################################################
######################### METHOD DEFINITIONS ########################
setMethod("summary", "prob.pool",  function(object, ...){
  cat("Probabilistic species pool \n\n")
  cat(paste("Pools:\t", paste(names(which(!sapply(object@pools, is.null))), collapse = ", "), sep = ""), "\n")
  cat(paste("Species:\t", length(object@species), "\n", sep = ""))
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
