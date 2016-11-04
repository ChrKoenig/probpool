library(methods)
library(raster)

# Constructor function
probpool = function(env.pool = NULL, disp.pool = NULL, bio.pool = NULL){
  pools = list(env.pool = env.pool, disp.pool = disp.pool, bio.pool = bio.pool)
  species = lapply(pools, names)
  new("probpool", 
      pools = pools,
      species = names(env.pool),
      PSI = list(PSI.env = sum(env.pool), PSI.disp = sum(disp.pool), PSI.bio = sum(bio.pool)),
      slots = c("pools", "species", "PSI")
  )
}

is.valid.probpool = function(object){
  errors = character()
  if(all(sapply(object@pools, is.null))){
    errors = c(errors,"no probabilities provided")
  }
  # Check types
  for(pool in object@pools){
    if(!(extends(class(pool), "Raster") | is.matrix(pool) | is.null(pool))){
      errors = c(errors, "Invalid pool type. Please provide a raster object or a matrix")
    }
  }
  
  # Check dimensions
  pool.dims = matrix(sapply(object@pools, dim))
  if(any(apply(pool.dims, 1, function(x) length(unique(x)) != 1))){
    errors = c(errors, "All pools need to have the same dimensions")
  }
  return(ifelse(length(errors) == 0, TRUE, errors))
}

# Class definition
setClass("probpool",
         slots = c(pools = "list",
                   species = "character",
                   PSI = "list",
                   slots = "character"),
         validity = is.valid.probpool)


plot(env.pool)

