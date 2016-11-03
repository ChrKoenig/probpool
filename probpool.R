library(methods)
library(raster)

# Constructor function
probpool = function(env.pool = NULL, disp.pool = NULL, bio.pool = NULL){
  new("probpool", 
      pools = list(env.pool = env.pool, disp.pool = disp.pool, bio.pool = bio.pool),
      PSI = list(PSI.env = sum(env.pool), PSI.disp = sum(disp.pool), 
                 PSI.bio = sum(bio.pool), PSI.total = sum(env.pool, disp.pool, bio.pool)),
      PDF = list(),
      CDF = list(),
      slots = c("pools", "PSI", "PDF", "CDF")
  )
}

is.valid.probpool = function(object){
  errors = character()
  if(all(sapply(object@pools, is.null))){
    errors = c(errors,"no probabilities provided")
  }
  # Check types
  sapply(object@pools, function(pool){
    if(!(extends(class(pool), "Raster") | is.matrix(pool) | is.null(pool))){
      errors = c(errors, "Invalid pool type. Please provide a raster object or a matrix")
    }
  })
  # Check dimensions
  pool.dims = matrix(sapply(object@pools, dim))
  if(length(unique(pool.dims[1,])) != 1 | length(unique(pool.dims[2,])) != 1){
    errors = c(errors, "All pools need to have the same dimensions")
  }
  
  return(ifelse(length(errors) == 0, TRUE, errors))
}

# Class definition
setClass("probpool",
         slots = c(pools = "list",
                   PSI = "list",
                   PDF = "list",
                   CDF = "list",
                   slots = "character"),
         validity = is.valid.probpool)
