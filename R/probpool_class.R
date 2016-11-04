library(methods)
library(raster)


###################################################################
######################### CLASS DEFINITION ########################
# Constructor function
probpool = function(env.pool = NULL, disp.pool = NULL, bio.pool = NULL){
  pools = list(env.pool = env.pool, disp.pool = disp.pool, bio.pool = bio.pool)
  species = lapply(pools, names)
  new("probpool", 
      pools = pools,
      species = names(pools[[which(!sapply(pools, is.null))]]),
      PSI = list(env.pool = sum(env.pool), disp.pool = sum(disp.pool), bio.pool = sum(bio.pool)),
      slots = c("pools", "species", "PSI")
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
  message(pool.dims)
  if(any(apply(pool.dims, 1, function(x) length(unique(x)) != 1))){
    errors = c(errors, "All pools need to have the same dimensions")
  }
  # check species names
  
  return(ifelse(length(errors) == 0, TRUE, errors))
}

# Class
setClass("probpool",
         slots = c(pools = "list",
                   species = "character",
                   PSI = "list",
                   slots = "character"),
         validity = is.valid.probpool)

###################################################################
######################### METHOD DEFINITIONS ########################
summary.probpool = function(object,...){
  cat("Probabilistic species pool \n\n")
  cat(paste("Pools present: ", names(which(!sapply(object@pools, is.null))), "\n", sep = ""))
  cat(paste("Species present: ", length(object@species), "\n", sep = ""))
}

print.probpool = function(object,...){
  summary(object)
}
