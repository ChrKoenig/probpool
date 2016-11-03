library(methods)
library(raster)

# Class definition
setClass("probpool",
         slots = c(env.pool = "Raster",
                   dist.pool = "Raster",
                   bio.pool = "Raster",
                   PSI.env = "Raster",
                   PSI.dist = "Raster",
                   PSI.bio = "Raster",
                   PSI.total = "Raster"))

# Constructor function
probpool = function(env.pool = NA, dist.pool = NA, bio.pool = NA){
  args = list(env.pool = env.pool, dist.pool = dist.pool, bio.pool = bio.pool)
  # Check arguments
  if(length(which(is.na(args))) == 3){stop("no probabilities provided")}
  
  # Check dimensions
  arg.dims = data.frame(lapply(args[!is.na(args)], function(x) dim(x)[1:2]))
  if(!all(apply(arg.dims, 1, function(x){length(unique(x)) == 1}))){stop("Dimension mismatch")}
  
  # fill 
  args[is.na(args)] = 1
  
  #
  PSI.env = sum(args$env.pool)
  PSI.dist = sum(args$dist.pool)
  PSI.bio = sum(args$bio.pool)
  PSI.total = sum(PSI.env, PSI.dist, PSI.bio)
  new("probpool", 
      env.pool = args$env.pool, 
      dist.pool = args$dist.pool, 
      bio.pool = args$bio.pool,
      PSI.env = PSI.env,
      PSI.dist = PSI.dist,
      PSI.bio = PSI.bio,
      PSI.total = PSI.total
  )
}