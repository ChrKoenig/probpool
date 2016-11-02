library(methods)
library(raster)

# Class definition
setClass("probpool",
         slots = c(env.pool = "Raster",
                   dist.pool = "Raster",
                   bio.pool = "Raster",
                   PSI.env = "numeric",
                   PSI.dist = "numeric",
                   PSI.bio = "numeric",
                   PSI.total = "numeric"))

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
  
  new("probpool", 
      env.pool = args$env.pool, 
      dist.pool = NA, 
      bio.pool = args$bio.pool,
      PSI.env = 1,
      PSI.dist = 1,
      PSI.bio = 1,
      PSI.total = 1
  )
}

#PSI = Reduce(f = "*", c(sum(values(env.pool), na.rm = T), sum(values(dist.pool), na.rm = T), sum(values(bio.pool), na.rm = T))