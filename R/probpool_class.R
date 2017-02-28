library(methods)
library(raster)

###################################################################
######################### CLASS DEFINITION ########################
# Class
setClass("prob.pool",
         slots = c(pools = "list",
                   species = "character",
                   interaction.method = "character",
                   PSI = "list",
                   slots = "character"),
         validity = is.valid.prob.pool)

# Function for validity check
is.valid.prob.pool = function(object){
  errors = character()
  if(all(sapply(object@pools, is.null))){ 
    errors = c(errors,"Please provide at least one of the following arguments: env.pool, disp.pool, occurences")
  }
  if(is.null(object@interaction.matrix) & is.null(object@pools$env.pool) & is.null(object@pools$disp.pool)){
    errors = c(errors,"Interaction.matrix is missing.")
  }
  if(!all(sapply(object@pools, function(pool) {extends(class(pool), "Raster") | is.null(pool)}))){ # check types
    errors = c(errors, "Invalid argument. Please provide a raster object.")
  }
  if(Reduce(all.equal, lapply(object@pools[!sapply(object@pools, is.null)], names))){ # check species names
    errors = c(errors, "Species names do not match.")
  }
  pool.dims = matrix(sapply(object@pools, dim)) # Check dimensions
  if(any(apply(pool.dims, 1, function(x) length(unique(x)) != 1))){
    errors = c(errors, "All raster objects need to have the same dimensions.")
  }
  if(is.na(interaction.method)){
    errors = c(errors, "Unknown interaction.method. Choose '1' for modification or '2' for multiplication")
  }
  if(interaction.method == 2 & !is.null(interaction.matrix) & is.null(occurences)){
    errors = c(errors, "Multiplication approach (interaction.method = 2) requires species occurences.")
  }
  return(ifelse(length(errors) == 0, TRUE, errors))
}

# Constructor function
prob.pool = function(env.pool = NULL, disp.pool = NULL, occurences = NULL,
                     interaction.matrix = NULL, interaction.method = 1){
  
  pp = new("prob.pool", 
           pools = list(env.pool = env.pool, disp.pool = disp.pool, occurences = occurences, prob.pool = NULL),
           species = names(pools[[min(which(!sapply(pools, is.null)))]]),
           interaction.method = c("Modification","Multiplication")[interaction.method],
           PSI = list(env.pool = sum(pools$env.pool), disp.pool = sum(pools$disp.pool), 
                      occurences = sum(pools$occurences), prob.pool = sum(pools$prob.pool)),
           slots = c("pools", "species", "interaction.method", "PSI"))
  
  valid_pool = c(!is.null(env.pool), !is.null(disp.pool))
  if(is.null(interaction.matrix)){ # Easy case: no interactions, ignore occurences
    pools$prob.pool = ifelse(valid_pool[1], env.pool, 1) * ifelse(valid_pool[1], disp.pool, 1)
  } else { # Interactions present
    if(interaction.method == 1){ # Modification approach
      if(any(valid_pool)){ # Base probabilities from env/disp layer
        prob.pool = ifelse(valid_pool[1], env.pool, 1) * ifelse(valid_pool[1], disp.pool, 1)
        multiply.prob(occurences, interaction.matrix, abundance = T) 
      } else{ # Estimate base probabilities from occurence layer
        n_total = dim(occurences)[3]
        n_mean = mean(raster::values(sum(occurences)), na.rm = T)
        b = n_mean/n_total
        
        interaction_values = values(sum(example_bio.pool))
        interaction_values_modified = sapply(interaction_values, FUN = modify_prob, p = N_mean)
        richness_modified = richness
        values(richness_modified) = interaction_values_modified
        N_mean_modified = mean(values(richness_modified), na.rm = T)
      }
    } else { # Multiplication approach
      warning("Caution: Results are not interpretable as probabilities using this interaction.method.")
      interactions = multiply.prob(occurences, interaction_matrix, abundance = T)
      if(any(valid_pool)){ # Multiply interaction values with probabilities from env/disp layer
        prob.pool = ifelse(valid_pool[1], pp@env.pool, 1) * ifelse(valid_pool[1], pp@disp.pool, 1)
        prob.pool.final = multiply.prob(prob.pool, interaction.matrix, abundance = T)
      } else{ # Multiply interaction values with occurence layer and rescale to 0-1
        prob.pool = multiply.prob(occurences, interaction.matrix, abundance = T) 
        prob.pool.final = (prob.pool.final + 1) / 2# rescale
      }
    }
  } 
}


modify.prob = function(probabilities, interaction.matrix){ 
  modify = function(p,b){
    if(is.na(p) | is.na(b)){return(NA)}
    if(b > 0){return(p + (1-p)*b)}
    if(b <= 0){return(p + (p*b))}
  }
  probabilities <- values(probabilities)
  probabilities <- probabilities[complete.cases(probabilities),]
}

multiply.prob <- function(probabilities, interaction.matrix, abundance = TRUE){
  probabilities <- values(probabilities)
  probabilities <- probabilities[complete.cases(probabilities),]
  
  if(abundance){ 
    probabilities <- probabilities/max(probabilities)
  } else {
    probabilities[probabilities > 0] <- 1 
  }
  
  # multiply the incoming interactions of each species x (columns in int.matrix)
  # with the occurrence of all other species for the given site y
  interactions <- lapply(1:dim(probabilities)[3], function(x) {
    interactions.x <- t(sapply(1:nrow(probabilities), function(y) probabilities[y,] * interaction.matrix[,x]))
    interactions.x <- rowMeans(interactions.x)
    interactions.x.rst <- probabilities[[x]]
    interactions.x.rst[!is.na(interactions.x.rst)] <- interactions.x
    return(interactions.x.rst)
  })
  interactions <- stack(interactions)     
  return(interactions)
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
