library(sp);library(raster); library(gdistance)

##############################################
# calcD is a function that calculates the 
# dispersal probabilities for a species to 
# all locations on the landscale from a focal
# cell. It takes a vector containing the species 
# presence/absence information, a vector containing
# the distances from a focal cell to other others
# and a constant describing the species' dispersal
# ability. It returns a single probability value
# describing the chance of dispersing into that cell

# Also, allow two methods of calculating the 
# dispersal kernel (negative exponential 
# and fat-tailed)

calcD = function(occupancy, distance, k, method = "negexp")
{
  index = which(occupancy > 0)
  if(method == "negexp") {distFunction = function(d,k){1 - prod(1-exp(-1*d/k))}}
  if(method == "fattail") {distFunction = function(d,k){1 - prod(1/d^k)}}
  return(distFunction(distance[index],k))
}


# occurence.surfaces needs to be a raster stack including rasters of species occurances or abundances with values of 0 for absences and values > 0 for occurances. Abundances are not considered though.

# disp.ability needs to be a vector of dispersal abilities. Scaling range? Larger values = higher dispersal abilities

# cond.surfaces needs to be a raster stack including rasters of conductance values between 0 and one
# NAs will be replaced by 0 and 0s will be replaced by small value; set variable for this?


disp_pool <- function(occurrence.surfaces, disp.ability, method=c("negexp","fattail"),cond.surfaces=NULL, longlat=TRUE) {
  
  # TODO: check extent of rasters, type of data etc.
  
  occurrences <- rasterToPoints(occurrence.surfaces)
  
  if (is.null(cond.surfaces)){
    
    distances <- spDists(occurrences[,c(1:2)], occurrences[,c(1:2)], longlat = longlat)
    
    dispersal <- lapply(1:dim(occurrence.surfaces)[3], function(x) {
      
      dispersal.x <- sapply(1:nrow(occurrences), function(y) calcD(occurrences[,x+2],distances[,y], ifelse(length(disp.ability)==1,disp.ability,disp.ability[x]), method = method[1]))
      dispersal.x.rst <- occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] <- dispersal.x
      return(dispersal.x.rst)
    })
    
  } else {
    
    # TODO add extra margin around conductance layers?
    
    # replace NAs by 0
    cond.surfaces[is.na(cond.surfaces)] <- 1/100
    # truncation: replace <1/100 by 0                                   
    cond.surfaces[cond.surfaces<1/100] <- 1/100
    
    if (class(cond.surfaces)=="RasterLayer") {
      # create transition raster
      spec.trans <- transition(cond.surfaces, mean, 8)
      # geocorrection
      spec.trans <- geoCorrection(spec.trans, type="c")
      
      # calculate commute distances
      distances <- commuteDistance(spec.trans, occurrences[,c(1:2)])/2
      distances <- as.matrix(distances)
    }
    
    
    dispersal <- lapply(1:dim(occurrence.surfaces)[3], function(x) {
      
      if (class(cond.surfaces) %in% c("RasterStack","RasterBrick")) {
        # create transition raster
        spec.trans <- transition(cond.surfaces[[x]], mean, 8)
        # geocorrection
        spec.trans <- geoCorrection(spec.trans, type="c")
      
        # calculate commute distances
        distances <- commuteDistance(spec.trans, occurrences[,c(1:2)])/2
        distances <- as.matrix(distances)
      }
      
      dispersal.x <- sapply(1:nrow(occurrences), function(y) calcD(occurrences[,x+2],distances[,y], ifelse(length(disp.ability)==1,disp.ability,disp.ability[x]), method = method[1]))  # check 5000
      dispersal.x.rst <- occurrence.surfaces[[x]]
      dispersal.x.rst[!is.na(dispersal.x.rst)] <- dispersal.x
      return(dispersal.x.rst)
    })
  }
  
  dispersal <- stack(dispersal)      
  
  return(dispersal)
}



# occurence.surfaces needs to be a raster stack including rasters of species occurences or abundances with values of 0 for absences and values > 0 for occurances. Values wil be scaled to range from 0 to 1
# or the disp.pool, the env.pool or the product of both

# int.matrix is a species by species matrix (may be asymmetric) with interactions assumed to be directed from the species in the row to the species in the colums

bio_pool <- function(occurrence.surfaces, interaction.matrix, abundance=TRUE) {
  occurrences <- values(occurrence.surfaces)
  occurrences <- occurrences[complete.cases(occurrences),]
  
  # Null probability
  richness.mean = mean(rowSums(occurrences), na.rm = TRUE)
  richness.total = dim(occurrences)[2]
  null_prob = richness.mean/richness.total 

  if(abundance){ 
    occurrences <- occurrences/max(occurrences)
  } else {
    occurrences[occurrences > 0] <- 1 
  }
  
  # multiply the incoming interactions of each species x (columns in int.matrix)
  # with the occurrence/probability of all other species for the given site y
  interactions <- lapply(1:dim(occurrence.surfaces)[3], function(x) {
    interactions.x <- t(sapply(1:nrow(occurrences), function(y) occurrences[y,] * interaction.matrix[,x]))
    interactions.x <- rowMeans(interactions.x)
    interactions.x.rst <- occurrence.surfaces[[x]]
    interactions.x.rst[!is.na(interactions.x.rst)] <- interactions.x
    return(interactions.x.rst)
  })

  interactions <- stack(interactions)     
  return(interactions)
}

prob.surface = example_env.pool
prob.surface.inverse = 1-example_env.pool
cells.positive = example_bio.pool > 0
cells.negative = example_bio.pool < 0
prob.surface.facilitated = prob.surface + (prob.surface.inverse[cells.positive] * example_bio.pool[cells.positive])
prob.surface.competed = prob.surface + (prob.surface * cells.negative)
# allD is a wrapper function for calcD that applies the calculation
# of calcD to all cells and all species. It takes the full occupancy
# matrix, the full pairwise distance matrix and a vector of each
# species' dispersal ability. As with calcD, you specify the dispersal
# function by specifying "method"
allD = function(occupancy,distance,k,method = "negexp")
{
  Pd = sapply(1:ncol(occupancy),function(i){sapply(1:nrow(occupancy), function(j){calcD(occupancy[,i],distance[j,],k[i],method)})})
  return(Pd)
}

# allE is a function that calculates the environmental suitability
# for all species across all sites. There are two methods available:
# 1. consider the distance between a focal site and the closest (in 
# environmental space) occupied site, and rank it against all distances
# 2. use beals smoothing
calcE = function(occupancy,envdist,site)
{
  index = which(occupancy > 0)
  es = 1-length(which(envdist[,site] < min(envdist[index,site])))/256
  return(es)
}

allE = function(occupancy,environment=NULL, method)
{
  if(method == "mindist")
  {
    dists = as.matrix(dist(environment))
    Pe = sapply(1:ncol(occupancy),function(i){sapply(1:nrow(occupancy), function(j){calcE(d[,i],dists,j)})})
  }
  if(method == "beals")
  {
    require(vegan)
    Pe = beals(occupancy)	
  }
  return(Pe)
}


# Function to plot pool for a single focal unit defined by focalunit
plotPoolProbs<-function(probpool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(probpool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(rev(sort(probpool@pools[[pool]][loc[1]][1,])), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = pool)
}

# Function to plot pool probability density function
plotPoolPDF<-function(probpool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(probpool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(dpoibin(seq(1,length(probpool@pools[[pool]][loc[1]][1,]),1), probpool@pools[[pool]][loc[1]][1,]), names.arg=1:length(probpool@pools[[pool]][loc[1]][1,]), ylab = "Probabilities", xlab= "Species", main = pool)
}

# Function to plot pool cumulative density function
plotPoolCDF<-function(probpool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(probpool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(ppoibin(seq(1,length(probpool@pools[[pool]][loc[1]][1,]),1), probpool@pools[[pool]][loc[1]][1,]), ylim=c(0,1), names.arg=1:length(probpool@pools[[pool]][loc[1]][1,]), ylab = "Probabilities", xlab= "Species", main = pool)
}

plotFacComp<-function(probpool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(probpool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(rev(sort(probpool@pools[[pool]][loc[1]][1,])), ylim=c(-1,1), ylab = "Probabilities", xlab= "Species", main = pool, col="red")
}

# Function to plot pool as a raster
plotRasterPool<-function(probpool,pool)
{
  raster::plot(probpool@PSI[[pool]], main=pool)
}


#TODO: Function to convert species by sites to raster stack


      