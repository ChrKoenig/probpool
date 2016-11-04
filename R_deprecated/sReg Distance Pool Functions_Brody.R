library(vegan)

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





####################################
# See how the function works with Adam's LBR data
setwd("C:\\Users\\Patrick\\Projects\\P24_SpeciesPool\\code_Brody")
d = read.csv("serpentineGrassland2005_speciesBySites.csv",as.is=T)
column = (d$X-1)%/%16 + 1
row = (d$X-1) - column*16 + 17
d$X = NULL

# Euclidean distance matrix
dists = as.matrix(dist(cbind(row,column)))


#####################
# Three examples
# For species i and cell j:
i = 1
j = 120
calcD(d[,i],dists[j,],1.5)

##########
# Apply to 28 species and 256 sites using sapply to calculate across
# all species and all sites
Pd = sapply(1:28,function(i){sapply(1:256, function(j){calcD(d[,i],dists[j,],1.5)})})

##########
# Calcualte Pd using the wrapper function
Pd2 = allD(d,dists,rep(1.5,ncol(d)))

identical(Pd,Pd2) #Good, the wrapper works.

# So, when using this, if you have the full occupancy marix 
# and the full pairwise distance matrix, you can just use
# allD. Otherwise, you can build it up using calcD calls 
# inside for loops.



################################################
# Visualize some distribution surfaces
require(raster)
r = raster(extent(c(0,1,0,1)),nrow = 16,ncol = 16)

par(mfrow = c(3,5),mar = c(1,1,4,1))
Pd3 = allD(d,dists,rep(0.75,ncol(d)))
plot(setValues(r,Pd3[,1]),col = rainbow(50)[25:50],main = paste(names(d)[1],0.75))
plot(setValues(r,Pd3[,5]),col = rainbow(50)[25:50],main = paste(names(d)[5],0.75))
plot(setValues(r,Pd3[,10]),col = rainbow(50)[25:50],main = paste(names(d)[10],0.75))
plot(setValues(r,Pd3[,16]),col = rainbow(50)[25:50],main = paste(names(d)[16],0.75))
plot(setValues(r,Pd3[,24]),col = rainbow(50)[25:50],main = paste(names(d)[24],0.75))

plot(setValues(r,Pd[,1]),col = rainbow(50)[25:50],main = paste(names(d)[1],1.5))
plot(setValues(r,Pd[,5]),col = rainbow(50)[25:50],main = paste(names(d)[5],1.5))
plot(setValues(r,Pd[,10]),col = rainbow(50)[25:50],main = paste(names(d)[10],1.5))
plot(setValues(r,Pd[,16]),col = rainbow(50)[25:50],main = paste(names(d)[16],1.5))
plot(setValues(r,Pd[,24]),col = rainbow(50)[25:50],main = paste(names(d)[24],1.5))

#Give dispersal abilities based on traits
traits = traits[match(names(d),traits$species),]
traits$seedMass_mg[which(is.na(traits$seedMass_mg))] = exp(mean(log(traits$seedMass_mg),na.rm=T))
dispAbility = traits$heightClass
dispAbility[which(traits$habit == "annual")] = dispAbility[which(traits$habit == "annual")]+1
dispAbility = dispAbility - decostand(log(traits$seedMass_mg),method = "range")-0.5
dispAbility = decostand(dispAbility,method = "range") * 3 + 0.25

Pd2 = allD(d,dists,dispAbility)
plot(setValues(r,Pd2[,1]),col = rainbow(50)[25:50],main = paste(names(d)[1],dispAbility[1]))
plot(setValues(r,Pd2[,5]),col = rainbow(50)[25:50],main = paste(names(d)[5],dispAbility[5]))
plot(setValues(r,Pd2[,10]),col = rainbow(50)[25:50],main = paste(names(d)[10],dispAbility[10]))
plot(setValues(r,Pd2[,16]),col = rainbow(50)[25:50],main = paste(names(d)[16],dispAbility[16]))
plot(setValues(r,Pd2[,24]),col = rainbow(50)[25:50],main = paste(names(d)[24],dispAbility[24]))

