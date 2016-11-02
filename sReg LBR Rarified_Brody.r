###############################################################
# for a site, compare the pool values across different methods
# using the jaccard similarity
require(vegan)
setwd("C:/Documents/Leipzig trip/sDiv/serpentineGrassland_adamSmith")
d = read.csv("serpentineGrassland2005_speciesBySites.csv",as.is=T,row.names = 1)

column = (1:256-1)%/%16 + 1
row = (1:256-1) - column*16 + 17


# Give dispersal abilities based on traits
traits = read.csv("Traits.csv",as.is=T)
traits = traits[match(names(d),traits$species),]
traits$seedMass_mg[which(is.na(traits$seedMass_mg))] = exp(mean(log(traits$seedMass_mg),na.rm=T))
dispAbility = traits$heightClass
dispAbility[which(traits$habit == "annual")] = dispAbility[which(traits$habit == "annual")]+1
dispAbility = dispAbility - decostand(log(traits$seedMass_mg),method = "range")-0.5
dispAbility = decostand(dispAbility,method = "range")

dispAbility0 = rep(0.45,ncol(d))

# Calculate original pool values
Pd = allD(d, as.matrix(dist(cbind(row,column))),dispAbility0)
Pe = allE(d, env,method = "beals")
Ps = Pd * Pe

# Throw away plot data
ssize = (1:16)*16
runs = 100
repSpMeans = NULL #Stores the mean correlation of species scores per run
spMeans = NULL #Stores the mean correlations of species scores per species
repSiteMeans = NULL #Stores the mean correlation of site scores per run
siteMeans = NULL #Stores the mean correlation of site scores per site
for(j in 1:length(ssize))
	{
	spCorMat = matrix(NA,nrow = runs, ncol = ncol(d))
	siteCorMat = matrix(NA,nrow = runs, ncol = ssize[j])
	for(i in 1:runs)
		{
		samp = sample(1:256,ssize[j],replace = F)
		PdS = allD(d[samp,], as.matrix(dist(cbind(row,column)))[samp,samp],dispAbility0)
		PeS = allE(d[samp,], env[samp,],method = "beals")
		PsS = PdS * PeS

		spCorMat[i,] = sapply(1:ncol(Ps),function(i){cor(PsS[,i],Ps[samp,i])})
		siteCorMat[i,] = sapply(1:nrow(PsS),function(i){cor(PsS[i,],Ps[samp[i],])})
		}
	repSpMeans[j] = mean(apply(spCorMat,1,mean,	na.rm=T),na.rm=T)
	spMeans[j] = mean(apply(spCorMat,2,mean,na.rm=T),na.rm=T) #somewhat related to species abundance

	repSiteMeans[j] = mean(apply(siteCorMat,1,mean,na.rm=T),na.rm=T)
	siteMeans[j] = mean(apply(siteCorMat,2,mean,na.rm=T),na.rm=T)#These two have the same mean because there are no NAs in siteCorMat
	print(j)
	}
#This looks really good, but it's a bit unfair at the moment
#because species with no estimated D or E are just excluded
#giving them a small
plot(ssize,repSpMeans)
plot(ssize,repSiteMeans)

