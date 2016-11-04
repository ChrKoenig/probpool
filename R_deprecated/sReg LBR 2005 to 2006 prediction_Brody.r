
###########################################################
# Get occurence data in 2005 and 2006, and make them match
setwd("C:/Documents/Leipzig trip/sDiv/serpentineGrassland_adamSmith")
d = read.csv("serpentineGrassland2005_speciesBySites.csv",as.is=T,row.names = 1)
d2 = read.csv("serpentineGrassland2006_speciesBySites.csv",as.is=T)
d2 = d2[,which(names(d2)%in%names(d))]
d2$uS = 0
d2$CL = 0
d2 = d2[,match(names(d),names(d2))]

# Convert plot numbers to rows and columns
column = (d$X-1)%/%16 + 1
row = (d$X-1) - column*16 + 17

d$X = NULL
d2$X = NULL

# Give dispersal abilities based on traits
traits = read.csv("Traits.csv",as.is=T)

traits = traits[match(names(d),traits$species),]
traits$seedMass_mg[which(is.na(traits$seedMass_mg))] = exp(mean(log(traits$seedMass_mg),na.rm=T))
dispAbility = traits$heightClass
dispAbility[which(traits$habit == "annual")] = dispAbility[which(traits$habit == "annual")]+1
dispAbility = dispAbility - decostand(log(traits$seedMass_mg),method = "range")-0.5
dispAbility = decostand(dispAbility,method = "range")

# Prepare environmental data
envTable = read.csv('serpentineGrassland_sitesByEnv.csv')
envStack = stack()

# convert each env layer to raster
for (thisVar in c('fineSoilVol_gPerCm3', 'pH', 'Ca_mgPerG', 'K_mgPerG', 'Mg_gPerG', 'Mn_mgPerG', 'Na_mgPerG', 'Ni_mgPerG')) 
	{
	thisEnv = raster(matrix(envTable[ , thisVar], ncol=8))
	newEnv = disaggregate(x=thisEnv, fact=2, method='bilinear')
	envStack = stack(envStack, newEnv)
	}
names(envStack) = c('fineSoilVol_gPerCm3', 'pH', 'Ca_mgPerG', 'K_mgPerG', 'Mg_gPerG', 'Mn_mgPerG', 'Na_mgPerG', 'Ni_mgPerG')

# add Ca/Mg
envStack = stack(envStack, subset(envStack, 'Ca_mgPerG') / subset(envStack, 'Mg_gPerG'))
names(envStack)[nlayers(envStack)] = 'Ca_to_Mg_ratio'
x = (rep(1:16, each=16) - 0.5) / 16
y = rep(seq(15.5, 0.5, by=-1), 16) / 16

# get env data for each cell and scale
env = extract(envStack, cbind(x, y))
env = scale(env)

#############################################################

#############################################################
# Calculate Pd and Pe for 2005
Pd = allD(d, as.matrix(dist(cbind(row,column))),rep(0.45,ncol(d)))
Pe = allE(d, env,method = "mindist")
Ps = Pd * Pe
#############################################################

#############################################################
# Assess predictions at the site level
AS1 = AssessSpecies(d,d2,Pd,T)
median(AS1$AUC,na.rm=T)

#Optimize dispersal ability
iAUC = NULL
for(i in seq(0.05,2,by = 0.05))
	{
	Pd = allD(d, as.matrix(dist(cbind(row,column))),rep(i,ncol(d)))
	AS1 = AssessSpecies(d,d2,Pd,F)
	iAUC[which((seq(0.05,2,by = 0.05)==i))] = mean(AS1$AUC,na.rm=T)
	}
plot(seq(0.05,2,by = 0.05),iAUC)
#Hmm. the best one is when dispersal = 0.45. That seems surprisingly small


#Try out all combinations of pool definitions, environment and dispersal
Pd = allD(d, as.matrix(dist(cbind(row,column))),dispAbility)
Pe = allE(d, env,method = "beals")
Ps = Pd * Pe

AS2 = AssessSpecies(d,d2,Ps,F)
mean(AS2$AUC,na.rm=T)

AS3 = AssessSites(d,d2,Ps)
mean(AS3$AUC,na.rm=T)


AssessSpecies = function(d,d2,P,plot = T)
	{
	if(plot){par(mfrow = c(5,5),mar = c(1,1,1,1))}

	AUC = NULL
	med00 = NULL
	med01 = NULL
	p = NULL
	for(i in 1:ncol(d))
		{
		s00 = which(d[,i] == 0 & d2[,i]==0)
		s01 = which(d[,i] == 0 & d2[,i]>=1)

		if(length(s00) >= 2 & length(s01) >= 2)
			{
			test = t.test(P[s00,i],P[s01,i])
			p[i] = test$p.value
			med00[i] = median(P[s00,i])
			med01[i] = median(P[s01,i])
			AUC[i] = mean(sample(P[s01,i],1000000,replace=T) > sample(P[s00,i],1000000,replace=T))

			if(plot)
				{
				den1 = density(P[s00,i],bw = 0.06)
				den2 = density(P[s01,i],bw = 0.06)

				plot(den2,main = paste(names(d)[i],round(AUC[i],3)),col = "red",lwd = 2,xlim = c(0,1),axes = F)
				points(den1,type = "l",lwd = 2,col = "grey")
				points(median(P[s00,i]),0,pch = 17,cex = 3,col = "grey")
				points(median(P[s01,i]),0,pch = 17,col = "red",cex = 2)
				box()
				axis(1,labels = c(0,1),at = c(0,1))
				}
			}
		}
	return(data.frame(med00,med01,p,AUC))
	}

AssessSites = function(d,d2,P)
	{
	p = NULL
	med00 = NULL
	med01 = NULL
	AUC = NULL
	for(i in 1:nrow(d))
		{
		s00 = which(d[i,] == 0 & d2[i,] == 0)
		s01 = which(d[i,] == 0 & d2[i,] >= 1)

		if(length(s01) > 1)
			{
			test = t.test(P[i,s00],P[i,s01])
			p[i] = test$p.value
			med00[i] = median(P[i,s00])
			med01[i] = median(P[i,s01])
			AUC[i] = mean(sample(P[i,s01],1000000,replace=T) > sample(P[i,s00],1000000,replace=T))
			}
		}	
	return(data.frame(med00,med01,p,AUC))
	}


