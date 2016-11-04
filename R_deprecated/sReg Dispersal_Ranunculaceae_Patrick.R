setwd("D:\\Projects\\P24_SpeciesPool\\statistics\\")

library(maptools);library(gdistance);library(raster);library(rgdal);library(RColorBrewer);library(classInt)


### read in shapefile:
mtb.shp <- readShapeSpatial("Germany_data/Data_Anna_final/MTB_species_final_abbr.shp", proj4string=CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +datum=potsdam +units=m +no_defs +ellps=bessel +towgs84=606.0,23.0,413.0"))
mtb.wgs84 <- spTransform(mtb.shp,CRS("+proj=longlat +datum=WGS84"))
save(mtb.wgs84, file="mtb.wgs84.RData")


head(mtb.wgs84@data)
dim(mtb.wgs84@data)

# create a raster
extent.rst <- extent(mtb.wgs84)
extent.rst@xmin <- 35*(1/6)
extent.rst@xmax <- 91*(1/6)
extent.rst@ymin <- 47.2
extent.rst@ymax <- 55.1

#(extent.rst@xmax - extent.rst@xmin)/(1/6)
#(extent.rst@ymax - extent.rst@ymin)/(1/10)

extentraster <- raster(ext=extent.rst,nrows=79,ncols=56)

# plot distribution of exemplary species
mtb.rst <- rasterize(mtb.wgs84,extentraster,field="Acon_nape")

plot(mtb.rst)
plot(mtb.wgs84, add=TRUE)


# get species names
species.abbr <- names(mtb.wgs84@data)[9:length(names(mtb.wgs84@data))]

### load species of interest table for full species names
species <- read.csv("Germany_data/Data_Anna_final/Species_of_interest_names_final.csv")

### create species table for further steps
species <- species[species$Abbreviation %in% species.abbr,]
species <- unique(species[,c(4,5,6)])
dim(species)

### load suitability scores
suit <- read.table("Germany_data/Z_Current.glm.2011.txt",header=TRUE,sep="\t")
suit$WNAME <- gsub("\\.x","",suit$WNAME)

# subset suit for Ranunculaceae species of interest
species[!species$New_species_name %in% suit$WNAME,]

species <- species[species$New_species_name %in% suit$WNAME,]
suit <- suit[suit$WNAME %in% species$New_species_name,]
dim(suit)

rownames(suit) <- suit$WNAME
suit <- suit[,-c(1:7)]

suit <- t(suit)
suit <- as.data.frame(suit)

suit$TK_NR <- gsub("X","",rownames(suit))

head(suit)
suit[1:5,1:10]


### merge suitability to grid
mtb.wgs84.suit <- mtb.wgs84
mtb.wgs84.suit@data <- mtb.wgs84.suit@data[,c(1:2)]

mtb.wgs84.suit@data$orig.order <- c(1:nrow(mtb.wgs84.suit@data))
mtb.wgs84.suit@data  <- merge(mtb.wgs84.suit@data,suit,by="TK_NR",all.x=TRUE,sort=FALSE) 
mtb.wgs84.suit@data <- mtb.wgs84.suit@data[order(mtb.wgs84.suit@data$orig.order),]
head(mtb.wgs84.suit@data)

dim(mtb.wgs84.suit@data)
dim(mtb.wgs84@data)



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

calcD = function(d, k, method = "negexp")
	{
	if(method == "negexp") { return (1 - prod(1-exp(-1*d/k)))}
	if(method == "fattail") { return (1 - prod(1/d^k))}
	}








### gridcells that are not NA in distribution data
gridcells <- mtb.wgs84@data$TK_NR[!is.na(mtb.wgs84@data$Acon_nape)]

# get coordinates of all gridcells 
gridcells.coords <- coordinates(mtb.wgs84)[!is.na(mtb.wgs84@data$Acon_nape),]
rownames(gridcells.coords) <- gridcells

# load dispersal ability
disp.cap <- read.csv("data_Karsten\\RanunculaceaePFTs.csv")
disp.cap$Species <- gsub(" ","_",disp.cap$Species)
disp.cap$Species[!disp.cap$Species %in% species$New_species_name]
species$New_species_name[!species$New_species_name %in% disp.cap$Species]
disp.cap$Species[disp.cap$Species=="Aconitum_variegatum"] <- "Aconitum_degenii"
species <- merge(species,disp.cap,by.x="New_species_name",by.y="Species",all.x=TRUE,sort=FALSE)
species$SumRdisp <- species$SumRdisp/max(species$SumRdisp)+0.5
species$LDD <- species$LDD+0.5

range(species$SumRdisp)
range(species$LDD)

dim(species)
dim(disp.cap)


# prepare results data.frame
dispersal.rst <- as.data.frame(matrix(NA,length(gridcells),nrow(species)))
rownames(dispersal.rst) <- gridcells
colnames(dispersal.rst) <- species$Abbreviation

dispersal.geo <- dispersal.rst
dispersal.rst.cap <- dispersal.rst
dispersal.geo.cap <- dispersal.rst
gamma <- dispersal.rst

# calculate geographic distances
distances.geo <- spDists(gridcells.coords, gridcells.coords, longlat = TRUE)
rownames(distances.geo) <- gridcells
colnames(distances.geo) <- gridcells

suit.rst.stack <- lapply(species$New_species_name, function(x) rasterize(mtb.wgs84.suit,extentraster,field=as.character(x)))
suit.rst.stack <- stack(suit.rst.stack, layers= as.character(species$Abbreviation))
save(suit.rst.stack, file="suit_rst_stack.RData")

occ.rst.stack <- lapply(species$Abbreviation, function(x) rasterize(mtb.wgs84,extentraster,field=as.character(x)))
occ.rst.stack <- stack(occ.rst.stack, layers=as.character(species$Abbreviation))
save(occ.rst.stack, file="occ_rst_stack.RData")

dispersal.ability <- species$LDD
names(dispersal.ability) <- species$Abbreviation
save(dispersal.ability, file = "dispersal_ability.RData")


#for (i in 1:nrow(species)){
for (i in 1:5){
i=1

spec.rst <- rasterize(mtb.wgs84.suit,extentraster,field=as.character(species$New_species_name[i]))
#plot(spec.rst)

# replace NAs by 0
spec.rst[is.na(spec.rst)] <- 1/100

# truncation: replace NAs by 0
spec.rst[spec.rst<1/100] <- 1/100

# create transition raster
spec.trans <- transition(spec.rst, mean, 8)

# geocorrection
spec.trans <- geoCorrection(spec.trans, type="c")

# calculate commute distances
distances.resist <- commuteDistance(spec.trans, gridcells.coords)
distances.resist <- as.matrix(distances.resist)

# get gridcells the species occurs in
spec.range <- mtb.wgs84@data$TK_NR[which(mtb.wgs84@data[,as.character(species$Abbreviation[i])]==1)]


# loop trough gridcells
for (k in 1:length(gridcells)){
#k=5
print(c(i,k))

#extract distances for focal gridcell to occurrance gridcells
distance.rst <- distances.resist[as.character(gridcells[k]),as.character(spec.range)]
distance.geo <- distances.geo[as.character(gridcells[k]),as.character(spec.range)]

dispersal.rst[k,i] <- calcD(distance.rst/5000, 1, method = "negexp")
dispersal.geo[k,i] <- calcD(distance.geo/50, 1, method = "negexp")
dispersal.rst.cap[k,i] <- calcD(distance.rst/5000, species$LDD[i], method = "negexp")
dispersal.geo.cap[k,i] <- calcD(distance.geo/50, species$LDD[i], method = "negexp")
gamma[k,i] <- ifelse(min(distances.geo[as.character(gridcells[k]),as.character(spec.range)])<=100,1,0)
}
}

save(dispersal.rst, file="dispersal.rst.RData")
save(dispersal.geo, file="dispersal.geo.RData")
save(dispersal.rst.cap, file="dispersal.rst.cap.RData")
save(dispersal.geo.cap, file="dispersal.geo.cap.RData")
save(gamma, file="gamma.RData")




#plot
plot(spec.rst)
plot(mtb.wgs84, add=TRUE)

points(spec.range,cex=0.5, pch=16, col="blue")
points(coord[1,1],coord[1,2],cex=0.5, pch=16, col="red")

AtoB <- shortestPath(spec.trans, coord[1,], coord[2,], output="SpatialLines")
AtoC <- shortestPath(spec.trans, coord[1,], coord[ceiling(nrow(coord)/2),], output="SpatialLines")
AtoD <- shortestPath(spec.trans, coord[1,], coord[nrow(coord),], output="SpatialLines")

plot(AtoB,add=TRUE)
plot(AtoC,add=TRUE,col="red")
plot(AtoD,add=TRUE,col="blue")





load("gamma.RData")

### alpha 
alpha <- mtb.wgs84@data[!is.na(mtb.wgs84@data$Acon_nape),]
rownames(alpha) <- alpha$TK_NR
alpha <- alpha[,as.character(species$Abbreviation)]
head(alpha)
head(gamma)
dim(alpha)
dim(gamma)
save(alpha, file="alpha.RData")


### model suit
env.modelsuit <- suit
rownames(env.modelsuit) <- gsub("X","",rownames(env.modelsuit))
env.modelsuit <- env.modelsuit[,as.character(species$New_species_name)]
env.modelsuit <- env.modelsuit[as.character(gridcells),]
colnames(env.modelsuit) <- as.character(species$Abbreviation)
head(env.modelsuit)
dim(env.modelsuit)
head(alpha)
dim(alpha)
save(env.modelsuit, file="env.modelsuit.RData")
                           

### Beals smoothing
env.beals <- read.csv("data_Karsten\\RanunculBeals.csv")
rownames(env.beals) <- env.beals$TK_NR
env.beals <- env.beals[as.character(gridcells),as.character(species$Abbreviation)]
head(env.beals)
dim(env.beals)
head(alpha)
dim(alpha)
save(env.beals, file="env.beals.RData")




load("alpha.RData")
load("gamma.RData")
load("dispersal.rst.RData")
load("dispersal.geo.RData")
load("dispersal.rst.cap.RData")
load("dispersal.geo.cap.RData")
load("env.modelsuit.RData")
load("env.beals.RData")

# combinations
comb.rst.modelsuit     <- dispersal.rst * env.modelsuit
comb.rst.beals         <- dispersal.rst * env.beals
comb.geo.modelsuit     <- dispersal.geo * env.modelsuit
comb.geo.beals         <- dispersal.geo * env.beals
comb.rst.cap.modelsuit <- dispersal.rst.cap * env.modelsuit
comb.rst.cap.beals     <- dispersal.rst.cap * env.beals
comb.geo.cap.modelsuit <- dispersal.geo.cap * env.modelsuit
comb.geo.cap.beals     <- dispersal.geo.cap * env.beals

head(comb.geo.cap.beals)


pool.list <- list(alpha,gamma,dispersal.rst,dispersal.geo,dispersal.rst.cap,dispersal.geo.cap,
                  env.modelsuit,env.beals,comb.rst.modelsuit,comb.rst.beals,comb.geo.modelsuit,
                  comb.geo.beals,comb.rst.cap.modelsuit,comb.rst.cap.beals,comb.geo.cap.modelsuit,comb.geo.cap.beals)

pool.names <- c("alpha","gamma","dispersal.rst","dispersal.geo","dispersal.rst.cap","dispersal.geo.cap",
                "env.modelsuit","env.beals","comb.rst.modelsuit","comb.rst.beals","comb.geo.modelsuit",
                "comb.geo.beals","comb.rst.cap.modelsuit","comb.rst.cap.beals","comb.geo.cap.modelsuit",
                "comb.geo.cap.beals")

save(pool.list, pool.names, file="species.pool.RData") 





### plot
dev.off()
windows(width=20,height=30)
#for (i in (1:nrow(species))){
for (i in (1:5)){

par(mfrow=c(4,4),mar=c(2,2,2,0))

for (k in (1:length(pool.list))){

### merge probabilities to grid
pool.values <- pool.list[[k]]
pool.values$TK_NR <- rownames(pool.values)

mtb.wgs84.pool <- mtb.wgs84
mtb.wgs84.pool@data <- mtb.wgs84.pool@data[,c(1:2)]

mtb.wgs84.pool@data$orig.order <- c(1:nrow(mtb.wgs84.pool@data))
mtb.wgs84.pool@data  <- merge(mtb.wgs84.pool@data,pool.values,by="TK_NR",all.x=TRUE,sort=FALSE) 
mtb.wgs84.pool@data <- mtb.wgs84.pool@data[order(mtb.wgs84.pool@data$orig.order),]
head(mtb.wgs84.pool@data)

my.pal <- brewer.pal(9,"Greens")

legendraster <- extentraster
legendraster[] <- rep(c(0,1),ncell(legendraster)/2) 



spec.rst <- rasterize(mtb.wgs84.pool,extentraster,field=species$Abbreviation[i])
plot(spec.rst, breaks=seq(0,1,length.out=10), col=my.pal, legend=F, box=F)
#plot(mtb.wgs84, add=TRUE)
points(coordinates(mtb.wgs84)[mtb.wgs84@data[,as.character(species$Abbreviation[i])]==1,],cex=0.25, pch=16, col="blue")
title(paste(pool.names[k]," ",species$Abbreviation[i],sep=""),cex.main=1)


}
savePlot(paste("plots//",species$Abbreviation[i],".pdf",sep=""),type="pdf")

}             










# plot sums 

dev.off()
windows(width=20,height=30)
par(mfrow=c(4,4),mar=c(2,2,2,0))
load("mtb.wgs84.window.RData")

for (k in (1:length(pool.list))){

### merge Dispersal probability to grid
pool.values <- pool.list[[k]]
pool.values$sums <- rowSums(pool.values)
pool.values$TK_NR <- rownames(pool.values)

mtb.wgs84.pool <- mtb.wgs84
mtb.wgs84.pool@data <- mtb.wgs84.pool@data[,c(1:2)]

mtb.wgs84.pool@data$orig.order <- c(1:nrow(mtb.wgs84.pool@data))
mtb.wgs84.pool@data  <- merge(mtb.wgs84.pool@data,pool.values,by="TK_NR",all.x=TRUE,sort=FALSE) 
mtb.wgs84.pool@data <- mtb.wgs84.pool@data[order(mtb.wgs84.pool@data$orig.order),]
head(mtb.wgs84.pool@data)

my.pal <- c("#F2F2F2",brewer.pal(7,"BuGn")[3:7],"#1A1A1A")
my.pal <- colorRampPalette(my.pal)
my.pal <- my.pal(20)

library(colorRamps)
my.pal <- matlab.like(20)

legendraster <- extentraster
legendraster[] <- rep(c(0,1),ncell(legendraster)/2) 

spec.rst <- rasterize(mtb.wgs84.pool,extentraster,field="sums")
plot(spec.rst, breaks=seq(0,50,length.out=21), col=my.pal, legend=F, box=F)
plot(mtb.wgs84.window,add=TRUE)
title(pool.names[k],cex.main=1)
}
savePlot("plots//sums.pdf",type="pdf")






# Figure paper

dev.off()
windows(width=18/2.54,height=4/2.54)
par(mfrow=c(1,6),mar=c(0,0,1,0.5),oma=c(0,0,0,2.5))
load("mtb.wgs84.window.RData")

selection <- c("alpha","gamma","dispersal.geo","dispersal.rst.cap","comb.geo.beals","comb.rst.cap.modelsuit")
                                                                                            


for (k in (1:length(selection))){

### merge Dispersal probability to grid
pool.values <- pool.list[[which(pool.names==selection[k])]]
pool.values$sums <- rowSums(pool.values)
pool.values$TK_NR <- rownames(pool.values)

mtb.wgs84.pool <- mtb.wgs84
mtb.wgs84.pool@data <- mtb.wgs84.pool@data[,c(1:2)]

mtb.wgs84.pool@data$orig.order <- c(1:nrow(mtb.wgs84.pool@data))
mtb.wgs84.pool@data  <- merge(mtb.wgs84.pool@data,pool.values,by="TK_NR",all.x=TRUE,sort=FALSE) 
mtb.wgs84.pool@data <- mtb.wgs84.pool@data[order(mtb.wgs84.pool@data$orig.order),]
head(mtb.wgs84.pool@data)

my.pal <- c("#F2F2F2",brewer.pal(7,"BuGn")[3:7],"#1A1A1A")
my.pal <- colorRampPalette(my.pal)
my.pal <- my.pal(10)

library(colorRamps)
my.pal <- matlab.like(10)

legendraster <- extentraster
legendraster[] <- rep(c(0,1),ncell(legendraster)/2) 

spec.rst <- rasterize(mtb.wgs84.pool,extentraster,field="sums")
plot(spec.rst, breaks=seq(0,50,length.out=11), col=my.pal, legend=F, axes=FALSE, box=FALSE)
plot(mtb.wgs84.window,add=TRUE)
title(c("Alpha","Gamma",expression("P"[D-D]),expression("P"[D-R-T]),expression("P"[D-D_E-B]),expression("P"[D-R-T_E-R]))[k],cex.main=1)
}
legendraster <- extentraster
legendraster[] <- rep(c(0,1),ncell(legendraster)/2) 
plot(legendraster,breaks=seq(0,50,length.out=11), smallplot=c(0.93, 0.99, 0.1, 0.8),col=my.pal, legend.only=TRUE, line=1, legend.width=1,axis.args=list(cex.axis=0.7),legend.args=list(text="Pool size", font=2, side=4, line=1.5, cex=0.7))
savePlot("plots//Fig.pdf",type="pdf")

                                            




# plot histograms

dev.off()
windows(width=20,height=30)
par(mfrow=c(4,4),mar=c(2,2,2,0))

for (k in (1:length(pool.list))){
pool.values <- pool.list[[k]]
hist(as.matrix(pool.values),main=pool.names[k])
}
savePlot("plots//histograms.pdf",type="pdf")



# plot sum histograms

dev.off()
windows(width=20,height=30)
par(mfrow=c(4,4),mar=c(2,2,2,0))

for (k in (1:length(pool.list))){
pool.values <- pool.list[[k]]
pool.values$sums <- rowSums(pool.values)
hist(pool.values$sums,main=pool.names[k])
}
savePlot("plots//histograms_sums.pdf",type="pdf")











### identify cells at border

### gridcells that are not NA in distribution data
gridcells <- mtb.wgs84@data$TK_NR[!is.na(mtb.wgs84@data$Acon_nape)]

# get coordinates of all gridcells 
gridcells.coords <- coordinates(mtb.wgs84)[!is.na(mtb.wgs84@data$Acon_nape),]
rownames(gridcells.coords) <- gridcells

# calculate geographic distances
distances.geo <- spDists(gridcells.coords, gridcells.coords, longlat = TRUE)
rownames(distances.geo) <- gridcells
colnames(distances.geo) <- gridcells

nneigh <- NA

for (k in 1:length(gridcells)){
nneigh[k] <- length(which(distances.geo[as.character(gridcells[k]),]<=100))
}

save(nneigh, file="nneigh.RData")

nneigh.frame <- cbind(gridcells,nneigh)


mtb.wgs84.nneigh <- mtb.wgs84
mtb.wgs84.nneigh@data <- mtb.wgs84.nneigh@data[,c(1:2)]

mtb.wgs84.nneigh@data$orig.order <- c(1:nrow(mtb.wgs84.nneigh@data))
mtb.wgs84.nneigh@data  <- merge(mtb.wgs84.nneigh@data,nneigh.frame,by.x="TK_NR",by.y="gridcells",all.x=TRUE,sort=FALSE) 
mtb.wgs84.nneigh@data <- mtb.wgs84.nneigh@data[order(mtb.wgs84.nneigh@data$orig.order),]
head(mtb.wgs84.nneigh@data)


plot(mtb.wgs84)
plot(mtb.wgs84[which(mtb.wgs84.nneigh@data$nneigh>200),],col="red",add=TRUE)
#plot(mtb.wgs84[which(mtb.wgs84.nneigh@data$nneigh>210),],col="orange",add=TRUE)
#plot(mtb.wgs84[which(mtb.wgs84.nneigh@data$nneigh>220),],col="green",add=TRUE)
plot(mtb.wgs84[which(mtb.wgs84.nneigh@data$nneigh>229),],col="darkgreen",add=TRUE)
#plot(mtb.wgs84[which(mtb.wgs84.nneigh@data$nneigh>240),],col="blue",add=TRUE)
#plot(mtb.wgs84[which(mtb.wgs84.nneigh@data$nneigh>250),],col="black",add=TRUE)

### use >229


### create polygon for correlation window
mtb.wgs84.window <- mtb.wgs84[which(mtb.wgs84.nneigh@data$nneigh>229),]
mtb.wgs84.window <- unionSpatialPolygons(mtb.wgs84.window,mtb.wgs84.window@data$COUNT)
plot(mtb.wgs84.window)

save(mtb.wgs84.window, file="mtb.wgs84.window.RData")




################################################################################
###  correlations
load("species.pool.RData") 
pool.size <- as.data.frame(sapply(pool.list, rowSums))
colnames(pool.size) <- pool.names

#exclude border gridcells
load("nneigh.RData")
pool.size.n <- pool.size[nneigh>229,]
dim(pool.size.n)

cor.matrix <- cor(pool.size.n)  ### correlate poolsize values of gridcells among pool methods




### Mantel tests
library(vegan)
pool.list.n <- pool.list

#exclude border gridcells
for (i in 1:length(pool.list.n)){ 
    pool.list.n[[i]] <- pool.list.n[[i]][nneigh>229,]
}

pool.dissim <- lapply(pool.list.n, function(x) (vegdist(x,"euclidean")))   # euclidean distance matrices for each pool method

for (i in 1:(nrow(cor.matrix)-1)){
   for (k in (1+i):nrow(cor.matrix)){ 
      cor.matrix[i,k] <- mantel(pool.dissim[[i]], pool.dissim[[k]],permutations=1)$statistic
}}

save(cor.matrix,file="cor.matrix.RData")
write.csv(cor.matrix,file="cor.matrix.csv")
################################################################################



cor.matrix[]<-1

for (i in 1:(nrow(cor.matrix)-1)){
  for (k in (1+i):nrow(cor.matrix)){ 
    cor.matrix[i,k] <- mantel(pool.dissim[[i]], pool.dissim[[k]],permutations=999)$signif
  }}

for (i in 2:nrow(cor.matrix)){
  for (k in 1:(i-1)){ 
    cor.matrix[i,k] <- cor.test(pool.size.n[,i], pool.size.n[,k])$p.value
  }}



cor.matrix
save(cor.matrix,file="cor.matrix.sign.RData")
write.csv(cor.matrix,file="cor.matrix.sign.csv")
################################################################################

library(ade4)
system.time(mantel(pool.dissim[[i]], pool.dissim[[k]],permutations=50)$signif)
system.time(mantel.rtest(pool.dissim[[i]], pool.dissim[[k]], nrepet = 50)$pvalue)

mantel(pool.dissim[[2]], pool.dissim[[4]],permutations=10)
mantel(pool.dissim[[2]], pool.dissim[[3]],permutations=10)



################################################################################
###  Community structure, Probabilty Mass Function (PMF) and Cumulative Density
###  Functions (CDF) for 1 gridcell (See Figure XPMF)
setwd("C:\\Users\\Patrick\\Projects_PD\\P24_SpeciesPool\\statistics\\")


load("species.pool.RData") 
pool.names

# choose species pool concept
pool <- pool.list[[13]] # "comb.rst.cap.modelsuit"
pool.disp <- pool.list[[4]] # "dispersal.geo"
pool.env <- pool.list[[7]] # "env.modelsuit"
pool.env.beals <- pool.list[[8]] # "env.beals"


# choose species pool grid cell: middle? TK 2520
pool <- t(pool["2520",]) # G?ttingen
pool <- pool[order(pool,decreasing=TRUE)]
names(pool) <- c(1:length(pool))

dev.off()
windows(width=18/2.54,height=5.5/2.54)
par(mfrow=c(1,3),mar=c(4,2.5,2.5,0.5),oma=c(0,2,0,0))

barplot(pool,ylim=c(0,1),main="Community structure",ylab="P",xpd=NA,xlab="Species")

pool.table <- table(rowSums(sapply(pool, function(x) rbinom(1000000,1,x))))/1000000

barplot(pool.table,ylim=c(0,1),main="PMF",xlab="Species",xpd=NA)

barplot(rbind(cumsum(pool.table),1-cumsum(pool.table)),ylim=c(0,1),main="CDF",xlab="Species",xpd=NA)

sum(pool)
savePlot("plots//PMF_CDF_comb.rst.cap.modelsuit.pdf",type="pdf")





# choose species pool grid cell: middle? TK 2520
pool.disp <- t(pool.disp["2520",]) # G?ttingen
pool.disp <- pool.disp[order(pool.disp,decreasing=TRUE)]
names(pool.disp) <- c(1:length(pool.disp))

dev.off()
windows(width=18/2.54,height=5.5/2.54)
par(mfrow=c(1,3),mar=c(4,2.5,2.5,0.5),oma=c(0,2,0,0))

barplot(pool.disp,ylim=c(0,1),main="Community structure",ylab="P",xpd=NA,xlab="Species")

pool.disp.table <- table(rowSums(sapply(pool.disp, function(x) rbinom(1000000,1,x))))/1000000

barplot(pool.disp.table,ylim=c(0,1),main="PMF",xlab="Species",xpd=NA)

barplot(rbind(cumsum(pool.disp.table),1-cumsum(pool.disp.table)),ylim=c(0,1),main="CDF",xlab="Species",xpd=NA)

sum(pool.disp)
savePlot("plots//PMF_CDF_dispersal.geo.pdf",type="pdf")




# choose species pool grid cell: middle? TK 2520
pool.env <- t(pool.env["2520",]) # G?ttingen
pool.env <- pool.env[order(pool.env,decreasing=TRUE)]
names(pool.env) <- c(1:length(pool.env))

dev.off()
windows(width=18/2.54,height=5.5/2.54)
par(mfrow=c(1,3),mar=c(4,2.5,2.5,0.5),oma=c(0,2,0,0))

barplot(pool.env,ylim=c(0,1),main="Community structure",ylab="P",xpd=NA,xlab="Species")

pool.env.table <- table(rowSums(sapply(pool.env, function(x) rbinom(1000000,1,x))))/1000000

barplot(pool.env.table,ylim=c(0,1),main="PMF",xlab="Species",xpd=NA)

barplot(rbind(cumsum(pool.env.table),1-cumsum(pool.env.table)),ylim=c(0,1),main="CDF",xlab="Species",xpd=NA)

sum(pool.env)
savePlot("plots//PMF_CDF_env.modelsuit.pdf",type="pdf")



# choose species pool grid cell: middle? TK 2520
pool.env.beals <- t(pool.env.beals["2520",]) # G?ttingen
pool.env.beals <- pool.env.beals[order(pool.env.beals,decreasing=TRUE)]
names(pool.env.beals) <- c(1:length(pool.env.beals))

dev.off()
windows(width=18/2.54,height=5.5/2.54)
par(mfrow=c(1,3),mar=c(4,2.5,2.5,0.5),oma=c(0,2,0,0))

barplot(pool.env.beals,ylim=c(0,1),main="Community structure",ylab="P",xpd=NA,xlab="Species")

pool.env.beals.table <- table(rowSums(sapply(pool.env.beals, function(x) rbinom(1000000,1,x))))/1000000

barplot(pool.env.beals.table,ylim=c(0,1),main="PMF",xlab="Species",xpd=NA)

barplot(rbind(cumsum(pool.env.beals.table),1-cumsum(pool.env.beals.table)),ylim=c(0,1),main="CDF",xlab="Species",xpd=NA)

sum(pool.env.beals)
savePlot("plots//PMF_CDF_env.beals.pdf",type="pdf")








dev.off()
windows(width=18/2.54,height=16.5/2.54)
par(mfrow=c(3,3),mar=c(4,2.5,2.5,0.5),oma=c(0,2,0,0))


# choose species pool grid cell: middle? TK 2520
pool <- t(pool["2520",]) # G?ttingen
pool <- pool[order(pool,decreasing=TRUE)]
names(pool) <- c(1:length(pool))


barplot(pool,ylim=c(0,1),main="Community structure",ylab="P",xpd=NA,xlab="Species")

pool.table <- table(rowSums(sapply(pool, function(x) rbinom(1000000,1,x))))/1000000

barplot(pool.table,ylim=c(0,1),main="PMF",xlab="Species",xpd=NA)

barplot(rbind(cumsum(pool.table),1-cumsum(pool.table)),ylim=c(0,1),main="CDF",xlab="Species",xpd=NA)

sum(pool)





# choose species pool grid cell: middle? TK 2520
pool.disp <- t(pool.disp["2520",]) # G?ttingen
pool.disp <- pool.disp[order(pool.disp,decreasing=TRUE)]
names(pool.disp) <- c(1:length(pool.disp))


barplot(pool.disp,ylim=c(0,1),main="Community structure",ylab="P",xpd=NA,xlab="Species")

pool.disp.table <- table(rowSums(sapply(pool.disp, function(x) rbinom(1000000,1,x))))/1000000

barplot(pool.disp.table,ylim=c(0,1),main="PMF",xlab="Species",xpd=NA)

barplot(rbind(cumsum(pool.disp.table),1-cumsum(pool.disp.table)),ylim=c(0,1),main="CDF",xlab="Species",xpd=NA)

sum(pool.disp)



# choose species pool grid cell: middle? TK 2520
pool.env <- t(pool.env["2520",]) # G?ttingen
pool.env <- pool.env[order(pool.env,decreasing=TRUE)]
names(pool.env) <- c(1:length(pool.env))


barplot(pool.env,ylim=c(0,1),main="Community structure",ylab="P",xpd=NA,xlab="Species")

pool.env.table <- table(rowSums(sapply(pool.env, function(x) rbinom(1000000,1,x))))/1000000

barplot(pool.env.table,ylim=c(0,1),main="PMF",xlab="Species",xpd=NA)

barplot(rbind(cumsum(pool.env.table),1-cumsum(pool.env.table)),ylim=c(0,1),main="CDF",xlab="Species",xpd=NA)

sum(pool.env)

savePlot("plots//PMF_CDF.pdf",type="pdf")







par(mfrow=c(3,4),mar=c(4,2.5,2.5,0.5),oma=c(0,2,0,0))
pool1 <- rep(0.1,10)
pool2 <- rep(0.5,10)
pool3 <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.5,0.5,0.9,0.9)
#pool3 <- pool3[order(pool3,decreasing = TRUE)]
pool4 <- c(rep(0.1,8),0.5,0.5)

barplot(pool1,ylim=c(0,1))
barplot(pool2,ylim=c(0,1))
barplot(pool3,ylim=c(0,1))
barplot(pool4,ylim=c(0,1))

dbinom(c(1:10),10,0.1)
dbinom(c(1:10),10,0.5)

apply(sapply(pool3, function(x) dbinom(c(1:length(pool3)),length(pool3),x)),1,prod) #bullshit
colSums(sapply(pool3, function(x) dbinom(c(1:length(pool3)),length(pool3),x))/length(pool3))   #bullshit

#check: http://en.wikipedia.org/wiki/Poisson_binomial_distribution
#poibin-package {poibin}
#?dpoibin

barplot(dbinom(c(1:length(pool1)),10,pool1),ylim=c(0,1))  # http://en.wikibooks.org/wiki/R_Programming/Probability_Functions/Binomial#Probability_Mass_Function
barplot(dbinom(c(1:length(pool2)),length(pool2),pool2),ylim=c(0,1))
barplot(dbinom(c(1:length(pool3)),length(pool3),pool3),ylim=c(0,1))
barplot(dbinom(c(1:length(pool4)),length(pool4),pool4),ylim=c(0,1))

barplot(pbinom(c(1:length(pool1)),length(pool1),pool1),ylim=c(0,1))
barplot(pbinom(c(1:length(pool2)),length(pool2),pool2),ylim=c(0,1))
barplot(pbinom(c(1:length(pool3)),length(pool3),pool3),ylim=c(0,1))
barplot(pbinom(c(1:length(pool4)),length(pool4),pool4),ylim=c(0,1))





par(mfrow=c(3,4))
pool1 <- rep(0.1,10)
pool2 <- rep(0.5,10)
pool3 <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.5,0.5,0.9,0.9)
#pool3 <- pool3[order(pool3,decreasing = TRUE)]
pool4 <- c(rep(0.1,8),0.5,0.5)

barplot(pool1,ylim=c(0,1))
barplot(pool2,ylim=c(0,1))
barplot(pool3,ylim=c(0,1))
barplot(pool4,ylim=c(0,1))


dbinom(c(0:length(pool1)),10,pool1)
table(rowSums(sapply(pool1, function(x) rbinom(1000000,1,x))))/1000000 #similar

barplot(table(rowSums(sapply(pool1, function(x) rbinom(1000000,1,x))))/1000000,xlim=c(1,length(pool1)),ylim=c(0,1))
barplot(table(rowSums(sapply(pool2, function(x) rbinom(1000000,1,x))))/1000000,xlim=c(1,length(pool2)),ylim=c(0,1))
barplot(table(rowSums(sapply(pool3, function(x) rbinom(1000000,1,x))))/1000000,xlim=c(1,length(pool3)),ylim=c(0,1))
barplot(table(rowSums(sapply(pool4, function(x) rbinom(1000000,1,x))))/1000000,xlim=c(1,length(pool4)),ylim=c(0,1))


barplot(rbind(cumsum(table(rowSums(sapply(pool1, function(x) rbinom(1000000,1,x))))/1000000),1-cumsum(table(rowSums(sapply(pool1, function(x) rbinom(1000000,1,x))))/1000000)),xlim=c(1,length(pool1)),ylim=c(0,1))
barplot(rbind(cumsum(table(rowSums(sapply(pool2, function(x) rbinom(1000000,1,x))))/1000000),1-cumsum(table(rowSums(sapply(pool2, function(x) rbinom(1000000,1,x))))/1000000)),xlim=c(1,length(pool2)),ylim=c(0,1))
barplot(rbind(cumsum(table(rowSums(sapply(pool3, function(x) rbinom(1000000,1,x))))/1000000),1-cumsum(table(rowSums(sapply(pool3, function(x) rbinom(1000000,1,x))))/1000000)),xlim=c(1,length(pool3)),ylim=c(0,1))
barplot(rbind(cumsum(table(rowSums(sapply(pool4, function(x) rbinom(1000000,1,x))))/1000000),1-cumsum(table(rowSums(sapply(pool4, function(x) rbinom(1000000,1,x))))/1000000)),xlim=c(1,length(pool4)),ylim=c(0,1))








