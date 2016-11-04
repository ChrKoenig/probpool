source("R/probpool_functions.R")


load("data/Ranunculaceae_dispersal_ability.RData")
load("data/Ranunculaceae_occurrences.RData")
load("data/Ranunculaceae_env_prob.RData")


par(mfrow=c(1,1))

disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack, disp.ability = 50)
plot(disp.rst.stack[[5]])

save(disp.rst.stack,file="data/Ranunculaceae_disp_prob.RData")

disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack[[4:5]], disp.ability = 2500, cond.surfaces=suit.rst.stack[[4:5]])
plot(disp.rst.stack[[2]])


occurrence.surfaces = occ.rst.stack[[4:5]]
disp.ability = 2500
method=c("negexp","fattail")
cond.surfaces=suit.rst.stack[[4:5]]
longlat=TRUE

class(suit.rst.stack[[4:5]])
class(cond.surfaces)


disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack[[1:5]], disp.ability=dispersal.ability[1:5])
plot(disp.rst.stack[[5]])

disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack[[1:5]], disp.ability=dispersal.ability[1:5], cond.surfaces=suit.rst.stack[[1:5]])
plot(disp.rst.stack[[5]])


plot(suit.rst.stack[[5]])



load("data/Simu_occurrences.RData")
load("data/Simu_env_prob.RData")

disp.simu.stack <- disp_pool(occurrence.surfaces = simu_occ_rst_stack, disp.ability = 1, longlat = FALSE)
plot(disp.simu.stack[[5]])

values(disp.simu.stack[[5]])

plot(simu_occ_rst_stack[[2]])

###Checking bio pool for hte simulated data
load("data/Simu_InteractionMatrix.RData")

int.matrix <- Simu_InteractionMatrix
mean(Simu_InteractionMatrix)

occurrence.surfaces <- simu_occ_rst_stack
Simu_bio.rst.stack <- bio_pool(occurrence.surfaces, int.matrix, abundance=FALSE)

head(values(occurrence.surfaces))
head(values(Simu_bio.rst.stack))
range(values(Simu_bio.rst.stack[[5]]))

plot(Simu_bio.rst.stack[[5]])

par(mfrow=c(2,3))

plot(simu_occ_rst_stack[[5]])
plot(Simu_env_prob[[5]])
plot(disp.simu.stack [[5]])
plot(Simu_bio.rst.stack[[5]])
plot(Simu_env_prob[[5]]*disp.simu.stack[[5]])


plot(Simu_env_prob[[5]]*disp.simu.stack[[5]]*Simu_bio.rst.stack[[5]])



river <- raster(nrows=15, ncols=15, vals = t(matrix(c(rep(1,6*15),rep(0.5,15),rep(0,15),rep(0.5,15),rep(1,6*15)),nrow=15,ncol=15)),xmn=-1, xmx=16, ymn=-1, ymx=16)
river <- suit.rst.stack[[1]]
extent(river)

# replace NAs by 0
river[is.na(river)] <- 1/100

# truncation: replace <1/100 by 0                                   
river[river<1/100] <- 1/100
plot(river)

occurrences <- rasterToPoints(simu_occ_rst_stack)
summary(occurrences)

spec.trans <- transition(river, mean, 8)
image(transitionMatrix(spec.trans))

# geocorrection
spec.trans <- geoCorrection(spec.trans, type="c")

# calculate commute distances
distances <- commuteDistance(spec.trans, occurrences[,c(1:2)])/2
distances <- as.matrix(distances)

commuteDistance(spec.trans, matrix(c(14.5,0.5,6,7),2,2))
commuteDistance(spec.trans, matrix(c(13.9,6,0,7),2,2))

commuteDistance(spec.trans, matrix(c(47.21,47.3,5.8333331,5.9),2,2))
commuteDistance(spec.trans, matrix(c(5.84, 6, 47.21, 48),2,2))

extent(river)














plot(simu_occ_rst_stack[[1]])


disp.simu.stack <- disp_pool(occurrence.surfaces = simu_occ_rst_stack, cond.surfaces=river, disp.ability = 100, longlat = FALSE)
plot(simu_occ_rst_stack[[5]])
plot(river)
plot(disp.simu.stack[[5]])







### Bio pool

load("data/Ranunculaceae_dispersal_ability.RData")
load("data/Ranunculaceae_occurrences.RData")
load("data/Ranunculaceae_env_prob.RData")
load("data/Ranunculaceae_disp_prob.RData")

    
names(suit.rst.stack[[1]])


# create interaction matrix(
interactions <- matrix(runif(min = -1, max = 1, n = 51^2),nrow = 51, ncol = 51)
colnames(interactions) <- names(dispersal.ability)
rownames(interactions) <- names(dispersal.ability)
diag(interactions) <- 0
int.matrix <- interactions

occurrence.surfaces <- suit.rst.stack*disp.rst.stack

bio.rst.stack <- bio_pool(occurrence.surfaces, int.matrix)

save(bio.rst.stack, file="data/Ranunculaceae_bio_prob.RData")



plot(occ.rst.stack[[1]])
plot(suit.rst.stack[[1]])
plot(disp.rst.stack[[1]])
plot(bio.rst.stack[[1]])

plot(suit.rst.stack[[1]]*disp.rst.stack[[1]]*interactions[[1]])


names(occ.rst.stack) <- names(dispersal.ability)
names(bio.rst.stack)

