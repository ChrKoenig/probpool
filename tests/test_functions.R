source("R/probpool_functions.R")


load("data/Ranunculaceae_dispersal_ability.RData")
load("data/Ranunculaceae_occurrences.RData")
load("data/Ranunculaceae_env_prob.RData")


par(mfrow=c(1,1))

disp.rst.stack <- disp_pool(occurrence.surfaces = occ.rst.stack[[1:5]], disp.ability = 50)
plot(disp.rst.stack[[5]])

# save(disp.rst.stack,file="Ranunculaceae_disp_prob")

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

disp.simu.stack <- disp_pool(occurrence.surfaces = simu_occ_rst_stack, disp.ability = 1, longlat = FALSE)
plot(disp.simu.stack[[5]])

values(disp.simu.stack[[5]])

plot(simu_occ_rst_stack[[2]])

###Checking bio pool for hte simulated data
load("data/Simu_InteractionMatrix.RData")
cooccur.probs2<-as.matrix(cooccur.probs)
summary(cooccur.probs2)
range(cooccur.probs2)
dim( cooccur.probs2)

int.matrix <- cooccur.probs2

occurrence.surfaces <- simu_occ_rst_stack
Simu_bio.rst.stack <- bio_pool(occurrence.surfaces, int.matrix, abundances=)

occurrence.surfaces2<-occurrence.surfaces[occurrence.surfaces,]

head(values(occurrence.surfaces))
head(values(Simu_bio.rst.stack))

plot(Simu_bio.rst.stack[[5]])


river <- raster(nrows=15, ncols=15, vals = t(matrix(c(rep(1,6*15),rep(0.5,15),rep(0,15),rep(0.5,15),rep(1,6*15)),nrow=15,ncol=15)),xmn=0, xmx=15, ymn=0, ymx=15)

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


plot(simu_occ_rst_stack[[1]])


disp.simu.stack <- disp_pool(occurrence.surfaces = simu_occ_rst_stack, cond.surfaces=river, disp.ability = 10, longlat = FALSE)
plot(simu_occ_rst_stack[[5]])
plot(river)
plot(disp.simu.stack[[5]])




load("occ_rst_stack.RData")
load("suit_rst_stack.RData")
load("disp_rst_stack.RData")
load("dispersal_ability.RData")


names(suit.rst.stack[[1]])

interactions <- matrix(runif(min = -1, max = 1, n = 51^2),nrow = 51, ncol = 51)
colnames(interactions) <- names(dispersal.ability)
rownames(interactions) <- names(dispersal.ability)
diag(interactions) <- 0
int.matrix <- interactions

occurrence.surfaces <- suit.rst.stack*disp.rst.stack





bio.rst.stack <- bio_pool(occurrence.surfaces, int.matrix)

save(bio.rst.stack, file="bio_rst_stack.RData")


plot(occ.rst.stack[[1]])
plot(suit.rst.stack[[1]])
plot(disp.rst.stack[[1]])
plot(bio.rst.stack[[1]])

plot(suit.rst.stack[[1]]*disp.rst.stack[[1]]*interactions[[1]])


names(occ.rst.stack) <- names(dispersal.ability)
names(bio.rst.stack)

load("occ_rst_stack.RData")
