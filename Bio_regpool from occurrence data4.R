
setwd("C:/Users/Juliano/Documents/SpeciesPool")

library(foreign)

#Species Traits
#This files should contain the species IDs in the first column (numeric)
#other columns give the traits
SpeciesTraits<-read.table(file="SpeciesTraits.TXT",sep = "", header=T)
summary(SpeciesTraits)

#Species Co-occurrences in community simulations 
#First two columns are X and Y coordinates
#Third column is the species ID
#The other columns are abundances
#TO-DO: presence as well?
Abundances_in_Community_perCell_perSp_perStage<-read.table(file="Abundances_in_Community_perCell_perSp_perStage.TXT",sep = "", header=T)
summary(Abundances_in_Community_perCell_perSp_perStage)

#Species potential co-occurrences obtained from simulating each species alone
#Similar structure than previous 
#Only possible for simulated or experimental data
#It gives the total number of species in the pool (in this data 400)
Abundances_Alone_perCell_perSp_perStage<-read.table(file="Abundances_Alone_perCell_perSp_perStage.TXT",sep = "", header=T)
summary(Abundances_Alone_perCell_perSp_perStage)
#GammaRichness is the total number of species considered, the total pool:
GammaRichness<-length(levels(as.factor(as.character(Abundances_Alone_perCell_perSp_perStage[,3]))))

#####
#Merging trait data into the abundance matrices via species ID
Abundances_in_Community_Traits_Merged<-merge(Abundances_in_Community_perCell_perSp_perStage, SpeciesTraits, by = intersect(names(Abundances_in_Community_perCell_perSp_perStage),names(SpeciesTraits)))
Abundances_Alone_Traits_Merged<-merge(Abundances_Alone_perCell_perSp_perStage, SpeciesTraits, by = intersect(names(Abundances_Alone_perCell_perSp_perStage),names(SpeciesTraits)))

#Adding columns for the probabilities:
Abundances_in_Community_Traits_Merged$DemoProbTemp<-0
Abundances_in_Community_Traits_Merged$DemoInteractionProbTemp<-0
Abundances_in_Community_Traits_Merged$InteractionTemp<-0
#
Abundances_Alone_Traits_Merged$DemoProbTemp<-0
Abundances_Alone_Traits_Merged$DemoInteractionProbTemp<-0
Abundances_Alone_Traits_Merged$InteractionTemp<-0

#calculating the probability of surviving for each species
#at each cell
MaxNumber_UniqueEnvironments<-15
ReplicateNumber_perUniqueEnvironment<-15

for (i in 1:GammaRichness){ #GammaRichness is the total number of species considered, the total pool
  SpeciesTraits2<-SpeciesTraits[SpeciesTraits[,1]==i,]
  Abundances_Alone_Traits_Merged2<-Abundances_Alone_Traits_Merged[Abundances_Alone_Traits_Merged$SpeciesID==i,] 
  Abundances_in_Community_Traits_Merged2<-Abundances_in_Community_Traits_Merged[Abundances_in_Community_Traits_Merged$SpeciesID==i,]
  
  #if the species survived in the community
  if(nrow(Abundances_in_Community_Traits_Merged2)>0){
    #for each y value (altitudinal level), there are x replicates (e.g. y is a unique environmental state, in this case, 15)
    #To-DO: make y being the environmental ID!!!
    for (y in 1:MaxNumber_UniqueEnvironments){
      Abundances_in_Community_Traits_Merged3<-Abundances_in_Community_Traits_Merged2[Abundances_in_Community_Traits_Merged2$Y==y,]
      Abundances_Alone_Traits_Merged3<-Abundances_Alone_Traits_Merged2[Abundances_Alone_Traits_Merged2$Y==y,]
      #TO-DO: make only one relevant abundance column
      #TO-DO: use only presence as well?
      #for perennials:
      if(Abundances_in_Community_Traits_Merged3$Annual[1]==0) {  
        #this is the probability of surviving in this environment when the species is alone (env + demographic + dist filter):
        ProbAlone<-nrow(Abundances_Alone_Traits_Merged3[which(Abundances_Alone_Traits_Merged3$AdultAbundance>0),])/ReplicateNumber_perUniqueEnvironment
        #this is the probability of surviving in this environment when within diffuse resource competition with other co-occurring species (env + demographic + dist +interaction filter)
        ProbComm<-nrow(Abundances_in_Community_Traits_Merged3[which(Abundances_in_Community_Traits_Merged3$AdultAbundance>0),])/ReplicateNumber_perUniqueEnvironment           
        
        Abundances_in_Community_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)-2]<-ProbAlone
        Abundances_in_Community_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)-1]<-ProbComm
        
        #this is rule making sure about the 
        #the probability of surviving in this environment when within diffuse resource competition with other co-occurring species (env + demographic + dist +interaction filter)
        if((ProbComm!=0)&&(ProbComm>=ProbAlone)) Abundances_in_Community_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)]<-1 #here the species is facilitated (if the species lived alone and in even more in community, its survival rate will be always 100% in community)
        if((ProbComm!=0)&&(ProbAlone>ProbComm)) Abundances_in_Community_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)]<-ProbComm #here the species is outcompeted (the maximum survival prob is then given by the community prob)
      }
      #the same for annuals
      if(Abundances_in_Community_Traits_Merged3$Annual[1]==1){  #annual
        ProbAlone<-nrow(Abundances_Alone_Traits_Merged3[which(Abundances_Alone_Traits_Merged3$SeedAbundance>0),])/ReplicateNumber_perUniqueEnvironment
        ProbComm<-nrow(Abundances_in_Community_Traits_Merged3[which(Abundances_in_Community_Traits_Merged3$SeedAbundance>0),])/ReplicateNumber_perUniqueEnvironment
        Abundances_in_Community_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)-2]<-ProbAlone
        Abundances_in_Community_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)-1]<-ProbComm
        if((ProbComm!=0)&&(ProbComm>=ProbAlone)) Abundances_in_Community_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)]<-1
        if((ProbComm!=0)&&(ProbAlone>ProbComm)) Abundances_in_Community_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)]<-ProbComm
      }
      #building a final maxtrix
      #TO-DO: This could be an array (to be simmilar as a RasterStack? i.e. one prob matrix per sp?)
      if((i==1)&&(y==1)){
        Abundances_Traits_Merged0<-Abundances_in_Community_Traits_Merged3
      }
      else{
        Abundances_Traits_Merged0<-rbind(Abundances_Traits_Merged0,Abundances_in_Community_Traits_Merged3)
      }
    }
  }
  #if the species did not survive in the community
  else{
    for (y in 1:MaxNumber_UniqueEnvironments){
      Abundances_Alone_Traits_Merged3<-Abundances_Alone_Traits_Merged2[Abundances_Alone_Traits_Merged2$Y==y,]
      #for perennials:
      if(Abundances_Alone_Traits_Merged3$Annual[1]==0) ProbAlone<-nrow(Abundances_Alone_Traits_Merged3[which(Abundances_Alone_Traits_Merged3$AdultAbundance>0),])/ReplicateNumber_perUniqueEnvironment
      #for annuals:
      if(Abundances_Alone_Traits_Merged3$Annual[1]>0) ProbAlone<-nrow(Abundances_Alone_Traits_Merged3[which(Abundances_Alone_Traits_Merged3$SeedAbundance>0),])/ReplicateNumber_perUniqueEnvironment
      
      Abundances_Alone_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)-2]<-ProbAlone
      Abundances_Alone_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)-1]<-0
      Abundances_Alone_Traits_Merged3[,ncol(Abundances_in_Community_Traits_Merged3)]<-0
      
      if((i==1)&&(y==1)){
        Abundances_Traits_Merged0<-Abundances_Alone_Traits_Merged3
      }
      else{
        Abundances_Traits_Merged0<-rbind(Abundances_Traits_Merged0,Abundances_Alone_Traits_Merged3)
      }
    } 
  }                  
}
summary(Abundances_Traits_Merged0)

length(levels(as.factor(as.character(Abundances_Traits_Merged0[which(Abundances_Traits_Merged0[,ncol(Abundances_in_Community_Traits_Merged3)-2]>0),1]))))
length(levels(as.factor(as.character(Abundances_Traits_Merged0[which(Abundances_Traits_Merged0[,ncol(Abundances_in_Community_Traits_Merged3)-1]>0),1]))))
length(levels(as.factor(as.character(Abundances_Traits_Merged0[which(Abundances_Traits_Merged0[,ncol(Abundances_in_Community_Traits_Merged3)]>0),1]))))
length(levels(as.factor(as.character(Abundances_Traits_Merged0[,1]))))
#levels(as.factor(as.character(AlphaInCommunity[,7])))
#AlphaInCommunity[which(AlphaInCommunity[,7]==1),]
#AlphaAlone0[which(AlphaAlone0[,7]==1),]


#comparing with co-occurr package:
install.packages("cooccur")
library("cooccur")
 
#making a spp_site matrix:
MaxNumber_UniqueEnvironments<-15
ReplicateNumber_perUniqueEnvironment<-15

spp_site_Matrix<-matrix(0,nrow=GammaRichness, ncol=MaxNumber_UniqueEnvironments*ReplicateNumber_perUniqueEnvironment)
spp_site_Matrix_AbsPres<-matrix(0,nrow=GammaRichness, ncol=MaxNumber_UniqueEnvironments*ReplicateNumber_perUniqueEnvironment)
CountingcellID<-1
for (i in 1:GammaRichness){ #GammaRichness is the total number of species considered, the total pool
  Abundances_Traits_Merged2<-Abundances_Traits_Merged0[Abundances_Traits_Merged0$SpeciesID==i,]
  CountingcellID<-1
  for(y in 1:MaxNumber_UniqueEnvironments) {
    Abundances_Traits_Merged3<-Abundances_Traits_Merged2[Abundances_Traits_Merged2$Y==y,]
    for(x in 1:ReplicateNumber_perUniqueEnvironment) {
      Abundances_Traits_Merged4<-Abundances_Traits_Merged3[Abundances_Traits_Merged3$X==x,]

      if(Abundances_Traits_Merged2$Annual[1]==0)  {
          spp_site_Matrix[i,CountingcellID]<-Abundances_Traits_Merged4$AdultAbundance
          if(Abundances_Traits_Merged4$AdultAbundance>0) spp_site_Matrix_AbsPres[i,CountingcellID]<-1 
      }
      
      if(Abundances_Traits_Merged2$Annual[1]==1) {
        spp_site_Matrix[i,CountingcellID]<-Abundances_Traits_Merged4$SeedAbundance
        if(Abundances_Traits_Merged4$SeedAbundance>0) spp_site_Matrix_AbsPres[i,CountingcellID]<-1 
        
      } 
      
      CountingcellID<-CountingcellID+1
      
    }
  }
}  
install.packages("raster")
library(raster)

specslist <- lapply(1:nrow(spp_site_Matrix), function(x) raster(nrows=15, ncols=15,vals=spp_site_Matrix[x,],
                                                                xmn=1, xmx=15, ymn=1, ymx=15))
Simu_specslist.stack <- stack(specslist)
plot(Simu_specslist.stack[[5]])

save(Simu_specslist.stack, file = "Simu_specslist.stack.RData")


#spp_site_Matrix
cooccur.probs<-cooccur(spp_site_Matrix_AbsPres, type = "spp_site", thresh = TRUE, spp_names = FALSE)
cooccur.probs<-cooccur(spp_site_Matrix_AbsPres, type = "spp_site", thresh = TRUE, spp_names = FALSE, eff_matrix = TRUE)
cooccur.probs<-cooccur(spp_site_Matrix_AbsPres, type = "spp_site", thresh = FALSE,spp_names = FALSE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)

summary(cooccur.probs)

cooccur.probs

head(cooccur.probs)

#head(prob.table(cooccur.probs))
plot(cooccur.probs)

save(cooccur.probs, file = "Simu_InteractionMatrix.RData")




#abundance-based cooccur is very computational demanding
#cooccur.probs<-cooccur(spp_site_Matrix, type = "spp_site", thresh = TRUE, spp_names = FALSE)
 















# ploting means:
install.packages("poibin") 
library(poibin)                                                     

AlphaPotentialRichness<-matrix(0,nrow=15*15, ncol=3)
AlphaPotentialNeedleRichness<-matrix(0,nrow=15*15, ncol=3) 
AlphaPotentialRichnessPoiBin<-matrix(0,nrow=15*15, ncol=3)
AlphaPotentialNeedleRichnessPoiBin<-matrix(0,nrow=15*15, ncol=3) 


SpeciesNeedle<-Abundances_Traits_Merged0[Abundances_Traits_Merged0$Y==1,] #this is the geographic needle of the paper (25 C=altitude 1)
SpeciesNeedle2<-SpeciesNeedle[SpeciesNeedle$DemoInteractionProbTemp>0,] #demographic needle
speciesinneedle<-levels(as.factor(as.character(SpeciesNeedle2$SpeciesID)))
for(z in 1:length(speciesinneedle)){
  if(z==1)SpeciesNeedle3<-Abundances_Traits_Merged0[Abundances_Traits_Merged0$SpeciesID==as.numeric(as.character(speciesinneedle))[z],] 
  else SpeciesNeedle3<-rbind(SpeciesNeedle3,Abundances_Traits_Merged0[Abundances_Traits_Merged0$SpeciesID==as.numeric(as.character(speciesinneedle))[z],] )
} 
countingRow<-1
for (g in 1:15){
  for (f in 1:15){  
    Abundances_Traits_Merged2<-Abundances_Traits_Merged0[Abundances_Traits_Merged0$X==g,]          
    Abundances_Traits_Merged3<-Abundances_Traits_Merged2[Abundances_Traits_Merged2$Y==f,]             
    Abundances_Traits_Merged4<-Abundances_Traits_Merged3[Abundances_Traits_Merged3$DemoInteractionProbTemp>0,]
    
    AlphaPotentialRichness[countingRow,1]<-g
    AlphaPotentialRichness[countingRow,2]<-f
    AlphaPotentialRichness[countingRow,3]<-length(levels(as.factor(Abundances_Traits_Merged4$SpeciesID)))
    #####              
    SpeciesNeedle4<-SpeciesNeedle3[SpeciesNeedle3$X==g,]          
    SpeciesNeedle5<-SpeciesNeedle4[SpeciesNeedle4$Y==f,]             
    SpeciesNeedle6<-SpeciesNeedle5[SpeciesNeedle5[,ncol(Abundances_in_Community_Traits_Merged3)-1]>0,] #env + disp+ demographic +interaction needle, where survival in community is !=0
    
 
    AlphaPotentialNeedleRichness[countingRow,1]<-g
    AlphaPotentialNeedleRichness[countingRow,2]<-f
    AlphaPotentialNeedleRichness[countingRow,3]<-length(levels(as.factor(SpeciesNeedle6$SpeciesID)))
    ############
    
    AlphaPotentialRichnessPoiBin[countingRow,1]<-g
    AlphaPotentialRichnessPoiBin[countingRow,2]<-f
    nOfSpp<-0
    for(nsp in 1: nrow(Abundances_Traits_Merged3)){                
      nOfSpp<-nOfSpp+rpoibin(1, Abundances_Traits_Merged3[nsp,ncol(Abundances_in_Community_Traits_Merged3)-1],wts=NULL) #To-do: give the column name 8 = $DemoInteractionProbTemp
    }
    AlphaPotentialRichnessPoiBin[countingRow,3]<-nOfSpp
    #####              
   
    AlphaPotentialNeedleRichnessPoiBin[countingRow,1]<-g
    AlphaPotentialNeedleRichnessPoiBin[countingRow,2]<-f
    nOfSpp<-0
    for(nsp in 1: nrow(SpeciesNeedle5)){                
      nOfSpp<-nOfSpp+rpoibin(1, SpeciesNeedle5[nsp, ncol(Abundances_in_Community_Traits_Merged3)-1],wts=NULL) #To-do: give the column name 8 = $DemoInteractionProbTemp
    }
    AlphaPotentialNeedleRichnessPoiBin[countingRow,3]<-nOfSpp
    
    ##
    countingRow<-countingRow+1
    
  }
}

AlphaPotentialRichness<-as.data.frame(AlphaPotentialRichness)
colnames(AlphaPotentialRichness)<-c("X", "Y","Richness")
AlphaPotentialNeedleRichness<-as.data.frame(AlphaPotentialNeedleRichness)
colnames(AlphaPotentialNeedleRichness)<-c("X", "Y","Richness") 
AlphaPotentialRichnessPoiBin<-as.data.frame(AlphaPotentialRichnessPoiBin)
colnames(AlphaPotentialRichnessPoiBin)<-c("X", "Y","Richness") 
AlphaPotentialNeedleRichnessPoiBin<-as.data.frame(AlphaPotentialNeedleRichnessPoiBin)
colnames(AlphaPotentialNeedleRichnessPoiBin)<-c("X", "Y","Richness") 
summary(AlphaPotentialNeedleRichness)
summary(AlphaPotentialRichness)
summary(AlphaPotentialRichnessPoiBin)
summary(AlphaPotentialNeedleRichnessPoiBin)


#ploting means
install.packages("fields")

library(fields)

AlphamatrixPotential<-matrix(0,nrow=15, ncol=15)
AlphamatrixPotentialNeedle<-matrix(0,nrow=15, ncol=15)
AlphamatrixPotentialPoiBin<-matrix(0,nrow=15, ncol=15)
AlphamatrixPotentialNeedlePoiBin<-matrix(0,nrow=15, ncol=15)
for(j in 1:nrow(AlphamatrixPotential)) {
  for(i in 1:ncol(AlphamatrixPotential)) {
    AlphaPotentialRichness2<-AlphaPotentialRichness[AlphaPotentialRichness[,1]==i,]
    AlphaPotentialNeedleRichness2<-AlphaPotentialNeedleRichness[AlphaPotentialNeedleRichness[,1]==i,]
    AlphaPotentialRichnessPoiBin2<-AlphaPotentialRichnessPoiBin[AlphaPotentialRichnessPoiBin[,1]==i,]
    AlphaPotentialNeedleRichnessPoiBin2<-AlphaPotentialNeedleRichnessPoiBin[AlphaPotentialNeedleRichnessPoiBin[,1]==i,]
    
    AlphaPotentialRichness3<-AlphaPotentialRichness2[AlphaPotentialRichness2[,2]==j,]
    AlphaPotentialNeedleRichness3<-AlphaPotentialNeedleRichness2[AlphaPotentialNeedleRichness2[,2]==j,]
    AlphaPotentialRichnessPoiBin3<-AlphaPotentialRichnessPoiBin2[AlphaPotentialRichnessPoiBin2[,2]==j,]
    AlphaPotentialNeedleRichnessPoiBin3<-AlphaPotentialNeedleRichnessPoiBin2[AlphaPotentialNeedleRichnessPoiBin2[,2]==j,]
    
      AlphamatrixPotential[i,j]<-AlphaPotentialRichness3[,3]
      AlphamatrixPotentialNeedle[i,j]<-AlphaPotentialNeedleRichness3[,3]
      AlphamatrixPotentialPoiBin[i,j]<-AlphaPotentialRichnessPoiBin3[,3]
      AlphamatrixPotentialNeedlePoiBin[i,j]<- AlphaPotentialNeedleRichnessPoiBin3[,3]
  }
}


# standardize over filters:    
install.packages("colorRamps")
library(colorRamps)
Maxvalue<- 210
Minvalue<- 10
deltavalue<-Maxvalue-Minvalue
par(las=1,cex=1.5)

image.plot(1:15,1:15,AlphamatrixPotential,col=matlab.like(40) ,  breaks=c(Minvalue,
                                                                          Minvalue+1:39*(deltavalue/40),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")
image.plot(1:15,1:15,AlphamatrixPotentialPoiBin,col=matlab.like(40) ,  breaks=c(Minvalue,
                                                                                Minvalue+1:39*(deltavalue/40),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")


image.plot(1:15,1:15,AlphamatrixPotentialNeedlePoiBin,col=matlab.like(40) ,  breaks=c(Minvalue,
                                                                                Minvalue+1:39*(deltavalue/40),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")

#Old:
Maxvalue<-max(AlphamatrixPotentialPoiBin)
deltavalue<-Maxvalue-min(AlphamatrixPotential[AlphamatrixPotential!=0])
par(las=1,cex=1.5)

image.plot(1:15,1:15,AlphamatrixPotential, col=two.colors(n=20, start="grey5", end="grey95", middle="grey50"), breaks=c(min(AlphamatrixPotential[AlphamatrixPotential!=0]),
                                                                                                                        min(AlphamatrixPotential[AlphamatrixPotential!=0])+1:19*(deltavalue/20),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")

Maxvalue<-max(AlphamatrixPotentialPoiBin)
deltavalue<-Maxvalue-min(AlphamatrixPotentialPoiBin[AlphamatrixPotentialPoiBin!=0])
par(las=1,cex=1.5)
image.plot(1:15,1:15,AlphamatrixPotentialPoiBin, col=two.colors(n=20, start="grey5", end="grey95", middle="grey50"), breaks=c(min(AlphamatrixPotentialPoiBin[AlphamatrixPotentialPoiBin!=0])-1,
                                                                                                                              min(AlphamatrixPotentialPoiBin[AlphamatrixPotentialPoiBin!=0])+1:19*(deltavalue/20),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")


Maxvalue<-max(AlphamatrixPotentialNeedlePoiBin)
deltavalue<-Maxvalue-min(AlphamatrixPotentialNeedle[AlphamatrixPotentialNeedle!=0])
par(las=1,cex=1.5)
image.plot(1:15,1:15,AlphamatrixPotentialNeedle, col=two.colors(n=20, start="grey5", end="grey95", middle="grey50"), breaks=c(min(AlphamatrixPotentialNeedle[AlphamatrixPotentialNeedle!=0])-1,
                                                                                                                              min(AlphamatrixPotentialNeedle[AlphamatrixPotentialNeedle!=0])+1:19*(deltavalue/20),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")
image.plot(1:15,1:15,AlphamatrixPotentialNeedle, col=heat.colors(n=20), breaks=c(min(AlphamatrixPotentialNeedle[AlphamatrixPotentialNeedle!=0]),
                                                                                 min(AlphamatrixPotentialNeedle[AlphamatrixPotentialNeedle!=0])+1:19*(deltavalue/20),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")


Maxvalue<-max(AlphamatrixPotentialNeedlePoiBin)
deltavalue<-Maxvalue-min(AlphamatrixPotentialNeedlePoiBin[AlphamatrixPotentialNeedlePoiBin!=0])
par(las=1,cex=1.5)
image.plot(1:15,1:15,AlphamatrixPotentialNeedlePoiBin, col=two.colors(n=20, start="grey5", end="grey95", middle="grey50"), breaks=c(min(AlphamatrixPotentialNeedlePoiBin[AlphamatrixPotentialNeedlePoiBin!=0])-1,
                                                                                                                                    min(AlphamatrixPotentialNeedlePoiBin[AlphamatrixPotentialNeedlePoiBin!=0])+1:19*(deltavalue/20),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")
image.plot(1:15,1:15,AlphamatrixPotentialNeedlePoiBin, col=heat.colors(n=20), breaks=c(min(AlphamatrixPotentialNeedle[AlphamatrixPotentialNeedle!=0])-1,
                                                                                       min(AlphamatrixPotentialNeedle[AlphamatrixPotentialNeedle!=0])+1:19*(deltavalue/20),Maxvalue+1), xaxt="n", yaxt="n", xlab="Mean Potential Total Richness")



