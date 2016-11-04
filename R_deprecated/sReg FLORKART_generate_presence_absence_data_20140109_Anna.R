#############################################################
# Script for the preparation of the FLORKART database to be used as case study
# sDiv workshop sRegPool
# Anna Cord

setwd("C:/Daten/Konferenzen/2013/Workshop_iDiv/Germany_data")

require(raster)        
require(maps)
require(maptools)
require(rgdal)
require(lattice)

gis.directory <- paste(getwd(),"/GIS/", sep="")
florkart.directory <- paste(getwd(),"/FLORKART/", sep="")

#####################################################################
#################### IMPORTING AND DISPLAYING DATA ##################
#####################################################################

# ---------------- Import and display GIS layers (for maps only)-----------------
mtb.shp <- readShapeSpatial("GIS/mtb_brd_GK_Zone3.shp", proj4string=CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +datum=potsdam +units=m +no_defs +ellps=bessel +towgs84=606.0,23.0,413.0"))
plot(mtb.shp)            

#########################################################################
################### LOAD FLORKART DATABASE ##############################
#########################################################################

# Load overview of MTBs from FLORKART
florkart.atlmtb <-read.table(paste(florkart.directory,"ATLMTB.txt",sep=""), sep=",", head=T) 

# Load file showing transformation between TK_NR (used in MTB shape file) and MTBID (used in ATLMTB within FLORKART)
tknr.mtbid <- read.table(paste(florkart.directory,"TKNR_MTBID.csv",sep=""), sep=",", head=T) 

# Add column to MTB shape file which includes the MTBID
mtb.shp$MTBID <- NaN    #add new column to attribute table

for (i in 1:length(mtb.shp$MTBID)) {          # Note! TK 2407 (=MTBID 361) does not include terrestrial area and is therefore not included in mtb.shp
l <- mtb.shp$TK_NR[i]                         # This results in only 2,994 MTBs.
m <- tknr.mtbid$MTBID[tknr.mtbid$TK_NR == l]  # However, FLORKART (ATLMTB.txt) does record 1 species for this MTB.
mtb.shp$MTBID[i] <- m                         # Since all analyses here are based on mtb.shp, TK 2407 was excluded from the analysis
}  

# Create list with MTBs     
list.mtbid <- unique(mtb.shp$MTBID)      # List with all MTBs based on mtb.shp, MTBs for which no information from FLORKART is available are flagged with -999
list.mtbid.2994 <- subset(list.mtbid,list.mtbid>0)
length(list.mtbid.2994)

# Create list with all taxa (Species1,Species2, etc.)
list.species <- unique(florkart.atlmtb$TAXNR)     
total.number.species <- length(list.species)

#########################################################################
### CREATE MATRIX WITH MTBs AND SPECIES PRESENCE/ABSENCE ################
#########################################################################

mtb.species.matrix <- data.frame(matrix(NaN,ncol=total.number.species, nrow=length(list.mtbid.2994)))
dim(mtb.species.matrix)
rownames(mtb.species.matrix) <- paste("MTB",list.mtbid.2994,sep="")
colnames(mtb.species.matrix) <- list.species
head(mtb.species.matrix)
     
for (k in 1:length(list.mtbid.2994)) {
focal.mtbid <- list.mtbid.2994[k]
print(focal.mtbid)
subset.florkart.atlmtb <- florkart.atlmtb[florkart.atlmtb$MTBID==focal.mtbid,]
subset2.florkart.atlmtb <- subset.florkart.atlmtb[subset.florkart.atlmtb$SYM %in% c("O","Q","X","Z"),]      # Include only data after 1950 and with specific SYM code
temp <- subset2.florkart.atlmtb$TAXNR    # Extract all species that were reported per MTB
mtb.species.matrix[k,]  <- (list.species %in% temp) + 0
}

mtb.species.matrix <- cbind(list.mtbid.2994,mtb.species.matrix)
names(mtb.species.matrix)[1] <- "MTBID"
     
# Save species list per MTB
write.csv(mtb.species.matrix,"MTB_species_matrix.csv",row.names = FALSE)   

#########################################################################
################### SPECIES OF INTEREST #################################
#########################################################################

# Define taxa of interest
family <- "Ranunculaceae"
genus.list <- c("Aconitum","Actaea","Adonis","Anemone","Aquilegia","Caltha","Ceratocephala","Clematis","Consolida","Delphinium","Eranthis","Helleborus","Hepatica","Myosurus","Nigella","Pulsatilla","Ranunculus","Thalictrum","Trollius")

# Identify genera of interest to later compile the final species list
florkart.species <- read.csv(paste(florkart.directory,"SPECIES.csv",sep=""), sep=",", head=T) 
interest.species <- florkart.species[florkart.species$GATTUNG %in% genus.list,]
interest.species.number <- interest.species$NR_NAME

species.names.table <- data.frame(matrix(NaN,ncol=2, nrow=length(interest.species.number)))
names(species.names.table) <- c("Species_number","Species_name")

# Select only species of interest from MTB species matrix 
selection <- colnames(mtb.species.matrix) %in% interest.species.number
mtb.interest.species.matrix <- mtb.species.matrix[,selection]

write.csv(mtb.interest.species.matrix,"MTB_Species_of_interest_matrix.csv",row.names = FALSE)

# Extract species number and name of the species of interest
for (i in as.integer(names(mtb.interest.species.matrix)))
{
  florkart.species.subset <- florkart.species[florkart.species$NR_NAME == i,]
  j <- which(as.integer(names(mtb.interest.species.matrix)) == i)
  species.names.table[j,1] <- i
  species.names.table[j,2] <- paste(unique(florkart.species.subset$GATTUNG),"_",unique(florkart.species.subset$ART_EPIT),"_",unique(florkart.species.subset$RANG),"_",unique(florkart.species.subset$SUB_EPIT),sep="")
}

write.csv(species.names.table,"Species_of_interest_names.csv",row.names = FALSE) 

# Merge columns for species (subspecies etc.) based on expert knowledge (outside of R)
species.table.final <- read.csv("Species_of_interest_names_final.csv", head=TRUE)
names(species.table.final)
rclmat <- cbind(species.table.final$Species_number, species.table.final$New_species_number)

# Replace old species names with new names in mtb.interest.species.matrix -> only 119 species remain
index <- match(rclmat[,1],as.integer(names(mtb.interest.species.matrix)))
names(mtb.interest.species.matrix)[index] = rclmat[,2]

# Exclude all species which were not modelled on the European scale -> only 52 species remain
list.new.species.numbers <- unique(names(mtb.interest.species.matrix))
list.modelled.species <- unique(as.integer(species.table.final$New_species_number[which(species.table.final$Modelled == 1)]))
length(list.modelled.species)

# Create empty data frame
MTB.species <- data.frame(matrix(NaN,ncol=length(list.modelled.species), nrow=length(list.mtbid.2994)))
names(MTB.species) <- list.modelled.species

for (i in list.modelled.species)
{
  temp <- mtb.interest.species.matrix[,names(mtb.interest.species.matrix)==i]
  j <- which(list.modelled.species == i)
  if (length(dim(temp)) > 0)  { 
      presence <- apply(temp[,],1,max)
      MTB.species[,j] <- presence    
      } 
      else 
      {
      MTB.species[,j] <- temp 
      }
}

# Replace species numbers by names
#names(MTB.species) <- list.modelled.species
rclmat2 <- data.frame(cbind(species.table.final$New_species_number, levels(species.table.final$New_species_name)[species.table.final$New_species_name]))
rclmat2 <- rclmat2[1:69,]

# Remove duplicates
dupl <- which(duplicated(rclmat2))
rclmat2 <- rclmat2[-dupl, ]

index2 <- match(rclmat2[,1],as.integer(names(MTB.species)))
names(MTB.species)[index2] <- levels(rclmat2[,2])[rclmat2[,2]]

# Add TK numbers
tknr.mtbid <- read.table(paste(florkart.directory,"TKNR_MTBID.csv",sep=""), sep=",", head=T) 

MTB.species.final <- cbind(list.mtbid.2994,MTB.species)
names(MTB.species.final) <- c("MTBID",names(MTB.species))
head(MTB.species.final)

MTB.species.final <- merge(MTB.species.final,tknr.mtbid,by="MTBID",all.x=TRUE, sort=FALSE)

write.csv(MTB.species.final,"MTB_species_final.csv",row.names = FALSE)

# Merge results to mtb.shp
mtb.shp@data$orig.order <- c(1:nrow(mtb.shp@data))
mtb.shp@data  <- merge(mtb.shp@data,MTB.species.final,by="TK_NR",all.x=TRUE,sort=FALSE)
mtb.shp@data <- mtb.shp@data[order(mtb.shp@data$orig.order),]

# Plot distribution maps of all species
#mtb.shp@data <- merge(mtb.shp@data,MTB.species.final,by="MTBID",all.x=TRUE, sort=FALSE)  

for (i in colnames(mtb.shp@data)[8:59])
{
trellis.device(color=T, device="png", width = 700, height = 1000, file=paste(getwd(),"/Species_distributions/",i,".png",sep=""))
plot(mtb.shp,col=mtb.shp@data[,i])
dev.off()
}

# Replace full species names by abbreviations
MTB.species.final.abbr <- read.csv("MTB_species_final_abbreviations.csv", head=TRUE)
head(MTB.species.final.abbr)

# Write shape file with species information using abbreviations for species names
mtb.shp.abbr <- readShapeSpatial("GIS/mtb_brd_GK_Zone3.shp", proj4string=CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +datum=potsdam +units=m +no_defs +ellps=bessel +towgs84=606.0,23.0,413.0"))
mtb.shp.abbr@data$orig.order <- c(1:nrow(mtb.shp.abbr@data))
mtb.shp.abbr@data  <- merge(mtb.shp.abbr@data,MTB.species.final.abbr,by="TK_NR",all.x=TRUE,sort=FALSE)
mtb.shp.abbr@data <- mtb.shp.abbr@data[order(mtb.shp.abbr@data$orig.order),]

writeSpatialShape(mtb.shp.abbr, "MTB_species_final_abbr.shp")  
distribution.shp <- readShapeSpatial("MTB_species_final_abbr.shp",proj4string=CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +datum=potsdam +units=m +no_defs +ellps=bessel +towgs84=606.0,23.0,413.0"))
spplot(distribution.shp,"Adon_aest") # test plot

save.image("FLORKART_20140109.RData")
     
