#SA2 - Spatial Ecology Script

#Student no.: 100966533

#Project title: Using citizen science species observation data to evaluate the potential for Castor fiber (Eurasian Beaver) rewilding in the Dark Peak SSSI, Peak District

#Aims:
# To use citizen science species observation data (iNaturalist) to assess Castor fiber (Eurasian Beaver) resource use in central Europe
# To characterise mean, minimum and maximum foraging distance from the river network for the species
# To generalise foraging distance information to assess potential resource use by Castor fiber within the Dark Peak SSSI

#Data used:
# iNaturalist Research Grade Castor fiber (Eurasian Beaver) species presence data from Central Europe
# Corine Land Cover 2018 
# European river network data 
# Dark Peak River network data
# Dark Peak SSSI shapefile 

#---------------------------------

#Sections:
#0 - Workspace Set-up
#1 - Data Processing
#2 - Beaver frequency from river network
#3 - Beaver Resource Selection Frequency
#4 - Beaver Resource Selection Functions
#5 - Calculating conductivity in Europe and Dark Peak
#6 - Estimating foraging distance from the Dark Peak river network 

#------------------------------------------------------------------------
#0 - Workspace Set-up
setwd("E:/0.Masters/1.5_GIS_Workspace/4.SE/SA2")

#Load packages
library(raster) #To manage spatial raster data
library(sp) #To manage spatial vector data
library(tidyverse) #To manage data columns
library(reshape2) #For data formatting into matrices
library(rgdal) #to produce coordinate reference systems
library(RColorBrewer) #To generate a landcover legend
library(ggplot2) #To visualise the data
library(maptools) #To provide the global vector basemap
library(adehabitatHS) #To compute resource selection
library(rgeos) #For clipping data layers
library(gdistance) #To calculate distance layer from Dark Peak river
library(sf) #To generate a raster bounding box 
library(fitdistrplus) #To generate kernel distributions of estimated beaver distance from river network

#Load the data in:
Beaver <- read.csv("CastorEU.csv") #iNaturalist Castor fiber species observation data
LandcoverLabels <- read.csv("clc_legend.csv") #Read in the Corine 2018 Landcover Labels

#Spatial data - Europe
EULandcoverRas <- raster("CLC100md.tif") #Read in the EU landcover layer (raster)
EURiversRas <- raster("rivCastor.tif") #Read in European river layer (raster)
EURiversVec <- readOGR("E:/0.Masters/1.5_GIS_Workspace/4.SE/SA2",layer="EURivers") #Read in the European river layer (vector)

#Spatial data - Dark Peak
DarkPeak<-readOGR("E:/0.Masters/1.5_GIS_Workspace/4.SE/SA2/DarkPeakSSSI",layer="DarkPeakMCP") #Read in the Dark Peak SSSI basemap (vector)
PeakRiversVec<-readOGR("E:/0.Masters/1.5_GIS_Workspace/4.SE/SA2/PeakRivers",layer="rivers_Peak") #Read in the Dark Peak SSSI rivers layer (vector)

#-----------------------------------------------------------------------------------------------------------------------------------------
#1 - Data Processing

#A: Beaver Data Processing
#Clean the Beaver data
#1 - Remove unwanted columns
Beaver <- Beaver[,-c(2,4:5,7:9,10:21,25:31,34:35)]

#2 - Relocate columns for ease of use 
Beaver <- Beaver %>% relocate(common_name, .before = id)
Beaver <- Beaver %>% relocate(scientific_name, .before = id)
Beaver <- Beaver %>% relocate(user_id, .after = positional_accuracy)

#3 - Positional accuracy
Beaver <- subset(Beaver,!is.na(Beaver$positional_accuracy)) #Remove observations with NA positional accuracy
Beaver <- Beaver[Beaver$positional_accuracy  <=  100,] #Remove beaver observations of >100m (0.1km) positional accuracy

#4 - Beaver data exploration
length(unique(Beaver[,"id"])) #1152 unique Beaver IDs
length(unique(Beaver[,"user_id"])) #553 unique users contributed to the dataset

#5 - Producing a shapefile from Beaver data
BeaverCoords <- data.frame(x = Beaver$longitude, y = Beaver$latitude) #Create an object containing the Beaver coordinate data
BeaverAttributes <- data.frame(Scientific_Name = Beaver$scientific_name, ID = Beaver$id, Observation_Date = Beaver$observed_on) #Create an object containing the relevant Beaver data
CRS <- crs(EULandcoverRas) #Set the Central European coordinates system
CRS <-  CRS("+proj=longlat +datum=WGS84") 
BeaverPoints <- SpatialPointsDataFrame(BeaverCoords, data = BeaverAttributes, proj4string = CRS) #Generate the SpatialPointsDataFrame

#B: Spatial processing - Europe
#Set the Beaver points CRS to be the same as CORINE landcover raster
crs(EULandcoverRas)
crs(BeaverPoints)
BeaverPoints <- spTransform(BeaverPoints, crs(EULandcoverRas))

#Manage the European Rivers data
#Vector - Check coordinates system and transform
crs(EURiversVec) #EU river vector data
crs(EULandcoverRas) #EU landcover
EURiversVec <- spTransform(EURiversVec,crs(EULandcoverRas)) #Make the rivers vector layer CRS The same as the landcover dataset

#Raster - Check coordinates system and transform
crs(EURiversRas) #Check the coordinates system
crs(EULandcoverRas) #EU landcover
EURiversRas <- spTransform(EURiversRas,crs(EULandcoverRas)) #Make the basemap's CRS The same as the landcover dataset

#Crop the landcover map to the rivers raster
st_bbox(EURiversRas) #find the boundaries of the rivers raster
x.min <-  4283321  -100 #Add 100m to support buffer processing
y.min <- 2744061 + 100
x.max <-  4825021 - 100
y.max <- 2888461 + 100
extent.new <- extent(x.min, x.max, y.min, y.max) #Create bounding box

EURiverCrop <- crop(EULandcoverRas,extent.new) #Generate the Central Europe Landcover data by cropping the CORINE landcover layer by the bounding box
EURiversVec <- crop(EURiversVec,EURiverCrop) #Crop rivers to the Central Europe landcover layer
EUBeaverRivers <- crop(BeaverPoints,EURiverCrop) #Crop Beavers to the Central Europe landcover layer

#Create a landcover legend - Central Europe
LandcoverLabels <- read.csv("clc_legend.csv") #Read in the Corine 2018 Landcover Labels
LandcoverLabels[21,] <- c(21, "Agriculture") #Rename long landcover label
LandcoverLabels[16,] <- c(4, "Fruit trees & berry plantations") #Rename long landcover label
LandcoverLabels[4,] <- c(4, "Roads & Rail") #Rename long landcover label
EURiverCropFac<-as.factor(EURiverCrop) #Set the landcover raster data as factors
levelsEU <-levels(EURiverCropFac) #Create an object to hold the levels
levelsEU <-data.frame(levelsEU) #Convert to a data frame to access the levels as data
namesEU<-LandcoverLabels[LandcoverLabels$Code %in% levelsEU$ID,] #Create an object containing only the legend items (landcover) present 
coulEU <- colorRampPalette(brewer.pal(12, "Paired"))(25) #Create a new palette object with 12 combinations of the 8 colours in the "Paired" palette

#Visualise the study area
plot(EURiverCrop,legend = FALSE,axes=T,box=T,col = coulEU) #Plot the raster with the colour palette

par(xpd = TRUE)
legend(4830000,3065000,paste(namesEU[,2]), fill = coulEU, #Add the legend
       cex = 0.7,bty="n") 
points(EUBeaverRivers, col="black", pch=12) #Plot the beaver observation data

#Plot to check European data extent 
#1 - River vector
plot(EURiversVec)
points(EUBeaverRivers, pch=16, col="blue")

#2 - River raster
plot(EURiversRas)
points(EUBeaverRivers, pch=16, col="blue")

#3 - Landcover
plot(EURiverCrop) #Plot river landcover raster crop
points(EUBeaverRivers, pch=16, col="blue") #Plot cropped beaver points

#C: Spatial processing -  Dark Peak
#Use the Central Europe CRS for the processing in Dark Peak
crs(DarkPeak) #Check coordinate systems
crs(PeakRiversVec)
crs(EULandcoverRas)
DarkPeak <- spTransform(DarkPeak, crs(EULandcoverRas)) #Project the dark peak basemap to the CORINE landcover map CRS
DarkPeakVec <- spTransform(PeakRiversVec, crs(DarkPeak)) #Project the dark peak rivers layer to the CORINE landcover map CRS

#Generate landcover layer for the Dark Peak SSSI
PeakLandcover <- crop(EULandcoverRas,DarkPeak) #Crop the Landcover values to the Dark Peak basemap
DarkPeakVec <- crop(DarkPeakVec,DarkPeak) #Crop the Dark Peak rivers layer to the Dark Peak basemap 

#Create a landcover legend - Dark Peak
PeakLandcoverFac<-as.factor(PeakLandcover) #Set the landcover raster data as factors
levelsPeak <-levels(PeakLandcoverFac) #Create an object to hold the levels
levelsPeak <-data.frame(levelsPeak) #Convert to a data frame to access the levels as data
namesPeak<-LandcoverLabels[LandcoverLabels$Code %in% levelsPeak$ID,] #Create an object containing only the legend items (landcover) present 
coulPeak <- colorRampPalette(brewer.pal(8, "Paired"))(12) #Create a new palette object with 12 combinations of the 8 colours in the "Set3" palette

#Visualise the study area
plot(PeakLandcover,legend = FALSE,axes=T,box=T,col = coulPeak) #Plot the raster with the colour palette

par(xpd = TRUE)
legend(3559000,3453500,paste(namesPeak[,2]), fill = coulPeak, #Add the legend
       cex = 0.7,bty="n") 
lines(DarkPeak)  #Plot the Dark Peak basemap
lines(DarkPeakVec, col="grey50") #Plot the rivers

#-----------------------------------------------------------------------------
#2 - Beaver frequency from river network

#What frequency is there of Beaver observations with distance from the Central European river network?
#Step 1: Merge the river vector dataset to reduce processing time
length(EURiversVec) #46711 individual rivers
RiverMerge<- gLineMerge(EURiversVec) #Merge the lines dataset
length(RiverMerge) #Only one shape

#Step 2: Buffer around the rivers dataset at set widths:
Buff50 <- raster::buffer(RiverMerge, width=50) #50m
Buff100 <- raster::buffer(RiverMerge, width=100) #100m
Buff150 <- raster::buffer(RiverMerge, width=150) #150m
Buff200 <- raster::buffer(RiverMerge, width=200) #200m
Buff500 <- raster::buffer(RiverMerge, width=500) #500m
Buff750 <- raster::buffer(RiverMerge, width=750) #750m
Buff1000 <- raster::buffer(RiverMerge, width=1000) #1000m

#Step 3: Find the cumulative frequency of beavers with distance from the river network
Buf1 <- length(EUBeaverRivers[Buff50, ]) #There were 29/130 beavers within 50m of the river network
Buf2 <- length(EUBeaverRivers[Buff100, ]) #There were 51/130 beavers within 100m of the river network
Buf3 <- length(EUBeaverRivers[Buff150, ]) #There were 72/130 beavers within 150m of the river network
Buf4 <- length(EUBeaverRivers[Buff200, ]) #There were 79/130 beavers within 200m of the river network
Buf5 <- length(EUBeaverRivers[Buff500, ]) #There were 116/130 beavers within 500m of the river network
Buf6 <- length(EUBeaverRivers[Buff750, ]) #There were 125/130 beavers within 750m of the river network
Buf7 <- length(EUBeaverRivers[Buff1000, ]) #There were 128/130 beavers within 1000m of the river network

#Plot the cumulative frequency of Beavers from the river network
BeaverCumul <- c(Buf1,Buf2,Buf3,Buf4,Buf5,Buf6,Buf7)
Distance <- c(50, 100, 150, 200, 500, 750, 1000)
CumulBeavers <- data.frame(cbind(Distance, BeaverCumul))
BeavCumul<- plot(CumulBeavers$Distance, CumulBeavers$BeaverCumul, frame = FALSE, pch = 19,  #Plot the cumulative frequency of beavers from the river network
     xlab = "Distance from river (m)", ylab = "Cumulative number of beavers", title = "Beaver cumulative frequency with distance (m) from the river network")
plot(BeavCumul)

#Step 4: Find the number of beavers situated within each buffer zone
#Calculate number of beavers within distance of rivers
Beav50 <- length(EUBeaverRivers[Buff50,]) #There were 29/130 beavers between 0-50m of the river network
Beav100 <- length(EUBeaverRivers[Buff100-Buff50,]) #There were 22/130 beavers between 0-50m of the river network
Beav150 <- length(EUBeaverRivers[Buff150-Buff100,]) #There were 21/130 beavers between 100-150m of the river network
Beav200 <- length(EUBeaverRivers[Buff200-Buff150,]) #There were 7/130 beavers between 150-200m of the river network
Beav500 <- length(EUBeaverRivers[Buff500-Buff200,]) #There were 37/130 beavers between 200-500m of the river network
Beav750 <- length(EUBeaverRivers[Buff750-Buff500,]) #There were 9/130 beavers between 500-750m of the river network
Beav1000 <- length(EUBeaverRivers[Buff1000-Buff750,]) #There were 3/130 beavers between 750-1000m of the river network

#Create new dataframe
BeaverNum <- c(Beav50,Beav100,Beav150,Beav200,Beav500,Beav750,Beav1000)
Distance <- c(50, 100, 150, 200, 500, 750, 1000)
NumBeavers <- data.frame(cbind(Distance, BeaverNum))

#Plot the frequency of beavers from the river network per buffer zone
BeavFreq<- plot(NumBeavers$Distance, NumBeavers$BeaverNum, frame = FALSE, pch = 19, 
      xlab = "Distance from river (m)", ylab = "Number of beavers", title = "Beaver frequency with distance (m) from the river network")

 ggplot(NumBeavers, aes(x=Distance, y=BeaverNum))+
   geom_point()+
   geom_line()+
   ggtitle("Observation frequency with distance (m) from the river network") +
        xlab("Distance (m)")+
    ylab("Number of Observations")

#Within which buffer were beavers most likely to be found?
opt<-NumBeavers[which(NumBeavers$BeaverNum==max(NumBeavers$BeaverNum)),]
opt #Most beavers at 200-500m = 37

#Select only the beaver observations within 500m of the river network
Beavers500m <- crop(EUBeaverRivers,Buff500) #116 beavers

plot(Beavers500m)

#rm(Buff50, Buff100, Buff150, Buff200, Buff750, Buff1000, Buf1,Buf2,Buf3,Buf4,Buf5,Buf6,Buf7,Beav50,Beav100,Beav150,Beav200,Beav500,BEav750,Beav1000,BeavNum, NumBeavers,BeaverCumul, CumulBeavers) #Tidy workspace
#-----------------------------------
#3 - Beaver Resource Selection Frequency
#A: What landcover types did beavers select most frequently in Central Europe at 500m from the river network?
set.seed(11) #Set random number generator
back.xy <- dismo::randomPoints(EURiverCrop, p=EUBeaverRivers, n=1000,ext = extent(EUBeaverRivers)) #Generate 1000 pseudoabsence points from the beaver points extent and landcover values
back.xy<-SpatialPoints(back.xy) #Convert into a shapefile
length(back.xy) #1000 pseudoabsence samples

#Crop the pseudoabsence points with the 500m buffer to ensure pseudoabsence points extract only landcover values within 500m of the river network
back.xy <- crop(back.xy,Buff500)
length(back.xy) #654 pseudoabsence samples
back.xy <- data.frame(back.xy)

#Visualise the presence and pseudoabsence points on the map
plot(EURiverCrop)
plot(Beavers500m,add=T, pch=16)
points(back.xy,add=TRUE, col='red', pch=16) #Plot presence and absence points

#Convert background points to the European Landcover CRS and extract landcover values for  presence and absence points
A<-SpatialPoints(back.xy,crs(EURiverCrop)) #Set to SpatialPointsDataFrame
P<-EUBeaverRivers
eA<-raster::extract(EURiverCrop,A) #Extract landcover values at absence points
eP<-raster::extract(EURiverCrop,P) #Extract landcover values at presence points
table(eA) 
table(eP) 

#Plot the proportion of points (presence and pseudoabsence) occurring in the region
par(mfrow=c(1,2))
hist(eA,freq=FALSE,breaks=c(0:40),xlim=c(0,40),ylim=c(0,1))
hist(eP,freq=FALSE, breaks=c(0:40),xlim=c(0,40),ylim=c(0,1))

dev.off()

#Visualise the frequency of beavers present per landcover type 
eP<-table(eP)
eP<- as.data.frame(eP)
colnames(eP)<-c("Code","Freq")
eP<- inner_join(eP, LandcoverLabels , by = "Code")
BeavPres<- ggplot(eP, aes(x=Freq, y=LandCover)) +
  geom_bar(stat="identity") #Plot beaver presence frequency
BeavPres

#Visualise the available landcover types within 500m of the river network
eA<-table(eA)
eA<- as.data.frame(eA)
colnames(eA)<-c("LandCover","Freq")
eA<- inner_join(eA, LandcoverLabels , by = "Code") #Remove irrelevant landcover types
eA<- eA[!(eA$LandCover=="Vineyards"),]
eA<- eA[!(eA$LandCover=="Moorland and Marsh"),]
eA <- eA[!(eA$LandCover=="Fruit trees and berry plantations"),]
BeavAbs <- ggplot(eA, aes(x=Freq, y=LandCover)) +
  geom_bar(stat="identity") #Plot landcover pseudoabsence frequency
BeavAbs

#B: What landcover types were available within 500m of the river network in Dark Peak?
back.xy <- dismo::randomPoints(PeakLandcover, n=1000,ext = extent(PeakLandcover)) #Generate 1000 pseudoabsence points from the beaver points extent and landcover values
back.xy<-SpatialPoints(back.xy) #Convert into a shapefile
length(back.xy) #1000 pseudoabsence samples

#Crop the background points with merged and buffered Dark Peak river data
PeakRivMerge <- gLineMerge(DarkPeakVec) #Merge Dark Peak Rivers
PeakBuff <- raster::buffer(PeakRivMerge, width=500) #Create a 500m buffer around the rivers
back.xy <-  crop(back.xy,PeakBuff) #Crop the background points with the buffered 500m
length(back.xy) #499 points

#Visualise the output
plot(PeakLandcover)  #Dark Peak Landcover
lines(DarkPeakVec) #River vector
lines(DarkPeak) #Basemap
points(back.xy, pch=16, col="red") #Background points

#Convert background points to PeakLandcover CRS and extract landcover values
pA<-SpatialPoints(back.xy,crs(PeakLandcover))
pA<-raster::extract(PeakLandcover,pA) 
table(pA) 

#Convert to dataframe, add landcover labels and visualise the available landcover types within 500m of the
pA<-table(pA)
peA<- as.data.frame(pA)
colnames(pA)<-c("Code","Freq")
pA<- inner_join(pA, LandcoverLabels , by = "Code")
PeakLand<- ggplot(pA, aes(x=Freq, y=LandCover)) +
  geom_bar(stat="identity") #Plot beaver presence frequency
PeakLand

#Create a bar plot including all landcover frequencies (Presence vs Central Europe availability vs Dark Peak availability)
eA[,3] <- "Europe"
colnames(eA) <-c("Code","Freq", "Availability") #Central Europe presence
eP[,3] <- "Beavers" #Central Europe absence
colnames(eP) <-c("Code","Freq", "Availability") #Central Europe presence
pA[3] <- "Peak" #Dark Peak
colnames(pA) <-c("Code","Freq", "Availability") #Central Europe presence
All_Landcover<- full_join(x=eA,y=eP, by = c("Code", "Availability", "Freq"))
All_Landcover<- full_join(x=All_Landcover,y=pA, by = c("Code", "Availability", "Freq"))
All_Landcover<- inner_join(All_Landcover, LandcoverLabels , by = "Code") #Add landcover labels to the dataset

#Create the all landcover bar plot
ggLandcover<- ggplot(All_Landcover, aes(fill=Availability, y=LandCover, x=Freq)) +
  geom_bar(position='dodge', stat='identity', color="Black")  +
  theme(panel.grid.major=element_line(colour="grey75"),panel.grid.minor=element_line(colour="grey75"))
   
ggLandcover
#----------------------------------------------------------------------------
#4 -  Beaver Resource Selection Functions

#What foraging behaviours did the Central European beavers display in the European study area within 500m of the river network?

#A: Resource Selection Functions (foraging behaviour) from the European Beaver presence data
#Take a systematic sample from the raster
set.seed(11);back.xy <- sampleRegular(EURiverCrop,size=50000,useGDAL=TRUE, sp = TRUE) #50000 samples selected from the landcover raster and set to produce a spatialpoints layer

#Select only the background points within 500m of the river network
back.xy <- crop(back.xy,Buff500) #Select only the landcover values within 500m of the river network
back.xy <- data.frame(back.xy) #Set to dataframe
length(back.xy) #31746 samples
head(back.xy) #Unwanted columns
back.xy <- as.matrix(back.xy[,-c(2:4)]) #Remove unecessary columns and return to matrix format
 
#Extract raster values from the presence data (within 500m of river network) by their locations
pres.cov <- raster::extract(EURiverCrop, Beavers500m) 
pres.cov<-data.frame(pres.cov,pres=1) #Convert the extracted values into a dataframe
colnames(pres.cov)<-c("Cover","Pres") #Step 4: Rename the columns 
head(pres.cov) #inspect the raster values within the dataset

#Calculate the 'length' of the vector for each habitat (i.e. the count)
avail <- tapply(back.xy, back.xy, length) 
avail #inspect the data
length(avail)
class(avail) #data type = array 
avail<-as.numeric(avail) #Convert the array data into numeric data

#Create columns in pres.cov for analysis
pres.cov$ID<-as.numeric(Beavers500m$Scientific_Name) #Design I analysis all animals into one group
USEAll<-dcast(pres.cov, ID ~ Cover, length, value.var = "ID") #Format the data
namesAll<-as.data.frame(LandcoverLabels[LandcoverLabels$Code %in% pres.cov$Cover,]) #Set the landcover names
colnames(USEAll)<-c("ID",as.character(namesAll$LandCover))

#Clean the data  
length(avail) #inspect the length of avail
avail
length(USEAll)
avail<-avail[-c(4:8, 12:13,18,20:23)] #remove irrelevant habitat landcovers 
avail

#Reformat land use to matrix
USEAll<-as.matrix(USEAll) 

#Analyses I (Europe): 
#Compute Selection Ratios for Habitat Selection Studies - Type I
sel.ratioI <- widesI(u = USEAll[,c(2:ncol(USEAll))],a = avail, avknown = T, alpha = 0.05)

summary(sel.ratioI) #Type I list of statistics, View the Type I Analyses statistics

par(mai=c(2,0.82,0.82,0.42)) #Set the plotting scale
sel.ratioI #Type I - matrix of Wi values and significance
plot(sel.ratioI) #Plot the Type I Analyses - graphical visualisation, Manly selectivity measure, Scaled selection ratios, Used and available proportions

selRatioEU <- c(sel.ratioI$wi) #Extract the selection ratios per habitat within 500m of the river network for the Beavers in Central Europe

#B - Predicted Resource Selection Functions (foraging behaviour) of Beavers in the Dark Peak SSSI, Peak District

#Aims:
#To select only land use (USEAll) values for the landcover types available in the Dark Peak SSSI 
#To estimate the selection ratios of Beavers from the landcover available in the Dark Peak SSSI to predict beaver foraging behaviour in the Dark Peak study area

#Select only the landcover types required
levels(as.factor(EURiverCrop)) #25 IDs, Landcover layers = 1, 2, 3, 4 ,5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 18, 20, 21, 22, 26, 29, 30, 31, 32, 35, 40
levels(as.factor(PeakLandcover)) #12 IDs, #Landcover layers = 2, 3, 7, 11, 12, 18, 21, 22, 26, 29, 35, 40
                  #Cannot use the following classes due to absence from pres and abs data: 7, 26, 35
                  #Resulting layers: 2, 3, 11, 12, 18, 21, 22, 29, 40

set.seed(11); back.xy <- sampleRegular(PeakLandcover,size=50000,useGDAL=TRUE, sp=TRUE) #50000 samples selected from the Dark Peak landcover raster

#Crop Dark Peak background points to within 500m of the river network
back.xy <- crop(back.xy,PeakBuff) #Select only the landcover values within 500m of the river network
length(back.xy) #24595 background points
back.xy <- data.frame(back.xy)
head(back.xy) #Identify unwanted columns
unique(back.xy$CLC100md) #12 landcover types, (2, 3,  11,  12, 18, 22,  26, 29, 35, 40)
back.xy <- as.matrix(back.xy[,-c(2:4)]) #Remove unecessary columns and return to matrix format

avail <- tapply(back.xy, back.xy, length) 
avail #inspect the data
length(avail)
USEAll
length(USEAll)
avail<-as.numeric(avail) 

#Remove landcover classes not present at Dark Peak from USEAll and avail: 26 and 35
avail<-avail[-c(7, 9)] #remove irrelevant habitat landcovers

USEAll<-data.frame(USEAll) #Convert the USEAll data back to dataframe
USEAll<-USEAll[-c(2, 5, 9,10)] #Remove unwanted landocver classes
length(USEAll) #Check length
length(avail)

USEAll<-as.matrix(USEAll) #Reformat European land use data back to matrix

#Analyses I (Dark Peak): 
sel.ratioI <- widesI(u = USEAll[,c(2:ncol(USEAll))],a = avail, avknown = T, alpha = 0.05)
plot(sel.ratioI) #Plot the Type I Analyses - graphical visualisation, Manly selectivity measure, Scaled selection ratios, Used and available proportions
selRatioPeak <- c(sel.ratioI$wi) #Extract the selection ratios per habitat for the Beavers in Central Europe

dev.off()

#---------------------------------------------
#5 - Calculating conductivity

#A: Calculate the cost layer from the Central European River Network
selRatioEU #European landcover selection ratios

#Find the relevant landcover categories
EURiverCrop<-as.factor(EURiverCrop) #raster categories as factors
levels(EURiverCrop) #inspect levels (25 landcover classes)

#Select only the 12 landcover types in selRatioEU, adding NAs to the rest
selRatioEU <- c(4.1912716, 3.1918298 ,1.6159119,NA,NA,NA,NA,NA,NA,99.3486590 ,10.8746505,0.1549334 ,NA,NA, 0.6696468 ,1.0859975,1.2072069 , 0.7478077 ,NA, 1.1366160 ,NA,NA,NA,NA, 9.6182168)
length(selRatioEU) #25 classes
#Calculate the cost layer
  vectorOfOnes <- rep(1,25) #produce a dataset of X repetitions of 1s
Cost<-vectorOfOnes/selRatioEU #Calculate cost by dividing 1 with the selection ratios (1/25)
RCmat<- cbind(levels(EURiverCrop)[[1]],Cost) #cbind (bind by column) the raster values and the Cost values
RCmat #Inspect the output

#Reclassify the EU raster by the cost layer produced
EURiverCropCost<-reclassify(EURiverCrop, RCmat)
raster::plot(EURiverCropCost) #Plot the new layer
points(EUBeaverRivers, pch=16) #EU Beaver Points
#Interpretation: Beavers in Central Europe were able to exploit a wide range of landcover types

#B - Foraging distance from the River Network - Peak District
#Aim: to generate a resistance layer representing habitat selection values from the Dark Peak river network

PeakLandcover<-as.factor(PeakLandcover) #raster categories as factors
levels(PeakLandcover) #inspect levels (list of categories)

#Select only those landcover values relevant to the Dark Peak SSSI and beaver selRatios
#Landcover layers = 2, 3, 7, 11, 12, 18, 21, 22, 26, 29, 35, 40 (12 landcover types)
selRatioPeak       #2, 3,    11, 12, 18,     22,     29,     40 (8 landcover types)
#Landcover types of PeakLandcover not in selRatio: 7 (minerals), 26 (natural grasslands), 35 (moorland, marsh)

#Input empty values and remove landcover types not in for these areas in selRatioPeak according to levels(PeakLandcover)
selRatioPeak <- c(5.6989070, 5.8854167, NA, 3.4073465, 2.7945144, 0.2957110, NA, 0.9317531 ,NA,  0.5138062 , NA ,1.8038893 )
length(selRatioPeak)

#Calculate the cost layer
vectorOfOnes <- rep(1,12) #produce a dataset of X repetitions of 1s
Cost<-vectorOfOnes/selRatioPeak #Calculate cost by dividing 1 with the selection ratios (1/18)
RCmat<- cbind(levels(PeakLandcover)[[1]],Cost) #cbind (bind by column) the raster values and the Cost values
RCmat #Inspect the output

#Reclassify the raster by the cost layer produced
PeakLandcoverCost<-reclassify(PeakLandcover, RCmat)
raster::plot(PeakLandcoverCost) #Plot the new layer
lines(DarkPeakVec) #Plot with the Dark Peak rivers
lines(DarkPeak)

#Conductance maps 

#Calculate the mean, maximum and minimum conductance
#Europe
mean_cond <- transition(1/EURiverCropCost, transitionFunction=mean, 8)
max_cond <- transition(1/EURiverCropCost, transitionFunction=max, 8)
min_cond <- transition(1/EURiverCropCost, transitionFunction=min, 8)
 
#Rasterise the conductance transition layers
mean_rasEU <- raster(mean_cond)
max_rasEU <- raster(max_cond)
min_rasEU <- raster(min_cond)

#Plot the conductance layers
plot(mean_rasEU)
plot(max_rasEU)
plot(min_rasEU)

#Calculate the mean, maximum and minimum conductance
#Peak
mean_cond <- transition(1/PeakLandcoverCost, transitionFunction=mean, 8)
max_cond <- transition(1/PeakLandcoverCost, transitionFunction=max, 8)
min_cond <- transition(1/PeakLandcoverCost, transitionFunction=min, 8)

#Rasterise the conductance transition layers
mean_rasPeak <- raster(mean_cond)
max_rasPeak <- raster(max_cond)
min_rasPeak <- raster(min_cond)

#Plot the conductance layers
plot(mean_rasPeak)
lines(DarkPeakVec, col="grey55")
lines(DarkPeak, col="grey40")

#Plot the conductance layers
plot(max_rasPeak)
lines(DarkPeakVec, col="grey55")
lines(DarkPeak, col="grey40")

#Plot the conductance layers
plot(min_rasPeak)
lines(DarkPeakVec, col="grey55")
lines(DarkPeak, col="grey40")

#---------------------------------------------

#6 - Estimating foraging distance from the Dark Peak river network 

#Calculate the predicted dispersal of Beavers from the river network using the cost layer
#Distance matrix
#A.obs = DistPeak = Peaklandcovercost
#nodes = RiverNodes = DarkPeakVec

DistPeak<- as.matrix(PeakLandcoverCost) #Convert the cost layer to a distance matrix
head(DistPeak)
RiverNodes<- data.frame(DarkPeakVec) #Convert the Dark Peak river data to data frame
head(RiverNodes)

#Set up the cost  matrix 
colnames(DistPeak) <- 1:217 #Rename columns and rows
rownames(DistPeak) <- 1:285
diag(DistPeak) <- 0  #Set the diagonal to 0
DistPeak[is.na(DistPeak)] <- 0 #Set NA values to 0

View(DistPeak)

#Calculate the distance matrix
coords <- cbind(RiverNodes$FROM_NODE, RiverNodes$TO_NODE) #Select the river node coordinates
length(coords) #612 coordinates
distmat <- pointDistance(coords, lonlat=F) #Calculate the geographical distance between river points
distmat <- distmat/100 #Convert distances into metres

View(distmat) #View the distance matrix

matLocations <- which(DistPeak>0, arr.ind=T) #Return array indices greater than 0

within_disp <- cbind(distmat[matLocations], DistPeak[matLocations]) #Bind the distance matrix with the cost matrix
head(within_disp) #Peek at the resistance matrix head

within_disp <- rep(within_disp[,1], within_disp[,2])

#Plot the distance matrix
hist(within_disp,col="grey50",border="white", xlab="Distance from the river network (m)", ylab="Frequency", main="") #plot

#Get summary stats for the dispersal distances
mean.dist <- mean(within_disp)  #Mean dispersal distance       
max.dist <- max(within_disp)  #Maximum dispersal distance
min.dist <- min(within_disp) #Minimum dispersal distance

mean.dist #158.1016m
max.dist #427.4081m
min.dist #0.01m

#Kernel distributions of predicted beaver distance from the river network
disp.lnorm <- fitdistrplus::fitdist(data=within_disp, distr="lnorm", method="mle") #Fit a log normal curve to the dispersal data
disp.exp <- fitdistrplus::fitdist(data=within_disp, distr="exp", method="mle") #Fit an exponential curve to the dispersal data
disp.AIC <- gofstat(list(disp.exp, disp.lnorm),
                    fitnames=c("exponential", "lognormal")) #Item to list the goodness-of-fit (AIC) values 

#Check the model goodness-of-fit statistics
disp.AIC$aic #exponential =597558.4,  lognormal = 581858.1 
summary(disp.exp) #Summary of the exponential curve
summary(disp.lnorm) #Summary of the log-normal curve, a better fit (lower AIC value)

#Plot the kernel distributions of predicted beaver distance from the river network
denscomp(disp.exp, addlegend=FALSE) #Plot the exponential distribution of dispersal
denscomp(disp.lnorm,addlegend = FALSE, xlab="Distance from the river network (m)") #Plot the log-normal distribution of dispersal

