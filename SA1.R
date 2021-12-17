#SA1 - Spatial Ecology Code
#Student no.: 10966533

#Sections:
#Section 1 - Scale
#Section 2 - Resource Selection & Home Range
#Section 3 - Species Distribution Models

#======================================================================================================================
#Section 1 - Scale

#Aim: To investigate the characteristic scales of Meles meles (Eurasian badger) to deciduous tree cover

#Tasks:
#Calculate the characteristic scale(s) of the response by the Eurasian badger (Meles meles) to:
#A - deciduous tree cover, B - buffer distances of between 100 and 2000m in 100m intervals, C - the same extent (study area) as the practical
#Calculate the optimum scale for modelling badger response to broadleaf tree cover using the log likelihood statistic for model selection
#Produce plots showing the relationship between buffer radii and model performance (loglikelihood statistic)
#Data processing: A - Remove any records that contain unconfirmed (under the Identification verification status field), B - Remove any points with coordinate uncertainty > 1000m. 

#Sections:
#0- Setting up the working space
#1 - Data Preparation
#2 - Data Processing 
#3 - Generating pseudo-absence points
#4 - Raster Reclassification
#5 - Buffer Analysis
#6 - Log-likelihood Statistic
#7 - Determining the optimum buffer size

#----------------------------------------------------------------------------------------------------------------------
#0 - Setting up the working space

setwd("E:/0.Masters/1.5_GIS_Workspace/4.SE/SA1")

#Load in packages
library(raster) #for working with raster data
library(sp) #for working with vector data
library(rgdal) #for working with coordinate systems
library(dismo) #for species distribution modelling
library(mapview) #to provide interactive visualisations of spatial data
library(rgeos) #for providing spatial geometry operations (e.g. buffers)

Melesmeles <- read.csv("Melesmeles.csv")

colnames(Melesmeles) #check data
head(Melesmeles)

#----------------------------------------------------------------------------------------------------------------------
#1 - Data Preparation
Melesmeles <- Melesmeles[, -c(1,4, 9)] #Remove irrelevant columns

#Processing observation data 
unique(Melesmeles$Identification.verification.status) #remove unconfirmed observations
Melesmeles <- Melesmeles[!Melesmeles$Identification.verification.status == "Unconfirmed",] 
Melesmeles <- Melesmeles[!Melesmeles$Identification.verification.status == "Unconfirmed - plausible",] 
Melesmeles <- Melesmeles[!Melesmeles$Identification.verification.status == "Unconfirmed - not reviewed",] 
unique(Melesmeles$Identification.verification.status) #Only accepted observations left

#Processing coordinate certainty
Melesmeles <- Melesmeles[Melesmeles$Coordinate.uncertainty_m  <=  1000,] #Remove coordinate uncertainty >1000m

sum(is.na(Melesmeles)) #Check for NA values, None

#Rename column names
colnames(Melesmeles)[3] <- "Latitude"
colnames(Melesmeles)[4] <- "Longitude" 

#Spatial processing
SpatialMelesmeles<-data.frame(x=Melesmeles$Longitude,y=Melesmeles$Latitude) #Extract coordinate data
crs.latlong<-crs("+init=epsg:4326") #Set desired coordinates system
SpatialMelesmeles <-SpatialPoints(SpatialMelesmeles,proj4string = crs.latlong)
crs(SpatialMelesmeles) #Check the points have the desired crs

#Exploring the badger data
plot(SpatialMelesmeles) #Outputs the Melesmeles observations for the whole UK
mapview(SpatialMelesmeles) #Interactive map of UK-wide badger distribution

#Crop the spatial extent to study area
studyExtent<-c(-4.2,-2.7,56.5,57.5) #list coordinates in the order: min Long, max Long, min Lat, max Lat
SpatialMelesmeles<-raster::crop(SpatialMelesmeles,studyExtent) #crop(x, y), x = spatial object, y = extent object

mapview(SpatialMelesmeles) #Interactive map of study area badgers (Scotland)

#Raster set-up
LCM<-raster("LCMUK.tif") #read in the raster data
crs(LCM) #compare coordinates systems - LCM is set to OSGB, meaning both datasets have different projections
crs(SpatialMelesmeles)
SpatialMelesmeles<-spTransform(SpatialMelesmeles,crs(LCM)) #Project badger crs to LCM crs
MelesmelesCoords <-coordinates(SpatialMelesmeles) #Extract badger coordinates

#Use projected badger coordinates to set up a new map extent
#Set geographical limits
x.min <- min(MelesmelesCoords[,1]) - 5000
x.max <- max(MelesmelesCoords[,1]) + 5000
y.min <- min(MelesmelesCoords[,2]) - 5000
y.max <- max(MelesmelesCoords[,2]) + 5000

extent.new <- extent(x.min, x.max, y.min, y.max) #Produce an object containing these limits

#Crop the raster to the desired parameters
LCM <- crop(LCM, extent.new)

#Visualise the mapping data
plot(LCM)
plot(SpatialMelesmeles,add=TRUE)

#Data Clean-up
rm(x.max, x.min, y.max, y.min, studyExtent, crs.latlong, extent.new, MelesmelesCoords)
#----------------------------------------------------------------------------------------------------------------------
#3 - Generating pseudo-absence points
set.seed(11) #Initialise the pseudorandom number generator

back.xy <- randomPoints(LCM, p=SpatialMelesmeles, n=1000,ext = extent(SpatialMelesmeles)) #Generate random points
back.xy<-SpatialPoints(back.xy) #Convert the random points into a SpatialPoints dataframe

#Plot the background points
plot(LCM) #Plot the map
plot(SpatialMelesmeles,add=T) #Plot the presence points on the map with the "add=" argument
plot(back.xy,add=TRUE, col='red') #Plot the background points as red 

#Calibrating coordinate systems between species data & background points data
A<-SpatialPoints(back.xy,crs(LCM)) #create a spatial points layer from the background points using the crs of the LCM object.
crs(A) #check the crs
P<-SpatialMelesmeles

#Characterising point locations
eA<-extract(LCM,A) #Extract values from the LCM layer to the background points
eP<-extract(LCM,P)  #Extract values from the LCM layer to the presence points

table(eA) #Peek at the background points data
table(eP) #Peek at the species points data

par(mfrow=c(1,2)) #Set up the plotting pane for data presentation (rows=1, columns =2) 

#Producing histograms of landcover at background points
hist(eA,freq=FALSE,breaks=c(0:21),xlim=c(0,21),ylim=c(0,1)) #plot a histogram of land cover at background points (landcover type on x-axis, proportions on y-axis.) 
hist(eP,freq=FALSE,breaks=c(0:21),xlim=c(0,21),ylim=c(0,1)) #plot a histogram of land cover at species presence points (landcover type on x-axis, proportions on y-axis.) 

hist(eA,freq=FALSE,breaks=c(0:21),xlim=c(0,21),ylim=c(0,0.4))
hist(eP,freq=FALSE,breaks=c(0:21),xlim=c(0,21),ylim=c(0,0.5))

#Produce a dataframe to store the presence/absence data
Abs<-data.frame(A,Pres=0) # create a data frame from the extracted background values
Pres<-data.frame(P,Pres=1) # create a data frame from the extracted presence values 
MelesmelesData<-rbind(Pres,Abs) # bind the two data frames by row

head(MelesmelesData) #display the first several rows
tail(MelesmelesData) #display the last several rows
class(MelesmelesData) #check the object's class, data.frame

#4 - Raster Reclassification
#Reclassifying the raster to extract broadleaf woodland landcover
LCM<-as.factor(LCM) #Redefine the raster as a categorical variable
levels(LCM) #Shows the level attributes of the redefined raster 
reclass <- c(0,1,rep(0,19)) #Set the reclassification criteria
RCmatrix<- cbind(levels(LCM)[[1]],reclass) #Combine the reclassification criteria and LCM levels to produce a matrix containing only broadleaf landcover data
broadleaf <- reclassify(LCM, RCmatrix) #Reclassify the raster using the matrix to include only the broadleaf (deciduous) landcover data
plot(broadleaf) #Plot the broadleaf raster
plot(SpatialMelesmeles,add=TRUE) #Add the badger data

#Data Clean-up
rm(reclass, back.xy, RCmatrix)

#----------------------------------------------------------------------------------------------------------------------
#5 - Buffer analysis 

#Objects required for buffer analysis:
#bufferlandcover2km - raster object for the broadleaf woodland within each buffer
#grainarea - numeric holding area described by each cell in the raster 
#bufferArea - numeric to hold the buffer area in hectares (ha)
#landcover2km - Numeric to hold the total area of woodland within the buffer
#percentlandcover2km - Object containing the final result as a percentage of the buffered area

buf2km<-2000 #Create an object to hold the distance parameter for the 2km buffer
buffer.site1.2km<-buffer(SpatialMelesmeles[1,],width=buf2km) #Apply the buffer() function to the first item of SpatialMelesmeles
zoom(broadleaf,buffer.site1.2km) #Close-up of the result
plot(buffer.site1.2km,border="red",lwd=5,add=T) #Add the buffer
dev.off() #Clear the active graphics pane

#Crop the broadleaf layer to the extent of the buffer 
buffer2km <- crop(broadleaf, buffer.site1.2km)
bufferlandcover2km <- mask(broadleaf, buffer.site1.2km) #Clip the cropped layer to the area defined by the buffer
bufferArea <- (3.14159*buf2km^2) #calculate the area of the buffer according to the buffer width
landcover2km <- cellStats(bufferlandcover2km, 'sum')*625 #Calculate the total area of woodland inside the buffer, 625m = the area of each cell of the 25m landcover raster) 
percentlandcover2km <- landcover2km/bufferArea*100 #calculate the percentage of landcover within the buffer
percentlandcover2km #31.11482m of broadleaf woodland in the buffer zone

#Create functions to automate the producing of buffers
radii<-seq(100,2000,by=100) #Set up an object containing instructions to produce buffers between 100-2000m at 100m intervals
res <- matrix(NA, nrow = nrow(MelesmelesData), ncol = length(radii)) #Apply the buffer instructions to the badger data by a matrix

#Function to automate the buffer calculations
landBuffer <- function(i){         
  
  MelesmelesPoints<-SpatialPoints(cbind(MelesmelesData[i,1],MelesmelesData[i,2]))  #Select each badger point to draw buffers around
  MelesmelesBuffer <- buffer(MelesmelesPoints, width=radii[j])           #buffer each point
  bufferlandcover <- crop(broadleaf, MelesmelesBuffer)                 #crop the landcover layer to the buffer extent
  masklandcover <- mask(bufferlandcover, MelesmelesBuffer)            # mask the above to the buffer
  landcoverArea <- cellStats(masklandcover, 'sum')*625      #use cellStats() to sum landcover area
  percentcover <- landcoverArea/gArea(MelesmelesBuffer)*100           # convert to precentage cover
  
  return(percentcover)} #Return the result

#For loop to apply each buffer size to each coordinate in MelesmelesData (presence/absence data)
for(i in 1:nrow(MelesmelesData)){
  for(j in seq_along(radii)){
    res[i,j]<-landBuffer(i)}
}

#Convert the output into a dataframe to prepare it for data analysis
glmData<-data.frame(res)
colnames(glmData)<-c("w100","w200","w300","w400","w500","w600","w700","w800","w900","w1000","w1100","w1200","w1300","w1400","w1500","w1600","w1700","w1800","w1900","w2000")
head(glmData)

#Data Clean-up
rm(landcover2km, res, landBuffer, bufferlandcover2km, bufferArea, buf2km, buffer.site1.2km, percentlandcover2km, j, buffer2km)

#----------------------------------------------------------------------------------------------------------------------
#6 - Log-likelihood Statistic
glmData$Pres<-MelesmelesData$Pres #Add the presences data
head(glmData) #Peek at the GLM dataset

#Generate General Linear Models (GLMs) to test for the relationship between broadleaf cover & presence per buffer size:
glm100<-glm(Pres~w100,family = "binomial",data = glmData)
glm200<-glm(Pres~w200,family = "binomial",data = glmData)
glm300<-glm(Pres~w300,family = "binomial",data = glmData)
glm400<-glm(Pres~w400,family = "binomial",data = glmData)
glm500<-glm(Pres~w500,family = "binomial",data = glmData)
glm600<-glm(Pres~w600,family = "binomial",data = glmData)
glm700<-glm(Pres~w700,family = "binomial",data = glmData)
glm800<-glm(Pres~w800,family = "binomial",data = glmData)
glm900<-glm(Pres~w900,family = "binomial",data = glmData)
glm1000<-glm(Pres~w1000,family = "binomial",data = glmData)
glm1100<-glm(Pres~w1100,family = "binomial",data = glmData)
glm1200<-glm(Pres~w1200,family = "binomial",data = glmData)
glm1300<-glm(Pres~w1300,family = "binomial",data = glmData)
glm1400<-glm(Pres~w1400,family = "binomial",data = glmData)
glm1500<-glm(Pres~w1500,family = "binomial",data = glmData)
glm1600<-glm(Pres~w1600,family = "binomial",data = glmData)
glm1700<-glm(Pres~w1700,family = "binomial",data = glmData)
glm1800<-glm(Pres~w1800,family = "binomial",data = glmData)
glm1900<-glm(Pres~w1900,family = "binomial",data = glmData)
glm2000<-glm(Pres~w2000,family = "binomial",data = glmData)

#Generate a list of loglikelihood values per buffer zone
ll100<-as.numeric(logLik(glm100))
ll200<-as.numeric(logLik(glm200))
ll300<-as.numeric(logLik(glm300))
ll400<-as.numeric(logLik(glm400))
ll500<-as.numeric(logLik(glm500))
ll600<-as.numeric(logLik(glm600))
ll700<-as.numeric(logLik(glm700))
ll800<-as.numeric(logLik(glm800))
ll900<-as.numeric(logLik(glm900))
ll1000<-as.numeric(logLik(glm1000))
ll1100<-as.numeric(logLik(glm1100))
ll1200<-as.numeric(logLik(glm1200))
ll1300<-as.numeric(logLik(glm1300))
ll1300<-as.numeric(logLik(glm1300))
ll1400<-as.numeric(logLik(glm1400))
ll1500<-as.numeric(logLik(glm1500))
ll1600<-as.numeric(logLik(glm1600))
ll1700<-as.numeric(logLik(glm1700))
ll1800<-as.numeric(logLik(glm1800))
ll1900<-as.numeric(logLik(glm1900))
ll2000<-as.numeric(logLik(glm2000))

#Compile the loglikelihood values into a dataframe
ll<-c(ll100,ll200,ll300,ll400,ll500,ll600,ll700,ll800,ll900,ll1000,ll1100,ll1200,ll1300,ll1400,ll1500,ll1600,ll1700,ll1800,ll1900,ll2000)
Dist<-c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000) #Generate distance values
glmRes<-data.frame(cbind(Dist,ll)) #Combine loglikelihood and distance values into a dataset

#Plot the log-likelihood values
plot(glmRes$Dist, glmRes$ll, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "buffer", ylab = "logLik")

#Data Cleanup
rm(glm100,glm200,glm300,glm400,glm500,glm600,glm700,glm800,glm900,glm1000,glm1100,glm1200,glm1300,glm1400,glm1500,glm1600,glm1700,glm1800,glm1900,glm2000)
rm(ll100,ll200,ll300,ll400,ll500,ll600,ll700,ll800,ll900,ll1000,ll1100,ll1200,ll1300,ll1400,ll1500,ll1600,ll1700,ll1800,ll1900,ll2000)
rm(Dist, ll, glmData)
#----------------------------------------------------------------------------------------------------------------------
#7 - Determining the optimum buffer size
opt<-glmRes[which(glmRes$ll==max(glmRes$ll)),] #Finding the optimum buffer size, based upon the log likelihood statistics
opt #print value, 1400m as the optimum buffer size

fw.test <- focalWeight(broadleaf, c(1,50),type= 'Gauss') #Apply a Gaussian 
sigma<-c(100,200,300,400,500) #Set up a list of sigma values to test
loglikFocal<- rep(NA,length=length(sigma)) #Create empty vector to hold results

#Evaluating the effects of the selected optimum buffer
#Function to loop through the sigma values
for(i in sigma){
  print (i)
  fw.i <- focalWeight(broadleaf, c(i,1400),type= 'Gauss') #Buffer around each cell of the broadleaf layer and apply a Gaussian filter
  focalBroadleaf<-focal(broadleaf,fw.i,sum) #Focal() function creates the raster
  
  #Producing a presence/absence dataframe
  eA<-extract(focalBroadleaf,A) #Extract absence values from the raster
  eP<-extract(focalBroadleaf,P) #Extract presence values from the raster
  Abs<-data.frame(eA,Pres=0) #Produce absence dataframe
  Pres<-data.frame(eP,Pres=1) #Produce presence dataframe
  names(Abs)<-c("area","Pres") #Rename columns
  names(Pres)<-c("area","Pres") #Rename columns
  head(Pres) #Peek at the data
  MelesmelesData<-rbind(Pres,Abs) #Combine the presence and absence datasets together
  head(MelesmelesData) #Peek at the data
  
  pres.i <- glm(Pres ~ area, family = "binomial", data = MelesmelesData) #Test for the relationship between species presence data and focalbroadleaf using a GLM
  LLi<-as.vector(logLik(pres.i)) #Calculate the log-likelihood statistic
  loglikFocal[i]<-as.numeric(LLi) #Set to as numeric for data compilation
  
}

#Plot the logLikelihood statistic
loglikFocal1400<-loglikFocal[!is.na(loglikFocal)] #remove NA values
result<-data.frame(sigma,loglikFocal1400) #Combine sigma and loglikelihood values into a dataframe
plot(result$sigma, result$loglikFocal1400, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "sigma", ylab = "logLik") #Plot the log Likelihood statistic

#Final Data Cleanup
rm(Melesmeles, SpatialMelesmeles, MelesmelesData, LCM)
rm(broadleaf, i, loglikFocal1400, radii, sigma, fw.test, P, opt, result, LLi, pres.i,fw.i, loglikFocal,focalBroadleaf, eA, eP, A, Pres, Abs,glmRes)
dev.off()

#======================================================================================================================
#Section 2 - Resource Selection & Home Range

#Aim: to investigate the resource selection & home ranges of three specific Canis lupus individuals in Canada using Kernel Density Estimations (KDEs)

#Tasks:
#Use the kernel Use Density approach to map the core home ranges (polygon) of animals 1 to 3 in the dataset (animal IDs 3,6 and 10)
#Produce spatial plots (polygons) of core ranges 
#Calculate the core range in km2 of each animal separately and all three animals combined.

#Sections:
#0- Setting up the working space
#1 - Data Preparation
#2 - Data Mapping
#3 - Calculating Core Home Ranges 
#4 - Home range size calculations (km²)

#----------------------------------------------------------------------------------------------------------------------
#0 - Setting up the working space

setwd("E:/0.Masters/1.5_GIS_Workspace/4.SE/SA1")

#Load packages
library(sp) # for spatial (vector) objects
library(dismo) # for generating background data
library(adehabitatHR) # for computing home ranges
library(adehabitatHS) # for computing resource selection
library(rgdal) #for projection coordinate systems
library(raster) #for spatial conversions
library(maptools) #for mapping reference data
library(reshape2) #for data management
library(rgeos) #for buffers & shape geometry
library(RColorBrewer) #for generating map legends

#Read in files
Canislupus<-read.csv("Wolves.csv") #telemetry data
CA<-raster("AlbertaCanis.tif") #bring in raster data
legendCA<-read.csv("LegendCA.csv") #Landcover legend data

#----------------------------------------------------------------------------------------------------------------------
#1 - Data Preparation
Canis <- Canislupus[Canislupus$animal_ID %in% c(3, 6, 10),] #Select data from the three selected individuals (Wolves 3, 6, 10)
Canis<-Canis[,4:9] #select only the desired data columns (columns 6 to 9 of the .csv file)
Canis.latlong<-data.frame(x=Canis$Lon,y=Canis$Lat) #Create an object containing the coordinate data (Canis telemetry data)
crs.latlong<-crs("+init=epsg:4326") #Create an object for the desired projection: (WGS84)
SpatialWolves<-SpatialPointsDataFrame(coords=Canis.latlong,data=Canis,proj4string = crs.latlong) #Convert the Canis telemetry data into a spatial points dataframe
crs(SpatialWolves) #check crs: +proj=longlat +datum=WGS84 +no_defs 
plot(SpatialWolves) #plot the data for all 3 individuals

#Setting up the raster coordinate system
crs(SpatialWolves) #Two different coordinate systems: +proj=longlat +datum=WGS84 +no_defs 
crs(CA) #+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs 
SpatialWolves <- spTransform(SpatialWolves,crs(CA)) #Set the Wolf telemetry data to the raster coordinates system

SpatialWolvesCol <- factor(SpatialWolves$animal_ID) #Colour code per animal_ID in preparation for plotting

#----------------------------------------------------------------------------------------------------------------------
#2 - Data Mapping
#Plotting telemetry point data with landcover
#Generating the landcover legend
CA<-as.factor(CA) #Convert the raster to a factor
levelsCA<-levels(CA) #Create an object to store the landcover levels
levelsCA<-data.frame(levelsCA) #Convert the levels object to a data frame to be accessible
namesCA<-legendCA[legendCA$Code %in% levelsCA$ID,] #Extract only the legend items (landcover types) present in the raster
coul <- colorRampPalette(brewer.pal(8, "Accent"))(14) #Create a new palette object with 14 combinations of the 8 colours in the "Accent" palette

plot(CA,legend = FALSE,axes=T,box=T,col = coul) #Plot again with the new palette
par(xpd = TRUE)# set "xpd" to TRUE to allow legend plotting outside of the map frame

#Add a legend using the labels from namesCA and the colours from the new palette object
legend(-1000000,975000,paste(namesCA[,2]), fill = coul,
       cex = 0.7,bty="n", title = "Landcover Types" ) # cex = character size, bty ="n" removes the box around the legend
palette("Polychrome 36") #set a new default palette, this time for plotting the Canis points
plot(SpatialWolves, pch = 16,col=SpatialWolvesCol,add=T)#plot them using the factor canisCol created earlier but with the more distinct colour palette

#use legend function to fit a legend on the opposite side of the plot
legend(-1000000,1030000,legend = paste(levels(SpatialWolvesCol)), col = 1:length(SpatialWolvesCol), cex = 1, pch = 19, title = "Canis lupis IDs", text.font  =2)

#Data Clean-up
rm(coul, SpatialWolvesCol, coul, levelsCA, namesCA, legendCA, crs.latlong, Canis.latlong)

#----------------------------------------------------------------------------------------------------------------------
#3 - Calculating Core Home Ranges
#Kernel density estimations - individual
IndividualWolves <- kernelUD(SpatialWolves[,5], h = "href", kern ="bivnorm") #Kernel Density Estimation (KDE) per individual: Estimating the kernel home-range (KDE) per individual

#Kernel density estimations - group
AllWolves<-kernelUD(SpatialWolves[,3], h = "href", kern ="bivnorm") #Kernel Density Estimation (KDE) for all animals: Estimating the kernel home-range (KDE) for all animals

#Visualising the KDE outputs
class(IndividualWolves) #data type = estUDm, estimated UseDensity
image(IndividualWolves) #Plotting the Kernel Density Estimation (KDE) per individual 

class(AllWolves) #data type = estUDm, estimated UseDensity
image(AllWolves) #Plotting the Kernel Density Estimation (KDE) grouped

#Setting up the KDEs for mapping home ranges
#Converting the all-animals KDE output to raster format
r<-estUDm2spixdf(AllWolves) #Conversion of estUDm result to a Spatial Pixels Dataframe
class(r) #data type = SpatialPixelsDataFrame
rCanis<-raster(r) #Rasterise the KDE for all individuals
rCanis <- reclassify(rCanis, cbind(-Inf, 0, NA), right=TRUE) #remove zeros values from the  map by setting all values <=0 to NA. 
polyCanisAll<-getverticeshr(AllWolves, percent = 50) #Extract the general-range contours of the rasterised KDEs 

#Calculating core home ranges for all Wolves
plot(rCanis) #Plot the KDE output
plot(polyCanisAll,add=TRUE) #Add the contour

dev.off() #Reset the plot environment

#For loop to generate home ranges per individual:
for (i in unique(SpatialWolves$animal_ID)){ #for each individual wolf
  
  SpatialWolves.i<-SpatialWolves[which(SpatialWolves$animal_ID==i),]
  
  SpatialWolvesID.i<-data.frame(SpatialWolves.i$animal_ID) #Specify i as representing each individual wolf
  xy.i<-cbind(SpatialWolves.i$Lon,SpatialWolves.i$Lat) #Create a new object combining individuals with coordinate data (long/lat) 
  coordinates(SpatialWolvesID.i)<-xy.i #Extract coordinate data
  kud.i <- kernelUD(SpatialWolves.i[,3], h="href", kern ="bivnorm") #Create kernels per individual
  
  r.i<-estUDm2spixdf(kud.i) #Convert the Kernel Density Estimation (KDE) to a spatial points dataframe
  Ras.i <-raster(r.i) #Rasterise the spatialpixelsdataframe
  rCanis.i <- reclassify(Ras.i, cbind(-Inf, 0, NA), right=TRUE) #Remove zero & NA values
  polyCanis.i<-getverticeshr(kud.i, percent = 50) #Generate polygon home-range contours per individual KDE
  clipCanis.i<-mask(rCanis.i,polyCanis.i) #Clip overlapping polygon home ranges
  
  #Plot the individual Kernel Density Estimations (KDEs) for home ranges:
  plot(clipCanis.i,axes=TRUE,box=FALSE,legend=FALSE,add=T)
  plot(polyCanis.i,main=i,add=TRUE)
  
}

plot(polyCanisAll,add=TRUE) #Add the contour to the map for comparison

#Data Clean-up
rm(IndividualWolves, AllWolves, r, rCanis, SpatialWolves.i, i, SpatialWolvesID.i, polyCanis.i, polyCanisAll, r.i, Ras.i, rCanis.i, polyCanis.i, clipCanis.i, kud.i, xy.i)

#----------------------------------------------------------------------------------------------------------------------
#4 - Home range size calculations (km²)
#Calculating the size of each individual's home range (km²):
IndividualRangeWolves<-mcp.area(SpatialWolves[,5],percent = seq(20,100, by = 5),
                      unin = "m",
                      unout = "km2") #Creates a graphical output of 3 graphs representing home ranges per individual - #x = home-range level (km2), y = home-range size (km)

#Calculating the size of all individual's home range (km²):
AllRangeWolves<-mcp.area(SpatialWolves[,3],percent = seq(20,100, by = 5),
                      unin = "m",
                      unout = "km2") #Creates a graphical output of 3 graphs representing home ranges per individual - #x = home-range level (km2), y = home-range size (km)


#Final Data Clean-up
rm(CA, Canis, Canislupus, SpatialWolves)
rm(IndividualRangeWolves, AllRangeWolves)

dev.off() #Reset the plot environment

#======================================================================================================================
#Section 3 - Species Distribution Modelling

#Aim: to investigate the species distribution of Solanum acaule using maxent & GLM approaches

#Tasks:
#Use the gbif() function to download records for Solanum acaule
#Generate species distribution models using the maxnet() and glm() approaches
#Use the same environmental data and approach to background sampling as used in the practical
#Evaluate the models using a k-fold approach based on 5 folds
#Select the best performing model and create 2 mapped outputs:
  #1. a probability plot of the predictions
  #2. a presence/absence map based on a sensible threshold

#Sections:
#0- Setting up the working space
#1 - Data Preparation & Mapping
#2 - Model Set-up
#3 - GLM Model
#4 - Maxnet Model
#5 - Model Evaluation
#6 - Presence/Absence Map

#----------------------------------------------------------------------------------------------------------------------
#0- Setting up the working space

setwd("E:/0.Masters/1.5_GIS_Workspace/4.SE/SA1")

#Load in packages
library(dismo)
library(maptools)
library(glmnet)
library(maxnet)
library(raster)
library(sp)
library(randomForest)

#----------------------------------------------------------------------------------------------------------------------
#1 - Data Preparation & Mapping
#Locate the environmental data in dismo
files <- list.files(path=paste(system.file(package="dismo"),
                               '/ex', sep=''), pattern='grd', full.names=TRUE )

files #Peek at the file location of the environmental data

#A - Load the raster data into the workspace
env<-stack(files[1:8]) #8 raster layers
#Label each covariate layer
names(env) <- c("MeanTemp.","AnnualPrecip.","Precip.WettestQuarter","PrecipDriestQuarter","Max.Temp.","Min.Temp.","Temp.Range","MeanTemp.WettestQuarter")

plot(env) #Plot data

#B - Load in the vector basemap
data(wrld_simpl) 

#C - Species presence data
Solanum <- gbif("Solanum", "acaule*",sp=TRUE) #6934 records found

#Plot the Solanum data with MeanTemp
plot(env, 1)
plot(wrld_simpl,add=TRUE)
plot(Solanum,add=TRUE,col='red') #Distributed across Latin & South America

#Select subregions relevant to Bradypus species distributions
View(wrld_simpl@data)
SA<-subset(wrld_simpl,wrld_simpl$SUBREGION==5|wrld_simpl$SUBREGION==13) #Select subregions 5 & 13
plot(SA) #Just the species distribution data now

crs(Solanum)<-crs(SA) #set the solanum crs to that of the SA object created above
Solanum<-Solanum[SA, ] #subset presence points to only those inside SA
plot(SA) #Plot the background map
plot(Solanum, add=TRUE) #Plot the species data

#----------------------------------------------------------------------------------------------------------------------
#2 - Model Set-up

#Extracting background environmental data
#Create random points on cells of the 'env' object within the extent of solanum and avoiding cells containing points from Solanum
set.seed(11)
back.xy <- randomPoints(env, n=3000,p= Solanum,ext = extent(Solanum)) #create matrix (coordinates) of background points 
back<-SpatialPoints(back.xy,crs(Solanum)) #create spatial points object from back.xy
crs(back)<-crs(SA) #make sure crs matches other data
back<-back[SA,] #subset to SA object

#Visualise the background points
plot(SA) #Plot the background map 
plot(Solanum,add=TRUE) #Plot the species data
plot(back,add=TRUE, col='red') #Plot the background points

#Extract information from the 'env' raster stack to the presence & background points
eA<-extract(env,back) #Extract (pseudo)absence values
eP<-extract(env,Solanum) #Extract presence values

Pres.cov<-data.frame(eP,Pres=1) #Create data frame from presence values, species data
Back.cov<-data.frame(eA,Pres=0) #Create data frame from absence values, background data
all.cov<-rbind(Pres.cov,Back.cov) #combine the presence and absence data

#Investigate the datasets
head(all.cov)
tail(all.cov)
summary(all.cov)
nrow(all.cov)

#----------------------------------------------------------------------------------------------------------------------
#3 - GLM Model
GLM.Solanum<-glm(Pres~.,binomial(link='logit'), data=all.cov) #Generate a GLM testing the species distribution relationship between presence and absence values
summary(GLM.Solanum) #List the GLM outputs

glm.map <- predict(env, GLM.Solanum, type = "response") #Using the GLM to estimate probabilities of Solanum presence based on the environmental data
plot(glm.map,main="GLM Model of Solanum acaule species distribution") #Plot the GLM Species Distribution Model

#----------------------------------------------------------------------------------------------------------------------
#4 - Maxnet Model
#Run a presence-only Maxnet model selecting feature "linear" and "quadratic" using the "lq" argument.  
SolanumMaxent <- maxnet(all.cov$Pres, all.cov[,1:8],maxnet.formula(all.cov$pres,all.cov[,c("MeanTemp.", "AnnualPrecip.", "Precip.WettestQuarter", "PrecipDriestQuarter", "Max.Temp.", "Min.Temp.", "Temp.Range")],classes="lq")) 
maxnet.cloglog.map <- predict(env, SolanumMaxent, clamp=F, type="cloglog") #Estimate the probabilities of Solanum presence based upon the env data
plot(maxnet.cloglog.map, axes=T, box=T, main="Maxnet-clog of of Solanum acaule species distribution") #Plot the Maxnet Species Distribution Model

#----------------------------------------------------------------------------------------------------------------------
#5 - Model Evaluation
#Selecting the best model using the k-fold approach

set.seed(5) #Initialise the pseudorandom number generator
folds=5 #set the number of folds to use, = 5 folds

#Separate presence and absence data according to the five folds, using the kfold() function:
kfold_pres <- kfold(Pres.cov, folds) #Separating presence data
kfold_back <- kfold(Back.cov, folds) #Separating absence data

#A - GLM k-fold
#Create an empty list to hold the GLM k-fold results (creates five sets)
eGLM<-list() #Empty list
par(mfrow=c(2,3)) #Set up the plot environment

#For loop function to iterate over folds
for (i in 1:folds) {
  train <- Pres.cov[kfold_pres!= i,] #For presence values, select all folds which are not 'i' to train the model
  test <- Pres.cov[kfold_pres == i,] #Use the remaining fold to test the model
  backTrain<-Back.cov[kfold_back!=i,] #For background values, select training folds
  backTest<-Back.cov[kfold_back==i,] #Use the remaining fold to test the model
  dataTrain<-rbind(train,backTrain) #Bind the presence and background training data together
  dataTest<-rbind(test,backTest) #Bind the presence and background test data together
  glm_eval <- glm(Pres~.,binomial(link = "logit"), data=dataTrain)#GLM model trained on presence and absence points
  eGLM[[i]] <- evaluate(p=dataTest[ which(dataTest$Pres==1),],a=dataTest[which(dataTest$Pres==0),], glm_eval) #Use the testing data (kfold==i) for model evaluation
  
  #check the AUC by plotting ROC values
  plot(eGLM[[i]],'ROC') #Plot the model evaluation outputs (Produces a Receiver Operator’s Characteristic Curve)
  #Relatively high AUC values per K-fold: 0.973, 0.97, 0.976, 0.982, 0.979
}

#Inspect the model evaluation output
eGLM #Display model outputs for each K-fold (n presences, n absences, AUC, core, max TPR+TNR)
aucGLM <- sapply( eGLM, function(x){slot(x, 'auc')} ) #Extract the AUC values from the GLM K-fold

#calculate the mean AUC value for comparison with other model
mean(aucGLM) #Mean = 0.9759049

#B - Maxnet k-fold
#Create an empty list to hold the Maxnet k-fold results (creates five sets)
par(mfrow=c(2,3)) #Set up the plot environment
eMAX<-list() #Empty list
folds=5 #Set the number of folds to use, = 5 folds

#Separate presence and absence data according to the five folds, using the kfold() function:
kfold_pres <- kfold(Pres.cov, folds) #Separating presence data
kfold_back <- kfold(Back.cov, folds) #Separating absence data

#For loop function to iterate over folds
for (i in 1:folds) {
  train <- Pres.cov[kfold_pres!= i,] #For presence values, select all folds which are not 'i' to train the model
  test <- Pres.cov[kfold_pres == i,] #Use the remaining fold to test the model
  backTrain<-Back.cov[kfold_back!=i,] #For background values, select training folds
  backTest<-Back.cov[kfold_back==i,] #Use the remaining fold to test the model
  dataTrain<-rbind(train,backTrain) #Bind the presence and background training data together
  dataTest<-rbind(test,backTest) #Bind the presence and background test data together
  maxnet_eval <- maxnet(dataTrain$Pres, dataTrain[,1:8]) #Maxnet model trained on presence data
  eMAX[[i]] <- evaluate(p=dataTest[which(dataTest$Pres==1),],a=dataTest[which(dataTest$Pres==0),], maxnet_eval) #Use the testing data (kfold==i) for model evaluation
 
  #check the AUC by plotting ROC values
  plot(eMAX[[i]],'ROC') #Plot the model evaluation outputs (Produces a Receiver Operator’s Characteristic Curve)
  #Relatively high AUC values per K-fold: 0.989, 0.994, 0.987, 0.982, 0.988
}

#Inspect the model evaluation output
eMAX #Display model outputs for each K-fold (n presences, n absences, AUC, cor, max TPR+TNR)
aucMAX <- sapply(eMAX, function(x){slot(x, 'auc')} ) #Extract the AUC values from the Maxnet K-fold

#calculate the mean values for comparison with other model
mean(aucMAX) #Mean = 0.9878774

#Selected the Maxnet model as the best performing model, according to the mean AUC values calculated using the k-fold approach

#6 - Presence/Absence Map
#Get the maxTPR+TNR for the maxnet model
Opt_MAX<-sapply( eMAX, function(x){ x@t[which.max(x@TPR + x@TNR)] } )
Opt_MAX # -8.707803 -8.796520 -8.832066 -8.799503 -8.743151
Mean_OptMAX<-mean(Opt_MAX)
Mean_OptMAX #-8.775809

#Plot the Maxnet presence/absence map
prMAX <- predict(env, SolanumMaxent) #Using Maxnet model to estimate probabilities of Solanum presence based upon the environmental data
par(mfrow=c(1,2)) #Set up the plot environment
plot(prMAX, main='Maxent Prediction of Solanum acaule')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(prMAX > Mean_OptMAX, main='Presence/absence Map of Solanum acaule')
plot(wrld_simpl,add=TRUE, border ='dark grey')
