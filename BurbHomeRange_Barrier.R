################################################################################
#################### HOME RANGE ESTIMATES OF BURBOT ############################
################################################################################
rm(list=ls(all=TRUE))

# NEED THESE PACKAGES TO RUN ALL OF THE CODE IN THIS SCRIPT
library(readr)
library(adehabitatHR)
library(tidyverse)
library(anytime)
library(PBSmapping)

# Merge all depth data to one file - file includes only unique detections #
Burb_Location = list.files(path = "E:/BSUGradBurbot/DATA/FishData/fourfish/",
                           pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                                                      # Store all files in list
  bind_rows                                                                 # Combine data sets into one data set 
Burb_Location 

#filter by desired HPE
Burb_Location1 <- Burb_Location %>% filter(HPE < 17)

#refomat text date from depth data as a factor to POSIXct
Burb_Location1$DATETIME = anytime(as.factor(Burb_Location1$DATETIME))

#convert lat long to utm
tagll <- data.frame(X = Burb_Location1$LON, Y = Burb_Location1$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
Burb_Location1$EASTING<- xyrec[,1]
Burb_Location1$NORTHING<- xyrec[,2]

x=Burb_Location1$EASTING
y=Burb_Location1$NORTHING

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Burb_Location1[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Burb_Location1[ , c(1,3,10,12,13)])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb <- SpatialPointsDataFrame(coords      = coords,
                               data        = data, 
                               proj4string = crs)

#create grid to calculate volume and area of home range
#https://stackoverflow.com/questions/41683905/grid-too-small-for-kernelud-getverticeshr-adehabitathr-home-range-estimation
xmin <- 314519.03
xmax <- 321519.03
ymin <- 5217540.97
ymax <- 5225540.97

x <- seq(xmin, xmax, by = 15.)
y <- seq(ymin, ymax, by = 15.)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

#calculate kernel density estimate for all fish
kud <- kernelUD(burb[,1], h="href",  grid = xy)

#print all fish estimates
image(kud)

#vectorized home range estimate 95 kud
homerange <- getverticeshr(kud)
class(homerange)
plot(homerange,col = 1:4)

#rasterized home range estimate 95 kud
vud <- getvolumeUD(kud)
vud

## Set up graphical parameters
par(mfrow=c(1,2))
par(mar=c(0,0,2,0))
## The output of kernelUD for the animal
image(kud[[2]])
title("Output of kernelUD")
## Convert into a suitable data structure for the use of contour
xyz <- as.image.SpatialGridDataFrame(kud[[2]])
contour(xyz, add=TRUE)
## and similarly for the output of getvolumeUD
par(mar=c(0,0,2,0))
image(vud[[2]])
title("Output of getvolumeUD")
xyzv <- as.image.SpatialGridDataFrame(vud[[2]])
contour(xyzv, add=TRUE)

## store the volume under the UD (as computed by getvolumeUD)
## of the animal in fud
fud <- vud[[2]]
## store the value of the volume under UD in a vector hr95
hr95 <- as.data.frame(fud)[,1]
## if hr95 is <= 95 then the pixel belongs to the home range
## (takes the value 1, 0 otherwise)
hr95 <- as.numeric(hr95 <= 95)
## Converts into a data frame
hr95 <- data.frame(hr95)
## Converts to a SpatialPixelsDataFrame
coordinates(hr95) <- coordinates(vud[[2]])
gridded(hr95) <- TRUE
## display the results
image(hr95)

#write 95 kud as data frame
HR <- as.data.frame(homerange)

#calculate kernel area estimates
ii <- kernel.area(kud, percent = seq(50,95, by = 5), unout = "km2")
ii <- as.data.frame(t(ii))
colnames(ii) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_FishInfo.csv")
fish <- fish[c(1,12,24,44),]

#condense all data in to file with HR estimates and fish info
HR_Burbs <- cbind(fish,ii)
HR_Burbs <- HR_Burbs %>% select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN, HR_50,HR_75,HR_80,HR_85,HR_90,HR_95)

#################################################################################
## INCORPORATE LAKE AS A BARRIER
#################################################################################

################################################################################
################# CODE FOR BARRIER JUST ON EAST SIDE OF LAKE ###################
################################################################################
library(rgdal)
library(maptools)

# Load KML coordinates
coords = getKMLcoordinates('E:\\BSUGradBurbot\\DATA\\BML Data\\BadMedEastBorderNEW.kml')
coords = SpatialPoints(coords, CRS('+proj=longlat'))

coords.proj = spTransform(coords, CRS=CRS('+init=epsg:32615 +units=m'))

a <- coords.proj@coords
b <- as.data.frame(a)

colnames(b) <- c("e","n","elev")

setwd("E:/BSUGradBurbot/DATA/BML Data")

#write.csv(b, "BML_east.csv")

Lake = read_csv("E:/BSUGradBurbot/DATA/BML Data/BML_east.csv")
str(Lake)

#convert read.table to list for bound layer
bound = structure(list(x = Lake$e, y = Lake$n), .Names = c("x","y"))
lines(bound, lwd=1)

str(bound)
## We convert bound to SpatialLines:
bound <- do.call("cbind",bound)
Slo1 <- Line(bound)
Sli1 <- Lines(list(Slo1), ID="frontier1")
barrier <- SpatialLines(list(Sli1), CRS("+init=epsg:32615"))
str(barrier)
## estimation of the UD with barrier
kud2 <- kernelUD(burb[,1], h="href", boundary = barrier)

################################################################################
############ CODE FOR WHOLE LAKE AS BARRIER ####################################
################################################################################
# Load KML coordinates
coords = getKMLcoordinates('E:\\BSUGradBurbot\\DATA\\BML Data\\BML_outline.kml')
coords = SpatialPoints(coords, CRS('+proj=longlat'))

coords.proj = spTransform(coords, CRS=CRS('+init=epsg:32615 +units=m'))

a <- coords.proj@coords
b <- as.data.frame(a)

colnames(b) <- c("e","n","elev")

setwd("E:/BSUGradBurbot/DATA/BML Data")

#write.csv(b, "BML_outline.csv")

Lake = read_csv("E:/BSUGradBurbot/DATA/BML Data/BML_outline.csv")

#convert read.table to list for bound layer
bound = structure(list(x = Lake$e, y = Lake$n), .Names = c("x","y"))
lines(bound, lwd=1)

## We convert bound to SpatialLines:
bound <- do.call("cbind",bound)
Slo1 <- Line(bound)
Sli1 <- Lines(list(Slo1), ID="frontier1")
barrier <- SpatialLines(list(Sli1), CRS("+init=epsg:32615"))

## estimation of the UD with barrier
kud3 <- kernelUD(burb[,1], h="href", boundary = barrier)
