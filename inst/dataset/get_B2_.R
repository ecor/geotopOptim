# TODO: Add comment
# 
# Author: ecor
###############################################################################

#!/usr/bin/env Rscript

rm(list=ls())

library(stringr)
library(zoo)

file_SWP  <- "/Volumes/AlpEnv/Projekte/HiResAlp/06_Workspace/BrJ/02_data/Station_data_Mazia/processed/SWP/SWP_h_B2.csv"
file_SWC  <- "/Volumes/AlpEnv/Projekte/HiResAlp/06_Workspace/BrJ/02_data/Station_data_Mazia/processed/SWC/SWC_h_B2.csv"
file_SoilTemp <- "/Volumes/AlpEnv/Projekte/HiResAlp/06_Workspace/BrJ/02_data/Station_data_Mazia/processed/SoilTemp/SoilTemp_h_B2.csv"

fileCopy <- TRUE
##datadir <- "/Users/ecor/Dropbox/R-activity-2/mazia/data"
datadir <- "/Users/ecor/ownCloud/activity/mazia/data"

if (fileCopy==TRUE) {
	
	
	
	file.copy(from=file_SWP,to=datadir)
	file.copy(from=file_SWC,to=datadir)
	file.copy(from=file_SoilTemp,to=datadir)
	
}

file_SWC <- list.files(datadir,pattern="SWC_h_B2.csv",full.name=TRUE)[1]
file_SWP <- list.files(datadir,pattern="SWP_h_B2.csv",full.name=TRUE)[1]
file_SoilTemp <- list.files(datadir,pattern="SoilTemp_h_B2.csv",full.name=TRUE)[1]

SWC <- read.table(file_SWC,sep=",",header=TRUE)
SWP <- read.table(file_SWP,sep=",",header=TRUE)
SoilTemp <- read.table(file_SoilTemp,sep=",",header=TRUE)
		
		
format_date <- "%m/%d/%Y %H:%M:"
tz <- "Etc/GMT+1"
SWC$Date <- as.POSIXct(SWC$Date,format=format_date,tz=tz)
SWP$Date <- as.POSIXct(SWP$Date,format=format_date,tz=tz)
SoilTemp$Date <- as.POSIXct(SoilTemp$Date,format=format_date,tz=tz)


Date <- SWC$Date
Date <- Date[Date %in% SWP$Date]
Date <- Date[Date %in% SoilTemp$Date]

SWC <- SWC[SWC$Date %in% Date,]
SWP <- SWP[SWP$Date %in% Date,]
SoilTemp <- SoilTemp[SoilTemp$Date %in% Date,]

Date <- SWC$Date
SWC <- SWC[,!(names(SWC) %in% "Date")]
SWC <- as.zoo(SWC)
index(SWC) <- Date

Date <- SoilTemp$Date
SoilTemp <- SoilTemp[,!(names(SoilTemp) %in% "Date")]
SoilTemp <- as.zoo(SoilTemp)
index(SoilTemp) <- Date

####


Date <- SWP$Date
SWP <- SWP[,!(names(SWP) %in% "Date")]
### Unit Conversion from hPa to mm
waterdensity <- 1000 ##  kg/m^3
gravity <- 9.81 # m/s^2
SWP <- (SWP*100)/(waterdensity*gravity)*1000

### END Unit Conversion from hPa to mm
## Water Sction as Water Pressure (minus added!)
SWP <- -SWP

## Exclusion of "very dry" values 
#SWP[SWP<=(-300)] <- NA

## zoo object 
SWP <- as.zoo(SWP)
index(SWP) <- Date




##########

datapath <- "/Users/ecor/ownCloud/activity/R-Packages/geotopOptim/data"


SoilTemp00z <- SoilTemp[,str_detect(names(SoilTemp),"z0")]


SWC05z <- SWC[,str_detect(names(SWC),"z5") & !str_detect(names(SWC),"z50")]
SWP05z <- SWP[,str_detect(names(SWP),"z5") & !str_detect(names(SWP),"z50")] 
SoilTemp05z <- SoilTemp[,str_detect(names(SoilTemp),"z5") & !str_detect(names(SoilTemp),"z50")]

SWC20z <- SWC[,str_detect(names(SWC),"z20")]
SWP20z <- SWP[,str_detect(names(SWP),"z20")]
SoilTemp20z <- SoilTemp[,str_detect(names(SoilTemp),"z20")]

###


SWC50z <- SWC[,str_detect(names(SWC),"z50")]
SoilTemp50z <- SoilTemp[,str_detect(names(SoilTemp),"z50")]


MonteciniB2 <- list(
		SoilTemp00z=SoilTemp00z,
		SWC05z=SWC05z,SWP05z=SWP05z,SoilTemp05z=SoilTemp05z,
		SWC20z=SWC20z,SWP20z=SWP20z,SoilTemp20z=SoilTemp20z,
		SWC50z=SWC50z,SoilTemp50z=SoilTemp50z
		)

save(MonteciniB2,file=paste(datapath,"MonteciniB2.rda",sep="/"))


stop()
##########

sensors_SWC  <-  unlist(lapply(X=str_split(names(SWC),pattern="_"),FUN=function(x){x[2]}))
sensors_SWP  <-   unlist(lapply(X=str_split(names(SWP),pattern="_"),FUN=function(x){x[2]}))
sensors_SoilTemp  <-   unlist(lapply(X=str_split(names(SoilTemp),pattern="_"),FUN=function(x){x[2]}))


depth_SWC  <-  unlist(lapply(X=str_split(names(SWC),pattern="_"),FUN=function(x){x[3]}))
depth_SWP  <-   unlist(lapply(X=str_split(names(SWP),pattern="_"),FUN=function(x){x[3]}))
depth_SoilTemp <-  unlist(lapply(X=str_split(names(SoilTemp),pattern="_"),FUN=function(x){x[3]}))


names(SWC) <- paste(sensors_SWC,depth_SWC,sep="_")
names(SWP) <- paste(sensors_SWP,depth_SWP,sep="_")
names(SoilTemp) <- paste(sensors_SoilTemp,depth_SWP,sep="_")
##### 

depth <- "z20"

SWCd <- SWC[,depth_SWC==depth]
SWPd <- SWP[,depth_SWP==depth]
SoilTempd <- SoilTemp[,depth_SoilTemp==depth]


sensors_SWCd <- sensors_SWC[depth_SWC==depth]
sensors_SWPd <- sensors_SWP[depth_SWP==depth]
sensors_SoilTempd <- sensors_SoilTemp[depth_SoilTemp==depth]


names(SWCd) <- sensors_SWCd
names(SWPd) <- sensors_SWPd
names(SoilTempd) <- sensors_SoilTempd 

#########


plot(SWCd$CSn,SWPd$CSn)
plot(SWCd$CSt,SWPd$CSt)
 