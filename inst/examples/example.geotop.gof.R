library(geotopOptim)
library(ggplot2)

data(MonteciniB2)
## OBSERVATION PROCESSING


obs_SWC <- MonteciniB2[str_detect(names(MonteciniB2),"SWC")]
zvalues <-  as.numeric(unlist(lapply(X=str_split(names(obs_SWC), pattern="", n = Inf),FUN=function (x) {
							out <- as.numeric(x)
							out <- out[!is.na(out)]
							out <- paste(out,collapse="")
							return(out)
						})))
zformatter = "z%04d"
names(obs_SWC) <- sprintf(zformatter,zvalues)
obs_SWC <- lapply(X=obs_SWC,FUN=function(x){
			
			if (length(dim(x))>1) {
				max <- apply(X=x,MARGIN=1,FUN=max,na.rm=TRUE)
				min <- apply(X=x,MARGIN=1,FUN=min,na.rm=TRUE)
				mean <- apply(X=x,MARGIN=1,FUN=mean,na.rm=TRUE)
				sd <- apply(X=x,MARGIN=1,FUN=sd,na.rm=TRUE)
			} else {
				
				mean <- as.vector(x)
				max <-  as.vector(x)
				min <-  as.vector(x)
				sd <- NA
				
			}
			out <- data.frame(min=min,mean=mean,max=max,sd=sd)
			
			out <- as.zoo(out)
			index(out) <- as.POSIXlt(index(x))
			
			return(out)
			
		})
###########
###########




simpath <- system.file("Montecini_pnt_1_225_B2_004",package="geotopOptim")
bin <-   "/home/ecor/local/geotop/GEOtop/bin/geotop-2.0.0"
runpath <- "/home/ecor/temp/geotopOptim_tests"

vars <- "SoilLiqContentProfileFile"

sim <- geotopZProfile(bin=bin,simpath=simpath,runpath=runpath,clean=TRUE,variable=vars,data.frame=TRUE,level=1,zformatter=zformatter)


modeled <- sim$SoilLiqContentProfileFile$z0020
observed <- obs_SWC$z0020
m <- merge(observed,modeled)
m <- m[!is.na(m$modeled),]
## Plotting standard deviation of observation vs observed range (max-min)

plot(m$max-m$min,m$modeled-m$min)
abline(0,1)
points(m$max-m$min,m$modeled-m$mean,pch=2)


plot(m$sd,m$mean-m$modeled)



## plot(