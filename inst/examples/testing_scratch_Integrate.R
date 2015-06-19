# TODO: Add comment
# 
# Author: ecor
###############################################################################
source('~/Dropbox/R-packages/geotopOptim/R/integratefunDataFrame.R', chdir = TRUE)



data(MuntatschiniB2)
## OBSERVATION PROCESSING


obs_SWC <- MuntatschiniB2[str_detect(names(MuntatschiniB2),"SWC")]
zvalues <-  as.numeric(unlist(lapply(X=str_split(names(obs_SWC), pattern="",
								n = Inf),FUN=function (x) {
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


simpath <- system.file("Muntatschini_pnt_1_225_B2_004",package="geotopOptim")
bin <-   "/Users/ecor/local/bin/geotop_zh"
runpath <- "/Users/ecor/ownCloud/job"

vars <- "SoilLiqContentProfileFile"

sim <- geotopZProfile(bin=bin,simpath=simpath,runpath=runpath,
		clean=TRUE,variable=vars,data.frame=TRUE,level=1,zformatter=zformatter,intern=TRUE)[[vars]]


sim_int <- integratefunDataFrame(df=sim,formatter=formatter)
