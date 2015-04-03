#!/usr/bin/env Rscript

rm(list=ls())



library(geotopOptim)
library(ggplot2)

set.seed(1223)

data(MuntatschiniB2)
## OBSERVATION PROCESSING

##http://stackoverflow.com/questions/2151212/how-can-i-read-command-line-parameters-from-an-r-script
args<-commandArgs(TRUE)
### WORKING OUTPUT 
print(args)
###wpath_pso <-  "/Users/ecor/ownCloud/job/trial001"
wpath_pso <- argsParser(option="-wpath_out",args=args,novalue_response = ".")

print(wpath_pso)
print(getwd())
setwd(wpath_pso)
print(getwd())


####
 

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
				len <- apply(X=x,MARGIN=1,FUN=function(x) {length(x[!is.na(x)])})
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
bin <-   "/Users/ecor/ownCloud/geotop_se27xx/GEOtop/bin/geotop-2.0.0" 
###/Users/ecor/local/bin/geotop_zh"
runpath <- "/Users/ecor/ownCloud/job"
#
t <- str_split(simpath,"/")[[1]]
simdir <- t[length(t)]

rundir <- paste(runpath,simdir,sep="/")


#
vars <- "SoilLiqContentProfileFile"

### Use geotopGOF with an internal GEOtop simulation

## create a list with arguments for geotopGOF
#
#x <- param <- c(N=1.4,Alpha=0.0021,ThetaSat=0.5,ThetaRes=0.05,LateralHydrConductivity=0.021,NormalHydrConductivity=0.021)
#upper <- x*3
#
#upper["LateralHydrConductivity"] <- 0.1
#upper["NormalHydrConductivity"] <- 0.1
#upper["ThetaSat"] <- 0.55
#
#lower <- x/3
#lower["N"] <- 1.1
#lower["LateralHydrConductivity"] <- 0
#lower["NormalHydrConductivity"] <- 0
#lower["ThetaSat"] <- 0.35
###### 
#param.template.file <- ""
#param.template <- read.table("")
#####
#####
geotop.soil.param.file <- system.file("examples/param/param.csv",package="geotopOptim")
geotop.soil.param <- read.table(geotop.soil.param.file,header=TRUE,sep=",",stringsAsFactors=FALSE)
lower <- geotop.soil.param$lower
upper <- geotop.soil.param$upper
names(lower) <- geotop.soil.param$name
names(upper) <- geotop.soil.param$name


geotop.model <- list(bin=bin,simpath=simpath,runpath=runpath,
		clean=TRUE,variable=vars,data.frame=TRUE,level=1,zformatter=zformatter,intern=TRUE,names_par=names(upper))
control <- list(maxit=10,npart=6) ## instead of 10  Maximim 20 iterations!!
#control <-   list(maxit=2,npart=1) #list(maxit=5,npart=3) ## instead of 10  Maximim 20 iterations!!
######


#####
dir.create(rundir)
layer=c("z0005","z0020","z0050")
#pso <- geotopPSO(obs=obs_SWC,geotop.model=geotop.model,layer=c("z0020"),gof.mes="KGE",lower=lower,upper=upper,control=control)
pso <- geotopPSO(obs=obs_SWC,geotop.model=geotop.model,layer=layer,gof.mes="KGE",lower=lower,upper=upper,control=control)
#####





sim_SWC <- pso$gof$sim
#sim_SWC <- get.geotop.inpts.keyword.value(vars,wpath=rundir,data.frame=TRUE,date_field="Date12.DDMMYYYYhhmm.",zlayer.formatter="z%04d")
		

time_sim <- index(obs_SWC[[layer[1]]])
###[1] "2009-10-01 23:00:00 GMT+1" "2009-10-02 00:00:00 GMT+1"
time_obs <- index(sim_SWC[,layer[1]])
#### "2012-10-04 00:00:00 GMT+1" "2012-10-04 01:00:00 GMT+1"
time <- intersect(time_sim,time_obs)

sim_SWC_t <- sim_SWC[time,]
obs_SWC_t <- lapply(X=obs_SWC,FUN=function(x,time){x[time,]},time=time)
df <- NULL
for (itl in layer) {
		dfl <- data.frame(time=time,sim=sim_SWC_t[,itl],obs=obs_SWC_t[[itl]][,"mean"],
				max=obs_SWC_t[[itl]][,"max"],min=obs_SWC_t[[itl]][,"min"])
		dfl$layer <- itl
		df <- rbind(df,dfl)	
	

}
save(pso,sim_SWC_t,obs_SWC_t,upper=upper,lower=lower,file="pso.rda")

g <- ggplot(df,aes(x=time,y=obs))+geom_line()+geom_line(aes(y=sim,color=3))+geom_ribbon(aes(ymax=max,ymin=min),alpha=0.4)+xlab("Time")+ylab("SWC")
g <- g+facet_grid(layer ~ ., scale = "fixed")
#hsc <- ggplot(df,aes(x=sim,y=obs))+geom_point()+geom_ribbon(aes(ymax=max,ymin=min),alpha=0.4)+xlab("Simulated")+ylab("Observed")

png("plot_SWC.png")
print(g)
dev.off()

######



