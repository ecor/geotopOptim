#!/usr/bin/env Rscript



## Launch example: 
##./example.geotop.pso.R  -wpath_out $GM_WPATH_OUT  -optim-soil-param $GM_OPTIM_PARAM_CSV_FILE -geotopbin $GM_GEOTOP_BIN -wpath_simpath $GM_GEOTOP_DATA 

## /Users/ecor/Dropbox/R-packages/geotopOptim/inst/examples/example.geotop.lhoat.obs.R -wpath_out /home/ecor/temp/geotopOptim_tests/trial020lhoat_KGE_variable_layers -optim-soil-param /home/ecor/temp/geotopOptim_tests/trial14_param/param.csv 
##
##
##
##




rm(list=ls())

# LOADING PACKAGES
library(geotopOptim)
library(ggplot2)



##http://stackoverflow.com/questions/2151212/how-can-i-read-command-line-parameters-from-an-r-script
args<- commandArgs(TRUE)
### WORKING OUTPUT 
help_flag <- "--help"


### EC 20150609 
###/Users/ecor/Dropbox/R-packages/geotopOptim/inst/examples/example.geotop.lhoat.obs.R 
wpath_out <- "/home/ecor/temp/geotopOptim_tests/trial020lhoat_KGE_variable_layers" 
param_csv <- "/Users/ecor/Dropbox/R-packages/geotopOptim/inst/examples/param/param_2nd.csv" 
args <- ""
### END EC 20150609 
Sys.setenv(GM_WPATH_OUT=wpath_out,GM_OPTIM_PARAM_CSV_FILE=param_csv)




if (length(args)==0) args <- help_flag

simpath_DEF <- Sys.getenv("GM_GEOTOP_DATA")
bin_DEF     <- Sys.getenv("GM_GEOTOP_BIN")
wpath_pso_DEF <- Sys.getenv("GM_WPATH_OUT")
geotop.soil.param.file_DEF <- Sys.getenv("GM_OPTIM_PARAM_CSV_FILE")
obs_rda_DEF  <- Sys.getenv("GM_OBS_RDAFILE")
obs_ts_DEF  <- Sys.getenv("GM_OBS_TS")


if (simpath_DEF=="")   simpath_DEF <- system.file("Muntatschini_pnt_1_225_B2_004",package="geotopOptim")
if (bin_DEF=="")        bin_DEF <-   "/Users/ecor/ownCloud/geotop_se27xx/GEOtop/bin/geotop-2.0.0"
if (wpath_pso_DEF=="")   wpath_pso_DEF <- "."
if (geotop.soil.param.file_DEF=="") geotop.soil.param.file_DEF <- system.file("examples/param/param.csv",package="geotopOptim")
if (obs_rda_DEF=="") obs_rda_DEF <- "UseInternalData"
if (obs_ts_DEF=="") obs_ts_DEF <- ""


wpath_pso <- argsParser(option="-wpath_out",args=args,novalue_response = wpath_pso_DEF)

print(wpath_pso)
print(getwd())
if (class(try(setwd(wpath_pso),silent=TRUE))=="try-error") {
	
	dir.create(wpath_pso,recursive=TRUE)
	
}
setwd(wpath_pso)
print(getwd())

 ###/home/ecor/local/geotop/GEOtop/bin/geotop-2.0.0"





simpath <- argsParser(option="-wpath_simpath",args=args,novalue_response=simpath_DEF)
runpath <- argsParser(option="-wpath_runpath",args=args,novalue_response=wpath_pso)
bin <- argsParser(option="-geotopbin",args=args,novalue_response=bin_DEF)
geotop.soil.param.file <- argsParser(option="-optim-soil-param",args=args,novalue_response=geotop.soil.param.file_DEF)
obs_rda <- argsParser(option="-obs_rda",args=args,novalue_response=obs_rda_DEF)
obs_ts <- argsParser(option="-obs_ts",args=args,novalue_response=obs_ts_DEF)		
		
		
needHelp <- argsParser(option=help_flag,args=args,novalue_response=FALSE)
print(needHelp)

if (needHelp==TRUE) {
	
	
	helpco <- readLines(system.file('examples/help/option.txt',package="geotopOptim"))
    vvout <- lapply(X=helpco,FUN=message)
	
	stop("Mandatory Running Arguments missing!! ")
	
	
}


set.seed(1223)
print(obs_rda)
## START OBS 
if (obs_rda=="UseInternalData") {
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


	} else {
		
		
		load(obs_rda)
		if (obs_ts=="") {
			
			obs_ts <- rev(str_replace(str_split(file,"/")[[1]],".txt",""))[1]
			
			
			
			
		}
		
		obs_SWC <- get(obs_ts)
		
		
	}
# END OBS


t <- str_split(simpath,"/")[[1]]
simdir <- t[length(t)]

rundir <- paste(runpath,simdir,sep="/")








vars <- "SoilLiqContentProfileFile"

### Use geotopGOF with an internal GEOtop simulation


geotop.soil.param <- read.table(geotop.soil.param.file,header=TRUE,sep=",",stringsAsFactors=FALSE)
lower <- geotop.soil.param$lower
upper <- geotop.soil.param$upper
names(lower) <- geotop.soil.param$name
names(upper) <- geotop.soil.param$name


geotop.model <- list(bin=bin,simpath=simpath,runpath=runpath,
		clean=TRUE,variable=vars,data.frame=TRUE,level=1,zformatter=zformatter,intern=TRUE,names_par=NULL,temporary.runpath=TRUE) #names(upper))
control <- list(parallel="parallel",f=0.3,N=10) 
#control <-   list(maxit=2,npart=1) #list(maxit=5,npart=3) ## instead of 10  Maximim 20 iterations!!
######


#####
dir.create(rundir)
layer=c("z0005","z0020","z0050")
#pso <- geotopPSO(obs=obs_SWC,geotop.model=geotop.model,layer=c("z0020"),gof.mes="KGE",lower=lower,upper=upper,control=control)
lhoat <- geotoplhoat(obs=obs_SWC,geotop.model=geotop.model,layer=layer,gof.mes="KGE",lower=lower,upper=upper,control=control)
#####

save(lhoat,file="lhoat.rda")
