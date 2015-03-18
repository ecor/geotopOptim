# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL 

#' GEOtop Goodnes of fit 
#' 
#' @param x vector with parameters to calibrate. See \code{param} in  \code{\link{geotopZProfile}} or \code{\link{geotopExec}}
#' @param geotop.model a list with arguments for \code{\link{geotopZProfile}}. It is used if \code{sim} is \code{NULL}.
#' @param sim simulated data as a  object returned by \code{\link{geotopZProfile}}
#' @param obs observed data 
#' @param layer layers corresponding to soil depth at whch GOF indices are calculated
#' @param obs_field obs field used in the observation data frame. Deafault is \code{"mean"}, it is used in case varaiables are measured with different sensors at the same depth and location. 
#' @param gof.mes string(s) containing adopted numerical goodness-of-fit measure. If it is \code{NULL} (Default), all mesasures returned by \code{\link{gof}} are calculated.
#' @param gof.expected.value.for.optim expected value for goodness-of-fit mesure, e.g. 0 or 1. It is used if this function is called by \code{link{geotopPSO}},\code{link{hydroPSO}} or \code{link{optim}}.
#' @param ... further aguments for \code{\link{gof}}
#' 
#' @export
#' @seealso \code{\link{geotopZProfile}},\code{\link{gof}}
#' 
#' @importFrom hydroGOF gof
#' @examples 
#' 
#' data(MuntatschiniB2)
#' ## OBSERVATION PROCESSING
#' 
#' obs_SWC <- MuntatschiniB2[str_detect(names(MuntatschiniB2),"SWC")]
#' zvalues <-  as.numeric(unlist(lapply(X=str_split(names(obs_SWC), pattern="", 
#' 				n = Inf),FUN=function (x) {
#' 					out <- as.numeric(x)
#' 				    out <- out[!is.na(out)]
#' 					out <- paste(out,collapse="")
#' 				return(out)
#' })))
#' zformatter = "z%04d"
#' names(obs_SWC) <- sprintf(zformatter,zvalues)
#' obs_SWC <- lapply(X=obs_SWC,FUN=function(x){
#' 
#' 				if (length(dim(x))>1) {
#' 					max <- apply(X=x,MARGIN=1,FUN=max,na.rm=TRUE)
#' 					min <- apply(X=x,MARGIN=1,FUN=min,na.rm=TRUE)
#' 					mean <- apply(X=x,MARGIN=1,FUN=mean,na.rm=TRUE)
#' 					sd <- apply(X=x,MARGIN=1,FUN=sd,na.rm=TRUE)
#' 				} else {
#' 				
#' 					mean <- as.vector(x)
#' 				 	max <-  as.vector(x)
#' 					min <-  as.vector(x)
#' 					sd <- NA
#' 
#' 				}
#' 				out <- data.frame(min=min,mean=mean,max=max,sd=sd)
#' 
#' 				out <- as.zoo(out)
#'              index(out) <- as.POSIXlt(index(x))
#' 
#' 				return(out)
#' 
#' })
#' ###########
#' ###########
#' 
#' 
#' simpath <- system.file("Muntatschini_pnt_1_225_B2_004",package="geotopOptim")
#' bin <-   "/Users/ecor/local/bin/geotop_zh"
#' runpath <- "/Users/ecor/ownCloud/job"
#' 
#' vars <- "SoilLiqContentProfileFile"
#' 
#' sim <- geotopZProfile(bin=bin,simpath=simpath,runpath=runpath,
#' clean=TRUE,variable=vars,data.frame=TRUE,level=1,zformatter=zformatter,intern=TRUE)
#' 
#' gof <- geotopGOF(obs=obs_SWC,sim=sim[[vars[1]]],layer=c("z0005","z0020"))
#' 
#' gof
#' 
#' ### Use geotopGOF with an internal GEOtop simulation
#' 
#' ## create a list with arguments for geotopZProfile
#' 
#' x <- param <- c(N=1.4,Alpha=0.0021,ThetaRes=0.05) 
#' 
#' geotop.model <- list(bin=bin,simpath=simpath,runpath=runpath,
#' clean=TRUE,variable=vars,data.frame=TRUE,level=1,zformatter=zformatter,intern=TRUE)
#' 
#' gof_geotop <- geotopGOF(x=x,obs=obs_SWC,geotop.model=geotop.model,layer=c("z0005","z0020"),gof.mes="KGE")
#' 
#' ## PLAY WITH THE PLOTS!!!!
#' 
#' ###
#' ###
#' ###
#' #modeled <- sim$SoilLiqContentProfileFile$z0020
#' #observed <- obs_SWC$z0020
#' #m <- merge(observed,modeled)
#' #m <- m[!is.na(m$modeled),]
#' ## Plotting standard deviation of observation vs observed range (max-min)
#' 
#' #plot(m$max-m$min,m$modeled-m$min)
#' #abline(0,1)
#' #points(m$mean-m$min,m$modeled-m$min,pch=2)
#' #
#' # 
#' #plot(m$sd,m$mean-m$modeled)
#' #
#' #
#' #
#' ## plot(
#' #

geotopGOF <- function(x=NULL,geotop.model=NULL,sim=NULL,obs,layer=c("z0005","z0020"),obs_field="mean",gof.mes=NULL,gof.expected.value.for.optim=NULL,...) {
	 
	
	if (!is.null(geotop.model)) {
		
		variables <- geotop.model[["variable"]]
		
		if (length(variables)>1) {
			
			variables <- variables[1]
			warning("No more than one variable name! Only first variable is considered in this analysis!!!")
			geotop.model[["variable"]] <- variables[1]
		}
		
		geotop.model[["param"]] <- x
		
		
		sim <- do.call(what=geotopZProfile,args=geotop.model)
		sim <- sim[[geotop.model[["variable"]]]]
		
		
		
	}

	layer <- intersect(layer,names(obs))
	layer <- intersect(layer,names(sim))
	
	if (length(layer)==0) {
		
		stop("No layers set!")
		
	}
	
	
	out <- NULL
	
	for (i in 1:length(layer)) {
	it <- layer[i]	
	modeled <- sim[,it]
	

	
	m <- merge(obs[[it]],modeled)
	m <- m[!is.na(m$modeled),]
	
	## da provare 
	
	val <- gof(obs=m[,obs_field],sim=m$modeled,...)

	if (i==1) {
		
		out <- array(NA,c(length(val),length(layer)))
		rownames(out) <- rownames(val)
		
	}
	out[,i] <- val
	
	}
	
	colnames(out) <- layer
	
	if (!is.null(gof.mes)) {
		
		out <- out[gof.mes,]
		
	}
	
	if (!is.null(gof.expected.value.for.optim)) {
		
		if (!is.na(gof.expected.value.for.optim)) {
			
			if (length(out)>1) {
				
				
				
			} else {
				
				out <- as.numeric(out[1])
			#	print(out)
				out <- abs(out-gof.expected.value.for.optim) ## heck this passage if the optimal is 1 or  0, we consider te minimizad distance Symmetrically!!!
			#
			##	stop()
			}
			
		}
		
		
		
		
	}
	return(out)
	
	
}


