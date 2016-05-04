# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL 

#' GEOtop Goodnes of fit 
#' 
#' @param x vector with parameters to calibrate. See \code{param} in  \code{\link{geotopZProfile}} or \code{\link{geotopExec}}
#' @param geotop.model a list with arguments for \code{\link{geotopZProfile}}. It is used if \code{sim} is \code{NULL}.
#' @param approx.list list of arguments for \code{\link{approxfunDataFrame}} (optional) or \code{\link{integratefunDataFrame}} (in this case rember to specify the \code{zrange} argument).
#' @param sim simulated data as a  object returned by \code{\link{geotopZProfile}}
#' @param obs observed data 
#' @param layer layers corresponding to soil depth at whch GOF indices are calculated
#' @param weights vector of weights to assing to each layer, in case of a weigeted-averaged goodness-of-fit measure over all layers. It is \code{NULL} (Default) the gooness-of-fit measures are separately calculated for each layer. If it is \code{"uniform"}, the weights are uniformly distributed in all layers. 
#' @param obs_field obs field used in the observation data frame. Deafault is \code{"mean"}, it is used in case varaiables are measured with different sensors at the same depth and location. 
#' @param gof.mes string(s) containing adopted numerical goodness-of-fit measure. If it is \code{NULL} (Default), all mesasures returned by \code{\link{gof}} are calculated.
#' @param gof.expected.value.for.optim expected value for goodness-of-fit mesure, e.g. 0 or 1. It is used if this function is called by \code{link{geotopPSO}},\code{link{hydroPSO}} or \code{link{optim}}.
#' @param output_simulation logical value. If it is \code{TRUE}, function returns a list with the GOF values and the simulated time series
#' @param names_par names of \code{x}
#' @param temporary.runpath logical value. If it \code{TRUE} GEOtop is running in a temporary sub-directory, see \code{\link{geotopExec}}.  Default is \code{FALSE}.
#' @param useSoilIntegratedValues logical values. Default is \code{FALSE}. If it is \code{TRUE} output is integrated with a soil thckness through \code{\link{integratefunDataFrame}}.
#' @param ... further aguments for \code{\link{gof}} .
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
#' bin <-   "/home/ecor/local/geotop/GEOtop/bin/geotop-2.0.0"
#' runpath <- "/home/ecor/temp/geotopOptim_tests"
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
#' gof_geotop_ <- geotopGOF(x=x,obs=obs_SWC,geotop.model=geotop.model,layer=c("z0005","z0020"),gof.mes="KGE",useSoilIntegratedValues=TRUE,approx.list=list(zrange=c(0,500),rescaleWithDz=TRUE))
#' 
# gof_geotop <- geotopGOF(x=x,obs=obs_SWC,geotop.model=geotop.model,layer=c("z0005","z0020"),gof.mes="KGE",useSoilIntegratedValues=FALSE,approx.list=list(zrange=c(0,500))
# ## PLAY WITH THE PLOTS!!!!
# 
# ###
# ###
# ###
# #modeled <- sim$SoilLiqContentProfileFile$z0020
# #observed <- obs_SWC$z0020
# #m <- merge(observed,modeled)
# #m <- m[!is.na(m$modeled),]
# ## Plotting standard deviation of observation vs observed range (max-min)
# 
# #plot(m$max-m$min,m$modeled-m$min)
# #abline(0,1)
# #points(m$mean-m$min,m$modeled-m$min,pch=2)
# #
# # 
# #plot(m$sd,m$mean-m$modeled)
# #
# #
# #
# ## plot(
#' #

geotopGOF <- function(x=NULL,geotop.model=NULL,approx.list=list(),sim=NULL,obs,layer=c("z0005","z0020"),obs_field="mean",gof.mes=NULL,gof.expected.value.for.optim=NULL,weights=NULL,output_simulation=FALSE,names_par=NULL,useSoilIntegratedValues=FALSE,temporary.runpath=FALSE,...) {
	
	##print("x:")
	##print(x)
	
	if (!is.null(geotop.model)) {
		
		variables <- geotop.model[["variable"]]
		
		if (length(variables)>1) {
			
			variables <- variables[1]
			warning("No more than one variable name! Only first variable is considered in this analysis!!!")
			geotop.model[["variable"]] <- variables[1]
		}
		
		geotop.model[["param"]] <- x
		
		if (!is.null(names_par)) {
			
			geotop.model[["names_par"]] <- names_par
			
		}
		
		geotop.model[["temporary.runpath"]] <- temporary.runpath
		
		sim <- do.call(what=geotopZProfile,args=geotop.model)
		sim <- sim[[geotop.model[["variable"]]]]
		
		
		
		
	}
	
	###
	
	
	
	###
	##layer0 <- layer
	
	if (useSoilIntegratedValues==TRUE) {
		
		### TO DO 
		##	integratefunDataFrame <- function(df,z=NULL,zrange=c(0,500),formatter="z%04d",factor=10,rescaleWithDz=FALSE,...) 
		approx.list[["df"]] <- sim
		##approx.list[["zrange"]] <- zrange
		
		sim <- do.call(what=integratefunDataFrame,args=approx.list)
		if (length(layer)==0) {
			
			layer <- "value"
		}
		layer <- layer[1]
		
		if (class(obs)!="list") {
			
			
			obs <- list(obs)
			names(obs) <- layer
			
		}
		
		
	} else {
		
		layer <- intersect(layer,names(obs))
		
		approx.list[["df"]] <- sim
		approx.list[["zout"]] <- layer
		
		sim <- do.call(what=approxfunDataFrame,args=approx.list)
		
		if (length(layer)==0) {
			
			stop("No layers set!")
			
		}
		
		
	} 
	## SEE 
	
	
	
	
	
	
	
	
	out <- NULL
	
	for (i in 1:length(layer)) {
		it <- layer[i]	
		modeled <- sim[,it]
		
		str(obs)
		str(sim)
		
		m <- merge(obs[[it]],modeled)
		m <- m[!is.na(m$modeled),]
		
		str(m)
		print(m[index(m)[1:10],])	
		
		## da provare 
		
		if (nrow(m)>2) {
			
			val <- gof(obs=m[,obs_field],sim=m$modeled,...)
			
		} else {
			
			val <- gof(1:10,1:10,...)
			val [,] <- NA 
			
		}	
		
		if (i==1) {
			
			out <- array(NA,c(length(val),length(layer)))
			rownames(out) <- rownames(val)
			
		}
		out[,i] <- val
		
	}
	
	colnames(out) <- layer
	outncol <- ncol(out)
	if (!is.null(gof.mes)) {
		
		out <- out[gof.mes,]
		
	}
	
	## INSERT HERE MARGING
	#print(out)
	
	
	
	
	if (!is.null(gof.expected.value.for.optim)) {
		
		if (!is.na(gof.expected.value.for.optim)) {
			
			#	if (length(out)>1) {
			
			
			
			#	} else {
			
			#	out <- as.numeric(out[1])
			#	print(out) weight
			out <- abs(out-gof.expected.value.for.optim) ## heck this passage if the optimal is 1 or  0, we consider te minimizad distance Symmetrically!!!
			#
			##	stop()
			#	}
			
		}
		
		
		
		
	}
	
	####print(weights)
	if (!is.null(weights)) {
		str(out)
		
		
		if (weights=="uniform") weights <- array(1/outncol,outncol)
		
		if (length(weights)==outncol) {
			
			weights <- array(weights,c(length(weights),1))
			out <- out %*% weights
			##	print("out:::")
			###	print(out)
			
		}
		
	}
	print(output_simulation)
	str(sim)
	if (output_simulation==TRUE) {
		
		l_out <- list()
		l_out$gof <- out
		l_out$sim <- sim
		
		out <- l_out 
		
		
		
	}
	print(out)
	return(out)
	
	
}