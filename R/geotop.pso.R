NULL
#' GEOtoop calibration through Particle Swam Optimization 
#'
#' @param par model parameters. See \code{\link{hydroPSO}},\code{\link{geotopGOF}} and \code{\link{geotopExec}}
#' @param fn function to optimize (minimize or maximize). Default is \code{\link{geotopGOF}}. See \code{\link{hydroPSO}}. 
#' @param gof.mes string(s) containing adopted numerical goodness-of-fit measure. If it is \code{NULL} (Default), all mesasures returned by \code{\link{gof}} are calculated.
#' @param gof.expected.value.for.optim expected value for Goodness-of-fit mesure, e.g. 0 or 1. It is used if this function is called by \code{link{geotopPSO}},\code{link{hydroPSO}} or \code{link{optim}}.
#' 
#' @param ... further arguments for \code{\link{hydroPSO}}.
#' 
#' @details The function \code{fn}, in case it is different from the default value \code{\link{geotopGOF}} , must always have the arguments \code{gof.mes} and \code{gof.expected.value.for.optim}.
#' 
#' 
#'@export
#'
#' 
#' @examples 
#' 
#' #' 
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

#' ### Use geotopGOF with an internal GEOtop simulation
#' 
#' ## create a list with arguments for geotopGOF
#' 
#' x <- param <- c(N=1.4,Alpha=0.0021,ThetaRes=0.05,LateralHydrConductivity=0.021,NormalHydrConductivity=0.021) 
#' upper <- x*3
#' 
#' upper["LateralHydrConductivity"] <- 0.1
#' upper["NormalHydrConductivity"] <- 0.1
#' 
#' lower <- x/3
#' lower["N"] <- 1.1
#' lower["LateralHydrConductivity"] <- 0
#' lower["NormalHydrConductivity"] <- 0
#' 
#' 
#' geotop.model <- list(bin=bin,simpath=simpath,runpath=runpath,
#' clean=TRUE,variable=vars,data.frame=TRUE,level=1,zformatter=zformatter,intern=TRUE)
#' control <- list(maxit=10) ## Maximim 20 iterations!! 
#' 
#' pso <- geotopPSO(par=x,obs=obs_SWC,geotop.model=geotop.model,layer=c("z0020"),gof.mes="KGE",lower=lower,upper=upper,control=control)
#' 
#' 
#' 
#' 
#' @seealso \code{\link{hydroPSO}},\code{\link{gof}}
#'

geotopPSO <- function(fn=geotopGOF,gof.expected.value.for.optim=NA,gof.mes="KGE",...) {

###		if (is.charecter(fn)) fn <- get(fn)
	   	if (is.null(gof.expected.value.for.optim))	gof.expected.value.for.optim <- NA
		if (is.na(gof.expected.value.for.optim)) {
			
			x <- 1:100
			gof.expected.value.for.optim <- gof(x,x)[gof.mes,1][1]
			
		}
		print(gof.expected.value.for.optim)
		
		out <- hydroPSO(fn=fn,gof.mes=gof.mes,gof.expected.value.for.optim=gof.expected.value.for.optim,...)
		print("out:")
		print(out)
		out$gof <- fn(x=out$par,gof.mes=gof.mes,gof.expected.value.for.optim=gof.expected.value.for.optim,...)
		
		
		return(out)
		



}