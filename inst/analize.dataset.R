#'
#' 

rm(list=ls())

library(geotopOptim)
library(fastICA)

data(MuntatschiniB2_unproc)

## x <- apply(X=MuntatschiniB2[[1]],MARGIN=1,FUN=mean)
MuntatschiniB2_  <- lapply(X=MuntatschiniB2,FUN=function(x) {
			
		
			out <- x
			print(dim(x))
			
			
			if (is.null(dim(x))) {
				
				str(x)
				out  <- as.data.frame(x)
				x0 <- as.data.frame(x)
				
				out$mean <- apply(X=x0,MARGIN=1,FUN=mean,na.rm=TRUE)
				out$sd <-  apply(X=x0,MARGIN=1,FUN=sd,na.rm=TRUE)
				out$len <- apply(X=x0,MARGIN=1,FUN=function(y){length(which(!is.na(y)))})
				out$range <- apply(X=x0,MARGIN=1,FUN=function(y,...){max(y,...)-min(y,...)},na.rm=TRUE)
				out <- as.zoo(out)
				index(out) <- index(x)
				
			} else {
				
				out <- x
				out$mean <-  apply(X=x,MARGIN=1,FUN=mean,na.rm=TRUE)
				out$sd <-  apply(X=x,MARGIN=1,FUN=sd,na.rm=TRUE)
				out$len <- apply(X=x,MARGIN=1,FUN=function(y){length(which(!is.na(y)))})
				out$range <- apply(X=x,MARGIN=1,FUN=function(y,...){max(y,...)-min(y,...)},na.rm=TRUE)
				
			}
			
			return(out)
		}) 


print(names(MuntatschiniB2))
datapath <- "/Users/ecor/ownCloud/activity/R-Packages/geotopOptim/data"

save(MuntatschiniB2,file=paste(datapath,"MuntatschiniB2.rda",sep="/"))

