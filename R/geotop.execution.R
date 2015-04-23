NULL 
#' Execution of GEOtop Hydrological Model
#' 
#' @param param vector of parameters to set for GEOtop simulation. Actually implemented only for soil type parameters. Deafault is \code{NULL}.
#' @param bin binary executable file of GEOtop with full path
#' @param inpts.file GEOtop configuretion \code{*.inpts} file. Default is \code{geotop.inpts}. 
#' @param simpath directory containg the \code{inpts.file}.
#' @param runpath directory where to run GEOtop
#' @param intern logical value, see \code{\link{system}}.
#' @param clean logical value. If it is \code{TRUE} previous simulations and other stuff in \code{runpath} were removed. See \code{\link{file.remove}} for functionality.
#' @param recovery logical value. If it is \code{TRUE}, GEOtop simulation are recovered or overwriiten according to the GEOtop settings in \code{inpts.file}. 
#' @param getKeywords character vector containing the keywords of \code{inpts.file} which can be imported after the GEOtop run. Default is \code{NULL}, nothing is imported.
#' @param data.frame logical vaue, see  \code{\link{get.geotop.inpts.keyword.value}}.  Default is \code{TRUE}, is enabled if \code{getKeywords} is not \code{NULL}.
#' @param date_field character vaue, see  \code{\link{get.geotop.inpts.keyword.value}}. Default is \code{"Date12.DDMMYYYYhhmm."}, is enabled if \code{getKeywords} is not \code{NULL}.
#' @param formatter character decimal formatter.  See \code{\link{get.geotop.inpts.keyword.value}}. Default is \code{"\%04d"}.
#' @param param.soil logical value If it \code{TRUE} (default) the vaues in \code{param} concerns soil parameters. 
#' @param paramPrefix character string. Default is \code{"Header"}. If \code{param.soil==TRUE}, the soil parameters in \code{param} are named with the corrosponding parameter keywords in \code{inpts.file} to which the  string \code{paramPrefix} is attached as a prefix. If \code{paramPrefix} is \code{NA} or \code{NULL} , the names of \code{param} elements must be matched with the corresponimg paramater name in the soil parameter data frame. 
#' @param ... further arguments for \code{\link{get.geotop.inpts.keyword.value}}
#' 
#' 
#' @details In this function, implementation, The parameters entered through \code{param} replace  only the ones of  the first (\code{"0001"}) soil type  declared in the GEOtop simulation directory.
#' 
#'  @export
#' 
#' @importFrom stringr str_split
#' @importFrom geotopbricks get.geotop.inpts.keyword.value
#' @importFrom soilwater swc
#' @examples
#' 
#' simpath <- system.file("Muntatschini_pnt_1_225_B2_004",package="geotopOptim")
#' bin <-   "/Users/ecor/local/bin/geotop_zh"
#' runpath <- "/Users/ecor/ownCloud/job"
#' 
#' vars <- c("SoilAveragedTempProfileFile",	"SoilLiqWaterPressProfileFile",
#' "SoilLiqContentProfileFile","SoilIceContentProfileFile")
#' 
#' 
#' out <- geotopExec(bin=bin,simpath=simpath,runpath=runpath,
#'       clean=TRUE,getKeywords=vars,data.frame=TRUE,level=1)
#' 
#' param <- c(N=1.4,Alpha=0.0021,ThetaRes=0.05) 
#' out1 <- geotopExec(param=param,bin=bin,simpath=simpath,
#' 			runpath=runpath,clean=TRUE,getKeywords=vars,
#' 			data.frame=TRUE,level=1,intern=TRUE)
#' 
#' 
#' 

#
#PointOutputFile	=	"tabs/point"
#PointAll	=	1
#
#!SnowProfileFile	=	"tabs/snow"
#SnowDepthLayersFile	="tabs/snowDepth"
#SnowTempProfileFile	=	"tabs/snowT"
#SnowLiqContentProfileFile	=	"tabs/snowLiq"
#SnowIceContentProfileFile	=	"tabs/snowIce"
#SnowAll	=	1
#
#BasinOutputFile	=	"tabs/basin"
#BasinAll	=	1
#
#SoilAveragedTempProfileFile		=	"tabs/soilTz"
#SoilLiqWaterPressProfileFile	=	"tabs/psiz"
#SoilLiqContentProfileFile		=	"tabs/thetaliq"
#SoilIceContentProfileFile		=	"tabs/thetaice"
#SoilAll=1

#
#
#
#
#





geotopExec <- function (param=NULL,bin="/Users/ecor/local/bin/geotop_zh",simpath,inpts.file="geotop.inpts",
		runpath="/Users/ecor/ownCloud/job",clean=TRUE,recovery=!clean,getKeywords=NULL,
		data.frame=TRUE,date_field="Date12.DDMMYYYYhhmm.",intern=FALSE,param.soil=TRUE,formatter = "%04d",paramPrefix="Header",names_par=NULL,...) {
	
	print(simpath)
	t <- str_split(simpath,"/")[[1]]
	simdir <- t[length(t)]
	
	rundir <- paste(runpath,simdir,sep="/")
	
	if (recovery==TRUE) clean <- FALSE
	
	
	if (clean==TRUE) {
		
		torv <- list.files(rundir,full.names=TRUE,include.dirs=TRUE,recursive=TRUE)
		file.remove(torv)
		torv <- list.files(rundir,full.names=TRUE,include.dirs=TRUE,recursive=TRUE)
		file.remove(torv)
		
	}
	
	file.copy(from=simpath,to=runpath,recursive=TRUE,overwrite=TRUE)
	
	## MOdify input file according to param:
	
	
	print(param)
	print(names(param))
	print(param.soil)
	print(names_par)
	if (!is.null(param)) {
		
		if (is.null(names(param))) {
			
			names(param) <- names_par
			
		}
		
		if (is.null(names(param))) {
			
			warning("param has NO NAMES and will be ignored!")
			
		} else if (param.soil==TRUE) {
			
			
#			HeaderLateralHydrConductivity	=	"Kh"
#			HeaderNormalHydrConductivity	=	"Kv"
#			HeaderThetaRes	=	"vwc_r"
#			HeaderWiltingPoint	=	"vwc_w"
#			HeaderFieldCapacity	=	"vwc_fc"
#			HeaderThetaSat	=	"vwc_s"
#			HeaderAlpha	=	"alpha"
#			HeaderN	=	"n"
#			HeaderSpecificStorativity	=	"stor"
#			HeaderKthSoilSolids	=	"Kth"
#			HeaderCthSoilSolids	=	"Cth"
			
			## 
			variable.soil.depth <- FALSE
			if (all(c("SoilDepth","NumberOfSoilLayers") %in% names(param))) {
				
				SoilDepth <- as.numeric(param["SoilDepth"])
				SurfaceSoilLayer <- as.numeric(param["SurfaceSoilLayer"])
				
				NumberOfSoilLayers  <- ceiling(as.numeric(param["NumberOfSoilLayers"]))
				NumberOfSoilLayers[NumberOfSoilLayers<4] <- 4
				
				
				variable.soil.depth <- TRUE
				param <- param[!(names(param) %in% c("SoilDepth","NumberOfSoilLayers","SurfaceSoilLayer"))]
			}
		
		
			if (c("PsiGamma") %in% names(param)) {
			
				psiGamma <- as.numeric(param["PsiGamma"])	
				
				
				
			}	else {
				
				psiGamma <- NA
				
			}
			
			######
			param <- param[!(names(param) %in% c("PsiGamma"))]
			### INSERT BOTTOM LAYER!!
			
			param_bottomlayer <- param[str_detect(names(param),"_bottomlayer")]
			param <- param[!(param %in% param_bottomlayer)]
			names(param_bottomlayer) <- str_replace_all(names(param_bottomlayer),"_bottomlayer","")
			
			
		    if (c("SoilInitPresL0001") %in% names(param))	{
				
				parSoilInitPresL <- which(str_detect(names(param),"SoilInitPresL"))
				param_soil <- param[parSoilInitPresL]
				param <- param[-parSoilInitPresL]
				len <- length(param_soil)
				
				names_param_soil <- sprintf("SoilInitPresL%04d",1:len)
				psisoil <- param_soil[names_param_soil]
				
				
				
			}	else {
				
				psisoil <- NULL
			}	
			
			
			
			
			
			param.soil.df.filename <- 	get.geotop.inpts.keyword.value("SoilParFile",wpath=rundir,inpts.file=inpts.file,add_wpath=TRUE)
			param.soil.df.filename <- paste(param.soil.df.filename,formatter,".txt",sep="")
			layer <- 1 
			param.soil.df.filename <- sprintf(param.soil.df.filename,layer)
		###	param.soil.df.filename <<- param.soil.df.filename
			param.soil.df <- read.table(param.soil.df.filename,header=TRUE,sep=",")
			
			
			if (is.null(paramPrefix)) ParamPrefix <- NA
			if (!is.na(paramPrefix)) {
				ids <- paste(paramPrefix,names(param),sep="")
				names(param) <- get.geotop.inpts.keyword.value(ids,wpath=rundir,inpts.file=inpts.file)
				### check Initial Condition
				if (length(param_bottomlayer)>0) {
					
					ids_b <- paste(paramPrefix,names(param_bottomlayer),sep="")
					names(param_bottomlayer) <- get.geotop.inpts.keyword.value(ids_b,wpath=rundir,inpts.file=inpts.file)
					
				}
				WiltingPoint <- get.geotop.inpts.keyword.value(paste(paramPrefix,"WiltingPoint",sep=""),wpath=rundir,inpts.file=inpts.file)
				FieldCapacity <- get.geotop.inpts.keyword.value(paste(paramPrefix,"FieldCapacity",sep=""),wpath=rundir,inpts.file=inpts.file)
				ThetaSat <- get.geotop.inpts.keyword.value(paste(paramPrefix,"ThetaSat",sep=""),wpath=rundir,inpts.file=inpts.file)
				ThetaRes <- get.geotop.inpts.keyword.value(paste(paramPrefix,"ThetaRes",sep=""),wpath=rundir,inpts.file=inpts.file)
				VG_Alpha <- get.geotop.inpts.keyword.value(paste(paramPrefix,"Alpha",sep=""),wpath=rundir,inpts.file=inpts.file)
				VG_N <- get.geotop.inpts.keyword.value(paste(paramPrefix,"N",sep=""),wpath=rundir,inpts.file=inpts.file)
				Dz <- get.geotop.inpts.keyword.value(paste(paramPrefix,"SoilDz",sep=""),wpath=rundir,inpts.file=inpts.file)
				SoilInitPres <- get.geotop.inpts.keyword.value(paste(paramPrefix,"SoilInitPres",sep=""),wpath=rundir,inpts.file=inpts.file)
				print(Dz)
				
			}
			
			if (variable.soil.depth==TRUE) {
				
				
				
				if(is.na(SurfaceSoilLayer)) SurfaceSoilLayer <- param.soil.df[1,Dz]
				
				print(SoilDepth)
				print(SurfaceSoilLayer)
				param.soil.df <- param.soil.df[1:NumberOfSoilLayers,]
				
				param.soil.df[1:NumberOfSoilLayers,] <- param.soil.df[1,]
				
				polycoeff <- array(1,NumberOfSoilLayers)
				polycoeff[1] <- 1-SoilDepth/SurfaceSoilLayer
				
				lambda <- polyroot(polycoeff)
				lambda <- lambda[Re(lambda)>=0]
				
				ail <- abs(Im(lambda))
				
				print(lambda)
				lambda <- Re(lambda[which.min(ail)])
				lambda <- lambda[lambda>=0][1]    ## Get positive or null solution!!
				param.soil.df[,Dz] <- SurfaceSoilLayer*lambda^(0:(NumberOfSoilLayers-1))
				
				### TO TEST!!!
				
				
				
			}
			
			
			
			
			
			
			### Adjust bottom layer 
			
			dz <- param.soil.df[,Dz]
			z <- dz/2
			for(i in 2:length(z)) {
				z[i] <- z[i-1]+dz[i]/2+dz[i-1]/2
				
			}
			
			zm <- (z-z[1])/(z[length(z)]-z[1])
			
			for (it in names(param)) {
				
				param.soil.df[,it] <- param[[it]]
				
				
				if (it %in% names(param_bottomlayer)) {
					
					#### DO INTERPOLATION 
					
					param.soil.df[,it] <- param[[it]]^(1-zm)*param_bottomlayer[[it]]^zm
					   
				}	
				
			}
		
			
			
			if (!is.na(psiGamma)) {
				print(Dz)
				print(param.soil.df)
				
				
				
				
				if (SoilInitPres %in% names(param)) {
					
					param.soil.df[,SoilInitPres] <- param[[SoilInitPres]]+psiGamma*z
				}
				
				if (!is.null(psisoil)) {
					
					param.soil.df[,SoilInitPres][1:length(psisoil)] <- psisoil
					
					zfront <- z[length(psisoil)]
					
					param.soil.df[,SoilInitPres][z>zfront] <- psisoil[length(psisoil)]+psiGamma*(z[z>zfront]-zfront)
			
					
				}	
				
				
			}
			
			
			### CORRECT SOIL WATER
			alpha <-  param.soil.df[,VG_Alpha]
			n <- param.soil.df[,VG_N]
			theta_sat <- param.soil.df[,ThetaSat]
			theta_res <- param.soil.df[,ThetaRes]
		###	=alpha,n=n,theta_sat=theta_sat,theta_res=theta_res)
			waterdensity <- 1000 ## kg/m^3
			gravity <- 9.81 ## m/s^2
			
			
			psi_WP <- -1500*1000 ## Pa
			psi_FC <- -33*1000  ##  Pa
			
			psi_WP <- psi_WP/(waterdensity*gravity)*1000 ## converted to water millimiters according to GEOtop
			psi_FC <- psi_FC/(waterdensity*gravity)*1000 ## converted to water millimiters according to GEOtop
			
			param.soil.df[,FieldCapacity] <- swc(psi=psi_FC,alpha=alpha,n=n,theta_sat=theta_sat,theta_res=theta_res,type_swc="VanGenuchten")
			param.soil.df[,WiltingPoint] <- swc(psi=psi_WP,alpha=alpha,n=n,theta_sat=theta_sat,theta_res=theta_res,type_swc="VanGenuchten")
			###
			print(param.soil.df)
			write.table(x=param.soil.df,file=param.soil.df.filename,sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)
		#####	param.soil.df <<- param.soil.df
			
					###data.frame=data.frame,date_field=date_field,...)
			
			
			
			
			
			
		} else {
			
			warning("param.soil is not TRUE: no other methods are implemented, param will be ignored")
			param <- NULL
		}
		
		
		
	}
	
	command.line <- paste(bin,rundir,sep=" ")
	cc <- system(command.line,intern=intern) 
	print(cc)
	

	if (length(getKeywords)>0) {
		
	####	date_field0 <<- date_field
		out <- lapply(X=getKeywords,FUN=get.geotop.inpts.keyword.value,wpath=rundir,inpts.file=inpts.file,data.frame=data.frame,date_field=date_field,formatter=formatter,...)
		names(out) <- getKeywords
		
	} else {
		
		out <- rundir
		
	}
	
	str(out)
	print(names(out[[1]]))
	
	
	
	
	return(out)
	
}