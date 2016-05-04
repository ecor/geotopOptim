# TODO: Add comment
# 
# Author: ecor
###############################################################################


source('~/Dropbox/R-packages/geotopOptim/R/geotop.execution.R', chdir = TRUE)


simpath <- system.file("Muntatschini_pnt_1_225_B2_004",package="geotopOptim")
bin <-   "/home/ecor/local/geotop/GEOtop/bin/geotop-2.0.0"
runpath <- "/home/ecor/temp/geotopOptim_tests"

vars <- c("AvailableSoilWaterContent","SoilLiqContentProfileFile")

param <- c(N=1.4,Alpha=0.0021,ThetaRes=0.05)
out1 <- geotopExec(param=param,bin=bin,simpath=simpath,
		runpath=runpath,clean=TRUE,getKeywords=vars,
		data.frame=TRUE,level=1,intern=TRUE)

