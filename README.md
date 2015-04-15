
# geotopOptimPSO Calibration/Optimization Wrapper R Package for the hydrological model GEOtop.

## Installation note on GEOtop 

GEOtop (Endrizzi et al, 2014 and references therein) is a distributed model of the mass and energy balance of the hydrological cycle. GEOtop is applicable to simulations in continuum in small or relatively large montain catchments. GEOtop deals with the effects of topography on the interaction between energy balance and hydrological cycle (water, glacier and snow) with peculiar solutions. The source code of GEOtop 2.0 with detailed documentation is available for researchers and environmental/software engineers through the following links:

* https://code.google.com/p/geotop/
* https://github.com/se27xx/GEOtop/
* https://github.com/skyglobe/geotop

For  Unix-like OS users, a C source code can be rapidly downloaded and built with the following instrunctions (https://github.com/se27xx/GEOtop/):


1. Open a console and go to the drectory where to clone GEOtop source code;
2. Clone the source code typing: "git clone https://github.com/se27xx/GEOtop";
3. Enter GEOtop directory typing: "cd GEOtop";
4. Create the subdirectory for the executable file  "mkdir bin";
5. Build and create GEOtop executale file: "make -f geotop.make"

GEOtop executable will be created in the subdirectory "bin". 
For other versions of GEOtop, please see the related URLs. 

## GEOtop calibration with "geotopOptim"



Rscript File 'example.geotop.pso.R'
HYDROPSO ALGORITM APPLIED TO GEOTOP HYDROLOGICAL MODEL

Options: 

-wpath_out            directory for output files
-geotopbin            GEOtop executable/binary file (full path)
-wpath_simpath        path to directory containing GEOtop simulation  (i.e. 'geotop.inpts' file)
					  See for example: 'system.file("Muntatschini_pnt_1_225_B2_004",package="geotopOptim") from an R console.

-wpath_runpath        directory where to run GEOtop (optional). By default is the same directory given by '-wpath-out'
-optim-soil-param     full name of the CSV file containing ranges of soil calibration parameter. 
					  See for example: 'system.file("examples/param/param.csv",package="geotopOptim")' from an R console.
				   
--help                help with  'example.geotop.pso.R' options and flags

Usage: ./example.geotop.pso.R  -wpath_out $GM_WPATH_OUT  -optim-soil-param $GM_OPTIM_PARAM_CSV_FILE -geotopbin $GM_GEOTOP_BIN -wpath_simpath $GM_GEOTOP_DATA 
Note: This scripts calibrates soil parameters for soil moisture modeling with GEOtop hydrological model. THe calibration algorithm is based on Particle Swarm Optimization. 
Actually this script works for a 1D (point/local scale) usage of the Hydrological model GEOtop. 
.... 


References: 

* Endrizzi, S., Gruber, S., Dall'Amico, M., and Rigon, R.: GEOtop 2.0: simulating the combined energy and water balance at and below the land surface accounting for soil freezing, snow cover and terrain effects, Geosci. Model Dev., 7, 2831-2857, doi:10.5194/gmd-7-2831-2014, 2014, http://www.geosci-model-dev.net/7/2831/2014/gmd-7-2831-2014.html

* Zambrano-Bigiarini, M.; R. Rojas (2013), A model-independent Particle Swarm Optimisation software for model
 calibration, Environmental Modelling & Software, 43, 5-25, doi:10.1016/j.envsoft.2013.01.004

*  Zambrano-Bigiarini, M., Rojas, R.(2014). hydroPSO: Particle Swarm Optimisation, with focus on Environmental Models. R
  package version 0.3-4.



