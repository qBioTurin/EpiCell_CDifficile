
#### This is with 1370 (without the Ex and Tr of Sheme) ####

#' 1) Let load the metabolic netwwork
#'
#' source("Rfunction/MatlabFileGeneration.R")
#' Matlab.file = readMATmod("Input/CD196_heme.mat")
#' save(Matlab.file, file = "Input/FBAmodel.RData")
#' source("Rfunction/Functions.R")
#' FBA.generation( model = "./Input/FBAmodel.RData", write = T)

#' 2) Let load the Epimod package

# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="master")
library(epimod)

#' 3) model generation

model.generation(net_fname = "./Net/EpitCellDifficile.PNPRO",
								 transitions_fname = "./Net/TransitionGeneral.cpp",
								 LP = T)

system("mv EpitCellDifficile.* ./Net")

#' 4) model analysis

model.analysis(solver_fname = "./Net/EpitCellDifficile.solver",
							 i_time = 0,
							 f_time = 96,
							 s_time = 0.25,
							 parameters_fname = "Input/ParametersList.csv",
							 functions_fname = "Rfunction/Functions.R",
							 debug = T)
