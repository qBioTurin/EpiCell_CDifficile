
# 1) Lets generate the RData storing the matlab FBA model
source("Rfunction/Functions.R")
FBA.generation( model = "./Input/CD196_heme.RData", write = T)

#### This is with 1370 (without the Ex and Tr of Sheme)
# source("Rfunction/MatlabFileGeneration.R")
# Matlab.file = readMATmod("Input/CD196_heme.mat")
# save(Matlab.file, file = "Input/FBAmodel.RData")
# source("Rfunction/Functions.R")
#FBA.generation( model = "./Input/FBAmodel.RData", write = T)

#2) Let load the Epimod package
#library(devtools)
#install_github("https://github.com/qBioTurin/epimod", ref="master")
library(epimod)

model.generation(net_fname = "./Net/SensitivityFBA.PNPRO",
								 transitions_fname = "./Net/TransitionGeneral.cpp",
								 LP = T)

#3)
paste0("unfolding2 ./home/SensitivityFBA -long-names ")
paste0("PN2ODE.sh ./home/SensitivityFBA -C ./home/TransitionGeneral.cpp -H -M")
model.analysis(solver_fname = "./Net/SensitivityFBA.solver",
							 i_time = 1,
							 f_time = 2,
							 s_time = 0.5,
							 parameters_fname = "Input/ParametersList.csv",
							 functions_fname = "Rfunction/Functions.R",debug = T)
