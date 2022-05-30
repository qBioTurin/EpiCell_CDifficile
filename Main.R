
# 1) Lets generate the RData storing the matlab FBA model
source("Rfunction/Functions.R")
FBA.generation(model = "./Input/CD196_heme.RData", write = T)
# 2) Lets generate the a flux distribution
load("./Input/CD196_heme.RData")
source("./Rfunction/flux_balance.R")
# 2) Lets loaded diets data
load("./Input/diets.RData")

FBA.allScenarios = list()

for (i in 1:length(names(diets))) {

	mtf_sol = flux_balance(model = model.mat,
												 typename = "CD196_heme",
												 diet = names(diets)[i])

	FBA.allScenarios[[i]] = list(mtf_sol = mtf_sol[[1]],
															 lb =  mtf_sol[[2]],
															 ub =  mtf_sol[[3]])

}

names(FBA.allScenarios) = names(diets)

ReactFlux = lapply(names(FBA.allScenarios), function(x){
	return(data.frame(Scenario = x,
										v = sybil::getFluxDist(FBA.allScenarios[[x]]$mtf_sol),
										ReactionName = model.mat@react_id,
										lb = FBA.allScenarios[[x]]$l,
										ub = FBA.allScenarios[[x]]$u))})

ReactFluxdf = do.call(rbind.data.frame, ReactFlux)

ReactFluxdf[which(ReactFluxdf$ReactionName %in% c("EX_biomass(e)",
																									"EX_pro_L(e)", "EX_leu_L(e)",
																									"EX_ile_L(e)", "EX_val_L(e)",
																									"EX_trp_L(e)", "EX_cys_L(e)"),), ]

#### This is with 1370 (without the Ex and Tr of Sheme)
# source("Rfunction/MatlabFileGeneration.R")
# Matlab.file = readMATmod("Input/CD196_heme.mat")
# save(Matlab.file, file = "Input/FBAmodel.RData")
# source("Rfunction/Functions.R")
#FBA.generation( model = "./Input/FBAmodel.RData", write = T)

#2) Let load the Epimod package

# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="master")

library(epimod)

model.generation(net_fname = "./Net/EpitCellDifficile.PNPRO",
								 transitions_fname = "./Net/TransitionGeneral.cpp",
								 LP = T)

#4) model analysis
model.analysis(solver_fname = "./Net/EpitCellDifficile.solver",
							 i_time = 1,
							 f_time = 2,
							 s_time = 0.5,
							 parameters_fname = "Input/ParametersList.csv",
							 functions_fname = "Rfunction/Functions.R",
							 debug = T)
