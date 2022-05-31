
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

