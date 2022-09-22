
setwd("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile")

# Lets compile the metabolic model from the RData storing the chemical reaction network

source("Rfunction/Functions.R")
# FBA.generation(model = "./Input/CDmodels/CD196HemeSink/CD196HemeSink.RData", write = T)

# Lets generate the a flux distribution to test the optimization of a mathematical programming object

source("./Rfunction/flux_balance.R")
load("./Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")

mtf_sol = flux_balance(model = model.mat, typename = model.mat@mod_desc, diet = F)
v = sybil::getFluxDist(mtf_sol[[1]])

# Through flux_balance method we can perform sensitivity analysis on target reactions

Sink = ("sink_pheme(c)")
EXs = c("EX_pro_L(e)", "EX_leu_L(e)", "EX_ile_L(e)", "EX_val_L(e)", "EX_trp_L(e)", "EX_cys_L(e)")

# acting on the permissible flux through target reactions
P = c(0.0, 0.001, 0.003, 0.005, 0.008, 0.015, 0.03, 0.05)
FBA.allScenarios = list()

EX_biomass = c(); EX_pro = c(); EX_leu = c();  EX_ile = c()
EX_val = c(); EX_trp = c(); EX_cys = c(); Sink_heme = c()

for (i in 1:length(P)) {
  
  # Lb
  model.mat@lowbnd[which(model.mat@react_id %in% EXs)] = 
    model.mat@lowbnd[which(model.mat@react_id %in% EXs)] * P[i]
  # Ub
  model.mat@uppbnd[which(model.mat@react_id %in% EXs)] = 
    model.mat@uppbnd[which(model.mat@react_id %in% EXs)] * P[i]
  
  cat("                                                               ", "\n")
  cat("lowbnd:", sybil::lowbnd(model.mat)[which(model.mat@react_id %in% EXs)])
  cat("                                                               ", "\n")
  cat("uppbnd:", sybil::uppbnd(model.mat)[which(model.mat@react_id %in% EXs)])
  cat("                                                               ", "\n")
  
  # FBA solutions
  mtf_sol = flux_balance(model = model.mat, typename = model.mat@mod_desc, diet = F)
  
  FBA.allScenarios[[i]] = list(mtf_sol = mtf_sol, 
                               lb = sybil::lowbnd(model.mat),
                               ub = sybil::uppbnd(model.mat))
  
  EX_biomass[i] = FBA.allScenarios[[i]][[1]][[1]]@alg_par[["wtobj"]]
  
  v = sybil::getFluxDist(mtf_sol[[1]])
  
  EX_pro[i] = v[which(model.mat@react_id %in% EXs)][1]
  EX_leu[i] = v[which(model.mat@react_id %in% EXs)][2]
  EX_ile[i] = v[which(model.mat@react_id %in% EXs)][3]
  EX_val[i] = v[which(model.mat@react_id %in% EXs)][4]
  EX_trp[i] = v[which(model.mat@react_id %in% EXs)][5]
  EX_cys[i] = v[which(model.mat@react_id %in% EXs)][6]
  
  Sink_heme[i] = v[which(model.mat@react_id == "sink_pheme(c)")]
  
  load("./Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")
  
}

objdata = data.frame(EX_biomass, EX_pro, EX_leu, EX_ile, 
                     EX_val, EX_trp, EX_cys, Sink_heme, 
                     BoundsActivity = P)

objdata = tidyr::gather(objdata, key, Value, -BoundsActivity)

ggplot2::ggplot(objdata, ggplot2::aes(x = BoundsActivity, y = Value, color = key)) + 
  ggplot2::geom_line() + ggplot2::ylim(-50, 150) + 
  ggplot2::geom_vline(xintercept = P, linetype="dotted") +
  ggplot2::geom_point(col = "black", shape = 21, size = 1) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 12), 
                 plot.title = ggplot2::element_text(size = 12),
                 axis.text = ggplot2::element_text(size = 10),
                 legend.key.size = ggplot2::unit(0.35, 'cm'),
                 legend.title=ggplot2::element_text(size = 10), 
                 legend.text=ggplot2::element_text(size = 10),
                 axis.title.x = ggplot2::element_text(color = "black", size = 12),
                 axis.title.y = ggplot2::element_text(color = "black", size = 12)) +
  ggplot2::scale_colour_discrete(name="Target Reactions") +
  ggplot2::ggtitle("Changing bounds Test (Sentivity analysis on essential AAs)") + 
  ggplot2::labs(y = "Flux [mmol/h*gDW]", x = "Bounds activity [x*lb; x*up]")
  
names(FBA.allScenarios) = c(paste(as.character(P*100), "%", sep = " "))

ReactFlux = lapply(names(FBA.allScenarios), function(x){
  return(data.frame(Scenario = x, v = sybil::getFluxDist(FBA.allScenarios[[x]]$mtf_sol[[1]]),
                    ReactionName = model.mat@react_id,
                    lb = FBA.allScenarios[[x]]$l, ub = FBA.allScenarios[[x]]$u))})

ReactFluxdf = do.call(rbind.data.frame, ReactFlux)

# 5) Lets loaded diets data
load("./Input/Diets/diets.RData")
load("./Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")

diets[["Template"]] = diets[["Unhealthy"]]

model.mat@react_id = gsub("\\(", replacement = "_", model.mat@react_id)
model.mat@react_id = gsub("\\)", replacement = "", model.mat@react_id)

diets[["Template"]]$`Flux Value`[which(diets[["Template"]]$Reaction %in% model.mat@react_id)] = 
  model.mat@uppbnd[which(model.mat@react_id %in% diets[["Template"]]$Reaction)]

diets[["Zero"]] = diets[["Unhealthy"]]
diets[["Zero"]]$`Flux Value` = (diets[["Zero"]]$`Flux Value`)*0

diets[["ZeroEssential"]] = diets[["Unhealthy"]]
diets[["ZeroEssential"]]$`Flux Value`[which(diets[["Unhealthy"]]$Reaction %in% 
                                              c("EX_pro_L_e", "EX_leu_L_e", "EX_ile_L_e", 
                                                "EX_val_L_e", "EX_trp_L_e", "EX_cys_L_e"))] = rep(0, 6)

FBA.allScenarios = list()

for (i in 1:length(names(diets))) {
  
  mtf_sol = flux_balance(model = model.mat,
                         typename = model.mat@mod_desc,
                         diet = names(diets)[i])
  
  FBA.allScenarios[[i]] = list(mtf_sol = mtf_sol[[1]],
                               lb =  mtf_sol[[2]],
                               ub =  mtf_sol[[3]])
  
}

names(FBA.allScenarios) = names(diets)

ReactFlux = lapply(names(FBA.allScenarios), function(x){
  return(data.frame(Scenario = x,
                    v = sybil::getFluxDist(FBA.allScenarios[[x]]$mtf_sol),
                    ReactionName = model.mat@react_id))})

ReactFluxdf = do.call(rbind.data.frame, ReactFlux)

df = ReactFluxdf[which(ReactFluxdf$ReactionName %in% 
                         c("EX_biomass(e)", "biomass205", "EX_pro_L(e)", "EX_leu_L(e)",
                           "EX_ile_L(e)", "EX_val_L(e)", "EX_trp_L(e)", "EX_cys_L(e)", "sink_pheme(c)"),), ]

FluxesvsDiet = ggplot2::ggplot(df, ggplot2::aes(x = ReactionName, y = v, fill = Scenario)) +
  ggplot2::geom_col(width = 0.33, position = ggplot2::position_dodge(0.5)) + 
  ggplot2::ggtitle("Fluxes Diets") + ggplot2::theme(axis.title = ggplot2::element_text(size = 12), 
                                                    plot.title = ggplot2::element_text(size=8),
                                                    axis.text = ggplot2::element_text(size = 10),
                                                    legend.key.size = ggplot2::unit(0.15, 'cm'),
                                                    legend.title=ggplot2::element_text(size = 6), 
                                                    legend.text=ggplot2::element_text(size=5))

library(patchwork)
p1 + p2 + p3 + p4 | p5 | p6 / p7

## Fluxes comparison
