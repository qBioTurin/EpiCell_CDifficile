
# 1) Lets generate the RData storing the matlab FBA model
# source("Rfunction/Functions.R")
# FBA.generation(model = "./Input/CD196_heme.RData", write = T)

# 2) Lets generate the a flux distribution
load("./Input/CD196_heme.RData")
source("./Rfunction/flux_balance.R")

EXs = c("EX_pro_L(e)", "EX_leu_L(e)", "EX_ile_L(e)", 
        "EX_val_L(e)", "EX_trp_L(e)", "EX_cys_L(e)")

# 4) Performing sensitivity analysis

P = c(0.008, 0.015, 0.05, 0.15, 1)

FBA.allScenarios = list()

source("./Rfunction/flux_balance.R")

for (i in 1:length(P)) {
  
  # Lb
  model.mat@lowbnd[which(model.mat@react_id %in% EXs)] = 
    model.mat@lowbnd[which(model.mat@react_id %in% EXs)] * P[i]
  # Ub
  model.mat@uppbnd[which(model.mat@react_id %in% EXs)] = 
    model.mat@uppbnd[which(model.mat@react_id %in% EXs)] * P[i]
  
  print(sybil::lowbnd(model.mat)[which(model.mat@react_id %in% EXs)])
  print(sybil::uppbnd(model.mat)[which(model.mat@react_id %in% EXs)])
  
  # FBA solutions
  mtf_sol = flux_balance(model = model.mat, typename = "CD196_heme", diet = F)
  
  FBA.allScenarios[[i]] = list(mtf_sol = mtf_sol, 
                               lb = sybil::lowbnd(model.mat),
                               ub = sybil::uppbnd(model.mat))
  
  load("./Input/CD196_heme.RData")
  
}

names(FBA.allScenarios) = c(paste("a ", P[1]*100, "%", sep = ""),
                            paste("b ", P[2]*100, "%", sep = ""),
                            paste("c ", P[3]*100, "%", sep = ""),
                            paste("d ", P[4]*100, "%", sep = ""),
                            paste("e ", P[5]*100, "%", sep = ""))

ReactFlux = lapply(names(FBA.allScenarios), function(x){ 
  return(data.frame(Scenario = x,
                    v = sybil::getFluxDist(FBA.allScenarios[[x]]$mtf_sol[[1]]),
                    ReactionName = model.mat@react_id,
                    lb = FBA.allScenarios[[x]]$l,
                    ub = FBA.allScenarios[[x]]$u))})

ReactFluxdf = do.call(rbind.data.frame, ReactFlux)

models = list()
model.diag.metrics = list()

for (s in 2:length(unique(ReactFluxdf$Scenario))) {
  
  Activity.1per = 
    ReactFluxdf$v[which(ReactFluxdf$Scenario == 
                          paste(unique(ReactFluxdf$Scenario)[1]))]
  
  models[[s-1]] =
    lm(ReactFluxdf$v[which(ReactFluxdf$Scenario == 
                             paste(unique(ReactFluxdf$Scenario)[s]))] ~ 
         Activity.1per, data = ReactFluxdf)
  
  model.diag.metrics[[s-1]] = broom::augment(models[[s-1]])
  
  names(model.diag.metrics[[s-1]])[1] = paste("v (bounds: ", P[s]*100, " %) [mmol gDW−1 h−1]", sep = "")
  names(model.diag.metrics[[s-1]])[2] = paste("v (bounds: ", P[1]*100, " %) [mmol gDW−1 h−1]", sep = "")
  
}

p1 = ggplot2::ggplot(model.diag.metrics[[1]], 
                     ggplot2::aes(`v (bounds: 0.8 %) [mmol gDW−1 h−1]`, 
                                  `v (bounds: 1.5 %) [mmol gDW−1 h−1]`)) + 
  ggplot2::geom_point() + ggplot2::stat_smooth(method = lm, se = FALSE) + 
  ggplot2::geom_segment(ggplot2::aes(xend = 
                                       `v (bounds: 1.5 %) [mmol gDW−1 h−1]`, 
                                     yend = .fitted), 
                        color = "red", size = 0.3) + 
  ggplot2::ggtitle(paste(unique(ReactFluxdf$Scenario)[2], 
                         sep = " ")) + 
  ggplot2::theme(axis.title=ggplot2::element_text(size=6), 
                 plot.title = ggplot2::element_text(size=8),
                 axis.text = ggplot2::element_text(size = 4))

p2 = ggplot2::ggplot(model.diag.metrics[[2]], 
                     ggplot2::aes(`v (bounds: 0.8 %) [mmol gDW−1 h−1]`, 
                                  `v (bounds: 5 %) [mmol gDW−1 h−1]`)) + 
  ggplot2::geom_point() + ggplot2::stat_smooth(method = lm, se = FALSE) + 
  ggplot2::geom_segment(ggplot2::aes(xend = 
                                       `v (bounds: 5 %) [mmol gDW−1 h−1]`, 
                                     yend = .fitted), 
                        color = "red", size = 0.3) + 
  ggplot2::ggtitle(paste(unique(ReactFluxdf$Scenario)[3], 
                         sep = " ")) + 
  ggplot2::theme(axis.title=ggplot2::element_text(size=6), 
                 plot.title = ggplot2::element_text(size=8),
                 axis.text = ggplot2::element_text(size = 4))

# Plot3: Dodged bar short
p3 = ggplot2::ggplot(model.diag.metrics[[3]], 
                     ggplot2::aes(`v (bounds: 0.8 %) [mmol gDW−1 h−1]`, 
                                  `v (bounds: 15 %) [mmol gDW−1 h−1]`)) + 
  ggplot2::geom_point() + ggplot2::stat_smooth(method = lm, se = FALSE) + 
  ggplot2::geom_segment(ggplot2::aes(xend = 
                                       `v (bounds: 15 %) [mmol gDW−1 h−1]`, 
                                     yend = .fitted), 
                        color = "red", size = 0.3) + 
  ggplot2::ggtitle(paste(unique(ReactFluxdf$Scenario)[4], 
                         sep = " ")) + 
  ggplot2::theme(axis.title=ggplot2::element_text(size=6), 
                 plot.title = ggplot2::element_text(size=8),
                 axis.text = ggplot2::element_text(size = 4))

p4 = ggplot2::ggplot(model.diag.metrics[[4]], 
                     ggplot2::aes(`v (bounds: 0.8 %) [mmol gDW−1 h−1]`, 
                                  `v (bounds: 100 %) [mmol gDW−1 h−1]`)) + 
  ggplot2::geom_point() + ggplot2::stat_smooth(method = lm, se = FALSE) + 
  ggplot2::geom_segment(ggplot2::aes(xend = 
                                       `v (bounds: 100 %) [mmol gDW−1 h−1]`, 
                                     yend = .fitted), 
                        color = "red", size = 0.3) + 
  ggplot2::ggtitle(paste(unique(ReactFluxdf$Scenario)[5], 
                         sep = " ")) + 
  ggplot2::theme(axis.title=ggplot2::element_text(size=6), 
                 plot.title = ggplot2::element_text(size=8),
                 axis.text = ggplot2::element_text(size = 4))

p5 = ggpubr::ggbarplot(cbind(
  ReactFluxdf[which(ReactFluxdf$ReactionName == "EX_biomass(e)"), ],
  r = ReactFluxdf[which(ReactFluxdf$ReactionName == "EX_biomass(e)"), ]$v - 
    ReactFluxdf[which(ReactFluxdf$ReactionName == "EX_biomass(e)"), ]$v[1]), 
  x = "Scenario", y = "r", fill = "Scenario", color = "white",
  palette = "Blues", # jco journal color palett. see ?ggpar
  sort.val = "none", sort.by.groups = TRUE, x.text.angle = 0, width=0.8) + 
  ggplot2::ggtitle("percentage bounds activity relative to defualt bounds") +
  ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                 plot.title = ggplot2::element_text(size=8),
                 axis.text = ggplot2::element_text(size = 7),
                 legend.key.size = ggplot2::unit(0.3, 'cm'),
                 legend.title=ggplot2::element_text(size = 6),
                 legend.text=ggplot2::element_text(size=5))

# source("./Rfunction/multiplot.R")
# multiplot(p1, p2, p3, p4, p5, layout = matrix(c(1, 3, 2 ,4, 5, 5), ncol=3, byrow=F))

df = ReactFluxdf[which(ReactFluxdf$ReactionName %in% 
                         c("EX_biomass(e)", "EX_pro_L(e)", "EX_leu_L(e)",
                           "EX_ile_L(e)", "EX_val_L(e)", "EX_trp_L(e)", "EX_cys_L(e)"),), ]

p6 = ggplot2::ggplot(df, ggplot2::aes(x = ReactionName, y = v, fill = Scenario)) +
  ggplot2::geom_col(width = 0.33, position = ggplot2::position_dodge(0.5)) + 
  ggplot2::ggtitle("Fluxes Essentials") + ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                                                         plot.title = ggplot2::element_text(size=8),
                                                         axis.text = ggplot2::element_text(size = 5),
                                                         legend.key.size = ggplot2::unit(0.15, 'cm'),
                                                         legend.title=ggplot2::element_text(size = 6), 
                                                         legend.text=ggplot2::element_text(size=5))

# 5) Lets loaded diets data
load("./Input/diets.RData")
diets[["Zero"]] = diets[["Unhealthy"]]
diets[["ZeroEssential"]] = diets[["Unhealthy"]]
diets[["Zero"]]$`Flux Value` = rep(0, length(length(diets[["Unhealthy"]]$`Flux Value`)))
diets[["ZeroEssential"]]$`Flux Value`[which(diets[["Zero"]]$Reaction %in% EXs)] = rep(0, 6)

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
                    ReactionName = model.mat@react_id))})

ReactFluxdf = do.call(rbind.data.frame, ReactFlux)

df = ReactFluxdf[which(ReactFluxdf$ReactionName %in% 
                         c("EX_biomass(e)", "EX_pro_L(e)", "EX_leu_L(e)",
                           "EX_ile_L(e)", "EX_val_L(e)", "EX_trp_L(e)", "EX_cys_L(e)"),), ]

p7 = ggplot2::ggplot(df, ggplot2::aes(x = ReactionName, y = v, fill = Scenario)) +
  ggplot2::geom_col(width = 0.33, position = ggplot2::position_dodge(0.5)) + 
  ggplot2::ggtitle("Fluxes Diets") + ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                                                    plot.title = ggplot2::element_text(size=8),
                                                    axis.text = ggplot2::element_text(size = 5),
                                                    legend.key.size = ggplot2::unit(0.15, 'cm'),
                                                    legend.title=ggplot2::element_text(size = 6), 
                                                    legend.text=ggplot2::element_text(size=5))

library(patchwork)
p1 + p2 + p3 + p4 | p5 | p6 / p7

## Fluxes comparison

saveRDS(subflux_total, "./glpkFluxes.rds")
glpkFluxes <- readRDS("./glpkFluxes.rds")
glpkFluxes["Time"] = NULL
glpkFluxes["Scenarion"] = c(rep("a 0.8%", 8), rep("b 1.5%", 8), 
                            rep("c 5%", 8), rep("d 15%", 8), rep("e 100%", 8))
names(glpkFluxes) = c("ReactionName", "v", "Scenario")
glpkFluxes = glpkFluxes[order(glpkFluxes$Scenario), ]
glpkFluxes = glpkFluxes[-c(which(glpkFluxes$ReactionName == "EX_pheme.e.")), ]

sybilFluxes <- readRDS("./sybilFluxes.rds")
sybilFluxes["lb"] = NULL
sybilFluxes["ub"] = NULL

p1 = ggplot2::ggplot(sybilFluxes, ggplot2::aes(x = ReactionName, y = v, fill = Scenario)) +
  ggplot2::geom_col(width = 0.33, position = ggplot2::position_dodge(0.5)) + 
  ggplot2::ggtitle("sybilFluxes") + ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                                                         plot.title = ggplot2::element_text(size=8),
                                                         axis.text = ggplot2::element_text(size = 5),
                                                         legend.key.size = ggplot2::unit(0.15, 'cm'),
                                                         legend.title=ggplot2::element_text(size = 6), 
                                                         legend.text=ggplot2::element_text(size=5))

p2 = ggplot2::ggplot(glpkFluxes, ggplot2::aes(x = ReactionName, y = v, fill = Scenario)) +
  ggplot2::geom_col(width = 0.33, position = ggplot2::position_dodge(0.5)) + 
  ggplot2::ggtitle("glpkFluxes") + ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                                                         plot.title = ggplot2::element_text(size=8),
                                                         axis.text = ggplot2::element_text(size = 5),
                                                         legend.key.size = ggplot2::unit(0.15, 'cm'),
                                                         legend.title=ggplot2::element_text(size = 6), 
                                                         legend.text=ggplot2::element_text(size=5))

p1 + p2