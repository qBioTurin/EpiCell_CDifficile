
source("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Rfunction/flux_balance_fba.R")
# FBA.generation(model = "./Input/CDmodels/CD196HemeSink/CD196HemeSink.RData", write = T)
load("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")

set.seed(1234)

P = c(1e-350, round(sort(runif(18, min = 0.001, max = 0.15)), 5))

FBA.allScenarios = list()

Sink = ("sink_pheme(c)")
EXs = c("EX_pro_L(e)", "EX_leu_L(e)", "EX_ile_L(e)", 
        "EX_val_L(e)", "EX_trp_L(e)", "EX_cys_L(e)")

EX_biomass = c(); EX_pro = c(); EX_leu = c();  EX_ile = c()
EX_val = c(); EX_trp = c(); EX_cys = c(); Sink_heme = c()

# for (i in 1:length(P)) {

# Lb
# model.mat@lowbnd[which(model.mat@react_id %in% EXs)] = 
# model.mat@lowbnd[which(model.mat@react_id %in% EXs)] * P[i]
# Ub
# model.mat@uppbnd[which(model.mat@react_id %in% EXs)] = 
# model.mat@uppbnd[which(model.mat@react_id %in% EXs)] * P[i]

# cat("                                                               ", "\n")
# cat("lowbnd:", sybil::lowbnd(model.mat)[which(model.mat@react_id %in% EXs)])
# cat("                                                               ", "\n")
# cat("uppbnd:", sybil::uppbnd(model.mat)[which(model.mat@react_id %in% EXs)])
# cat("                                                               ", "\n")

# FBA solutions
# mtf_sol = flux_balance(model = model.mat, typename = model.mat@mod_desc, diet = F)

# FBA.allScenarios[[i]] = list(mtf_sol = mtf_sol, 
# lb = sybil::lowbnd(model.mat),
# ub = sybil::uppbnd(model.mat))

# v = mtf_sol[[1]][["fluxes"]]

# EX_biomass[i] = v[which(model.mat@react_id == "EX_biomass(e)")]

# EX_pro[i] = v[which(model.mat@react_id == "EX_pro_L(e)")]
# EX_leu[i] = v[which(model.mat@react_id == "EX_leu_L(e)")]
# EX_ile[i] = v[which(model.mat@react_id == "EX_ile_L(e)")]
# EX_val[i] = v[which(model.mat@react_id == "EX_val_L(e)")]
# EX_trp[i] = v[which(model.mat@react_id == "EX_trp_L(e)")]
# EX_cys[i] = v[which(model.mat@react_id == "EX_cys_L(e)")]

# Sink_heme[i] = v[which(model.mat@react_id == "sink_pheme(c)")]

# load("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")

# }

set.seed(1234)

P = c(0.0, round(sort(runif(18, min = 0.001, max = 0.15)), 5))

names(FBA.allScenarios) = as.character(P*100)

ReactFlux = lapply(names(FBA.allScenarios), function(x){
  return(data.frame(Scenario = x, v = sybil::getFluxDist(FBA.allScenarios[[x]]$mtf_sol[[1]]),
                    ReactionName = model.mat@react_id,
                    lb = FBA.allScenarios[[x]]$l, ub = FBA.allScenarios[[x]]$u))})

ReactFluxdf = do.call(rbind.data.frame, ReactFlux)
# SubReactFluxdf = ReactFluxdf[which(ReactFluxdf$ReactionName %in% c("EX_leu_L(e)", "EX_pro_L(e)")), ]

df = data.frame(ReactionName = unique(ReactFluxdf$ReactionName))

for (i in unique(ReactFluxdf$Scenario)) {
  df = cbind(df, v = ReactFluxdf$v[which(ReactFluxdf$Scenario == i)])
} 

colnames(df) = c("ReactionName", names(FBA.allScenarios))

dfdiff = data.frame()

for (i in 2:length(colnames(df))) {
  ReactionName = unique(ReactFluxdf$ReactionName)
  Value = abs(df[, i]) - abs(df[, length(colnames(df))])
  key = rep(paste(LETTERS[i-1], paste(colnames(df)[i], colnames(df)[length(colnames(df))], 
                                      sep = " - "), sep = ": "), model.mat@react_num)
  
  data = data.frame(ReactionName, key, Value)
  
  dfdiff = rbind(dfdiff, data)
}

df = tidyr::gather(df, key, Value, -ReactionName)

p = df %>%
  ggplot2::ggplot(ggplot2::aes(x = ReactionName, y = Value, group = key, color = key)) +
  ggplot2::geom_line() + viridis::scale_color_viridis(discrete = TRUE) +
  ggplot2::ggtitle("Estimated fluxes differences between a bounds activity and max activity") +
  hrbrthemes::theme_ipsum() +
  ggplot2::ylab("DeltaFluxes") + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 12), 
                 plot.title = ggplot2::element_text(size = 12),
                 axis.text = ggplot2::element_text(size = 10),
                 legend.key.size = ggplot2::unit(0.35, 'cm'),
                 legend.title=ggplot2::element_text(size = 10), 
                 legend.text=ggplot2::element_text(size = 10),
                 axis.title.x = ggplot2::element_text(color = "black", size = 12),
                 axis.text.x = ggplot2::element_text(size = 1),
                 axis.title.y = ggplot2::element_text(color = "black", size = 12))

pdiff = dfdiff %>%
  ggplot2::ggplot(ggplot2::aes(x = ReactionName, y = Value, group = key, color = key)) +
  ggplot2::geom_line() + viridis::scale_color_viridis(discrete = TRUE) +
  ggplot2::ggtitle("Estimated fluxes differences between a bounds activity and max activity") +
  hrbrthemes::theme_ipsum() +
  ggplot2::ylab("DeltaFluxes") + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 12), 
                 plot.title = ggplot2::element_text(size = 12),
                 axis.text = ggplot2::element_text(size = 10),
                 legend.key.size = ggplot2::unit(0.35, 'cm'),
                 legend.title=ggplot2::element_text(size = 10), 
                 legend.text=ggplot2::element_text(size = 10),
                 axis.title.x = ggplot2::element_text(color = "black", size = 12),
                 axis.text.x = ggplot2::element_text(size = 1),
                 axis.title.y = ggplot2::element_text(color = "black", size = 12))

df$Value = df$Value

for (i in unique(df$key)) {
  print(mean(df[which(df$key == i), ]$Value)) 
}
