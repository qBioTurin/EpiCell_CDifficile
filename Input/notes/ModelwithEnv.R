
load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")
load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/Diets/diets.RData")

source("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Rfunction/flux_balance_fba.R")

d = diets[["EU_average"]]

ex <- sybil::findExchReact(model.mat)
ex = ex@react_id

value = d$`Flux Value`[which(d$Reaction %in% ex)]
name = d$Reaction[which(d$Reaction %in% ex)]
index = c()

for (i in 1:length(name)) {
  index[i] = which(model.mat@react_id == name[i])
}

react_id_ex = which(model.mat@react_id %in% unlist(ex))

lbNOTdiet = -1e-08

model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]] =
  rep(lbNOTdiet, length(model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]]))

model.mat@lowbnd[index] = value

save(model.mat, file = "~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSinkDiet.RData")

# check it out

# load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSinkDiet.RData")
# source("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Rfunction/flux_balance_fba.R")
# load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/Diets/diets.RData")
# 
# sol = flux_balance(model = model.mat, ShowExBounds = T,
#                    typename = model.mat@mod_desc,
#                    diet = "EU_average", lbNOTdiet = -1e-08, 
#                    symmetric = F, uppbndBiomass = "template")
# 
# v = sol[[1]][["fluxes"]]
# 
# EX_biomass = c()
# EX_ile = c(); EX_pro = c(); EX_leu = c()
# EX_val = c(); EX_cys = c(); EX_trp = c()
# Sink_heme = c()
# 
# EX_biomass = v[which(model.mat@react_id == 'EX_biomass(e)')]
# EX_leu = v[which(model.mat@react_id == 'EX_leu_L(e)')]
# EX_pro = v[which(model.mat@react_id == 'EX_pro_L(e)')]
# EX_ile = v[which(model.mat@react_id == 'EX_ile_L(e)')]
# EX_val = v[which(model.mat@react_id == 'EX_val_L(e)')]
# EX_cys = v[which(model.mat@react_id == 'EX_cys_L(e)')]
# EX_trp = v[which(model.mat@react_id == 'EX_trp_L(e)')]
# Sink_heme = v[which(model.mat@react_id == 'sink_pheme(c)')]
# 
# DietData = data.frame(EX_biomass,
#                       EX_pro, EX_leu, EX_val,
#                       EX_ile, EX_cys, EX_trp,
#                       Sink_heme, DietSet = "Template")
# 
# DietData = tidyr::gather(DietData, key, Value, -DietSet)
# 
# ggplot(DietData) + geom_bar(aes(x = DietSet, y = Value, fill = key),
#                             stat = "identity", position = "dodge") +
#   scale_y_continuous("Estimated Flux [mmol/gDW*h", expand = c(0, 0)) +
#   scale_x_discrete("Reaction") +
#   scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73",
#                                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
#   hrbrthemes::theme_ipsum() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0),
#         axis.line = element_blank(), axis.ticks.x = element_blank())
