
# Differential reaction fluxes between NoDrug and Therapy models. 

#' To compare the fluxes between the NoDrug and Therapy models, 
#' the comparison of individual fluxes could yield inaccurate results given the
#' large alternative optimal solution space. 
#'  
#'  We compared the flux distributions per reaction rather than the value provided 
#'  by one solution. To obtain the flux probability distributions, we performed a
#'  FVA of both the NoDrug and Therapy model.
#'   

library(dplyr)
library(ggplot2)

################################################################################
############## Bottom-Left image "Flux Variability Analysis" ###################
################################################################################

setwd("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Results")

conditions = c("NoDrug", "LowEff")

for (co in conditions) {
  
  load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")
  
  Vcd = ((pi*(0.7)^2)/4)*(((4*0.7)/6) + 9) # μm^3
  gDW_CDmax = (1/3)*Vcd*1.3 # [pg]
  
  Na = 6.022e20; c = 6.022e8
  
  load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/Diets/diets.RData")
  d = diets[["EU_average"]]
  
  ex <- sybil::findExchReact(model.mat)
  ex = ex@react_id
  value = d$`Flux Value`[which(d$Reaction %in% ex)]
  name = d$Reaction[which(d$Reaction %in% ex)]
  index = c()
  for (i in 1:length(name)) {index[i] = which(model.mat@react_id == name[i])}
  react_id_ex = which(model.mat@react_id %in% unlist(ex))
  lbNOTdiet = -10
  
  TracesPath = paste("./EpitCellDifficileHemeSink_analysis", co,
                     "/EpitCellDifficileHemeSink-analysis-1.trace", sep = "")
  
  trace = read.table(TracesPath, header = F)
  ReactionsNames = trace[1, ]
  ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
  ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
  
  trace = read.table(TracesPath, header = T)
  colnames(trace) = ReactionsNames
  
  subtrace = trace %>% 
    dplyr::select(Time, pheme_c, BiomassCD, CD, pro_L_e, leu_L_e, ile_L_e, 
                  val_L_e, trp_L_e, cys_L_e) %>% 
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  Biom = subtrace$Marking[which(subtrace$Places == "BiomassCD" & 
                                  subtrace$Time == max(subtrace$Time))]
  
  if ((gDW_CDmax - Biom) > 1e-12) {
    uppbndBiomass = (gDW_CDmax - Biom);
  } else {
    uppbndBiomass = 1e-12;
  }
  
  lowbndBiomass = lbNOTdiet
  
  model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]] =
    rep(lbNOTdiet, length(model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]]))
  
  nBac = subtrace$Marking[which(subtrace$Places == "CD" & 
                                  subtrace$Time == max(subtrace$Time))]
  
  model.mat@lowbnd[index] = (value*1e-04)/(nBac*Biom*1e-12) # [mmol/per community]
  
  Ub = 10; aaUb = 1e-12
  
  model.mat@uppbnd[react_id_ex[which(model.mat@uppbnd[react_id_ex] != 0)]] =
    rep(Ub, length(model.mat@uppbnd[react_id_ex[which(model.mat@uppbnd[react_id_ex] != 0)]]))
  
  model.mat@uppbnd[which(model.mat@react_id %in% 
                           c("EX_pro_L(e)", "EX_leu_L(e)",
                             "EX_ile_L(e)", "EX_val_L(e)", 
                             "EX_trp_L(e)", "EX_cys_L(e)"))] = rep(aaUb, 6)
  
  Met = subtrace$Marking[which(subtrace$Places %in% 
                                 c("pro_L_e", "leu_L_e", "ile_L_e", 
                                   "val_L_e", "trp_L_e", "cys_L_e") &
                                 subtrace$Time == max(subtrace$Time))]
  
  aaLb = model.mat@lowbnd[which(model.mat@react_id %in% 
                                  c("EX_pro_L(e)", "EX_leu_L(e)",
                                    "EX_ile_L(e)", "EX_val_L(e)", 
                                    "EX_trp_L(e)", "EX_cys_L(e)"))]
  
  if(nBac < 1) {
    aaLb = -((Met*c)/Na)*nBac
  } else {
    aaLb = -((Met*c)/Na)/(nBac*Biom*1e-12)
  }
  
  aaLb = -((Met*c)/Na)
  
  model.mat@lowbnd[which(model.mat@react_id %in% 
                           c("EX_pro_L(e)", "EX_leu_L(e)",
                             "EX_ile_L(e)", "EX_val_L(e)", 
                             "EX_trp_L(e)", "EX_cys_L(e)"))] = aaLb
  
  heme = subtrace$Marking[which(subtrace$Places == "pheme_c" & 
                                  subtrace$Time == max(subtrace$Time))]
  
  if(nBac < 1) {
    Lbheme = -(heme*1e+09)*nBac
  } else {
    Lbheme = -((heme*1e+09)/nBac*Biom*1e-12)
  }
  
  if(abs(Lbheme) > (heme*1e+09)){
    if(-(heme*1e+09) >= model.mat@lowbnd[which(model.mat@react_id == "sink_pheme(c)")]){
      Lbheme = -(heme*1e+09)
    } else {
      Lbheme = model.mat@lowbnd[which(model.mat@react_id == "sink_pheme(c)")]
    }
  }
  
  Ubheme = 1e-04;
  
  model.mat@uppbnd[which(model.mat@react_id == "sink_pheme(c)")] = Ubheme
  model.mat@lowbnd[which(model.mat@react_id == "sink_pheme(c)")] = Lbheme
  
  model.mat@uppbnd[which(model.mat@react_id == "EX_biomass(e)")] = uppbndBiomass
  model.mat@lowbnd[which(model.mat@react_id == "EX_biomass(e)")] = lowbndBiomass
  
  assign(paste("model", co, sep = ""), model.mat)
  Irr = sybil::mod2irrev(model.mat)
  assign(paste("model", co, "Irr", sep = ""), Irr)
  fv = sybil::fluxVar(Irr)
  
  assign(paste("lpMaxSol", co, sep = ""), sybil::maxSol(fv, "lp_obj"))
  assign(paste("lpMinSol", co, sep = ""), sybil::minSol(fv, "lp_obj"))
}

df = data.frame(lpMaxSolLowEff, lpMaxSolNoDrug, 
                Maxdiff = abs(lpMaxSolLowEff - lpMaxSolNoDrug),
                lpMinSolLowEff, lpMinSolNoDrug, 
                Mindiff = abs(lpMinSolLowEff - lpMinSolNoDrug),
                React_ID = fv@react@react_id, 
                ReactionType = c(rep("A", 985), rep("B", 985)))

a = dplyr::filter(df, grepl("EX_", React_ID))
a$ReactionType = rep("Exchange", length(a$ReactionType))
b = dplyr::filter(df, grepl("DM_", React_ID))
b$ReactionType = rep("Demand/Sink", length(b$ReactionType))
c = dplyr::filter(df, grepl("sink_", React_ID))
c$ReactionType = rep("Demand/Sink", length(c$ReactionType))
d = dplyr::filter(df, grepl("tra", React_ID))
d$ReactionType = rep("Transcription", length(d$ReactionType))
e = dplyr::filter(df, grepl("rep", React_ID))
e$ReactionType = rep("Replication", length(e$ReactionType))
f = dplyr::filter(df, !grepl("rep|tra|EX_|DM_|sink_", React_ID))
f$ReactionType = rep("Internals/Transporters", length(f$ReactionType))
h = dplyr::filter(df, grepl("biomass205", React_ID))
h$ReactionType = rep("Objective", length(h$ReactionType))

df = rbind(a, b, c, d, e, f, h)

df$lpMaxSolLowEff[which(df$lpMaxSolLowEff < 0)] = 
  numeric(length(df$lpMaxSolLowEff[which(df$lpMaxSolLowEff < 0)]))

df$lpMaxSolNoDrug[which(df$lpMaxSolNoDrug < 0)] = 
  numeric(length(df$lpMaxSolNoDrug[which(df$lpMaxSolNoDrug < 0)]))

df = filter(df, df$Maxdiff > 1e-12)
df <- df[order(df$Maxdiff), ]

# ex <- sybil::findExchReact(Irr)
# ex = ex@react_id
# 
# df = filter(df, df$React_ID %in% ex)

pdf("FVAplot.pdf", width = 6, height = 4.5, bg = "transparent")
ggplot(aes(x = lpMaxSolLowEff, y = lpMaxSolNoDrug), data = df) + 
  geom_point(aes(color = ReactionType, shape = ReactionType), size = 1.25) +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, size = 0.2,
              linetype="solid", color="black") +
  ggtitle("Flux Variability Analysis (mmol/gDW*h)") +
  xlab("log scale FV CDmodel Therapy (mmol/gDW*h)") + 
  ylab("log scale FV CDmodel No Drug (mmol/gDW*h)") +
  theme(plot.title = element_text(color="black", size = 12, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.80, 0.3),
        legend.background = element_rect(size = 0.25, fill = "white", color = "black")) +
  scale_x_log10() + scale_y_log10()
dev.off()

################################################################################
############## Bottom-Right image "Time course of aminoacid" ###################
################################################################################

setwd("~/EpiCell_CDifficile/Results/")

TracesPath = "EpitCellDifficileHemeSink_analysisNoDrug/EpitCellDifficileHemeSink-analysis-1.trace"
trace = read.table(TracesPath, header = F)

ReactionsNames = trace[1, ]
ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)

trace = read.table(TracesPath, header = T)
colnames(trace) = ReactionsNames

FluxVec = c("EX_biomass_e",
            "EX_pro_L_e", "EX_leu_L_e", "EX_ile_L_e",
            "EX_val_L_e", "EX_trp_L_e", "EX_cys_L_e",
            "sink_pheme_c")

subfluxNoDrug = trace[, c("Time", FluxVec)] %>% 
  tidyr::gather(key = "Reaction", value = "Flux", -Time)

subtraceNoDrug = trace %>% dplyr::select(Time, IECs, Drug, Damage,
                                         pheme_e, pheme_c, 
                                         BiomassCD, CD,
                                         pro_L_e, leu_L_e, ile_L_e, 
                                         val_L_e, trp_L_e, cys_L_e) %>% 
  tidyr::gather(key = "Places", value = "Marking", -Time)

subtraceNoDrug = cbind(subtraceNoDrug, 
                       Scenario = rep("No Drug", 
                                      length(subtraceNoDrug$Marking)))
subfluxNoDrug = cbind(subfluxNoDrug, 
                      Scenario = rep("No Drug", 
                                     length(subfluxNoDrug$Flux)))

TracesPath = "EpitCellDifficileHemeSink_analysisLowEff/EpitCellDifficileHemeSink-analysis-1.trace"
trace = read.table(TracesPath, header = F)

ReactionsNames = trace[1, ]
ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)

trace = read.table(TracesPath, header = T)
colnames(trace) = ReactionsNames

subfluxLowEff = trace[, c("Time", FluxVec)] %>% 
  tidyr::gather(key = "Reaction", value = "Flux", -Time)

subtraceLowEff = trace %>% dplyr::select(Time, IECs, Drug, Damage,
                                         pheme_e, pheme_c, 
                                         BiomassCD, CD,
                                         pro_L_e, leu_L_e, ile_L_e, 
                                         val_L_e, trp_L_e, cys_L_e) %>% 
  tidyr::gather(key = "Places", value = "Marking", -Time)

subtraceLowEff = cbind(subtraceLowEff, 
                       Scenario = rep("Low Efficacy", 
                                      length(subtraceLowEff$Marking)))
subfluxLowEff = cbind(subfluxLowEff, 
                      Scenario = rep("Low Efficacy", 
                                     length(subfluxLowEff$Flux)))

subtrace = rbind(subtraceNoDrug, subtraceLowEff)
subflux = rbind(subfluxNoDrug, subfluxLowEff)

amino = filter(subtrace, grepl("pro_L_e|leu_L_e|ile_L_e|val_L_e", Places))
for (i in 1:length(amino$Scenario)) {if(amino$Scenario[i] == "Low Efficacy") {amino$Scenario[i] = "Therapy"}}
for (i in 1:length(amino$Places)) {if(amino$Places[i] == "pro_L_e") {amino$Places[i] = "Proline"}}
for (i in 1:length(amino$Places)) {if(amino$Places[i] == "leu_L_e") {amino$Places[i] = "Leucine"}}
for (i in 1:length(amino$Places)) {if(amino$Places[i] == "ile_L_e") {amino$Places[i] = "Isoleucine"}}
for (i in 1:length(amino$Places)) {if(amino$Places[i] == "val_L_e") {amino$Places[i] = "Valine"}}

pdf("AAs.pdf", width = 6.5, height = 4.5, bg = "transparent")
ggplot(amino, aes(Time, (Marking/(6.022e20/6.022e08)),
                  colour = Scenario, linetype = Places)) +
  geom_line(size = 0.80) + theme_bw() +
  ggtitle("Changes in AAs abundance \n in FREE and THERAPY models") +
  xlab("Time (hour)") + ylab("Luminal aminoacids (mmol/mL)") +
  labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan')) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.12, 0.56),
        legend.background = element_rect(size = 0.25, 
                                         fill = "white", color = "black")) +
  labs(colour = "Condition") + 
  labs(linetype = "Aminoacid")
dev.off()

