
library(dplyr)
library(ggplot2)

################################################################################
############ Flux Variability Analysis (Multiple time-point) ###################
################################################################################

setwd("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Results")

conditions = c("NoDrug", "LowEff")
# times = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45); filter = 1e-12
times = 0

gDW_CDmax = (1/3)*((pi*(0.7)^2)/4)*(((4*0.7)/6) + 9)*1.3; Na = 6.022e20; c = 6.022e8

for (t in times) {
  for (co in conditions) {
    
    load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")
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
    
    subtrace = trace %>% dplyr::select(Time, 
                                       pheme_c, BiomassCD, CD, 
                                       pro_L_e, leu_L_e, ile_L_e, 
                                       val_L_e, trp_L_e, cys_L_e) %>% 
      tidyr::gather(key = "Places", value = "Marking", -Time)
    
    Biom = subtrace$Marking[which(subtrace$Places == "BiomassCD" & 
                                    subtrace$Time == t)]
    
    if ((gDW_CDmax - Biom) > 1e-12) {
      uppbndBiomass = (gDW_CDmax - Biom)
    } else {
      uppbndBiomass = 1e-12;
    }
    
    if ((gDW_CDmax - Biom) > 1e-12) {
      uppbndBiomass = (gDW_CDmax - Biom)
    } else {
      uppbndBiomass = 1e-12;
    }
    
    lowbndBiomass = lbNOTdiet
    
    model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]] =
      rep(lbNOTdiet, length(model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]]))
    
    nBac = subtrace$Marking[which(subtrace$Places == "CD" & 
                                    subtrace$Time == t)]
    
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
                                   subtrace$Time == t)]
    
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
                                    subtrace$Time == t)]
    
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
  
  if(max(df$Maxdiff) != 0) {
    df = filter(df, df$Maxdiff > filter)
  } else {}
  
  df <- df[order(df$Maxdiff), ]
  
  ggplot(aes(x = lpMaxSolLowEff, y = lpMaxSolNoDrug), data = df) + 
    geom_point(aes(color = ReactionType), size = 1.25) +
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
}
