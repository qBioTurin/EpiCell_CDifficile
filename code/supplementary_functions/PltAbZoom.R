
f.Time = 16
size = 0.25
Exper = "Model_Analysis"
colConNoDrug = "darkred"
colConTherapy = "darkcyan"
coltag = "darkgray"

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

tag = c("Ablated", "Unified")
Condition = c("NoDrug", "Therapy")

for (i in Condition) {
  for (j in tag) {
    
    subtrace = readr::read_csv(paste("./Results/CDiff_", j, "_", i, "_", 
                                     Exper, "/SubTrace", j, i, ".csv", sep = ""))[-1]
    
    subflux = readr::read_csv(paste("./Results/CDiff_", j, "_", i, "_", 
                                    Exper, "/SubFlux", j, i, ".csv", sep = ""))[-1]
    
    subbounds = readr::read_csv(paste("./Results/CDiff_", j, "_", i, "_", 
                                      Exper, "/SubBounds", j, i, ".csv", sep = ""))[-1]
    
    assign(paste("subtrace", j, i, sep = ""), subtrace)
    assign(paste("subflux", j, i, sep = ""), subflux)
    assign(paste("subbounds", j, i, sep = ""), subbounds)
  }
}

subtrace = rbind(subtraceAblatedNoDrug, subtraceAblatedTherapy,
                 subtraceUnifiedNoDrug, subtraceUnifiedTherapy)

subflux = rbind(subfluxAblatedNoDrug, subfluxAblatedTherapy,
                subfluxUnifiedNoDrug, subfluxUnifiedTherapy)

subbounds = rbind(subboundsAblatedNoDrug, subboundsAblatedTherapy,
                  subboundsUnifiedNoDrug, subboundsUnifiedTherapy)

#################################################################################

Cdiff = filter(subtrace, subtrace$Places == "CD" & subtrace$Time <= f.Time)
Biomass = filter(subtrace, subtrace$Places == "BiomassCD" & subtrace$Time <= f.Time)
Drug = filter(subtrace, subtrace$Places == "Drug" & subtrace$Time <= f.Time)
IECs = filter(subtrace, subtrace$Places == "IECs" & subtrace$Time <= f.Time)
pheme_e = filter(subtrace, subtrace$Places == "pheme_e" & subtrace$Time <= f.Time)

AA = filter(subtrace, subtrace$Places %in% 
              c("trp_L_e", "pro_L_e", "val_L_e", 
                "ile_L_e", "cys_L_e", "leu_L_e") & 
              subtrace$Time <= f.Time)

FluxBiom = filter(subflux, subflux$Reaction == "EX_biomass_e" & subflux$Time <= f.Time)
FluxPheme = filter(subflux, subflux$Reaction == "sink_pheme_c"& subflux$Time <= f.Time)

FluxAA = filter(subflux, subflux$Reaction %in% 
                  c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e", 
                    "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")
                & subflux$Time <= f.Time)

BoundsBiom = filter(subbounds, subbounds$Reaction %in% 
                      c("EX_biomass_e_Ub", "EX_biomass_e_Lb") 
                    & subbounds$Time <= f.Time)

BoundsAALb = filter(subbounds, subbounds$Reaction %in% 
                      c("EX_trp_L_e_Lb", "EX_pro_L_e_Lb","EX_val_L_e_Lb",
                        "EX_ile_L_e_Lb","EX_cys_L_e_Lb", "EX_leu_L_e_Lb")
                      & subbounds$Time <= f.Time)

BoundsAAUb = filter(subbounds, subbounds$Reaction %in% 
                      c("EX_trp_L_e_Ub", "EX_pro_L_e_Ub","EX_val_L_e_Ub",
                        "EX_ile_L_e_Ub","EX_cys_L_e_Ub", "EX_leu_L_e_Ub")
                    & subbounds$Time <= f.Time)

for (Con in Condition) {
  
  assign(paste("C", Con, sep = "."), 
       ggplot(filter(Cdiff, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
    geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
      colConNoDrug} else {colConTherapy})) +
      labs(x = "Time (h)", y = "C. difficile # cells",
           title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
           subtitle = paste("Place: ", Cdiff[1, 2], " (cells)", sep = "")) +
      theme(axis.text = element_text(size = 14), 
            axis.title=element_text(size = 14),
            plot.title = element_text(size = 16),
            plot.subtitle = element_text(size = 14)) +
      scale_y_continuous(label = scientific_10))
  
  assign(paste("B", Con, sep = "."), 
         ggplot(filter(Biomass, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
           geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "C. difficile average biomass (pg)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Place: ", Biomass[1, 2], " (pg)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)) +
           geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.3) +
           annotate("text", 7, 1.5, label = "biomass max", vjust = 0.5, color = "black", size = 3) +
           geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.3) +
           annotate("text", 7, 0.6, label = "biomass mean", vjust = 0.5, color = "black", size = 3) +
           geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.3) +
           annotate("text", 7, 0.15, label = "biomass min", vjust = 0.5, color = "black", size = 3))
  
  assign(paste("D", Con, sep = "."), 
         ggplot(filter(Drug, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
           geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "C. difficile # cells",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Place: ", Drug[1, 2], " (pmol)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)) +
           annotate("text", 3, 6500, label = "MIC value", vjust = 0.5, color = "black", size = 3))
  
  assign(paste("I", Con, sep = "."), 
         ggplot(filter(IECs, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
           geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "IECs # cells (cell)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Place: ", IECs[1, 2], " (cell)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)) +
           scale_y_continuous(label = scientific_10))
  
  assign(paste("He", Con, sep = "."), 
         ggplot(filter(pheme_e, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
           geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "Extr heme concentration (pmol)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Place: ", pheme_e[1, 2], " (pmol)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)) + 
           annotate("text", 10, 0.95e+06, label = "Hemorrhagic baseline", vjust = 0.5, color = "black", size = 3) +
           scale_y_continuous(label = scientific_10))
  
  assign(paste("A", Con, sep = "."), 
         ggplot(filter(AA, Condition == Con), aes(Time, Marking, colour = Scenario)) +
           theme_bw() + geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "Aminoacids concentration (pmol)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Place: aminoacids", " (molecule)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)) + 
           scale_y_continuous(label = scientific_10) +
           facet_wrap(~Places, nrow = 2))
  
  ################################################################################
  
  assign(paste("fB", Con, sep = "."), 
         ggplot(filter(FluxBiom, Condition == Con), aes(Time, Flux, colour = Scenario)) + theme_bw() +
           geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "EX_biomass(e) FBA flux (mmol/gDW*h)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Reaction: ", FluxBiom[1, 2], " (pg)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)))
  
  assign(paste("fH", Con, sep = "."), 
         ggplot(filter(FluxPheme, Condition == Con), aes(Time, Flux, colour = Scenario)) + theme_bw() +
           geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "sink_pheme(c) FBA flux (mmol/gDW*h)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Reaction: ", FluxPheme[1, 2], " (pmol)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)))
         
  assign(paste("fA", Con, sep = "."), 
         ggplot(filter(FluxAA, Condition == Con), aes(Time, Flux, colour = Scenario)) + 
           theme_bw() + geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "EX_Aminoacid(e) FBA flux (mmol/gDW*h)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Reaction: EX_Aminoacid_e", " (molecule)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)) +
           scale_y_continuous(label = scientific_10) +
           facet_wrap(~Reaction, nrow = 2))
  
  ################################################################################
  
  assign(paste("bB", Con, sep = "."), 
         ggplot(filter(BoundsBiom, Condition == Con), aes(Time, Flux, colour = Scenario, linetype = Reaction)) + 
           theme_bw() + geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "EX_biomass(e) Lb & Ub (mmol/gDW*h)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)))
         
  assign(paste("LbA", Con, sep = "."), 
         ggplot(filter(BoundsAALb, Condition == Con), aes(Time, Flux, colour = Scenario)) + 
           theme_bw() + geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "EX_Aminocid(e) Lb (mmol/gDW*h)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)))
  
  assign(paste("UbA", Con, sep = "."), 
         ggplot(filter(BoundsAAUb, Condition == Con), aes(Time, Flux, colour = Scenario)) + 
           theme_bw() + geom_line(size = size) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
             colConNoDrug} else {colConTherapy})) +
           labs(x = "Time (h)", y = "EX_Aminoacid(e) Ub (mmol/gDW*h)",
                title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
           theme(axis.text = element_text(size = 14), 
                 axis.title=element_text(size = 14),
                 plot.title = element_text(size = 16),
                 plot.subtitle = element_text(size = 14)) +
           facet_wrap(~Reaction, nrow = 2))
  
}

################################################################################

C.All = ggplot(Cdiff, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "C. difficile # cells",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", Cdiff[1, 2], " (cell)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) +
  scale_y_continuous(label = scientific_10)

B.All = ggplot(Biomass, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "C. difficile average biomass (pg)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", Biomass[1, 2], " (pg)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.3) +
  annotate("text", 7, 1.5, label = "biomass max", vjust = 0.5, color = "black", size = 3) +
  geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.3) +
  annotate("text", 7, 0.6, label = "biomass mean", vjust = 0.5, color = "black", size = 3) +
  geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.3) +
  annotate("text", 7, 0.15, label = "biomass min", vjust = 0.5, color = "black", size = 3)

D.All = ggplot(Drug, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "Drug concentration (pmol)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", Drug[1, 2], " (pmol)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  geom_hline(yintercept = 5800, linetype = "dashed", color = "darkgrey", lwd = 0.3) +
  annotate("text", 3, 6500, label = "MIC value", vjust = 0.5, color = "black", size = 3)

I.All = ggplot(IECs, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "IECs # cells (cell)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", IECs[1, 2], " (cell)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  scale_y_continuous(label = scientific_10)

He.All = ggplot(pheme_e, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "Extr heme concentration (pmol)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", pheme_e[1, 2], " (pmol)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) +
  scale_y_continuous(label = scientific_10)

A.All = ggplot(AA, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "Aminoacids concentration (pmol)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: aminoacids", " (molecule)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  scale_y_continuous(label = scientific_10) +
  facet_wrap(~Places, nrow = 2)

################################################################################

fB.All = ggplot(FluxBiom, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "EX_biomass(e) FBA flux (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Reaction: ", FluxBiom[1, 2], " (pg)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14))

fH.All = ggplot(FluxPheme, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "sink_pheme(c) FBA flux (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Reaction: ", FluxPheme[1, 2], " (pmol)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14))

fA.All = ggplot(FluxAA, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "EX_Aminoacid(e) FBA flux (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Reaction: EX_Aminoacid_e", " (molecule)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) +
  facet_wrap(~Reaction, nrow = 2)

################################################################################

bB.All = ggplot(BoundsBiom, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "EX_biomass(e) Bounds (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  facet_wrap(~Reaction, nrow = 1)

LbA.All = ggplot(BoundsAALb, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "EX_Aminoacid(e) Lb (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  facet_wrap(~Reaction, nrow = 2)

UbA.All = ggplot(BoundsAAUb, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
  geom_line(size = size) + scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Time (h)", y = "EX_Aminoacid(e) Ub (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  facet_wrap(~Reaction, nrow = 2)

################################################################################

C.Unified = ggplot(filter(Cdiff, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
  theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "C. difficile # cells",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", Cdiff[1, 2], " (cells)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  scale_y_continuous(label = scientific_10)

B.Unified = ggplot(filter(Biomass, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
  theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "C. difficile average biomass (pg)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", Biomass[1, 2], " (pg)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14))

D.Unified = ggplot(filter(Drug, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
  theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "Drug concentration (pmol)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", Drug[1, 2], " (pmol)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14))

I.Unified = ggplot(filter(IECs, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
  theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "IECs # cells (cell)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", IECs[1, 2], " (cell)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) +
  scale_y_continuous(label = scientific_10)

He.Unified = ggplot(filter(pheme_e, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
  theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "Extr heme concentration (pmol)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: ", pheme_e[1, 2], " (pmol)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) +
  scale_y_continuous(label = scientific_10)

A.Unified = ggplot(filter(AA, Scenario == "Unified"), 
                   aes(Time, (Marking/(6.022e20/6.022e08))*1e+09, colour = Condition)) + 
  theme_bw() + geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "Aminoacids concentration (pmol)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Place: aminoacids", " (molecule)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) +
  scale_y_continuous(label = scientific_10) +
  facet_wrap(~Places, nrow = 2)

################################################################################

fB.Unified = ggplot(filter(FluxBiom, Scenario == "Unified"), aes(Time, Flux, colour = Condition)) + 
  theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "EX_biomass(e) FBA flux (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Reaction: ", FluxBiom[1, 2], " (pg)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14))

fH.Unified = ggplot(filter(FluxPheme, Scenario == "Unified"), aes(Time, Flux, colour = Condition)) + 
  theme_bw() +
  geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "sink_pheme(c) FBA flux (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Reaction: ", FluxPheme[1, 2], " (pmol)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14))

fA.Unified = ggplot(filter(FluxAA, Scenario == "Unified"), 
                    aes(Time, Flux, colour = Condition)) + 
  theme_bw() + geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "EX_Aminoacid(e) FBA flux (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Reaction: EX_Aminoacid_e", " (molecule)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) + 
  facet_wrap(~Reaction, nrow = 2)

################################################################################

bB.Unified = ggplot(filter(BoundsBiom, Scenario == "Unified"), 
                    aes(Time, Flux, colour = Condition, linetype = Reaction)) + 
  theme_bw() + geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "EX_biomass(e) Lb & Ub (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14))

LbA.Unified = ggplot(filter(BoundsAALb, Scenario == "Unified"), 
                     aes(Time, Flux, colour = Condition)) + theme_bw() +
  geom_line(size = size) + scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "EX_Aminoacid(e) Lb (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) 
  facet_wrap(~Reaction, nrow = 2)

UbA.Unified = ggplot(filter(BoundsAAUb, Scenario == "Unified"), 
                     aes(Time, Flux, colour = Condition)) + 
  theme_bw() + geom_line(size = size) + 
  scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  labs(x = "Time (h)", y = "EX_Aminoacid(e) Ub (mmol/gDW*h)",
       title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
       subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
  theme(axis.text = element_text(size = 14), 
        axis.title=element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14)) +
  facet_wrap(~Reaction, nrow = 2)

###################################################################################
###################################################################################

ggsave(((C.NoDrug|C.Therapy)/(C.Unified|C.All)), 
       file = paste("AblationAnalysis_CD_zoom_", f.Time, ".pdf", sep = ""), 
       width = 15, height = 8)
ggsave(((B.NoDrug|B.Therapy)/(B.Unified|B.All)), 
       file = paste("AblationAnalysis_BiomassCD_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((D.NoDrug|D.Therapy)/(D.Unified|D.All)), 
       file = paste("AblationAnalysis_Drug_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((I.NoDrug|I.Therapy)/(I.Unified|I.All)), 
       file = paste("AblationAnalysis_IECs_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((He.NoDrug|He.Therapy)/(He.Unified|He.All)), 
       file = paste("AblationAnalysis_pheme_e_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((A.NoDrug|A.Therapy)/(A.Unified|A.All)), 
       file = paste("AblationAnalysis_AA_zoom_", f.Time, ".pdf", sep = ""), 
       width = 15, height = 8)
ggsave(((fB.NoDrug|fB.Therapy)/(fB.Unified|fB.All)), 
       file = paste("AblationAnalysis_EX_biomass_e_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((fH.NoDrug|fH.Therapy)/(fH.Unified|fH.All)), 
       file = paste("AblationAnalysis_sink_pheme_c.pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((fA.NoDrug|fA.Therapy)/(fA.Unified|fA.All)), 
       file = paste("AblationAnalysis_EX_aminocids_e_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((bB.NoDrug|bB.Therapy)/(bB.Unified|bB.All)), 
       file = paste("AblationAnalysis_EX_biomass_Lb_Ub_e_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((LbA.NoDrug|LbA.Therapy)/(LbA.Unified|LbA.All)), 
       file = paste("AblationAnalysis_EX_aminocids_e_Lb_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)
ggsave(((UbA.NoDrug|UbA.Therapy)/(UbA.Unified|UbA.All)), 
       file = paste("AblationAnalysis_EX_aminocids_e_Ub_zoom_", f.Time, ".pdf", sep = ""), 
       width = 14, height = 8)

system(paste("cp AblationAnalysis_*.pdf ", "./Results/images4paper/", sep = ""))
system(paste("mv AblationAnalysis_*.pdf ", "./Results/images4paper/", sep = ""))
