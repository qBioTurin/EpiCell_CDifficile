
library(patchwork)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("./Results")

AnalysisPlaces_Drug0_hemeHigh <- readRDS("./AnalysisPlaces_Drug0_hemeHigh.rds")
AnalysisPlaces_Drug0_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisPlaces_Drug0_hemeHigh$Time))
AnalysisPlaces_Drug0_hemeHigh$DrugCondition = rep("0", times = length(AnalysisPlaces_Drug0_hemeHigh$Time))

AnalysisPlaces_Drug0_hemeLow <- readRDS("./AnalysisPlaces_Drug0_hemeLow.rds")
AnalysisPlaces_Drug0_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisPlaces_Drug0_hemeLow$Time))
AnalysisPlaces_Drug0_hemeLow$DrugCondition = rep("0", times = length(AnalysisPlaces_Drug0_hemeLow$Time))

AnalysisPlaces_Drug005_hemeHigh <- readRDS("./AnalysisPlaces_Drug005_hemeHigh.rds")
AnalysisPlaces_Drug005_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisPlaces_Drug005_hemeHigh$Time))
AnalysisPlaces_Drug005_hemeHigh$DrugCondition = rep("0.05", times = length(AnalysisPlaces_Drug005_hemeHigh$Time))

AnalysisPlaces_Drug005_hemeLow <- readRDS("./AnalysisPlaces_Drug005_hemeLow.rds")
AnalysisPlaces_Drug005_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisPlaces_Drug005_hemeLow$Time))
AnalysisPlaces_Drug005_hemeLow$DrugCondition = rep("0.05", times = length(AnalysisPlaces_Drug005_hemeLow$Time))

AnalysisPlaces_Drug01_hemeHigh <- readRDS("./AnalysisPlaces_Drug01_hemeHigh.rds")
AnalysisPlaces_Drug01_hemeHigh$HemeCondition = rep("hemeHigh (g/day)", times = length(AnalysisPlaces_Drug01_hemeHigh$Time))
AnalysisPlaces_Drug01_hemeHigh$DrugCondition = rep("0.1", times = length(AnalysisPlaces_Drug01_hemeHigh$Time))

AnalysisPlaces_Drug01_hemeLow <- readRDS("./AnalysisPlaces_Drug01_hemeLow.rds")
AnalysisPlaces_Drug01_hemeLow$HemeCondition = rep("Heme Poor (1 g/day)", times = length(AnalysisPlaces_Drug01_hemeLow$Time))
AnalysisPlaces_Drug01_hemeLow$DrugCondition = rep("0.1", times = length(AnalysisPlaces_Drug01_hemeLow$Time))

AnalysisPlaces_Drug03_hemeHigh <- readRDS("./AnalysisPlaces_Drug03_hemeHigh.rds")
AnalysisPlaces_Drug03_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisPlaces_Drug03_hemeHigh$Time))
AnalysisPlaces_Drug03_hemeHigh$DrugCondition = rep("0.3", times = length(AnalysisPlaces_Drug03_hemeHigh$Time))

AnalysisPlaces_Drug03_hemeLow <- readRDS("./AnalysisPlaces_Drug03_hemeLow.rds")
AnalysisPlaces_Drug03_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisPlaces_Drug03_hemeLow$Time))
AnalysisPlaces_Drug03_hemeLow$DrugCondition = rep("0.3", times = length(AnalysisPlaces_Drug03_hemeLow$Time))

AnalysisPlaces_Drug05_hemeHigh <- readRDS("./AnalysisPlaces_Drug05_hemeHigh.rds")
AnalysisPlaces_Drug05_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisPlaces_Drug05_hemeHigh$Time))
AnalysisPlaces_Drug05_hemeHigh$DrugCondition = rep("0.5", times = length(AnalysisPlaces_Drug05_hemeHigh$Time))

AnalysisPlaces_Drug05_hemeLow <- readRDS("./AnalysisPlaces_Drug05_hemeLow.rds")
AnalysisPlaces_Drug05_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisPlaces_Drug05_hemeLow$Time))
AnalysisPlaces_Drug05_hemeLow$DrugCondition = rep("0.5", times = length(AnalysisPlaces_Drug05_hemeLow$Time))

AnalysisPlaces_Drug10_hemeHigh <- readRDS("./AnalysisPlaces_Drug10_hemeHigh.rds")
AnalysisPlaces_Drug10_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisPlaces_Drug10_hemeHigh$Time))
AnalysisPlaces_Drug10_hemeHigh$DrugCondition = rep("1", times = length(AnalysisPlaces_Drug10_hemeHigh$Time))

AnalysisPlaces_Drug10_hemeLow <- readRDS("./AnalysisPlaces_Drug10_hemeLow.rds")
AnalysisPlaces_Drug10_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisPlaces_Drug10_hemeLow$Time))
AnalysisPlaces_Drug10_hemeLow$DrugCondition = rep("1", times = length(AnalysisPlaces_Drug10_hemeLow$Time))

Places = rbind(AnalysisPlaces_Drug0_hemeLow, AnalysisPlaces_Drug03_hemeLow, AnalysisPlaces_Drug05_hemeLow, 
               AnalysisPlaces_Drug10_hemeLow, AnalysisPlaces_Drug0_hemeHigh, AnalysisPlaces_Drug03_hemeHigh,
               AnalysisPlaces_Drug05_hemeHigh, AnalysisPlaces_Drug10_hemeHigh)

Places$Time = as.numeric(as.character(Places$Time))
Places$CD = as.numeric(as.character(Places$CD))

colnames(Places) = c("Time [h]", "IECs [cell]", "Heme [molecule]", "BiomassCD", 
                     "pro_L_e", "leu_L_e", "ile_L_e", "val_L_e", "trp_L_e", "cys_L_e", 
                     "Drug [molecule]", "Damage", "pro_L_v", "leu_L_v", "ile_L_v", "val_L_v", 
                     "trp_L_v", "cys_L_v", "CD  [cell]", "HemeCondition", "DrugCondition")

dfPlaces = Places %>% select(`Time [h]`, `CD  [cell]`, `Drug [molecule]`, 
                             `Heme [molecule]`, `IECs [cell]`, DrugCondition, HemeCondition) %>%
  gather(`CD  [cell]`, `Drug [molecule]`, `Heme [molecule]`, `IECs [cell]`, key = "Places", value = "Marking")

colnames(dfPlaces) = c("Time [h]", "Drug Dose [g/day]", "HemeCondition", "Places", "Tokens")

log10pP = ggplot(dfPlaces, aes(x = `Time [h]`, y = log10(Tokens), group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Places~HemeCondition, scales = "free") + theme_bw() + ylab("log10 (Amount)") +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

pPCD = ggplot(dfPlaces[which(dfPlaces$Places == "CD  [cell]"), ], aes(x = `Time [h]`, y = Tokens, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Places~HemeCondition, scales = "free") + theme_bw() + ylab("Cells") +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

log10pPCD = ggplot(dfPlaces[which(dfPlaces$Places == "CD  [cell]"), ], aes(x = `Time [h]`, y = log10(Tokens), group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Places~HemeCondition, scales = "free") + theme_bw() + ylab("log10 (Cells)") +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

vec = c("Time", 
        "EX_biomass(e)", "EX_pheme(e)", "EX_cys_L(e)", 
        "EX_trp_L(e)", "EX_ile_L(e)", "EX_pro_L(e)", "EX_leu_L(e)", "EX_val_L(e)")
  
AnalysisFlux_Drug0_hemeHigh <- readRDS("./AnalysisFlux_Drug0_hemeHigh.rds")
AnalysisFlux_Drug0_hemeHigh = AnalysisFlux_Drug0_hemeHigh[, which(colnames(AnalysisFlux_Drug0_hemeHigh) %in% vec)]
AnalysisFlux_Drug0_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisFlux_Drug0_hemeHigh$Time))
AnalysisFlux_Drug0_hemeHigh$DrugCondition = rep("0", times = length(AnalysisFlux_Drug0_hemeHigh$Time))

AnalysisFlux_Drug0_hemeLow <- readRDS("./AnalysisFlux_Drug0_hemeLow.rds")
AnalysisFlux_Drug0_hemeLow = AnalysisFlux_Drug0_hemeLow[, which(colnames(AnalysisFlux_Drug0_hemeLow) %in% vec)]
AnalysisFlux_Drug0_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisFlux_Drug0_hemeLow$Time))
AnalysisFlux_Drug0_hemeLow$DrugCondition = rep("0", times = length(AnalysisFlux_Drug0_hemeLow$Time))

AnalysisFlux_Drug005_hemeHigh <- readRDS("./AnalysisFlux_Drug005_hemeHigh.rds")
AnalysisFlux_Drug005_hemeHigh = AnalysisFlux_Drug005_hemeHigh[, which(colnames(AnalysisFlux_Drug005_hemeHigh) %in% vec)]
AnalysisFlux_Drug005_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisFlux_Drug005_hemeHigh$Time))
AnalysisFlux_Drug005_hemeHigh$DrugCondition = rep("0.05", times = length(AnalysisFlux_Drug005_hemeHigh$Time))

AnalysisFlux_Drug005_hemeLow <- readRDS("./AnalysisFlux_Drug005_hemeLow.rds")
AnalysisFlux_Drug005_hemeLow = AnalysisFlux_Drug005_hemeLow[, which(colnames(AnalysisFlux_Drug005_hemeLow) %in% vec)]
AnalysisFlux_Drug005_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisFlux_Drug005_hemeLow$Time))
AnalysisFlux_Drug005_hemeLow$DrugCondition = rep("0.05", times = length(AnalysisFlux_Drug005_hemeLow$Time))

AnalysisFlux_Drug01_hemeHigh <- readRDS("./AnalysisFlux_Drug01_hemeHigh.rds")
AnalysisFlux_Drug01_hemeHigh = AnalysisFlux_Drug01_hemeHigh[, which(colnames(AnalysisFlux_Drug01_hemeHigh) %in% vec)]
AnalysisFlux_Drug01_hemeHigh$HemeCondition = rep("hemeHigh (g/day)", times = length(AnalysisFlux_Drug01_hemeHigh$Time))
AnalysisFlux_Drug01_hemeHigh$DrugCondition = rep("0.1", times = length(AnalysisFlux_Drug01_hemeHigh$Time))

AnalysisFlux_Drug01_hemeLow <- readRDS("./AnalysisFlux_Drug01_hemeLow.rds")
AnalysisFlux_Drug01_hemeLow = AnalysisFlux_Drug01_hemeLow[, which(colnames(AnalysisFlux_Drug01_hemeLow) %in% vec)]
AnalysisFlux_Drug01_hemeLow$HemeCondition = rep("Heme Poor (1 g/day)", times = length(AnalysisFlux_Drug01_hemeLow$Time))
AnalysisFlux_Drug01_hemeLow$DrugCondition = rep("0.1", times = length(AnalysisFlux_Drug01_hemeLow$Time))

AnalysisFlux_Drug03_hemeHigh <- readRDS("./AnalysisFlux_Drug03_hemeHigh.rds")
AnalysisFlux_Drug03_hemeHigh = AnalysisFlux_Drug03_hemeHigh[, which(colnames(AnalysisFlux_Drug03_hemeHigh) %in% vec)]
AnalysisFlux_Drug03_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisFlux_Drug03_hemeHigh$Time))
AnalysisFlux_Drug03_hemeHigh$DrugCondition = rep("0.3", times = length(AnalysisFlux_Drug03_hemeHigh$Time))

AnalysisFlux_Drug03_hemeLow <- readRDS("./AnalysisFlux_Drug03_hemeLow.rds")
AnalysisFlux_Drug03_hemeLow = AnalysisFlux_Drug03_hemeLow[, which(colnames(AnalysisFlux_Drug03_hemeLow) %in% vec)]
AnalysisFlux_Drug03_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisFlux_Drug03_hemeLow$Time))
AnalysisFlux_Drug03_hemeLow$DrugCondition = rep("0.3", times = length(AnalysisFlux_Drug03_hemeLow$Time))

AnalysisFlux_Drug05_hemeHigh <- readRDS("./AnalysisFlux_Drug05_hemeHigh.rds")
AnalysisFlux_Drug05_hemeHigh = AnalysisFlux_Drug05_hemeHigh[, which(colnames(AnalysisFlux_Drug05_hemeHigh) %in% vec)]
AnalysisFlux_Drug05_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisFlux_Drug05_hemeHigh$Time))
AnalysisFlux_Drug05_hemeHigh$DrugCondition = rep("0.5", times = length(AnalysisFlux_Drug05_hemeHigh$Time))

AnalysisFlux_Drug05_hemeLow <- readRDS("./AnalysisFlux_Drug05_hemeLow.rds")
AnalysisFlux_Drug05_hemeLow = AnalysisFlux_Drug05_hemeLow[, which(colnames(AnalysisFlux_Drug05_hemeLow) %in% vec)]
AnalysisFlux_Drug05_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisFlux_Drug05_hemeLow$Time))
AnalysisFlux_Drug05_hemeLow$DrugCondition = rep("0.5", times = length(AnalysisFlux_Drug05_hemeLow$Time))

AnalysisFlux_Drug10_hemeHigh <- readRDS("./AnalysisFlux_Drug10_hemeHigh.rds")
AnalysisFlux_Drug10_hemeHigh = AnalysisFlux_Drug10_hemeHigh[, which(colnames(AnalysisFlux_Drug10_hemeHigh) %in% vec)]
AnalysisFlux_Drug10_hemeHigh$HemeCondition = rep("Heme Rich (100 g/day)", times = length(AnalysisFlux_Drug10_hemeHigh$Time))
AnalysisFlux_Drug10_hemeHigh$DrugCondition = rep("1", times = length(AnalysisFlux_Drug10_hemeHigh$Time))

AnalysisFlux_Drug10_hemeLow <- readRDS("./AnalysisFlux_Drug10_hemeLow.rds")
AnalysisFlux_Drug10_hemeLow = AnalysisFlux_Drug10_hemeLow[, which(colnames(AnalysisFlux_Drug10_hemeLow) %in% vec)]
AnalysisFlux_Drug10_hemeLow$HemeCondition = rep("Heme Poor (1 g/day )", times = length(AnalysisFlux_Drug10_hemeLow$Time))
AnalysisFlux_Drug10_hemeLow$DrugCondition = rep("1", times = length(AnalysisFlux_Drug10_hemeLow$Time))

Fluxes = rbind(AnalysisFlux_Drug0_hemeHigh, AnalysisFlux_Drug0_hemeLow,
               AnalysisFlux_Drug03_hemeHigh, AnalysisFlux_Drug03_hemeLow, 
               AnalysisFlux_Drug05_hemeHigh, AnalysisFlux_Drug05_hemeLow, 
               AnalysisFlux_Drug10_hemeHigh, AnalysisFlux_Drug10_hemeLow)

Fluxes$Time = as.numeric(as.character(Fluxes$Time))
Fluxes$`EX_pheme(e)` = as.numeric(as.character(Fluxes$`EX_pheme(e)`))
Fluxes$`EX_cys_L(e)` = as.numeric(as.character(Fluxes$`EX_cys_L(e)`))
Fluxes$`EX_biomass(e)` = as.numeric(as.character(Fluxes$`EX_biomass(e)`))
Fluxes$`EX_ile_L(e)` = as.numeric(as.character(Fluxes$`EX_ile_L(e)`))
Fluxes$`EX_pro_L(e)` = as.numeric(as.character(Fluxes$`EX_pro_L(e)`))
Fluxes$`EX_trp_L(e)` = as.numeric(as.character(Fluxes$`EX_trp_L(e)`))
Fluxes$`EX_val_L(e)` = as.numeric(as.character(Fluxes$`EX_val_L(e)`))
Fluxes$`EX_leu_L(e)` = as.numeric(as.character(Fluxes$`EX_leu_L(e)`))

colnames(Fluxes) = c("Time [h]", "EX_pheme(e)", "EX_cys_L(e)", 
                     "EX_ile_L(e)", "EX_biomass(e)", "EX_pro_L(e)",
                     "EX_trp_L(e)", "EX_val_L(e)", "EX_leu_L(e)",
                     "HemeCondition", "DrugCondition")

dfFluxes = Fluxes %>% select(`Time [h]`, `EX_pheme(e)`, `EX_cys_L(e)`, 
                             `EX_biomass(e)`, `EX_ile_L(e)`, 
                             `EX_pro_L(e)`, `EX_trp_L(e)`, `EX_val_L(e)`,
                             `EX_leu_L(e)`,
                             HemeCondition, DrugCondition) %>%
  gather(`EX_pheme(e)`, `EX_cys_L(e)`, 
         `EX_biomass(e)`, `EX_ile_L(e)`, 
         `EX_pro_L(e)`, `EX_trp_L(e)`, `EX_val_L(e)`,
         `EX_leu_L(e)`, key = "Reactions", value = "Fluxes")

colnames(dfFluxes) = c("Time [h]", "HemeCondition", "Drug Dose [g/day]", "Reactions", "Fluxes [mmol/cell*h]")

dftime = dfFluxes %>% group_by(HemeCondition, `Drug Dose [g/day]`, Reactions) %>% 
  filter(`Time [h]` == max(`Time [h]`)) %>% 
  filter(`Time [h]` < 6) %>% 
  distinct() %>% 
  summarise(`Fluxes [mmol/cell*h]` = mean(`Fluxes [mmol/cell*h]`), `Time [h]` = 6)

dfFluxes = rbind(dfFluxes, dftime)

# pF = ggplot(dfFluxes, aes(x = `Time [h]`, y = Fluxes, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  #geom_line(size=0.75) + facet_grid(Fluxes~HemeCondition, scales = "free") + theme_bw() + ylab("Fluxes") +
  #theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  # geom_vline(xintercept = 0, linetype="dotted")

p1 = ggplot(dfFluxes[which(dfFluxes$Reactions == "EX_biomass(e)"), ], 
       aes(x = `Time [h]`, y = signif(`Fluxes [mmol/cell*h]`, digits = 2), group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Reactions~HemeCondition, scales = "free") + theme_bw() +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

p2 = ggplot(dfFluxes[which(dfFluxes$Reactions == "EX_pheme(e)"), ], 
       aes(x = `Time [h]`, y =  `Fluxes [mmol/cell*h]`, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.40) + facet_grid(Reactions~HemeCondition, scales = "free") + theme_bw() +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

ggplot2::ggsave(p2 + ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png", bg = "transparent",  width = 5.7, height =1.6, dpi = 1000, units = "in", device='png')

p3 = ggplot(dfFluxes[which(dfFluxes$Reactions == "EX_cys_L(e)"), ], 
       aes(x = `Time [h]`, y = `Fluxes [mmol/cell*h]`, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.40) + facet_grid(Reactions~HemeCondition, scales = "free") + theme_bw() +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

ggplot2::ggsave(p3 + ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png", bg = "transparent",  width = 5.5, height =1.6, dpi = 1000, units = "in", device='png')

p4 = ggplot(dfFluxes[which(dfFluxes$Reactions == "EX_ile_L(e)"), ], 
       aes(x = `Time [h]`, y = `Fluxes [mmol/cell*h]`, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Reactions~HemeCondition, scales = "free") + theme_bw() +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

p5 = ggplot(dfFluxes[which(dfFluxes$Reactions == "EX_pro_L(e)"), ], 
       aes(x = `Time [h]`, y = `Fluxes [mmol/cell*h]`, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Reactions~HemeCondition, scales = "free") + theme_bw() +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

ggplot2::ggsave(p3 + ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png", bg = "transparent",  width = 5.5, height =1.6, dpi = 1000, units = "in", device='png')

p6 = ggplot(dfFluxes[which(dfFluxes$Reactions == "EX_trp_L(e)"), ], 
       aes(x = `Time [h]`, y = `Fluxes [mmol/cell*h]`, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Reactions~HemeCondition, scales = "free") + theme_bw() +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

p7 = ggplot(dfFluxes[which(dfFluxes$Reactions == "EX_val_L(e)"), ], 
            aes(x = `Time [h]`, y = `Fluxes [mmol/cell*h]`, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Reactions~HemeCondition, scales = "free") + theme_bw() +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

p8 = ggplot(dfFluxes[which(dfFluxes$Reactions == "EX_leu_L(e)"), ], 
            aes(x = `Time [h]`, y = `Fluxes [mmol/cell*h]`, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.75) + facet_grid(Reactions~HemeCondition, scales = "free") + theme_bw() +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")

(pPCD | log10pP)

pPCD = ggplot(dfPlaces[which(dfPlaces$Places == "CD  [cell]"), ], aes(x = `Time [h]`, y = Tokens, group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.40) + facet_grid(Places~HemeCondition, scales = "free") + theme_bw() + ylab("Cells") +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position = 'right', 
        legend.direction = "horizontal", 
        legend.box = "horizontal") +
  guides(color = guide_legend(title.position = "top", 
                              # hjust = 0.5 centres the title horizontally
                              title.hjust = 0.5,
                              label.position = "bottom")) 

ggplot2::ggsave(pPCD + ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png", bg = "transparent",  width = 5.5, height =1.8, dpi = 1000, units = "in", device='png')

log10pPnoCD = ggplot(dfPlaces[which(!dfPlaces$Places == "CD  [cell]"), ], aes(x = `Time [h]`, y = log10(Tokens), group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.40) + facet_grid(Places~HemeCondition, scales = "free") + theme_bw() + ylab("log10 (Amount)") +
  theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted")+
  theme(legend.position = 'right', 
        legend.direction = "horizontal", 
        legend.box = "horizontal") +
  guides(color = guide_legend(title.position = "top", 
                              # hjust = 0.5 centres the title horizontally
                              title.hjust = 0.5,
                              label.position = "bottom")) 

ggplot2::ggsave(log10pPnoCD + ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png", bg = "transparent",  width = 5.3, height = 4.2, dpi = 1000, units = "in", device='png')

dfFluxes$`Fluxes [mmol/cell*h]` = signif(dfFluxes$`Fluxes [mmol/cell*h]`, digits=2)

pflux = ggplot(dfFluxes[which(dfFluxes$Reactions %in% c("EX_pheme(e)", "EX_cys_L(e)")), ], 
       aes(x = `Time [h]`, y = `Fluxes [mmol/cell*h]`, 
           group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.40) + facet_grid(Reactions~HemeCondition, scales = "free") + 
  theme_bw() + ylab("Fluxes") + theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted") + ylab("Fluxes [mmol/cell*h]") + theme(legend.position = "none")

ggplot2::ggsave(pflux + ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png", bg = "transparent",  width = 4.7, height =2.8, dpi = 1000, units = "in", device='png')

pflux = ggplot(dfFluxes[which(dfFluxes$Reactions %in% c("EX_cys_L(e)")), ], 
               aes(x = `Time [h]`, y = `Fluxes [mmol/cell*h]`, 
                   group = `Drug Dose [g/day]`, color = `Drug Dose [g/day]`))  +
  geom_line(size=0.40) + facet_grid(Reactions~HemeCondition, scales = "free") + 
  theme_bw() + ylab("Fluxes") + theme(panel.background = element_rect(fill = "#F8F8FF")) + 
  geom_vline(xintercept = 0, linetype="dotted") + ylab("Fluxes [mmol/cell*h]")

ggplot2::ggsave(pflux + ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png", bg = "transparent",  width = 6, height = 3.5, dpi = 1000, units = "in", device='png')


signif(2.60699382414341e-56, digits=3)
p1 + p2
  
((p1 + p2 + p3 + p4) / (p5 + p6 + p7 + p8))


ggplot2::ggsave(plot +ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "PlotTraces.png", 
                bg = "transparent", width = 25, height = 8, 
                dpi = 150, units = "in", device='png')
