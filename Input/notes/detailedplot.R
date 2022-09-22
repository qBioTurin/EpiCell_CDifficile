
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("~/EpiCell_CDifficile/Results/EpitCellDifficileHemeSink_sensitivityOnlyDeath/")
load("./EpitCellDifficileHemeSink-sensitivity.RData")

datalist = lapply(list.files(pattern = "\\.trace$") ,function(i){
  TracesPath = i
  n = as.numeric(gsub(pattern = "(EpitCellDifficileHemeSink-sensitivity-)|(.trace)",replacement = "",x = i) )
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
  
  subflux = trace[, c("Time", FluxVec)] %>% 
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  subtrace = trace %>% dplyr::select(Time, IECs, Drug, Damage,
                                     pheme_e, pheme_c, 
                                     BiomassCD, CD,
                                     pro_L_e, leu_L_e, ile_L_e, 
                                     val_L_e, trp_L_e, cys_L_e) %>% 
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  subtrace = cbind(subtrace, 
                   config = rep(as.character(i), length(subtrace$Time)),
                   Efflux = rep(config[[19]][[n]][[3]], length(subtrace$Time)),
                   Death4Treat = rep(config[[17]][[n]][[3]], length(subtrace$Time)))
  
  subflux = cbind(subflux, config = rep(as.character(i), length(subflux$Time)),
                  Efflux = rep(config[[19]][[n]][[3]], length(subflux$Time)),
                  Death4Treat = rep(config[[17]][[n]][[3]], length(subflux$Time)))
  
  df = rbind(subtrace, subflux)
  
  return(df)
})

subtrace = do.call(rbind, datalist)

setwd("~/EpiCell_CDifficile/Results")

grDevices::pdf(file = "./IntermediateDrugEfficacy.pdf",
               width = 16, height = 6)

hist = ggplot(subtrace %>% filter(Places %in% c("CD","Drug") & Time == max(Time)), aes(Marking)) +
  geom_histogram() + ggtitle("CD [cell] & Drug [pmol/mL]") + 
  facet_wrap(~Places, scales = "free") + geom_histogram(bins=30)

allCD = ggplot(subtrace[which(subtrace$Places == "CD"), ], 
               aes(x = Time, y = Marking, color = Death4Treat )) + 
  geom_line(aes(group = config)) + ggtitle("CD [cell]")

allDrug = ggplot(subtrace[which(subtrace$Places == "Drug"), ], 
                 aes(x = Time, y = Marking, color = Death4Treat)) + 
  geom_hline(yintercept = 5.8e-06*1e+09, linetype="dotted", color = "darkblue") +
  geom_line(aes(group = config)) + ggtitle("Drug [pmol]")

cd = filter(subtrace, subtrace$Places == "CD")
filtered = filter(cd, cd$Death4Treat < 0.0883 & cd$Death4Treat > 0.0878)

filtered$Death4Treat = as.character(filtered$Death4Treat)

for (i in 1:length(filtered$config)) {
  TracesPath = i
  n = gsub(pattern = "(EpitCellDifficileHemeSink-sensitivity-)|(.trace)", 
           replacement = letters[i], filtered$config[i])
  filtered$config[i] = n
}

# png("myplot.png", width=250, height=150, bg = "transparent")
# ggplot(filtered, aes(x = Time, y = Marking, color = Death4Treat)) +
#   geom_line(size=1.25) +
#   theme(legend.position = "none", 
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank()) +
#   geom_vline(xintercept = 8, linetype="dotted")
# dev.off()

filteredplot = ggplot(filtered, aes(x = Time, y = Marking, color = Death4Treat)) +
  geom_line(size=0.75) + theme_bw() + 
  theme(panel.background = element_rect(fill = "#F8F8FF")) +
  geom_vline(xintercept = 8, linetype="dotted") +
  ggtitle("Antibiotic Efficacy") +
  xlab("Time [h]") + ylab("CD [cell]")

pflux = ggplot(subtrace[which(subtrace$Places %in% 
                                c("BiomassCD", "EX_biomass_e",
                                  "EX_pro_L_e", "EX_leu_L_e", "EX_ile_L_e",
                                  "EX_val_L_e", "EX_trp_L_e", "EX_cys_L_e",
                                  "sink_pheme_c")), ],
               aes(x = Time, y = Marking, color = Death4Treat)) +
  geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
  ggtitle("Biomass [pg] & Fluxes [mmol/gDW*h]")

((allCD + allDrug) / filteredplot) | pflux

dev.off()

################################################################################
######################### Drug efficacy comparison #############################
################################################################################

grDevices::pdf(file = "Comparison.pdf",
               width = 12, height = 12)

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

heme = filter(subtrace, grepl("pheme_e|pheme_c", Places))
cell = filter(subtrace, grepl("IECs|Damage", Places))
for (i in 1:length(cell$Places)) {if(cell$Places[i] == "Damage") {cell$Places[i] = "IECs killed"}}
biom = filter(subtrace, subtrace$Places == "BiomassCD")
clo = filter(subtrace, subtrace$Places == "CD")
drug = filter(subtrace, subtrace$Places == "Drug")

cd = ggplot(clo, aes(x = Time, y = signif(Marking, digits = 3), color = Scenario)) +
  geom_line(size=0.75) + theme_bw() + 
  theme(panel.background = element_rect(fill = "#F8F8FF")) +
  geom_vline(xintercept = 8, linetype="dotted") + 
  ggtitle("C. difficile cells depending treatment") +
  xlab("Time [h]") + ylab("CD [cell]") + 
  scale_color_manual(values=c('darkred','darkblue'))

biom = ggplot(biom, aes(x = Time, y = signif(Marking, digits = 3), color = Scenario)) +
  geom_line(size=0.75) + theme_bw() + 
  geom_vline(xintercept = 8, linetype="dotted") + 
  ggtitle("C. difficile cells average biomass") + 
  geom_hline(yintercept = 1.57872, linetype="dotted", color = "magenta") + 
  geom_hline(yintercept = 0.496328, linetype="dotted", color = "magenta") +
  geom_hline(yintercept = 0.06738716, linetype="dotted",  color = "magenta") + 
  xlab("Time [h]") + ylab("Biomass [pg/cell]") +
  scale_color_manual(values=c('darkred','darkblue'))

drug = ggplot(drug, aes(x = Time, y = signif(Marking*1e-09, digits = 3), color = Scenario)) +
  geom_line(size=0.75) + theme_bw() + 
  geom_vline(xintercept = 8, linetype="dotted") + 
  geom_hline(yintercept = 5.8e-06, linetype="dotted", color = "darkblue") +
  ggtitle("Drug concentration") +
  xlab("Time [h]") + ylab("Drug [mmol/mL]") + 
  scale_color_manual(values=c('darkred','darkblue'))

cell = ggplot(cell, aes(Time, Marking, colour = Scenario, linetype = Places)) + 
  geom_line() + geom_point() + theme_bw() + 
  geom_vline(xintercept = 8, linetype="dotted") +
  ggtitle("Human Components") +
  xlab("Time [h]") + ylab("EICs [cell]") +
  labs(color  = "Conditions", linetype = "Places") +
  scale_color_manual(values=c('darkblue','darkred'))

heme = ggplot(heme, aes(Time, Marking, colour = Scenario, linetype = Places)) + 
  geom_line(size=1) + theme_bw() +
  geom_vline(xintercept = 8, linetype="dotted") +
  geom_hline(yintercept = 0.01*1e-03*1e+09, linetype="dotted", color = "darkred") + 
  ggtitle("Heme Trajectories") +
  xlab("Time [h]") + ylab("Heme [pmol/mL]") +
  labs(color  = "Scenario", linetype = "Places") +
  scale_color_manual(values=c('darkred', 'darkblue'))

f = ggplot(subflux, aes(Time, Flux, colour = Scenario, linetype = Reaction)) + 
  geom_line(size = 0.75) + 
  geom_vline(xintercept = 8, linetype="dotted") +
  ggtitle("FBA Estimated Fluxes") +
  xlab("Time [h]") + ylab("Flux [mmol/gDW*h]") +
  labs(color  = "Scenario", linetype = "Reaction") +
  scale_color_manual(values=c('darkred', 'darkblue'))

((cd + biom + drug) / (cell + heme)) / f

dev.off()
