
library(dplyr)
library(ggplot2)

################################################################################
# Drug metabolic indirect effects assessed by estimated reaction fluxes values #
################################################################################

setwd("~/EpiCell_CDifficile/Results/")

FluxVec = c("EX_biomass_e", "sink_pheme_c",
            "EX_pro_L_e", "EX_leu_L_e", "EX_ile_L_e",
            "EX_val_L_e", "EX_trp_L_e", "EX_cys_L_e")

conditions = c("NoDrug", "LowEff")

for (con in conditions) {
  
  trace = read.table(paste("EpitCellDifficileHemeSink_analysis", con, 
                           "/EpitCellDifficileHemeSink-analysis-1.trace", 
                           sep = ""), header = F)
  
  ReactionsNames = trace[1, ]
  ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
  ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
  
  trace = read.table(paste("EpitCellDifficileHemeSink_analysis", con, 
                           "/EpitCellDifficileHemeSink-analysis-1.trace", 
                           sep = ""), 
                     header = T)
  
  colnames(trace) = ReactionsNames
  
  subflux = trace[, c("Time", FluxVec)] %>%
    tidyr::gather(key = "Reaction", value = "Flux", -Time)
  
  subtrace = trace %>% 
    dplyr::select(Time, IECs, Drug, Damage, 
                  pheme_e, pheme_c, BiomassCD, CD,
                  pro_L_e, leu_L_e, ile_L_e, 
                  val_L_e, trp_L_e, cys_L_e) %>% 
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  assign(paste("subtrace", con, sep = ""),
         cbind(subtrace, Scenario = rep(con, length(subtrace$Marking))))
  
  assign(paste("subflux", con, sep = ""), 
         cbind(subflux, Scenario = rep(con, length(subflux$Flux))))
  
}

subtrace = rbind(subtraceNoDrug, subtraceLowEff)
subflux = rbind(subfluxNoDrug, subfluxLowEff)

fluxes = filter(subflux, grepl("EX_biomass_e|sink_pheme_c|EX_ile_L_e|EX_trp_L_e", Reaction))

for (i in 1:length(fluxes$Scenario)) {
  if(fluxes$Scenario[i] == "LowEff") {
    fluxes$Scenario[i] = "Therapy"
  } else {
    fluxes$Scenario[i] = "No drug"
  }
}

for (i in 1:length(fluxes$Reaction)) {
  if(fluxes$Reaction[i] == "EX_biomass_e") {
    fluxes$Reaction[i] = "EX_biomass(e)"}
}
for (i in 1:length(fluxes$Reaction)) {
  if(fluxes$Reaction[i] == "sink_pheme_c") {
    fluxes$Reaction[i] = "sink_pheme(c)"}
}
for (i in 1:length(fluxes$Reaction)) {
  if(fluxes$Reaction[i] == "EX_ile_L_e") {
    fluxes$Reaction[i] = "EX_ile_L(e)"}
}
for (i in 1:length(fluxes$Reaction)) {
  if(fluxes$Reaction[i] == "EX_trp_L_e") {
    fluxes$Reaction[i] = "EX_trp_L(e)"}
}

################################################################################
#################################### Legend ####################################
################################################################################

FbiomassLegend = ggplot(filter(fluxes, Reaction == "EX_biomass(e)"), 
       aes(Time, Flux*1e+03, colour = Scenario)) +
  geom_line(size = 0.80) + theme_bw() +
  ggtitle("Drug metabolic indirect effects assessed \n by estimated reaction fluxes values") +
  xlab("Time (hour)") + 
  ylab("Reaction Fluxes (× 10^3 mmol/mL per bacteria cell)") +
  labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan')) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.85))

legend = cowplot::get_legend(FbiomassLegend)

grid::grid.newpage()
grid::grid.draw(legend)

################################################################################
################################### Plotting ###################################
################################################################################

p1 = ggplot(filter(fluxes, Reaction == "EX_biomass(e)"), 
                        aes(Time, Flux*1e+03, colour = Scenario)) +
  geom_line(size = 0.80) + theme_bw() +
  # Growth homoeostasis
  ggtitle("Biomass production") +
  xlab("Time (hour)") + 
  labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan')) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(y = expression("Reaction Fluxes (×"~10^{-3}~"mmol per bacteria)"))

p2 = ggplot(filter(fluxes, Reaction == "sink_pheme(c)"), 
            aes(Time, Flux*1e+03, colour = Scenario)) +
  geom_line(size = 0.80) + theme_bw() +
  # Iron heme metabolism
  ggtitle("Intracellular heme consuption") +
  xlab("Time (hour)") + 
  labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan')) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(y = expression("Reaction Fluxes (×"~10^{-3}~"mmol per bacteria)"))

p3 = ggplot(filter(fluxes, Reaction == "EX_ile_L(e)"), 
            aes(Time, Flux*1e+03, colour = Scenario)) +
  geom_line(size = 0.80) + theme_bw() +
  # Amino acids homoeostasis
  ggtitle("Isoleucine uptake") +
  xlab("Time (hour)") + 
  labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan')) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(y = expression("Reaction Fluxes (×"~10^{-3}~"mmol per bacteria)"))

p4 = ggplot(filter(fluxes, Reaction == "EX_trp_L(e)"), 
            aes(Time, Flux*1e+03, colour = Scenario)) +
  geom_line(size = 0.80) + theme_bw() +
  # Amino acids homoeostasis
  ggtitle("Tryptophan uptake") +
  xlab("Time (hour)") + 
  labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan')) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(y = expression("Reaction Fluxes (×"~10^{-3}~"mmol per bacteria)"))

pdf("biomassFlux.pdf", width = 6, height = 3.5, bg = "transparent")
p1
dev.off()

pdf("hemeFlux.pdf", width = 6, height = 3.5, bg = "transparent")
p2
dev.off()

pdf("ileFlux.pdf", width = 6, height = 3.5, bg = "transparent")
p3
dev.off()

pdf("trpFlux.pdf", width = 6, height = 3.5, bg = "transparent")
p4
dev.off()