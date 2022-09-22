
library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="master", force = T)
library(epimod)
# downloadContainers()
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

################################################################################################
## Let's compile the FBA model from the RData storing the chemical reaction network with diet ##
################################################################################################

setwd("~/EpiCell_CDifficile")

source("./Input/notes/NewFBA.generation.R")

FBA.generation(RDataPresent = T, NutritionModeling = "VMH", ShowExBounds = T, write = T)

## FBAModelEUaverage main features ##

# ~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData
# ~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/Diets/diets.RData
# EU_average
# -10 [mmol/gDW*h]
# 1000 [mmol/gDW*h]
# -1e-12 [mmol/gDW*h]; biomass could be only produced and not consumed
# 10 [mmol/gDW*h]. We do not know exactly the maximun production rate for each metabolites.
# 1e-12 [mmol/gDW*h]; essential aa could be only produced and not consumed
# Reasonably, it could be in the magnitude of 10 [mmol/gDW*h]

## FBAModelZero main features ##

# ~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData
# ~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/Diets/diets.RData
# Zero
# -1e-50 [mmol/gDW*h]
# 1e-50 [mmol/gDW*h]
# -1e-50 [mmol/gDW*h]
# -1e-50 [mmol/gDW*h]
# -1e-50 [mmol/gDW*h]

################################################################################
######################### Generation (Dynamic Bounds) ##########################
################################################################################

# Sensitivity analysis on Bounds (sampled values)

model.generation(net_fname = "./Net/EpitCellDifficileHemeSink.PNPRO",
                 transitions_fname = "./Net/DynamicBoundsSinkHeme.cpp", 
                 fba_fname = "FBAModelEUaverage")

system("mv EpitCellDifficileHemeSink.* ./Net")

################################################################################
########################## Analysis (Dynamic Bounds) ###########################
################################################################################

model.analysis(solver_fname = "./Net/EpitCellDifficileHemeSink.solver",
               i_time = 0, f_time = 48, s_time = 0.5, 
               fba_fname = "FBAModelEUaverage", 
               atol = 1e-10, rtol = 1e-10, 
               debug = F, 
               event_times = c(8, 16, 24, 32, 40),
               event_function = "treat_generation",
               parameters_fname = "Input/csv/ParametersListSinkHeme.csv",
               functions_fname = "Rfunction/Functions.R")

################################################################################
########################## plotting (Dynamic Bounds) ###########################
################################################################################

grDevices::pdf(file = "./Results/noDrug.pdf",   
               # The directory you want to save the file in
               width = 16, height = 8)

TracesPath = "EpitCellDifficileHemeSink_analysis/EpitCellDifficileHemeSink-analysis-1.trace"
trace = read.table(TracesPath, header = F)
ReactionsNames = trace[1, ]
ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)

trace = read.table(TracesPath, header = T)
colnames(trace) = ReactionsNames

subtrace = trace %>% dplyr::select(Time, IECs, Drug, Damage,
                                   pheme_e, pheme_c, 
                                   BiomassCD, CD,
                                   pro_L_e, leu_L_e, ile_L_e, 
                                   val_L_e, trp_L_e, cys_L_e) %>% 
  tidyr::gather(key = "Places", value = "Marking", -Time)

amino = filter(subtrace, grepl("pro_L_e|leu_L_e|ile_L_e|val_L_e|trp_L_e|cys_L_e", Places))
heme = filter(subtrace, grepl("pheme_e|pheme_c", Places))
cell = filter(subtrace, grepl("IECs|Damage", Places))
biom = filter(subtrace, grepl("Biomass", Places))
clo = filter(subtrace, subtrace$Places == "CD")
drug = filter(subtrace, grepl("Drug", Places))

m = ggplot(subtrace, aes(x = Time, y = Marking)) + geom_line(aes(color = Places)) +
  facet_grid(Places ~ ., scales = "free_y") + theme(legend.position = "none")

FluxVec = c("EX_biomass_e",
            "EX_pro_L_e", "EX_leu_L_e", "EX_ile_L_e",
            "EX_val_L_e", "EX_trp_L_e", "EX_cys_L_e",
            "sink_pheme_c")

subflux = trace[, c("Time", FluxVec)] %>% 
  tidyr::gather(key = "Reaction", value = "Flux", -Time)

# colo = colorRampPalette(c("darkviolet", "magenta"))(length(unique(subflux$Reaction)))

g = data.frame(x = c(40, 40, 40, 10, 45, 5), 
               y = c(1.4, 0.73, 0.25, 1, 1.25e-04, 9000), 
               label = c("biomass max", "biomass mean", "biomass min",
                         "therapy begins",
                         "vigour threshold", 
                         "hemorrhagic heme"))

b = ggplot(biom) + 
  geom_line(aes(x = Time, y = Marking), color = 'darkred', size=0.75) + 
  geom_vline(xintercept = 8, linetype="dotted") + 
  geom_text(data = g[4, ], aes(x = x, y = y, label = label), size = 1.5, angle = 90) +
  ggtitle("CD biomass") + 
  geom_hline(yintercept = 1.57872, linetype="dotted", color = "red") + 
  geom_text(data = g[1, ], aes(x = x, y = y, label = label), size = 1.5) +
  geom_hline(yintercept = 0.496328, linetype="dotted", color = "red") +
  geom_text(data = g[2, ], aes(x = x, y = y, label = label), size = 1.5) +
  geom_hline(yintercept = 0.06738716, linetype="dotted",  color = "red") + 
  geom_text(data = g[3, ], aes(x = x, y = y, label = label), size = 1.5) +
  xlab("Time [h]") + ylab("Biomass [pg/cell]")

cd = ggplot(clo) + 
  geom_line(aes(x = Time, y = Marking), color = 'orange') + 
  geom_vline(xintercept = 8, linetype="dotted") + 
  ggtitle("CD") +
  xlab("Time [h]") + ylab("CD [cell]")

DrugMIC = 5.8e-06 # [mmol/mL]

d = ggplot(drug) + 
  geom_line(aes(x = Time, y = Marking*1e-09), color = '#56B4E7', size=0.75) + 
  geom_hline(yintercept = DrugMIC, linetype="dotted", color = "darkblue") +
  geom_vline(xintercept = 8, linetype="dotted") + 
  ggtitle("Drug conc") +
  xlab("Time [h]") + ylab("Drug [mmol/mL]")

c = ggplot(cell, aes(x=Time, y=Marking, color=Places, 
                     line=Places, group = Places))+ 
  scale_color_manual(values=c('#a5133f','#26B4E4')) +
  geom_point(size=1,shape=21,aes(color=Places)) +
  geom_vline(xintercept = 8, linetype="dotted") +
  ggtitle("Human Components") +
  xlab("Time [h]") + ylab("EICs [cell]")

# plot(subflux[subflux$Reaction == "sink_pheme_c", ]$Flux)
# subflux[which(subflux$Time == 0 | subflux$Time < 25), ]

minBiomFlux = 1e-04

f = ggplot(subflux, aes(x = Time, y = signif(Flux, digits = 3), color = Reaction)) +
  geom_line(size=0.75) + theme_bw() + 
  theme(panel.background = element_rect(fill = "#F8F8FF")) +
  geom_hline(yintercept = minBiomFlux, linetype="dotted", color = "darkgreen") +
  geom_vline(xintercept = 8, linetype="dotted") + ggtitle("Fluxes Trajectories") +
  xlab("Time [h]") + ylab("Estimated Flux [mmol/gDW*h]")

aa = ggplot(amino, aes(x = Time, y = signif((Marking/(6.022e20/6.022e08))*1e+09, digits = 3), color = Places)) +
  geom_line(size=0.75) + theme_bw() + theme(panel.background = element_rect(fill = "#F8F8FF")) +
  geom_vline(xintercept = 8, linetype="dotted") +
  ggtitle("Essential Aminoacids Trajectories") +
  xlab("Time [h]") + ylab("Concentration [pmol/ml]") +
  scale_color_manual(values=colorRampPalette(c("magenta", "gold", 
                                               "red", "blue", "green"))(length(unique(subflux$Reaction))))

sat = 0.01*1e-03*1e+09 # [pmol/mL] concentration easily found at hemorrhagic sites

h = ggplot(heme, aes(x = Time, y = signif(Marking, digits = 3), color = Places)) +
  geom_line(size=0.75) + theme_bw() + theme(panel.background = element_rect(fill = "#F8F8FF")) +
  geom_hline(yintercept = sat, linetype="dotted", color = "darkred") +
  geom_vline(xintercept = 8, linetype="dotted") + 
  ggtitle("Heme Trajectories") +
  xlab("Time [h]") + ylab("Concentration [pmol/mL]") + scale_color_manual(values=c('black','orange'))

print(((b | cd | d) / aa / h / c) | (f))

dev.off()
