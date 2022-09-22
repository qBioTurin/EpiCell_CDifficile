
library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="master", force = T)
library(epimod)
# downloadContainers()
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

setwd("~/EpiCell_CDifficile")

net_name = "EpitCellDifficileHemeSink10000"
  
model.generation(net_fname = paste("./Net/", net_name, ".PNPRO", sep = ""),
                 transitions_fname = "./Net/DynamicBoundsSinkHeme.cpp", 
                 fba_fname = "FBAModelEUaverage")

system(paste("mv ", net_name, ".* ./Net", sep = ""))

################################################################################
############## Analysis (Sensitivity on CD duplication rate) ###################
################################################################################

model.sensitivity(solver_fname = paste("./Net/", net_name, ".solver", sep = ""),
                  i_time = 0, f_time = 48, s_time = 0.5, 
                  fba_fname = "FBAModelEUaverage", 
                  parallel_processors = 20,
                  atol = 1e-10, rtol = 1e-10,
                  event_times = c(8, 16, 24, 32, 40), 
                  event_function = "treat_generation",
                  n_config = 20*50, target_value = "targetFN",
                  parameters_fname = "Input/csv/ParametersSensitivityAnalysis.csv",
                  functions_fname = "Rfunction/FunctionsSen.R")

setwd(paste("~/EpiCell_CDifficile/", net_name, "_sensitivity/", sep = ""))
load(paste("./", net_name, "-sensitivity.RData", sep = ""))

datalist = lapply( list.files(pattern = "\\.trace$") ,function(i){
  TracesPath = i
  n = as.numeric(gsub(pattern = paste("(", net_name, "-sensitivity-)|(.trace)", sep = ""), 
                      replacement = "", x = i) )
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

################################################################################
########################### Plotting (Sensitivity) #############################
################################################################################

setwd("~/EpiCell_CDifficile")

grDevices::pdf(file = "./Results/test.pdf",   
               # The directory you want to save the file in
               width = 20, height = 8)

g = data.frame(x = c(2, 10), y = c(5000, 1.52), 
               label = c("MIC conc", "MIC conc"))

DrugMIC = 5.8e-06*1e+09 # [pmol/mL]

# colo <- colorRampPalette(c("darkred", "orangered"))(length((subtrace[which(subtrace$Places == "CD"), ]$Time)))

pflux = ggplot(subtrace[which(subtrace$Places %in% 
                                c("BiomassCD", "EX_biomass_e",
                                  "EX_pro_L_e", "EX_leu_L_e", "EX_ile_L_e",
                                  "EX_val_L_e", "EX_trp_L_e", "EX_cys_L_e",
                                  "sink_pheme_c")), ],
               aes(x = Time, y = Marking, color = Death4Treat )) +
  geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
  ggtitle("Biomass [pg] & Fluxes [mmol/gDW*h]")

p1 = ggplot(subtrace[which(subtrace$Places == "CD"), ], 
            aes(x = Time, y = Marking, color = Death4Treat )) + 
  geom_line(aes(group = config)) + ggtitle("CD [cell]")

p2 = ggplot(subtrace[which(subtrace$Places == "Drug"), ], 
            aes(x = Time, y = Marking, color = Death4Treat)) + 
  geom_hline(yintercept = DrugMIC, linetype="dotted", color = "darkblue") +
  geom_line(aes(group = config)) + ggtitle("Drug [pmol]")

# colo <- colorRampPalette(c("darkblue", "cyan"))(length((subtrace[which(subtrace$Places == "CD"), ]$Time)))

p3 = ggplot(subtrace[which(subtrace$Places == "CD"), ], 
            aes(x = Time, y = Marking, color = Efflux)) + 
  geom_line(aes(group = config)) + 
  scale_colour_gradientn(colours = c("darkred", "red", "orangered", "orange")) + 
  ggtitle("CD [cell]")

p4 = ggplot(subtrace[which(subtrace$Places == "Drug"), ], 
            aes(x = Time, y = Marking, color = Efflux)) + 
  geom_line(aes(group = config)) +
  geom_hline(yintercept = DrugMIC, linetype="dotted", color = "darkblue") +
  scale_colour_gradientn(colours = c("darkred", "red", "orangered", "orange")) +
  ggtitle("Drug [pmol]")

((p1 + p2) / (p3 + p4)) | pflux

dev.off()
