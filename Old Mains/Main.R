
#### This is with 1370 (without the Ex and Tr of Sheme) ####

#' 1) Let load the metabolic netwwork
#'
#' source("Rfunction/MatlabFileGeneration.R")
#' Matlab.file = readMATmod("Input/CD196_heme.mat")
#' save(Matlab.file, file = "Input/FBAmodel.RData")
#' source("Rfunction/Functions.R")
#' FBA.generation( model = "./Input/FBAmodel.RData", write = T)

#' 2) Let load the Epimod package

setwd("~/EpiCell_CDifficile")

# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="master")
library(epimod)

################################################################################
################ Sensitivity analysis on Bounds (sampled values) ###############
################################################################################

model.generation(net_fname = "./Net/EpitCellDifficile.PNPRO",
								 transitions_fname = "./Net/TransitionGeneralFixedBounds.cpp",
								 LP = T)

system("mv EpitCellDifficile.* ./Net")

#' 4) model analysis

model.analysis(solver_fname = "./Net/EpitCellDifficile.solver",
               i_time = 0,
               f_time = 2,
               n_config = 12,
               s_time = 1,
               parameters_fname = "Input/ParametersListSensitivity.csv",
               functions_fname = "Rfunction/Functions.R")


FluxVec = c("EX_biomass.e.", "EX_pheme.e.", "EX_cys_L.e.", "EX_trp_L.e.", 
            "EX_val_L.e.", "EX_ile_L.e.", "EX_leu_L.e.", "EX_pro_L.e.")

library(dplyr)

FluxList = list.files(path = "./EpitCellDifficile_analysis/", 
                     pattern = "EpitCellDifficile-analysis-00")

subflux_total = data.frame()

for(i in FluxList) {
  subflux = read.table(paste("./EpitCellDifficile_analysis/", i, sep = ""), header = T, sep = ";") %>% 
    select(FluxVec, Time) %>% tidyr::gather(key = "Reaction", value = "Flux", -Time)
  subflux_total = rbind(subflux_total, subflux)
}

saveRDS(subflux_total, "./glpkFluxes.rds")
glpkFluxes <- readRDS("./glpkFluxes.rds")
glpkFluxes = glpkFluxes[glpkFluxes$Time == 0, ]
glpkFluxes = glpkFluxes[-c(which(glpkFluxes$Reaction == "EX_pheme.e.")), ]
glpkFluxes["Time"] = NULL

for (i in 1:length(FluxList)) {
  if(nchar(FluxList[i]) == 36) {
    str = substr(FluxList[i], 
                 nchar("EpitCellDifficile-analysis-00"), 
                 nchar(FluxList[i]))
    FluxList[i] = paste("EpitCellDifficile-analysis-00", str, sep = "")
  } else if (FluxList[i] == "EpitCellDifficile-analysis-000.trace") {
    FluxList[i] = "EpitCellDifficile-analysis-0000.trace"
    }
  }

glpkFluxes["Scenario"] = rep(substr(FluxList, 
                                nchar("EpitCellDifficile-analysis-00"), 
                                nchar(FluxList[1])-6), each = 7)

names(glpkFluxes) = c("ReactionName", "v", "Scenario")
glpkFluxes = glpkFluxes[order(glpkFluxes$Scenario), ]

ggplot2::ggplot(glpkFluxes, ggplot2::aes(x = ReactionName, y = v, fill = Scenario)) +
  ggplot2::geom_col(width = 0.33, position = ggplot2::position_dodge(0.5)) + 
  ggplot2::ggtitle("glpkFluxes") + ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                                                  plot.title = ggplot2::element_text(size=8),
                                                  axis.text = ggplot2::element_text(size = 5),
                                                  legend.key.size = ggplot2::unit(0.15, 'cm'),
                                                  legend.title=ggplot2::element_text(size = 6), 
                                                  legend.text=ggplot2::element_text(size=5))

################################################################################
############## Sensitivity analysis on Bounds (handpicked values) ##############
################################################################################

model.generation(net_fname = "./Net/EpitCellDifficile.PNPRO",
                 transitions_fname = "./Net/TransitionGeneralFixedBounds.cpp",
                 LP = T)

system("mv EpitCellDifficile.* ./Net")

input = c("ParametersListVeryLow.csv", "ParametersListLow.csv", 
          "ParametersListMid.csv", "ParametersListHigh.csv", 
          "ParametersListDefault.csv")

subflux_total = data.frame()

for(i in input) {
  model.analysis(solver_fname = "./Net/EpitCellDifficile.solver",
                 i_time = 0,
                 f_time = 2,
                 s_time = 1,
                 parameters_fname = paste("./Input/", i, sep = ""),
                 functions_fname = "Rfunction/Functions.R",
                 debug = T)
  
  FluxVec = c("EX_biomass.e.", "EX_pheme.e.", "EX_cys_L.e.", "EX_trp_L.e.", 
              "EX_val_L.e.", "EX_ile_L.e.", "EX_leu_L.e.", "EX_pro_L.e.")
  
  FluxPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-00.trace"
  
  subflux = read.table(FluxPath, header = T, sep = ";") %>% select(FluxVec, Time) %>%
    tidyr::gather(key = "Reaction", value = "Flux", -Time)
  subflux_total = rbind(subflux_total, subflux)
}

ggplot2::ggplot(subflux_total, ggplot2::aes(x = ReactionName, y = v, fill = Scenario)) +
  ggplot2::geom_col(width = 0.33, position = ggplot2::position_dodge(0.5)) + 
  ggplot2::ggtitle("Fluxes Essentials") + ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                                                         plot.title = ggplot2::element_text(size=8),
                                                         axis.text = ggplot2::element_text(size = 5),
                                                         legend.key.size = ggplot2::unit(0.15, 'cm'),
                                                         legend.title=ggplot2::element_text(size = 6), 
                                                         legend.text=ggplot2::element_text(size=5))

################################################################################
########################## Analysis (Dynamic Bounds) ###########################
################################################################################

model.generation(net_fname = "./Net/EpitCellDifficile.PNPRO",
                 transitions_fname = "./Net/TransitionGeneral.cpp",
                 LP = T)

system("mv EpitCellDifficile.* ./Net")

model.analysis(solver_fname = "./Net/EpitCellDifficile.solver",
							 i_time = 0,
							 f_time = 12,
							 s_time = 1,
							 parameters_fname = "Input/ParametersListTest.csv",
							 functions_fname = "Rfunction/Functions.R", 
							 ini_v = c(
							   # h # h = 0.5 half-life each hour
							   0.5,
							   # IECsDeath = 1e-10
							   1e-10,
							   # P
							   10000,
							   # nBacMax
							   4.57e05*100*100),
							 debug = T)

source("./Rfunction/ModelAnalysisPlot.R")

ModelAnalysisPlot(TracesPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-1.trace",
                  FluxPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-00.trace",
                  FluxVec = c("EX_biomass(e)", "EX_pheme(e)", "EX_cys_L(e)", "EX_trp_L(e)", 
                              "EX_val_L(e)", "EX_ile_L(e)", "EX_leu_L(e)", "EX_pro_L(e)"))
# + ggplot2::scale_y_log10()

TracesPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-1.trace"
FluxPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-00.trace"
trace = read.table(TracesPath, header = T)

ggplot2::ggsave(ModelAnalysisPlot(TracesPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-1.trace",
                                  FluxPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-00.trace",
                                  FluxVec = c("EX_biomass(e)", "EX_pheme(e)", "EX_cys_L(e)", "EX_trp_L(e)", 
                                              "EX_val_L(e)", "EX_ile_L(e)", "EX_leu_L(e)", "EX_pro_L(e)")) +
                  ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png",
                bg = "transparent",  width = 15, height = 8, dpi = 150, units = "in", device='png')
