
################################################################################
###################### Sensitivity analysis full model #########################
################################################################################

library(epimod)

model.generation(net_fname = "./Net/EpitCellDifficile.PNPRO",
                 transitions_fname = "./Net/TransitionGeneral.cpp",
                 LP = T)

system("mv EpitCellDifficile.* ./Net")

debug(model.analysis)

# drug 0
# parms$solver_type = "LSODA -atol 1e-13 -rtol 1e-13"

model.analysis(solver_fname = "./Net/EpitCellDifficile.solver",
               i_time = 0,
               f_time = 6,
               s_time = 0.01,
               parameters_fname = "Input/ParametersListSensitivity.csv",
               functions_fname = "SensitivityFBA_analysis/FunctionsSensitivity.R", 
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

TracesPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-1.trace"
trace = read.table(TracesPath, header = T)

FluxPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-00.trace"
flux = read.csv(FluxPath, header=FALSE, sep=";")
colnames(flux) = flux[1, ]
flux = flux[-1, ]

saveRDS(trace, "./Results/AnalysisPlaces_Drug005_hemeLow.rds")
saveRDS(flux, "./Results/AnalysisFlux_Drug005_hemeLow.rds")

# source("SensitivityFBA_analysis/SensitivityPlot.R")
# SensitivityPlot()

source("./Rfunction/ModelAnalysisPlot.R")

ModelAnalysisPlot(TracesPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-1.trace",
                  FluxPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-00.trace",
                  FluxVec = c("EX_biomass(e)", "EX_pheme(e)", "EX_cys_L(e)", "EX_trp_L(e)", 
                              "EX_val_L(e)", "EX_ile_L(e)", "EX_leu_L(e)", "EX_pro_L(e)"))

ggplot2::ggsave(ModelAnalysisPlot(TracesPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-1.trace",
                                  FluxPath = "EpitCellDifficile_analysis/EpitCellDifficile-analysis-00.trace",
                                  FluxVec = c("EX_biomass(e)", "EX_pheme(e)", "EX_cys_L(e)", "EX_trp_L(e)", 
                                              "EX_val_L(e)", "EX_ile_L(e)", "EX_leu_L(e)", "EX_pro_L(e)")) +
                  ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")),
                filename = "ggp_transparent1.png", bg = "transparent",  width = 10, height = 8, dpi = 150, units = "in", device='png')

