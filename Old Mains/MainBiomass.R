
############################################################
######## Dependencies on Biomass reaction Experiment #######
############################################################

setwd("./EpiCell_CDifficile")

# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="master")

library(epimod)
library(dplyr)
library(patchwork)

# Looking for Biomass components given the C. difficile model. 


load(file = "./Input/CD196_heme.RData")

BIGGdata_CD196_heme_react = readr::read_delim("./Input/BIGGdata_CD196_heme_react.tsv", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE)

react.index_Cdiff = sybil::checkReactId(model.mat, "biomass205")["biomass205"]@react_pos

BiomassData = list(React.ID = "biomass205", 
                   React.Index = react.index_Cdiff, 
                   React.Name = BIGGdata_CD196_heme_react$rxnNames[which(sybil::react_id(model.mat) %in% "biomass205")],
                   React.Lb = sybil::lowbnd(model.mat)[sybil::react_id(model.mat) %in% "biomass205"], 
                   React.Ub = sybil::uppbnd(model.mat)[sybil::react_id(model.mat) %in% "biomass205"],
                   Gene.Eq = model.mat@gpr[react.index_Cdiff],
                   React.Eq = BIGGdata_CD196_heme_react$equation[c(sybil::react_id(model.mat) %in% "biomass205")],
                   Met.biomass = model.mat@met_id[which(model.mat@S[, react.index_Cdiff] != 0)],
                   FirstOrderFromBiomass = model.mat@react_id[which(Matrix::colSums(model.mat@S[, react.index_Cdiff]*model.mat@S) != 0)])

FluxVec1 = c("EX_biomass(e)", "EX_cys_L(e)", "EX_trp_L(e)", 
             "EX_val_L(e)", "EX_ile_L(e)", "EX_leu_L(e)", "EX_pro_L(e)")

FluxVec2 = BiomassData[["FirstOrderFromBiomass"]]

################################################################################
########################## Experiment Testing Bounds ###########################
################################################################################

model.generation(net_fname = "./Net/EpitCellDifficile.PNPRO",
								 transitions_fname = "./Net/TransitionGeneralFixedBounds.cpp",
								 LP = T)

system("mv EpitCellDifficile.* ./Net")

model.analysis(solver_fname = "./Net/EpitCellDifficile.solver",
                 i_time = 0,
                 f_time = 2,
                 s_time = 1,
                 parameters_fname = "Input/ParametersListDefault.csv",
                 functions_fname = "Rfunction/Functions.R",
                 debug = T)

glpkFluxes.all <- read.csv("./EpitCellDifficile_analysis/EpitCellDifficile-analysis-0010000.trace", header=FALSE, sep=";")
colnames(glpkFluxes.all) = glpkFluxes.all[1, ]
glpkFluxes.all = glpkFluxes.all[-1, ]

glpkFluxes.aa = glpkFluxes.all %>% select(FluxVec1, Time) %>% tidyr::gather(key = "Reaction", value = "Flux", -Time)
glpkFluxes.aa = glpkFluxes.aa[glpkFluxes.aa$Time == 0, ]
glpkFluxes.aa["Time"] = NULL
glpkFluxes.aa["Scenario"] = rep("D", n = length(glpkFluxes.aa$Reaction))

glpkFluxes.all = glpkFluxes.all %>% select(FluxVec2, Time) %>% tidyr::gather(key = "Reaction", value = "Flux", -Time)
glpkFluxes.all = glpkFluxes.all[glpkFluxes.all$Time == 0, ]
glpkFluxes.all["Time"] = NULL
glpkFluxes.all["Scenario"] = rep("D", n = length(glpkFluxes.all$Reaction))

model.analysis(solver_fname = "./Net/EpitCellDifficile.solver",
               i_time = 0,
               f_time = 2,
               s_time = 1,
               parameters_fname = "Input/ParametersAlmostZero.csv",
               functions_fname = "Rfunction/Functions.R",
               debug = T)

glpkFluxesAZ.all <- read.csv("./EpitCellDifficile_analysis/EpitCellDifficile-analysis-000.trace", header=FALSE, sep=";")
colnames(glpkFluxesAZ.all) = glpkFluxesAZ.all[1, ]
glpkFluxesAZ.all = glpkFluxesAZ.all[-1, ]

glpkFluxesAZ.aa = glpkFluxesAZ.all %>% select(FluxVec1, Time) %>% tidyr::gather(key = "Reaction", value = "Flux", -Time)
glpkFluxesAZ.aa = glpkFluxesAZ.aa[glpkFluxesAZ.aa$Time == 0, ]
glpkFluxesAZ.aa["Time"] = NULL
glpkFluxesAZ.aa["Scenario"] = rep("AZ", n = length(glpkFluxesAZ.aa$Reaction))
glpkFluxes.aa = rbind(glpkFluxes.aa, glpkFluxesAZ.aa)
glpkFluxes.aa$Flux = as.double(glpkFluxes.aa$Flux)

Plot.aa = ggplot2::ggplot(glpkFluxes.aa, ggplot2::aes(x = Reaction, y = Flux, fill = Scenario)) +
  ggplot2::geom_col(width = 0.1, position = ggplot2::position_dodge(1)) + 
  ggplot2::ggtitle("Fluxes Essentials") + ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                                                         plot.title = ggplot2::element_text(size=8),
                                                         axis.text = ggplot2::element_text(size = 5),
                                                         legend.key.size = ggplot2::unit(0.15, 'cm'),
                                                         legend.title=ggplot2::element_text(size = 6), 
                                                         legend.text=ggplot2::element_text(size=5))

glpkFluxesAZ.all = glpkFluxesAZ.all %>% select(FluxVec2, Time) %>% tidyr::gather(key = "Reaction", value = "Flux", -Time)
glpkFluxesAZ.all = glpkFluxesAZ.all[glpkFluxesAZ.all$Time == 0, ]
glpkFluxesAZ.all["Time"] = NULL
glpkFluxesAZ.all["Scenario"] = rep("AZ", n = length(glpkFluxes.all$Reaction))
glpkFluxes.all = rbind(glpkFluxes.all, glpkFluxesAZ.all)
glpkFluxes.all$Flux = as.double(glpkFluxes.all$Flux)

Plot.all = ggplot2::ggplot(glpkFluxes.all, ggplot2::aes(x = Reaction, y = Flux, fill = Scenario)) +
  ggplot2::geom_col(width = 3, position = ggplot2::position_dodge(5)) + 
  ggplot2::ggtitle("Fluxes Comparison: D (P = 1) vs AZ (P = 1e-06)") + 
  ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                 plot.title = ggplot2::element_text(size=8),
                 axis.text.y = ggplot2::element_text(size = 5, angle = 0, vjust = 0.5, hjust=1),
                 axis.text.x = ggplot2::element_text(size = 1, angle = 90, vjust = 0.5, hjust=1),
                 legend.key.size = ggplot2::unit(0.15, 'cm'),
                 legend.title=ggplot2::element_text(size = 6), 
                 legend.text=ggplot2::element_text(size=5))

a = glpkFluxes.all[glpkFluxes.all$Reaction %in% glpkFluxes.all[which(abs(glpkFluxes.all$Flux[which(glpkFluxes.all$Scenario == "D")]) -
                       abs(glpkFluxes.all$Flux[which(glpkFluxes.all$Scenario == "AZ")]) > 50), ]$Reaction, ] 
  
Plot.allDiff = ggplot2::ggplot(a, ggplot2::aes(x = Reaction, y = Flux, fill = Scenario)) +
  ggplot2::geom_col(position = "dodge") + 
  ggplot2::geom_text(ggplot2::aes(label = round(Flux)), 
                     colour = "black", size = 2.5,
                     vjust = 0.5, position = ggplot2::position_dodge(.7)) + 
  ggplot2::theme(axis.title = ggplot2::element_text(size=7), 
                 plot.title = ggplot2::element_text(size=8),
                 axis.text.y = ggplot2::element_text(size = 5, angle = 0, vjust = 0.5, hjust=1),
                 axis.text.x = ggplot2::element_text(size = 5, angle = 90, vjust = 0.5, hjust=1),
                 legend.key.size = ggplot2::unit(0.15, 'cm'),
                 legend.title=ggplot2::element_text(size = 6), 
                 legend.text=ggplot2::element_text(size=5)) + 
  ggplot2::ggtitle("Differing fluxes (diff > 50 mmol/h): D (P = 1) vs AZ (P = 1e-06)")

Plot.aa | (Plot.all / Plot.allDiff)

