
setwd("~/EpiCell_CDifficile")

# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="master")
library(epimod)
library(dplyr)

# Lets compile the metabolic model from the RData storing the chemical reaction network

source("Rfunction/Functions.R")
# FBA.generation(model = "./Input/CDmodels/CD196HemeSink/CD196HemeSink.RData", write = T)

# Sensitivity analysis on Bounds (sampled values)

model.generation(net_fname = "./Net/EpitCellDifficileHemeSink.PNPRO",
                 transitions_fname = "./Net/ConstantBoundsSinkHeme.cpp",
                 LP = T)

system("mv EpitCellDifficileHemeSink.* ./Net")

################################################################################
########################## Analysis (Static Bounds) ############################
################################################################################

set.seed(1234)

FluxVec = c("EX_biomass.e.", 
            "sink_pheme.c.", 
            "EX_cys_L.e.", "EX_trp_L.e.", "EX_val_L.e.", 
            "EX_ile_L.e.", "EX_leu_L.e.", "EX_pro_L.e.")

P = c(1e-350, round(sort(runif(18, min = 0.001, max = 0.15)), 5))

SensitivityDataGLPK = data.frame()

for (i in 1:length(P)) {
  
  model.analysis(solver_fname = "./Net/EpitCellDifficileHemeSink.solver",
                 i_time = 0, f_time = 2, s_time = 1, 
                 parameters_fname = "Input/csv/ParametersListSinkHeme.csv",
                 functions_fname = "Rfunction/Functions.R", 
                 ini_v = c( 
                   # P
                   P[i],
                   # index
                   i,
                   # step time
                   s_time = 1), debug = T)
  
  FluxPath = paste("./EpitCellDifficileHemeSink_analysis/EpitCellDifficileHemeSink-analysis-00", 
              i,".trace", sep = "")
  
  subflux = read.table(FluxPath, header = T, sep = ";") %>% 
    select(FluxVec, Time) %>% 
    tidyr::gather(key = "Reaction", value = "Flux", -Time)
  
  scenario = rep(P[i], length(subflux$Time))
  subflux = cbind(subflux, scenario)
  
  SensitivityDataGLPK = rbind(SensitivityDataGLPK, subflux)
}

save(SensitivityDataGLPK, file = "./Input/CDmodels/CD196HemeSink/SensitivityDataGLPK.RData")
load("./Input/CDmodels/CD196HemeSink/SensitivityDataGLPK.RData")

SensitivityDataGLPK = SensitivityDataGLPK[, -1] %>% distinct()
colnames(SensitivityDataGLPK) = c("key", "Value", "BoundsActivity")

system("rm *.log")

ggplot2::ggplot(SensitivityDataGLPK, ggplot2::aes(x = BoundsActivity, y = Value, color = key)) +
  ggplot2::geom_line() + ggplot2::ylim(-150, 150) +
  ggplot2::geom_vline(xintercept = P, linetype="dotted") +
  ggplot2::geom_point(col = "black", shape = 21, size = 1) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 12), 
                 plot.title = ggplot2::element_text(size = 12),
                 axis.text = ggplot2::element_text(size = 10),
                 legend.key.size = ggplot2::unit(0.35, 'cm'),
                 legend.title=ggplot2::element_text(size = 10), 
                 legend.text=ggplot2::element_text(size = 10),
                 axis.title.x = ggplot2::element_text(color = "black", size = 12),
                 axis.title.y = ggplot2::element_text(color = "black", size = 12)) +
  ggplot2::scale_colour_discrete(name="Target Reactions") +
  ggplot2::scale_colour_brewer(palette = "Set1") +
  ggplot2::ggtitle("Changing bounds Test (essential AAs): algo = GLPK (from epimod)") + 
  ggplot2::labs(y = "Flux [mmol/h*gDW]", x = "Bounds activity [x*lb; x*up]")
