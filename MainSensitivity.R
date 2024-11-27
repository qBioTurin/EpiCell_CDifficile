
# loading library

# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="variabilityFBA", force = T)
library(epimod)
# downloadContainers(tag = "latest")
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(parallel)

# setting working directory
wd = "/home/raucello/EpiCell_CDifficile"
setwd(wd)
supp_function.dir = "/code/supplementary_functions/"

# FBA model compiling phase
source(paste0(wd, supp_function.dir, "Functions.R"))
source(paste0(wd, supp_function.dir, "FBAmodel.generation.R"))
# loading experiments executioner
source(paste0(wd, supp_function.dir, "Exe.R"))
source(paste0(wd, supp_function.dir, "PlottingAblation.R"))

# setting FBA model tags
model.type = "Clostridium"
model.name = "CD196HemeSink"

######

# # single cell C.diff mean observed biomass
gDW_CDmean = ComputeBiomassBac(diameter = 0.5, len = 5.5) # (pg)
# initial C. diff cells number
nBac = 45700000 # (cell)

# initial extracellular heme concentration:
# from -> (https://doi.org/10.1038/s41419-020-02950-8)
pheme_e.i = 100 # (μM)
# first conversion: 1 (μM) = 1.0e-06 (mmol/mL)
pheme_e.i = pheme_e.i*1.0e-06 # (mmol/mL)
# second conversion: 1 (mmol/mL) = 1.0e+09 (pmol/mL)
pheme_e.i = pheme_e.i*1e+09 # (pmol/mL)
pheme.lb = (-(pheme_e.i*1e-09)/(nBac*gDW_CDmean*1e-12)) # (mmol/gDW*h)

######

model_R.file = paste0(wd, "/input/", model.type, "/", model.name, "/", model.name, ".RData")

FBAmodel.gen(model.type = model.type,
           model.name = model.name,
           NutritionModeling = "VMH",
           diet_name = "EU_average",
           write = T, saveRData = T,
           change_manual = c("EX_pro_L_e", "EX_leu_L_e", 
                             "EX_ile_L_e", "EX_val_L_e", 
                             "EX_trp_L_e", "EX_cys_L_e", 
                             "sink_pheme_c", "EX_biomass_e"),
           lb_manual = c(rep(-1e-14, 6), pheme.lb, -1e-14),
           gene.exp = F, gene.exp_sample = 1, wd = wd) 

all.react = readRDS(paste0(wd, "/input/", model.type, "/", model.name, "/all.react_diet.rds"))
all.react = all.react[[1]]

Ib = dplyr::filter(all.react, ReactionDistrict == "boundary")

write.table(sort(Ib$ReactionPos, decreasing = FALSE), 
          file = paste0(wd, "/input/csv/react_index.txt"), 
          row.names = FALSE, col.names = FALSE)

system("rm -f dockerID *error.log *.log ExitStatusFile")

for (ablation in c("Unified", "Ablated", "SuperAblated")) {
load(model_R.file)
Exe.exp(model_cat = model.name,
        model_name = model.name,
        fba_fname = paste0(model.name, ".txt"),
        atol = 1e-06, rtol = 1e-06,
        time.step = 0.5, f_time = 72,
        event_times = c(8, 16, 24, 32, 40, 48),
        distance_measure = "ReferenceM",
        reference_data = paste0(wd, "/input/csv/ReferenceRanking.csv"),
        Exper = "Model_Sensitivity",
        Ablation = ablation,
        therapy = "Therapy",
        debug = F,
        FVA = F,
        FVAcomparison = F,
        react = "EX_biomass_e",
        debugFVA = F,
        replicates = 24,
        param_target = c("IECsDeath", "Death4Treat", "Detox"),
        net_fname = "EpitCellDifficileHemeSink",
        cores = detectCores(),
        wd = wd,
        supp_function.dir = supp_function.dir)
}

###### Model Simulation Phase

net_fname = "EpitCellDifficileHemeSink"

Exper = "Model_Sensitivity"
Ablation = "Unified"
therapy = "Therapy"

param_target = c("IECsDeath", "Death4Treat", "Detox")

setwd(paste0(wd, "/results/CDiff_", Ablation, "_", therapy, "_", Exper, "/", net_fname, "_sensitivity"))

id = as.integer(readLines("ConfigID.txt"))

setwd(paste0(wd, "/results/CDiff_", Ablation, "_", therapy, "_", Exper, "/"))
load(paste0(net_fname, "-analysis.RData"))

id_param <- sapply(param_target, function(x) {
which(sapply(config, function(y) y[[1]][[1]] == x))
})

for (p in seq_along(param_target)) {
assign(param_target[p], config[[id_param[p]]][[id]][[3]])
}

csv <- read.csv(file = paste0(wd, "/input/csv/ParametersListSinkHeme.csv"), header = F, quote = "")

par = c(paste0("c; IECsDeath; ", IECsDeath),
      paste0("g; Death4Treat; ", Death4Treat),
      paste0("c; Detox; ", Detox))

for(p in 1:length(id_param)) {csv[id_param[p], ] = par[p]}

# csv = csv %>% 
#   mutate_if(is.character, function(x) gsub("'", " ", x)) %>% 
#   mutate_if(is.character, function(x) paste0("'", x, "'"))

write.table(csv, file = paste0(wd, "/input/csv/ParametersListSinkHeme.csv"), col.names = FALSE, row.names = F, quote = F)

#####

load(model_R.file)
system("rm -f dockerID *error.log *.log ExitStatusFile")

ablation_types <- c("Unified", "Ablated", "SuperAblated")
therapy_types <- c("NoDrug", "Therapy")

for (ablation in ablation_types) {
for (therapy in therapy_types) {
  Exe.exp(model_cat = model.name,
          model_name = model.name,
          fba_fname = paste0(model.name, ".txt"),
          atol = 1e-06, rtol = 1e-06,
          time.step = 0.5, f_time = 72,
          event_times = c(8, 16, 24, 32, 40, 48),
          distance_measure = "ReferenceM",
          reference_data = paste0(wd, "/input/csv/ReferenceRanking.csv"),
          Exper = "Model_Analysis",
          Ablation = ablation,
          therapy = therapy,
          debug = F,
          FVA = F,
          FVAcomparison = F,
          react = "EX_biomass_e",
          debugFVA = F,
          replicates = 1,
          param_target = "IECsDeath",
          net_fname = "EpitCellDifficileHemeSink",
          cores = detectCores(),
          wd = wd,
          supp_function.dir = supp_function.dir)
}
}

################################ Plotting Phase ################################

PlottingAblation(Exper = "Model_Analysis", 
               colConNoDrug = "darkred", colConTherapy = "darkcyan", coltag = "darkgray",
               tag = c("Ablated", "Unified"),
               Condition = c("NoDrug", "Therapy"),
               wd = "/home/raucello/EpiCell_CDifficile")

# source("./Rfunction/PltAbZoom.R")
source(paste0(wd, supp_function.dir, "PlottingAblationViolin.R"))
source(paste0(wd, supp_function.dir, "PlottingAblationSensitivity.R"))
source(paste0(wd, supp_function.dir, "PlottingPanelPaper.R"))

# ################################ React2Path ####################################
# 
# source("./Rfunction/DatasetforPvalue.R")
# 
# ######################## Flux variability analysis #############################
# 
# setwd("~/EpiCell_CDifficile")
# 
# system("rm -f dockerID *error.log *.log ExitStatusFile")
# source("./Rfunction/Exe.R")
# 
# load(paste0("./Input/CDmodels/", "CD196HemeSink", "/", "CD196HemeSinkDiet", ".RData"))
# # react = "EX_biomass(e)"
# # react = as.vector(sybil::findExchReact(model.mat)@react_id)
# react = as.vector(model.mat@react_id)
# 
# Exe.exp(model_cat = "CD196HemeSink", model_name = "CD196HemeSinkDiet", 
#         fba_fname = "FBAModelEU_average", 
#         atol = 1e-06, rtol = 1e-06, f_time = 72, 
#         event_times = c(8, 16, 24, 32, 40, 48),
#         distance_measure = NULL, reference_data = NULL,
#         Exper = "Model_Analysis", Ablation = F, therapy = T, 
#         debug = F, FVA = T, FVAcomparison = T, react = react,
#         debugFVA = F, replicates = 1,
#         param_target = c("IECsDeath", "Death4Treat", "Detox"))
# 
# ################################################################################
# 
# load("~/EpiCell_CDifficile/subfvaMaxUnifiedTherapy.RData")
