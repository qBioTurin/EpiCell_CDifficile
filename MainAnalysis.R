
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

if(!file.exists(model_R.file)) {
  FBAmodel.gen(model.type = model.type,
               model.name = model.name,
               NutritionModeling = "VMH",
               diet_name = "EU_average",
               write = T, saveRData = T,
               change_manual = c("EX_pro_L_e", "EX_leu_L_e", 
                                 "EX_ile_L_e", "EX_val_L_e", 
                                 "EX_trp_L_e", "EX_cys_L_e", 
                                 "sink_pheme_c", "EX_biomass_e"),
               lb_manual = c(rep(1e-14, 6), pheme.lb, 1e-14),
               gene.exp = F, gene.exp_sample = 1, wd = wd)
}

all.react = readRDS(paste0(wd, "/input/", model.type, "/", model.name, "/all.react.rds"))
all.react = all.react[[1]]

Ib = dplyr::filter(all.react, ReactionDistrict == "boundary")

write.table(sort(Ib$ReactionPos, decreasing = FALSE), 
            file = paste0(wd, "/input/csv/react_index.txt"), 
            row.names = FALSE, col.names = FALSE)

######

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
            debug = T,
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

#######

PlottingAblation(Exper = "Model_Analysis", 
                 colConNoDrug = "darkred", colConTherapy = "darkcyan", coltag = "darkgray",
                 tag = c("Ablated", "Unified"),
                 Condition = c("NoDrug", "Therapy"),
                 wd = "/home/raucello/EpiCell_CDifficile")
