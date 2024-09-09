
# setting working directory
wd = getwd()

# setting FBA model tags
model.type = "Clostridium"
model.name = "CD196HemeSink"
net_fname = "EpitCellDifficileHemeSink"

supp_function_dir = paste0(wd, "/code/supplementary_functions/")
model_dir = paste0(wd, "/input/", model.type, "/", model.name, "/")

model_R.file = paste0(model_dir, model.name, ".RData")

# Workflow libraries
source(paste0(wd, "/code/scripts/libraries.R"))
# FBA model compiling phase
source(paste0(supp_function_dir, "Functions.R"))
source(paste0(supp_function_dir, "FBAmodel_generation.R"))
# loading experiments executioner
source(paste0(supp_function_dir, "Exe.R"))
# scripts references for biological quantities
source(paste0(wd, "/code/scripts/biological_numbers.R"))
# plotting
source(paste0(wd, "/code/scripts/plotting_ablation.R"))
source(paste0(wd, "/code/scripts/plotting_ablation_violin.R"))
source(paste0(wd, "/code/scripts/plotting_ablation_sensitivity.R"))
source(paste0(wd, "/code/scripts/plotting_panel_paper.R"))

manual_bounds = data.frame(
  reactions =
    c("EX_biomass_e", "sink_pheme_c",
      "EX_cys_L_e", "EX_leu_L_e", "EX_ile_L_e", "EX_val_L_e", "EX_trp_L_e", "EX_pro_L_e"),
  lb_manual = c(-1000, pheme_c.lb = ifelse(pheme_c.lb < P, P, pheme_c.lb), 
                iniAA.lb <- ifelse(iniAA.lb < P, P, iniAA.lb)),
  ub_manual = c((gDW_CDmax - gDW_CDmean), 1000, rep(0, 6)))

######

FBAmodel.gen(model_type = model.type,
             model_name = model.name,
             NutritionModeling = "VMH",
             diet_name = "EU_average",
             write = T, 
             saveRData = T,
             change_manual = manual_bounds$reactions,
             # (mmol/per community)
             lb_manual = rep(-1000, length(manual_bounds$reactions)),
             # lb_manual = manual_bounds$lb_manual,
             #  (mmol/gDW*h)
             ub_manual = rep(1000, length(manual_bounds$reactions)),
             # ub_manual = manual_bounds$ub_manual,
             lb_off = NA,
             gene.exp = F, 
             gene.exp_sample = 1,
             # (mmol/per community)
             GeneExpress_file = NA,
             wd = wd)

# extracting features for model boundary reactions
all_react = readRDS(paste0(model_dir, "param_all_react.rds"))
all_react = all_react[[1]]

Ib = dplyr::filter(all_react, ReactionDistrict == "boundary")

write.table(sort(Ib$ReactionPos, decreasing = FALSE), 
          file = paste0(wd, "/input/csv/react_index.txt"), 
          row.names = FALSE, col.names = FALSE)

system("rm -f dockerID *error.log *.log ExitStatusFile")

therapy_types <- c("NoDrug", "Therapy")
param_target = c("IECsDeath", "Death4Treat", "Detox")

load(model_R.file)

######

for ( ablation in c("Unified", "Ablated", "SuperAblated") ) {

  Exe.exp(model_name = model.name,
          fba_fname = paste0(model.name, ".txt"),
          atol = 1e-06, 
          rtol = 1e-06,
          time.step = 0.5, 
          f_time = 72,
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
          replicates = 8,
          param_target = param_target,
          net_fname = net_fname,
          cores = 4,
          wd = wd,
          supp_function.dir = supp_function_dir)
  
  setwd(wd)
}

######

Exper = "Model_Sensitivity"
Ablation = "Unified"
therapy = "Therapy"

setwd(paste0(
  wd, "/results/CDiff_", Ablation, "_", therapy, 
  "_", Exper, "/", net_fname, "_sensitivity"))

id = as.integer(readLines("ConfigID.txt"))

setwd(paste0(wd, "/results/CDiff_", Ablation, "_", therapy, "_", Exper, "/"))
load(paste0(net_fname, "-analysis.RData"))

id_param <- sapply(param_target, function(x) {
which(sapply(config, function(y) y[[1]][[1]] == x))
})

for (p in seq_along(param_target)) {
assign(param_target[p], config[[id_param[p]]][[id]][[3]])
}

csv <- read.csv(file = paste0(wd, "/input/csv/ParametersListSinkHeme.csv"), 
                header = F, quote = "")

par = c(paste0("c; IECsDeath; ", IECsDeath),
        paste0("g; Death4Treat; ", Death4Treat),
        paste0("c; Detox; ", Detox))

for(p in 1:length(id_param)) {csv[id_param[p], ] = par[p]}

write.table(csv, file = paste0(wd, "/input/csv/ParametersListSinkHeme.csv"), 
            col.names = FALSE, row.names = F, quote = F)

#####

load(model_R.file)

system("rm -f dockerID *error.log *.log ExitStatusFile")

for (ablation in c("Unified", "Ablated", "SuperAblated")) {
  for (therapy in c("Therapy", "NoDrug")) {
    
    setwd(wd)
    
    Exe.exp(model_name = model.name,
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
            net_fname = net_fname,
            cores = 1,
            wd = wd,
            supp_function.dir = paste0(wd, "/code/supplementary_functions"))
  }
}

#####

setwd(wd)

plotting_ablation(Exper = "Model_Analysis",
                  colConNoDrug = "darkred", 
                  colConTherapy = "darkcyan",
                  coltag = "darkgray",
                  tag = c("Ablated", "Unified"),
                  Condition = c("Therapy", "NoDrug"),
                  wd = wd)

plotting_panel_paper(Exper = "Model_Sensitivity",
                     Condition = "Therapy",
                     tag = c("SuperAblated", "Ablated", "Unified"),
                     param_target = "IECsDeath",
                     Tempi = c(0, 12, 24, 36, 48, 60),
                     colo1 = c("black", "magenta", "gold"),
                     colo2 = c("black", "#266867", "yellow"),
                     colo3 = c("#ffd166", "#ee6c4d", "#293241"),
                     wd)

# ######################## Flux variability analysis #############################
# 
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
# #####
# 
# load("subfvaMaxUnifiedTherapy.RData")
