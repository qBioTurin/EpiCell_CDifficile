
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
library(doParallel)
library(foreach)

# Function to set working directory to the closest folder with a specific name
set_closest_directory <- function(target_name) {
  wd = getwd()
  if( basename(wd) != target_name){
    # Get all directories relative to the current working directory
    all_dirs <- list.dirs(path = "~", full.names = TRUE, recursive = TRUE)
    
    # Filter directories with the specific name
    matching_dirs <- all_dirs[basename(all_dirs) == target_name]
    
    # Check if there are matches
    if (length(matching_dirs) == 0) {
      stop("No folder with the name '", target_name, 
           "' found in the current working directory or its subdirectories.")
    }
    
    # Find the directory with the shortest path (closest to the current working directory)
    closest_dir <- matching_dirs[which.min(nchar(matching_dirs))]
    
    # Set the working directory
    setwd(closest_dir)
  }
  # Return the new working directory
  return(getwd())
}

wd <- set_closest_directory("EpiCell_CDifficile")

setwd(wd)
supp_function.dir = "/code/supplementary_functions/"

# loading experiments executioner
source(paste0(wd, supp_function.dir, "Exe.R"))

# setting FBA model tags
model.type = "Clostridium"
model.name = "CD196HemeSink"

###### FBA model generation ######
model_R.file = paste0(wd, "/input/", model.type, "/", model.name, "/", model.name, ".RData")

if(!file.exists(model_R.file)) {
  # FBA model compiling phase
  source(paste0(wd, supp_function.dir, "Functions.R"))
  source(paste0(wd, supp_function.dir, "FBAmodel.generation.R"))
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
}

if(!file.exists(paste0(wd, "/input/csv/react_index.txt"))){
  all.react = readRDS(paste0(wd, "/input/", model.type, "/", model.name, "/all.react_diet.rds"))
  all.react = all.react[[1]]
  
  Ib = dplyr::filter(all.react, ReactionDistrict == "boundary")
  
  write.table(sort(Ib$ReactionPos, decreasing = FALSE), 
              file = paste0(wd, "/input/csv/react_index.txt"), 
              row.names = FALSE, col.names = FALSE)
  
}

######

system("rm -f dockerID *error.log *.log ExitStatusFile")

startTime = Sys.time()
for (ablation in c("Unified", "ParAblated", "Ablated")) {
  Exe.exp(model_cat = model.name,
          model_name = model.name,
          fba_fname = paste0(model.name, ".txt"),
          atol = 1e-06, rtol = 1e-06,
          time.step = 1, f_time = 72,
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
          replicates = 21,
          param_target = c("IECsDeath", "Death4Treat", "Detox"),
          net_fname = "EpitCellDifficileHemeSink",
          cores = detectCores(),
          wd = wd,
          supp_function.dir = supp_function.dir)
}
endTime = Sys.time()

# Parameters closest to the median trances in FIG2B ####
source("./code/supplementary_functions/MedianParamsConfiguration.R")

resParams = configuration.closeMedian(Condition = "Therapy", 
                                      tag = "Unified", numberOfSets = 3,
                                      subtrace = NULL,
                                      places2plot = c("CD", "IECs", "leu_L_e", "pheme_c"))

plot_config = resParams$plot
plot_config

saveRDS(resParams, file = "paramsProve.rds")

## Run the ConfSet analysis
numbConfig = 1:dim(resParams$Config)[1]
paramsConfig = resParams$Config
paramsgrid = rbind(
  expand.grid("1e-6", numbConfig, "Ablated"), 
  expand.grid("1e-6", numbConfig, "ParAblated"),
  expand.grid("1e-6", numbConfig, "Unified"))

MultipleAnalysis = lapply(seq_along(paramsgrid[, 1]),
                          function(c, paramsConfig, paramsgrid) {
                            # browser()
                            new_eps_value = paramsgrid[c, "Var1"]
                            nconfig = paramsgrid[c, "Var2"]
                            tag = paramsgrid[c, "Var3"]
                            
                            # Debug print 1: Print current configuration
                            cat(sprintf("\n\n=== Running configuration %d ===\n", c))
                            cat(sprintf("eps: %s, config: %d, tag: %s\n", new_eps_value, nconfig, tag))
                            cat("Parameters being set:\n")
                            print(paramsConfig[nconfig, c("Detox", "Death4Treat", "IECsDeath", "ConfParams")])
                            
                            model.generation(net_fname = paste0(wd, "/Net/", net_fname, ".PNPRO"),
                                             transitions_fname = paste0(wd, "/Net/", tag, ".cpp"),
                                             fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname))
                            
                            system(paste0("mv ", net_fname, ".* ./Net"))
                            
                            execution_start <- Sys.time()
                            model.analysis(solver_fname = paste0(wd, "/Net/", net_fname, ".solver"),
                                           i_time = 0,
                                           f_time = f_time,
                                           s_time = time.step,
                                           fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname),
                                           atol = atol,
                                           rtol = rtol,
                                           ini_v = as.numeric(paramsConfig[nconfig,c("Detox","Death4Treat","IECsDeath")]),
                                           parameters_fname = paste0(wd, "/", parameters_fname),
                                           functions_fname = paste0(wd, supp_function.dir, "Functions.R"),
                                           event_times = if (Condition != "NoDrug") event_times,
                                           event_function = if (Condition != "NoDrug") "treat_generation")
                            execution_end <- Sys.time()
                        
                            resFolder = paste0("./results/TimingEval/CDiff", tag, Condition, new_eps_value,
                                               gsub(pattern = " ", replacement = "", x = paramsConfig[nconfig, "ConfParams"]),
                                               collapse  = "_")
                            
                            if(dir.exists(resFolder)) system(paste("rm -r ", resFolder))
                            
                            system(paste0("mv ", net_fname, "_analysis* ", resFolder))
                            
                            traces = ModelAnalysisPlot(
                              paste0(resFolder,"/EpitCellDifficileHemeSink-analysis-1.trace"),
                              paste0(resFolder,"/EpitCellDifficileHemeSink-analysis-1-0.flux"), 
                              FluxVec = c("EX_biomass_e", "sink_pheme_c", "EX_cys_L_e", 
                                          "EX_trp_L_e", "EX_ile_L_e", "EX_pro_L_e", "EX_leu_L_e", "EX_val_L_e"),
                              new_eps_value,
                              paramsConfig[nconfig, "ConfParams"],
                              tag
                            )
                            
                            debug_plot = ggplot(
                              traces[["subtrace"]] %>% 
                                filter(Places %in% c("IECs", "pheme_e", "BiomassCD", "CD", "pheme_c"))) +
                              geom_line(aes(x = Time, y = Marking,
                                            linetype = tag),
                                        linewidth = 1) +
                              facet_wrap(~Places, scales = "free_y") +
                              labs(title = sprintf("Configuration %d: %s, eps=%s", 
                                                   nconfig, 
                                                   paramsConfig[nconfig, "ConfParams"],
                                                   new_eps_value),
                                   subtitle = sprintf("Detox=%.2e, Death4Treat=%.2e, IECsDeath=%.2e",
                                                      paramsConfig[nconfig, "Detox"],
                                                      paramsConfig[nconfig, "Death4Treat"],
                                                      paramsConfig[nconfig, "IECsDeath"])) +
                              theme_minimal() +
                              theme(legend.position = "bottom")
                            
                            # Save the debug plot
                            ggsave(filename = paste0(wd, "/results/", new_eps_value, "_",
                                                     "_Set_", nconfig, "_", tag, "_debug_plot.pdf"), 
                                   plot = debug_plot,
                                   width = 12, height = 8)
                            
                            return(list(traces = traces))
                          },
                          paramsConfig = resParams$Config, paramsgrid
)

flux = do.call(rbind, lapply(lapply(MultipleAnalysis,"[[", 1),"[[", 1))
trajectories = do.call(rbind, lapply(lapply(MultipleAnalysis, "[[",1), "[[", 2))
trajectories$new_eps_value = factor(trajectories$new_eps_value, levels = epstimes)
flux$new_eps_value = factor(flux$new_eps_value, levels = epstimes)
trajectories$config = factor(trajectories$config, levels = paste0("Set ", numbConfig))
flux$config = factor(flux$config, levels = paste0("Set ", numbConfig))
trajectories$ConfParams = trajectories$config

saveRDS(trajectories, file = "trajectories.rds")

##### plot FIG.2 B and C - Main paper ##### 
trajectories = readRDS(file = "trajectories.rds")

colors_new_confParams <- c("#3B9AB2", "#ff80e3ff", "#6e2076")
names(colors_new_confParams) <- unique(trajectories$ConfParams)

colors_new_eps_value <- rev(grey.colors(length(epstimes)))
names(colors_new_eps_value) <- unique(trajectories$new_eps_value)

source("./code/plot_functions/Fig2_plot.R")
Fig2 = plotting_Fig2_paper(Exper = "Model_Sensitivity",
                           Condition = "Therapy",
                           tag = c("Ablated", "ParAblated", "Unified"),
                           param_target = "IECsDeath",
                           Tempi = c(0, 12, 24, 36, 48, 60),
                           wd = wd,
                           trajectories = trajectories,
                           colConfigSets = colors_new_confParams,
                           variables_to_plot = c("CD", "IECs", "pheme_c", "leu_L_e", "trp_L_e"))

Fig2$pl2B
Fig2$pl2C

ggsave(plot = Fig2$pl2B,filename = "Figures/Fig2B.pdf", width = 21/2, height = 23/2)
ggsave(plot = Fig2$pl2C,filename = "Figures/Fig2C.pdf", width = 48/3, height = 30/3)

##### Extrapolate configID for the model analysis ##### 
##### considering the unified scenario from the sensitivity analysis

source(paste0(wd, "/code/supplementary_functions/ConfigID_generation.R"))
ListID = csvFileUpdating.minRankConfigID(wd, net_fname = "EpitCellDifficileHemeSink")

###### Run model Analysis with the new config ID ########

system("rm -f dockerID *error.log *.log ExitStatusFile")

ablation_types <- c("Unified", "ParAblated", "Ablated")
therapy_types <- c("NoDrug", "Therapy")

source(paste0(wd,"/code/supplementary_functions/write_tracefluxesbounds.R") )

for (ablation in ablation_types) {
  print(paste("#### " ,ablation)) 
  for (therapy in therapy_types) {
    print(paste("## " ,therapy)) 
    Exe.exp(model_cat = model.name,
            model_name = model.name,
            fba_fname = paste0(model.name, ".txt"),
            atol = 1e-06, rtol = 1e-06,
            time.step = 1, f_time = 72,
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

##### plot Fig S4 ##### 

source(paste0(wd,"/code/plot_functions/FigS4_plot.R"))

FigS4 = PlottingFigS4(Exper = "Model_Analysis", 
                      colConNoDrug = "darkred", 
                      colConTherapy = "darkcyan", 
                      coltag = "darkgray",
                      tag = c("Ablated", "Unified"),
                      Condition = c("NoDrug", "Therapy"),
                      wd = wd)

ggsave(plot = FigS4[[1]], filename = "Allplaces.pdf", path = "Figures/", width = 11, height = 12)
ggsave(plot = FigS4[[2]], filename = "Fluxes.pdf", path = "Figures/", width = 11, height = 6)

#### plotting biomass

Exper = "Model_Sensitivity"
Condition = "Therapy"
tag = c("Ablated", "ParAblated", "Unified")
param_target = "IECsDeath"
Tempi = c(0, 12, 24, 36, 48, 60)
colo1 = c("black", "magenta", "gold")
colo2 = c("black", "#266867", "yellow")
colo3 = c("#ffd166", "#ee6c4d", "#293241")

Na = 6.022e20
c = 6.022e08
pack = 1*(Na*(1/c))
fc = 1e-06

coloTag = data.frame(Ablated = "darkred", ParAblated = "#EC6D67", Unified = "#5351A2")
nameTag = data.frame(Ablated = "Ablated", ParAblated = "Partially Ablated", Unified = "Unified")

col_settings = c("darkred", "#EC6D67", "#5351A2")

aa_places <- c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
places <- c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", aa_places)
units = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))

subtrace <- do.call(rbind, lapply(tag, function(j) {
  subtrace = readRDS(
    file = paste0(wd, paste0("/results/CDiff", "_", j, "_", Condition, "_", Exper),"/subtrace_" , j, Condition, ".rds"))
  cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
}))

subtrace %>%
  dplyr::filter(Places == "BiomassCD") %>%
  dplyr::select(-c("Death4Treat", "Detox")) %>%
  dplyr::mutate(Time = as.numeric(Time))

df <- subtrace %>%
  dplyr::filter(Places == "BiomassCD") %>%
  dplyr::select(-c("Death4Treat", "Detox")) %>%
  dplyr::mutate(Time = as.numeric(Time))

p_bio = ggplot(df, aes(x = Time, y = Marking, color = Scenario)) +
  geom_line(aes(group = interaction(config, Scenario)), alpha = 0.08) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 9, face = "bold", color = "black"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
        plot.title.position = "plot",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values = col_settings)

ggsave(plot = p_bio, filename = "p_bio.pdf", path = "Figures/", width = 4, height = 2)
