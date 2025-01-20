
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

##### plot FIG.2 B and C - Main paper ##### 
# source(paste0(wd,"/code/plot_functions/Fig2_plot_old.R") )
# 
# pl = plotting_Fig2_paper(Exper = "Model_Sensitivity",
#                      Condition = "Therapy",
#                      tag = c("Ablated", "ParAblated", "Unified"),
#                      param_target = "IECsDeath",
#                      Tempi = c(0, 12, 24, 36, 48, 60),
#                      colo1 = c("black", "magenta", "gold"),
#                      colo2 = c("black", "#266867", "yellow"),
#                      colo3 = c("#ffd166", "#ee6c4d", "#293241"),
#                      wd )
# pl$pl2C
# pl$pl2B
# 
# ggsave(plot = pl$pl2C,filename = "Fig2C.pdf",path = "Figures/",width = 14,height = 10)
# ggsave(plot = pl$pl2B,filename = "Fig2B.pdf",path = "Figures/",width = 10,height = 14)

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
