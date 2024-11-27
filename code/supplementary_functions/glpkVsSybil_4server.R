
# reading react

# model_name = "CD196HemeSinkDiet"
# model_cat = "CD196HemeSink"
# react = as.vector(model.mat@react_id)
# tag = "Unified"
# Condition = "Therapy"
# Exper = "Model_Analysis"

load(paste0("./Input/CDmodels/", model_cat,"/", model_name, ".RData"))

react = gsub("\\(", replacement = "_", react)
react = gsub("\\)", replacement = "", react)
react = gsub("\\[", replacement = "_", react)
react = gsub("\\]", replacement = "", react)

# reading .trace

trace = read.table(
  "./EpitCellDifficileHemeSink_analysis/EpitCellDifficileHemeSink-analysis-1.trace",
  header = T)

# reading .flux

flux = read.table(
  "./EpitCellDifficileHemeSink_analysis/EpitCellDifficileHemeSink-analysis-1-0.flux",
  header = F)

ReactionsNames = flux[1, -c(1, 2)]
ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
ReactionsNames = gsub("\\[", replacement = "_", ReactionsNames)
ReactionsNames = gsub("\\]", replacement = "", ReactionsNames)

flux = read.table(
  "./EpitCellDifficileHemeSink_analysis/EpitCellDifficileHemeSink-analysis-1-0.flux",
  header = T)

colnames(flux) = c("Time", "Obj_0", ReactionsNames)

subflux = flux[, c("Time", ReactionsNames)] %>%
  tidyr::gather(key = "Reaction", value = "Flux", -Time)

obj = flux[, c("Time", "Obj_0")]
  
subbounds = dplyr::filter(subflux, grepl("_Ub|_Lb", Reaction))
obj_bounds = dplyr::filter(subbounds, grepl("biomass205", Reaction))

subflux = dplyr::filter(subflux, !grepl("_Ub|_Lb", Reaction))

assign("subflux", cbind(subflux, Scenario = rep(tag, length(subflux$Flux))))

react = gsub("\\(", replacement = "_", react)
react = gsub("\\)", replacement = "", react)
react = gsub("\\[", replacement = "_", react)
react = gsub("\\]", replacement = "", react)

subflux = dplyr::filter(subflux, Reaction %in% react)

Sb = dplyr::filter(subbounds, grepl(paste(react, collapse = "|"), Reaction))

SbUb = dplyr::filter(Sb, grepl("_Ub", Reaction))
SbLb = dplyr::filter(Sb, grepl("_Lb", Reaction))

assign("subbounds", cbind(subbounds, Scenario = rep(tag, length(subbounds$Flux))))

# reading FVA solutions

## min 

fvaMin = read.csv(paste("./Results/FVA_", tag, "_", Condition, "_", Exper, 
                        "/EpitCellDifficileHemeSink-fva-1-0-MIN.csv", 
                        sep = ""), header = F)

ReactionsNamesFVA = as.character(fvaMin[1, -c(1, 1)])

ReactionsNamesFVA = gsub("\\(", replacement = "_", ReactionsNamesFVA)
ReactionsNamesFVA = gsub("\\)", replacement = "", ReactionsNamesFVA)

fvaMin = read.csv(paste("./Results/FVA_", tag, "_", Condition, "_", Exper, 
                        "/EpitCellDifficileHemeSink-fva-1-0-MIN.csv", 
                        sep = ""), header = T)

colnames(fvaMin) = c("Time", ReactionsNamesFVA)

subfvaMin = fvaMin[, c("Time", ReactionsNamesFVA)] %>%
  tidyr::gather(key = "Reaction", value = "Flux", -Time)

save(subfvaMin, file = paste("subfvaMin", tag, Condition, ".RData", sep = ""))

assign("subfvaMin", cbind(subfvaMin, Scenario = rep(tag, length(subfvaMin$Flux))))

## max 

fvaMax = read.csv(paste("./Results/FVA_", tag, "_", Condition, "_", Exper, 
                        "/EpitCellDifficileHemeSink-fva-1-0-MAX.csv", 
                        sep = ""), header = F)

ReactionsNamesFVA = as.character(fvaMax[1, -c(1, 1)])

ReactionsNamesFVA = gsub("\\(", replacement = "_", ReactionsNamesFVA)
ReactionsNamesFVA = gsub("\\)", replacement = "", ReactionsNamesFVA)
ReactionsNamesFVA = gsub("\\[", replacement = "_", ReactionsNamesFVA)
ReactionsNamesFVA = gsub("\\]", replacement = "", ReactionsNamesFVA)

fvaMax = read.csv(paste("./Results/FVA_", tag, "_", Condition, "_", Exper, 
                        "/EpitCellDifficileHemeSink-fva-1-0-MAX.csv", 
                        sep = ""), header = T)

colnames(fvaMax) = c("Time", ReactionsNamesFVA)

subfvaMax = fvaMax[, c("Time", ReactionsNamesFVA)] %>%
  tidyr::gather(key = "Reaction", value = "Flux", -Time)

save(subfvaMax, file = paste("subfvaMax", tag, Condition, ".RData", sep = ""))

assign("subfvaMax", cbind(subfvaMax, Scenario = rep(tag, length(subfvaMax$Flux))))
