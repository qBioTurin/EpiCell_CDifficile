
library(readr)
library(dplyr)

wd  = "./Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/"
setwd(wd)

source("./Rfunction/epimod_FBAfunctions-main/R/class_generation.R")
source("./Rfunction/epimod_FBAfunctions-main/R/FBAgreatmodeClass.R")
source("./Rfunction/epimod_FBAfunctions-main/R/readMat.R")

model = FBA4Greatmod.generation(fba_mat = "./Input/CDmodels/CD196HemeSink/CD196HemeSink.mat")

diet = read_delim("./Input/Diets/vmh/EU_average.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

diet$Reaction = gsub("\\[", replacement = "\\(", diet$Reaction)
diet$Reaction = gsub("\\]", replacement = "\\)", diet$Reaction)

diet = diet %>% 
  rename(lwbnd = `Flux Value`) %>%
  mutate(lwbnd = -(lwbnd)/24, uppbwnd = 10)

model = setDiet(model, dietf = diet)
model = setDiet.name(model,diet_name = "EU_average")

for(r in model@react_id) model = setConstraints(model, reaction.name = r, newConstraints = c(-0.1, 10))

save(model, "./Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")
writeFBAfile(model, "./Input/CompiledModels/FBAModelEU_average")
