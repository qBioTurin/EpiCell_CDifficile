
diet_processing = function(metadata.diet) {
  
  # Reading diet formulation data
  diets_metadata = read_delim(metadata.diet, delim = "\t", escape_double = FALSE)
  
  # adding diet gut microbiota information to the diet metadata
  diets_metadata = rbind(diets_metadata, data.frame(
    ...1 = length(diets_metadata$...1) + 1,
    name = "Diet_default", 
    notes = "Table of all exchange reactions of the defined 
 essential metabolites, mucus glycans, and remaining
 metabolites with their respective concentrations",
    references = "Bauer, Eugen & Zimmermann, Johannes & 
 Baldini, Federico & Thiele, Ines & Kaleta, Christoph. 
 (2017). BacArena: Individual-based metabolic modeling of 
 heterogeneous microbes in complex communities. PLoS 
 computational biology. 13. e1005544. 
 10.1371/journal.pcbi.1005544. ",
    totalEnergyKcal = NA, lipids = NA, 
    lipids_kcal = NA, lipids_perc = NA,
    carbohydrates = NA, carbohydrates_kcal = NA,
    carbohydrates_perc = NA, proteins = NA, 
    proteins_kcal  = NA, proteins_perc = NA, 
    alcohol = NA, alcohol_kcal = NA,
    alcohol_perc = NA, dietFibers_mg = NA, 
    cholesterol_mg = NA, sucrose_mg = NA))
  
  diets = list()
  
  for (diet in 1:(length(diets_metadata$name) - 1)) {
    
    diets[[diet]] = read_delim(paste0("./diets/", diets_metadata$name[diet], ".tsv"), 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
    
    diets[[diet]][["Reaction"]] = gsub("\\(", replacement = "_", diets[[diet]][["Reaction"]])
    diets[[diet]][["Reaction"]] = gsub("\\)", replacement = "", diets[[diet]][["Reaction"]])
    
    diets[[diet]][["Reaction"]] = gsub("\\[", replacement = "_", diets[[diet]][["Reaction"]])
    diets[[diet]][["Reaction"]] = gsub("\\]", replacement = "", diets[[diet]][["Reaction"]])
    
  }
  
  # The exact diet definition with the predicted essential metabolites can be found in:
  diets[[length(diets_metadata$name)]] = read.csv(paste0("diets/Diet_default", ".tsv"), header = T)
  
  diets[[length(diets_metadata$name)]][["Exchange"]] = 
    gsub("\\(", replacement = "_", diets[[length(diets_metadata$name)]][["Exchange"]])
  diets[[length(diets_metadata$name)]][["Exchange"]] = 
    gsub("\\)", replacement = "", diets[[length(diets_metadata$name)]][["Exchange"]])
  
  diets[[length(diets_metadata$name)]][["Exchange"]] = 
    gsub("\\[", replacement = "_", diets[[length(diets_metadata$name)]][["Exchange"]])
  diets[[length(diets_metadata$name)]][["Exchange"]] = 
    gsub("\\]", replacement = "", diets[[length(diets_metadata$name)]][["Exchange"]])
  
  names(diets) = diets_metadata$name
  
  # saving diets data
  save(diets, file = "./diets/diets.RData")
  
}
