
flux_balance <- function(model.name,
                         model.type, 
                         diet, 
                         diet.name, 
                         uppbwnd.diet,
                         custom,
                         gene.exp,
                         GSE,
                         react.vec,
                         lb.vec,
                         ub.vec, 
                         wd) {
  
  model.file = paste0(wd, "/data/", model.type, "/", model.name, "/", model.name, ".RData")
  
  ## black model
  if(!file.exists(model.file)) {
    
    source(paste0(wd, "/code/epimod_FBAfunctions/R/class_generation.R"))
    source(paste0(wd, "/code/epimod_FBAfunctions/R/FBAgreatmodeClass.R"))
    source(paste0(wd, "/code/epimod_FBAfunctions/R/readMat.R"))
    
    model = FBA4Greatmod.generation(
      fba_mat = paste0(wd, "/data/", model.type, "/", model.name, "/", model.name, ".mat"))
    
    # moving files
    system(paste0("cd ", wd))
    
    system(paste0("mv ./geneAssociation.rds ./data/", model.type, "/", model.name, "/geneAssociation.rds"))
    system(paste0("mv ./geni.rds ./data/", model.type, "/", model.name, "/geni.rds"))
    system(paste0("mv ./officialName.rds ./data/", model.type, "/", model.name, "/officialName.rds"))
    system(paste0("mv ./subsystem.rds ./data/", model.type, "/", model.name, "/subsystem.rds"))
    system(paste0("mv ./genesFromGeneAss.rds ./data/", model.type, "/", model.name, "/genesFromGeneAss.rds"))
    
    source(paste0(wd, "/code/processing_functions/ExtractEx.R"))
    source(paste0(wd, "/code/processing_functions/createReactionString.R"))
    source(paste0(wd, "/code/processing_functions/modelorg2tsv.R"))
    
    metadataEQ(model = model, model.name = model.name,
               prefix = paste("BIGGdata", model.name, sep = "_"), 
               suffix = "tsv",
               extMetFlag = "b", 
               fielddelim = "\t", 
               entrydelim = ", ", 
               makeClosedNetwork = FALSE, 
               onlyReactionList = FALSE, 
               minimalSet = TRUE, 
               fpath = paste0(wd, "/data/", model.type,"/", model.name))
    
    bigg.path = paste0(wd, "/data/", model.type, "/", model.name, "/BIGGdata_", model.name, "_react.tsv")
    
    BIGGdata = read.delim2(bigg.path)
    saveRDS(BIGGdata$equation, file = paste0(wd, "/data/", model.type,"/", model.name, "/equation.rds"))
    
    ex = ExtractEx(model = model, model.name = model.name, bigg.path = bigg.path)
    
    saveRDS(ex, file = paste0(wd, "/data/", model.type, "/", model.name, "/boundary_react.rds"))
    save(model, file = paste0(wd, "/data/", model.type,"/", model.name, "/", model.name, ".RData"))
    
  }
  
  load(model.file)
  
  ## diet effect
  if (diet == T) {
    
    model@react_id = gsub("\\(", replacement = "_", model@react_id)
    model@react_id = gsub("\\[", replacement = "_", model@react_id)
    model@react_id = gsub("\\)", replacement = "", model@react_id)
    model@react_id = gsub("\\]", replacement = "", model@react_id)
    
    metadata.diet = paste0(wd, "/diets/diets_metadata.tsv")
    
    if(!file.exists(paste0(wd, "/diets/diets.RData"))) {
      
      source(paste0(wd, "/code/processing_functions/DietDataProcessing.R"))
      diet_processing(metadata.diet)
      load(paste0(wd, "/diets/diets.RData"))
      
    } 
    
    load(paste0(wd, "/diets/diets.RData"))
    
    equation = readRDS(paste0(wd, "/data/", model.type, "/", model.name, "/equation.rds"))
    officialName = readRDS(paste0(wd, "/data/", model.type, "/", model.name, "/officialName.rds"))
    geneAssociation = readRDS(paste0(wd, "/data/", model.type, "/", model.name, "/geneAssociation.rds"))
    subsystem = readRDS(paste0(wd, "/data/", model.type, "/", model.name, "/subsystem.rds"))
    
    diet = diets[[diet.name]] %>%
      dplyr::rename(lwbnd = `Flux Value`) %>%
      mutate(lwbnd = -((lwbnd)/24), uppbwnd = uppbwnd.diet)
    
    model = setDiet(model, dietf = diet)
    model = setDiet.name(model, diet_name = diet.name)
    
  }
  
  # gene espression integration
  if (gene.exp == T) {
    
    for (i in 1:length(GSE)) {
      
      if (GSE[i] >= 1) {
        
        model@lowbnd[i] = model@lowbnd[i] * (1 + log(GSE[i]))
        model@uppbnd[i] = model@uppbnd[i] * (1 + log(GSE[i]))
        
      } else {
        
        model@lowbnd[i] = model@lowbnd[i] / (1 + abs(log(GSE[i])))
        model@uppbnd[i] = model@uppbnd[i] / (1 + abs(log(GSE[i])))
        
      }
    }
  }
  
  ## custom bounds
  if (custom == T) {
    for(r in 1:length(react.vec)) {
      model = setConstraints(
        model,
        reaction.name = react.vec[r],
        newConstraints = c(lb.vec[r], ub.vec[r]))
    }
  }
  
  # fba solutions
  
  equation = readRDS(paste0(wd, "/data/", model.type, "/", model.name, "/equation.rds"))
  officialName = readRDS(paste0(wd, "/data/", model.type, "/", model.name, "/officialName.rds"))
  geneAssociation = readRDS(paste0(wd, "/data/", model.type,"/", model.name, "/geneAssociation.rds"))
  subsystem = readRDS(paste0(wd, "/data/", model.type, "/", model.name, "/subsystem.rds"))
  
  modelFBA = data.frame(abbreviation = model@react_id,
                        lowbnd = model@lowbnd,
                        uppbnd = model@uppbnd,
                        obj_coef = model@obj_coef,
                        equation = equation,
                        offcialName = officialName,
                        geneAssociation = geneAssociation,
                        subsystem = subsystem)
  
  fbasol = find_fluxes_df(modelFBA)
  
  return(fbasol)
  
}
