
HBInteraction <- function(bac.model,
                          human.model,
                          model.name.h,
                          model.name.b,
                          model.type.h,
                          model.type.b,
                          fbasol.hs) {
  
  # looking for human-bacterial shared exchange reactions
  source("./code/processing_functions/ExtractEx.R")
  
  bigg.path.h = paste0("./data/", model.type.h, "/", model.name.h, 
                       "/BIGGdata_", model.name.h, "_react.tsv")
  
  bigg.path.b = paste0("./data/", model.type.b, "/", model.name.b, 
                       "/BIGGdata_", model.name.b, "_react.tsv")
  
  BIGGdata.h = read.delim2(bigg.path.h)
  ex.h = ExtractEx(human.model, bigg.path = bigg.path.h)
  
  BIGGdata.b = read.delim2(bigg.path.b)
  ex.b = ExtractEx(bac.model, bigg.path = bigg.path.b)
  
  comm.ex = ex.h$react.id[which(ex.h$react.id %in% ex.b$react.id)]
  
  cat("number of common reactions = ", length(comm.ex))
  
  comm.EXfbasol.h = list()
  
  #  extracting solutions
  comm.EXfbasol.h$index = dplyr::filter(ex.h, react.id %in% comm.ex)$index
  comm.EXfbasol.h$name  = comm.ex
  comm.EXfbasol.h$v = fbasol.hs[comm.EXfbasol.h$index]
  # Showing the direction in which the reaction proceeds
  comm.EXfbasol.h$direction = sign(comm.EXfbasol.h$v)
  
  # Changing bacterial model's constraints according to the human FBA solutions

  reacts = which(bac.model@react_id %in% comm.EXfbasol.h$name)
  
  # ricordarsi di aggiungere dei weights basati sui valori dei flussi stimati dal modello umano
  
  for (i in 1:length(comm.EXfbasol.h$v)) {
    
    if (comm.EXfbasol.h$direction[i] == 1) {
      
      bac.model@lowbnd[reacts][i] = -1000 * (1 - comm.EXfbasol.h$v[i]/(bac.model@lowbnd[reacts][i]))
      bac.model@uppbnd[reacts][i] = 1000 * (1 - comm.EXfbasol.h$v[i]/(bac.model@uppbnd[reacts][i]))
      
      } else if (comm.EXfbasol.h$direction[i] == -1) {
      
        bac.model@lowbnd[reacts][i] = 1/(- 1000 * (1 - comm.EXfbasol.h$v[i]/(bac.model@lowbnd[reacts][i])))
        bac.model@uppbnd[reacts][i] = 1/(1000 * (1 - comm.EXfbasol.h$v[i]/(bac.model@uppbnd[reacts][i])))
        
        } else {
          bac.model@lowbnd[reacts][i] = bac.model@lowbnd[reacts][i]
          bac.model@uppbnd[reacts][i] = bac.model@uppbnd[reacts][i]
        }
    }
  
  save(bac.model, file = paste0("./data/", model.type.b, "/", model.name.b, "/", model.name.b, "_H_B.RData"))
  # investigating human metabolism metabolism performing FBA
  
  source("./code/processing_functions/flux_balance.R")
  
  modelFBA.bac = data.frame(abbreviation = bac.model@react_id,
                            lowbnd = bac.model@lowbnd,
                            uppbnd = bac.model@uppbnd,
                            equation = BIGGdata.b$equation,
                            obj_coef = bac.model@obj_coef)
                        
  fbasol.bac = find_fluxes_df(modelFBA.bac)
  
  comm.ReactFlux_bac = list()
  
  f = dplyr::filter(fbasol.bac, abbreviation %in% comm.EXfbasol.h$name)
  
  comm.ReactFlux_bac$name = f$abbreviation
  comm.ReactFlux_bac$direction = sign(f$flux)
  comm.ReactFlux_bac$v = f$flux
  
  return(comm.ReactFlux_bac)
  
}
 