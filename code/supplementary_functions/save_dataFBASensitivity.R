
save.dataFBASen = function(tag,
                           Condition, 
                           Exper,
                           r,
                           multi_param,
                           net_fname,
                           vecParam) {
  
  rep.row = function(x, n){
    matrix(rep(x, each = n), nrow = n)
    }
  
  if (multi_param) {
    setwd(paste0(wd, "/Results/CDiff", tag, Condition, Exper))
  } else {
    setwd(paste0(wd, "/Results/CDiff_", tag, "_", 
                Condition, "_", Exper))
  }
  
  numb = r
  
  load(paste("./", net_fname, "-analysis.RData", sep = ""))
  
  net_name = net_fname
  
  if (multi_param) {
    id_param = c()
    for (p in 1:length(param_target)) {
      for(i in 1:length(config)) {
        if(param_target[p] == config[[i]][[1]][[1]]) {
          id_param[p] = i
        }
      }
    }
  } else {
    for(i in 1:length(config)) {
      if(param_target == config[[i]][[1]][[1]]) {
        id_param = i
      }
    }
  }
  
  failed.tr <<- c()
  f.tr <<- 1
  
  if (vecParam) {
    
    list.Fl = list.files(pattern = "\\.flux$")
    
    subflux = lapply(sample(list.Fl, numb), function(i){
      
      TracesPath = i
      
      n = as.numeric(gsub(pattern = paste("(", net_name, "-analysis-)|(-0.flux)", sep = ""), 
                          replacement = "", x = i))
      
      flux = read.table(TracesPath, header = F)
      
      ReactionsNames = flux[1, -c(1, 2)]
      ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
      ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
      
      flux = read.table(TracesPath, header = T)
      
      colnames(flux) = c("Time", "Obj_0", ReactionsNames)
      
      subflux = flux[, c("Time", ReactionsNames)] %>%
        tidyr::gather(key = "Reaction", value = "Flux", -Time)
      
      obj = flux[, c("Time", "Obj_0")]
      
      subbounds = dplyr::filter(subflux, grepl("_Ub|_Lb", Reaction))
      obj_bounds = dplyr::filter(subbounds, grepl("biomass205", Reaction))
      
      subflux = dplyr::filter(subflux, !grepl("_Ub|_Lb", Reaction) & Time == 0.0)
      subflux = dplyr::filter(subflux, Time == 0.0)
      
      cf = rep.row(c(config[[id_param]][[n]][3][[1]], n), dim(subflux)[1])
      
      colnames(cf) = c(paste(param_target, 1:length(config[[id_param]][[n]][3][[1]]), sep = "."), "config")
      data.config = as.data.frame(cf)
      subflux = cbind(subflux, data.config)
      
      subflux = subflux[, -1]
      
      return(subflux)
      
    })
    
    subflux = do.call(rbind, subflux)
    
  } else {
    
    if (multi_param){
      
      list.Tr = list.files(pattern = "\\.trace$")
      subtrace = lapply(sample(list.Tr, numb), function(i){
        
        TracesPath = i
        
        n = as.numeric(gsub(pattern = paste("(", net_name, "-analysis-)|(.trace)", sep = ""), 
                            replacement = "", x = i))
        
        trace = read.table(TracesPath, header = T)
        
        if(!all(sapply(trace, is.numeric))) {
          print(i)
          failed.tr[f.tr] <<- n
          f.tr <<- f.tr + 1
          return()
        }
        
        subtrace = trace %>% 
          dplyr::select(Time, IECs, Drug, CD, Damage, pheme_e, 
                        pheme_c, BiomassCD, pro_L_e, leu_L_e, 
                        ile_L_e, val_L_e, trp_L_e, cys_L_e) %>% 
          tidyr::gather(key = "Places", value = "Marking", -Time)
        
        subtrace = cbind(subtrace, config = rep(as.character(i), length(subtrace$Time)))
        
        for (p in 1:length(param_target)) {
          subtrace = cbind(subtrace, p = rep(config[[id_param[p]]][[n]][[3]], length(subtrace$Time)))
        }
        
        return(subtrace)
        
      })
      
      list.Fl = list.files(pattern = "\\.flux$")
      subflux = lapply(sample(list.Fl, numb), function(i){
        
        TracesPath = i
        
        n = as.numeric(gsub(pattern = paste("(", net_name, "-analysis-)|(-0.flux)", sep = ""), 
                            replacement = "", x = i))
        
        flux = read.table(TracesPath, header = F)
        
        ReactionsNames = flux[1, -c(1, 2)]
        ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
        ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
        
        flux = read.table(TracesPath, header = T)
        
        if(TracesPath %in% paste(net_name,"-analysis-", failed.tr, "-0.flux", sep = "")) {
          print(TracesPath)
          return()
        }
        
        colnames(flux) = c("Time", "Obj_0", ReactionsNames)
        
        subflux = flux[, c("Time", ReactionsNames)] %>%
          tidyr::gather(key = "Reaction", value = "Flux", -Time)
        
        obj = flux[, c("Time", "Obj_0")]
        
        subbounds = dplyr::filter(subflux, grepl("_Ub|_Lb", Reaction))
        obj_bounds = dplyr::filter(subbounds, grepl("biomass205", Reaction))
        
        subflux = dplyr::filter(subflux, !grepl("_Ub|_Lb", Reaction))
        
        subflux = cbind(subflux, config = rep(as.character(i), length(subflux$Time)))
        
        for (p in 1:length(param_target)) {
          subflux = cbind(subflux, p = rep(config[[id_param[p]]][[n]][[3]], length(subflux$Time)))
        }
        
        return(subflux)
        
      })
      
    } else {
      
      list.Tr = list.files(pattern = "\\.trace$")
      subtrace = lapply(sample(list.Tr, numb), function(i){
        
        TracesPath = i
        
        n = as.numeric(gsub(pattern = paste("(", net_name, "-analysis-)|(.trace)", sep = ""), 
                            replacement = "", x = i))
        
        trace = read.table(TracesPath, header = T)
        
        subtrace = trace %>% 
          dplyr::select(Time, IECs, Drug, Damage, pheme_e, 
                        pheme_c, BiomassCD, CD, pro_L_e, leu_L_e, 
                        ile_L_e, val_L_e, trp_L_e, cys_L_e) %>% 
          tidyr::gather(key = "Places", value = "Marking", -Time)
        
        subtrace = cbind(subtrace, 
                         config = rep(as.character(i), length(subtrace$Time)),
                         param_sen = rep(config[[id_param]][[n]][[3]], length(subtrace$Time)))
        
        return(subtrace)
        
      })
      
      list.Fl = list.files(pattern = "\\.flux$")
      subflux = lapply(sample(list.Fl, numb), function(i){
        
        TracesPath = i
        
        n = as.numeric(gsub(pattern = paste("(", net_name, "-analysis-)|(-0.flux)", sep = ""), 
                            replacement = "", x = i))
        
        flux = read.table(TracesPath, header = F)
        
        ReactionsNames = flux[1, -c(1, 2)]
        ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
        ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
        
        flux = read.table(TracesPath, header = T)
        
        colnames(flux) = c("Time", "Obj_0", ReactionsNames)
        
        subflux = flux[, c("Time", ReactionsNames)] %>%
          tidyr::gather(key = "Reaction", value = "Flux", -Time)
        
        obj = flux[, c("Time", "Obj_0")]
        
        subbounds = dplyr::filter(subflux, grepl("_Ub|_Lb", Reaction))
        obj_bounds = dplyr::filter(subbounds, grepl("biomass205", Reaction))
        
        subflux = dplyr::filter(subflux, !grepl("_Ub|_Lb", Reaction))
        
        subflux = cbind(subflux, config = rep(as.character(i), length(subflux$Time)),
                        param_sen = rep(config[[id_param]][[n]][[3]], length(subflux$Time)))
        
        return(subflux)
        
      })
      
    }
    
    subtrace = do.call(rbind, subtrace)
    
    if(multi_param) {
      names(subtrace) = c("Time", "Places", "Marking", "config", param_target)
    }
    
    subflux = do.call(rbind, subflux)
    
    if(multi_param) {
      names(subflux) = c("Time", "Places", "Marking", "config", param_target)
    }
    
    setwd(paste0(wd, "/Results/CDiff_", tag, "_", Condition, "_", Exper))
    
    save(subtrace, file = paste("subtrace_" , tag, Condition, ".RData", sep = ""))
    
  }
  
  setwd(paste0(wd, "/Results/CDiff_", tag, "_", Condition, "_", Exper))
  
  save(subflux, file = paste("subflux_", tag, Condition, ".RData", sep = ""))
  
}