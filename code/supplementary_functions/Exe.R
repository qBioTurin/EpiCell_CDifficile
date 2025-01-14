
Exe.exp = function(model_cat, 
                   model_name,
                   fba_fname,
                   atol, 
                   rtol,
                   f_time, 
                   event_times,
                   distance_measure,
                   reference_data,
                   Exper, 
                   Ablation, 
                   therapy, 
                   debug, 
                   FVA, 
                   FVAcomparison, 
                   react, 
                   debugFVA,
                   replicates, 
                   param_target,
                   time.step,
                   net_fname,
                   wd,
                   supp_function.dir,
                   cores) {
  
  setwd(wd)
  
  system("rm -f dockerID *error.log *.log ExitStatusFile" )
  
  if(!dir.exists("./results")) dir.create("./results")
  
  tag <- Ablation
  Condition <- therapy
  
  Exper <- Exper
  r <- replicates*cores

  multi_param <- if (length(param_target) > 1) TRUE else FALSE
  
  transitions_fname = switch(tag,
                             "Ablated" = paste0(wd, "/Net/Ablated.cpp"),
                             "ParAblated" = paste0(wd, "/Net/ParAblated.cpp"),
                             "Unified" = paste0(wd, "/Net/Unified.cpp"))
  
  model.generation(net_fname = paste0(wd, "/Net/", net_fname, ".PNPRO"),
                   transitions_fname = transitions_fname,
                   fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname))
  
  system(paste0("mv ", net_fname, ".* ./Net"))
  
  switch(Exper, 
         Model_Analysis = {
           parallel_processors = 1
           n_config = 1
           parameters_fname <- "input/csv/ParametersListSinkHeme.csv"
         },
         Model_Sensitivity = {
           parallel_processors = cores
           n_config = parallel_processors*replicates
           parameters_fname <- "input/csv/ParametersSensitivityAnalysis.csv"
         },
         FBA_Sensitivity = {
           parallel_processors = cores
           n_config = parallel_processors*replicates
           parameters_fname <- "input/csv/ParametersSensitivityAnalysisBounds.csv"
         }
  )
  
  if(Exper == "Model_Analysis") {
    
    model.analysis(solver_fname = paste0(wd, "/Net/", net_fname, ".solver"),
                   i_time = 0, f_time = f_time, s_time = time.step,
                   parallel_processors = parallel_processors,
                   n_config = n_config,
                   FVA = T,
                   fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname),
                   atol = atol, rtol = rtol, 
                   debug = debug,
                   parameters_fname = paste0(wd, "/", parameters_fname),
                   functions_fname = paste0(wd, supp_function.dir, "Functions.R"),
                   event_times = if (Condition != "NoDrug") event_times,
                   event_function = if (Condition != "NoDrug") "treat_generation")
    
    if (FVA) {
      
      model.sensitivity(folder_trace = paste0(wd, "/", net_fname, "_analysis/"),
                        FVA = T,
                        fva_gamma = 0.9,
                        debug = debugFVA,
                        parallel_processors = (parallel::detectCores() - 1),
                        fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname),
                        flux_fname = react,
                        functions_fname = paste0(supp_function.dir, "Functions.R"))
      
      system(paste("rm -r ./results/FVA", tag, Condition, Exper, sep = "_"))
      system(paste("mv Model_sensitivity* ./results/FVA", tag, Condition, Exper, sep = "_"))
      
      if (FVAcomparison) {
        
        source(paste0(supp_function.dir, "glpkVsSybil_4server.R"))
        
        # p = ((deltaVFAmin|deltaVFAmax)/(deltaFVA|deltaFBA))/(pFBA|pFVA)
        # 
        # ggsave(p, file = paste("FVAComp", tag, Condition, ".pdf", sep = ""), 
        #        width = 10, height = 10)
        # 
        # system(paste("cp ./FVAComp", tag, Condition, ".pdf ", "./results/FVA", "_",
        #              tag, "_", Condition, "_", Exper, sep = ""))
        # 
        # system(paste("mv ./FVAComp", tag, Condition, ".pdf ", "./results/FVA", "_",
        #              tag, "_", Condition, "_", Exper, sep = ""))
        
      }
      
    }
    
    system(paste("rm -r ./results/CDiff", tag, Condition, Exper, sep = "_"))
    system(paste0("mv ", net_fname, "_analysis* ./results/CDiff", "_", tag, "_", Condition, "_", Exper))
    
    # source(paste0(wd, supp_function.dir, "tracefluxesbounds.R"))
    # source(paste0(wd, supp_function.dir, "PlottingAnalysisResults.R"))
    
    writedata.csv(wd, tag, Condition, Exper, net_fname)
    # PlottingPanorama(Exper, Condition, tag, pdf.name = "PanoramaPlot",wd)
  }
  
  if(Exper == "Model_Sensitivity") {
    
    model.sensitivity(
        solver_fname = paste0(wd, "/Net/", net_fname, ".solver"),
        i_time = 0, f_time = f_time, s_time = time.step,
        parallel_processors = parallel_processors,
        n_config = n_config,
        fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname),
        atol = atol, 
        rtol = rtol, 
        debug = debug,
        target_value = "CD",
        distance_measure = distance_measure, 
        reference_data = reference_data, 
        parameters_fname = paste0(wd, "/", parameters_fname),
        functions_fname = paste0(wd, supp_function.dir, "Functions.R"),
        event_times = if (Condition != "NoDrug") event_times,
        event_function = if (Condition != "NoDrug") "treat_generation")
    
    ### moving files and folders in results ####
    if(dir.exists(paste0(wd,"/results/CDiff", "_", tag, "_", Condition, "_", Exper))) 
      system(paste0("rm -r ",wd,"/results/CDiff","_", tag, "_", Condition, "_", Exper))
    system(paste0("mv ", net_fname, "_", "analysis* ",wd,"/results/CDiff", "_", tag, "_", Condition, "_", Exper))
    
    if(dir.exists(paste(wd,"/results/CDiff", "_", tag, "_", Condition, "_", Exper,"sensitivity"))) 
      system(paste("rm -r ",wd,"/results/CDiff", "_", tag, "_", Condition, "_", Exper, "sensitivity"))
    system(paste0("mv ", net_fname, "_", "sensitivity*  ",wd,"/results/CDiff", "_", tag, "_", Condition, "_", Exper))
    #####
    
    resfolder = paste0(wd,"/results/CDiff", "_", tag, "_", Condition, "_", Exper)
    load(paste0(resfolder,"/",net_fname, "-analysis.RData"))
    
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
    
    if (multi_param){
      failed.tr <- c()
      #f.tr <- 1
      
      registerDoParallel(cores = detectCores())
      subtrace <- foreach(i = sample(list.files(pattern = "\\.trace$",
                                                path = resfolder, full.names = T), r), .combine = rbind) %dopar% {
        
        TracesPath <- i
        
        n <- as.numeric(gsub(pattern = paste0("(", net_fname, "-analysis-)|(.trace)"), replacement = "", x = basename(i)) )
        
        trace <- read.table(TracesPath, header = TRUE)
        
        if (!all(sapply(trace, is.numeric))) {
          print(i)
          failed.tr[f.tr] <<- n
          f.tr <<- f.tr + 1
          return()
        }
        
        subtrace <- trace %>% 
          dplyr::select(Time, IECs, Drug, CD, Damage, pheme_e, 
                        pheme_c, BiomassCD, pro_L_e, leu_L_e, 
                        ile_L_e, val_L_e, trp_L_e, cys_L_e) %>% 
          tidyr::gather(key = "Places", value = "Marking", -Time)
        
        subtrace <- cbind(subtrace, config = rep(as.character(i), length(subtrace$Time)))
        
        for (p in 1:length(param_target)) {
          subtrace <- cbind(subtrace, p = rep(config[[id_param[p]]][[n]][[3]], length(subtrace$Time)))
        }
        
        return(subtrace)
      }
      stopImplicitCluster()
      
      registerDoParallel(cores = detectCores())
      subflux <- foreach(i = sample(list.files(pattern = "\\.flux$",path = resfolder,full.names = T), r), .combine = rbind) %dopar% {
        TracesPath <- i
        n <- as.numeric(gsub(pattern = paste("(", net_fname, "-analysis-)|(-0.flux)", sep = ""), replacement = "", x = basename(i)) )
        flux <- read.table(TracesPath, header = FALSE)
        ReactionsNames <- gsub("\\(|\\)", "", gsub("\\(", "_", flux[1, -c(1, 2)]))
        flux <- read.table(TracesPath, header = TRUE)
        if (TracesPath %in% paste0(net_fname, "-analysis-", failed.tr, "-0.flux")) {
          print(TracesPath)
          return(NULL)
        }
        colnames(flux) <- c("Time", "Obj_0", ReactionsNames)
        subflux <- flux[, c("Time", ReactionsNames)] %>%
          tidyr::gather(key = "Reaction", value = "Flux", -Time)
        obj <- flux[, c("Time", "Obj_0")]
        subbounds <- dplyr::filter(subflux, grepl("_Ub|_Lb", Reaction))
        obj_bounds <- dplyr::filter(subbounds, grepl("biomass205", Reaction))
        subflux <- dplyr::filter(subflux, !grepl("_Ub|_Lb", Reaction))
        subflux <- cbind(subflux, config = rep(as.character(i), length(subflux$Time)))
        for (p in 1:length(param_target)) {
          subflux <- cbind(subflux, p = rep(config[[id_param[p]]][[n]][[3]], length(subflux$Time)))
        }
        return(subflux)
      }
      stopImplicitCluster()
      
      names(subtrace) = c("Time", "Places", "Marking", "config", param_target)
      names(subflux) = c("Time", "Places", "Marking", "config", param_target)
    } else {
      subtrace <- lapply(sample(list.files(pattern = "\\.trace$",path = resfolder,full.names = T), r), function(i) {
        TracesPath <- i
        n <- as.numeric(gsub(pattern = paste0("(", net_name, "-analysis-)|(.trace)"), replacement = "", x = basename(i)) )
        trace <- read.table(TracesPath, header = TRUE)
        subtrace <- trace %>%
          dplyr::select(Time, IECs, Drug, Damage, 
                        pheme_e, pheme_c, BiomassCD, CD,
                        pro_L_e, leu_L_e, ile_L_e, val_L_e, trp_L_e, cys_L_e) %>%
          tidyr::gather(key = "Places", value = "Marking", -Time)
        subtrace <- cbind(subtrace, 
                          config = rep(as.character(i), length(subtrace$Time)), 
                          param_sen = rep(config[[id_param]][[n]][[3]], length(subtrace$Time)))
        
        return(subtrace)
      })
      subflux <- lapply(sample(list.files(pattern = "\\.flux$",path = resfolder,full.names = T), r), function(i) {
        TracesPath <- i
        n <- as.numeric(gsub(pattern = paste("(", net_name, "-analysis-)|(-0.flux)", sep = ""), replacement = "", x = basename(i)) )
        flux <- read.table(TracesPath, header = FALSE)
        ReactionsNames <- gsub("\\(|\\)", "", gsub("\\(", "_", flux[1, -c(1, 2)]))
        colnames(flux) <- c("Time", "Obj_0", ReactionsNames)
        subflux <- flux[, c("Time", ReactionsNames)] %>%
          tidyr::gather(key = "Reaction", value = "Flux", -Time)
        obj <- flux[, c("Time", "Obj_0")]
        subbounds <- dplyr::filter(subflux, grepl("_Ub|_Lb", Reaction))
        obj_bounds <- dplyr::filter(subbounds, grepl("biomass205", Reaction))
        subflux <- dplyr::filter(subflux, !grepl("_Ub|_Lb", Reaction))
        subflux <- cbind(subflux, 
                         config = rep(as.character(i), length(subflux$Time)), 
                         param_sen = rep(config[[id_param]][[n]][[3]], length(subflux$Time)))
        
        return(subflux)
      })
    }
    saveRDS(subtrace, file = paste0(resfolder,"/subtrace_" , tag, Condition, ".rds"))
    saveRDS(subflux, file = paste0(resfolder,"/subflux_" , tag, Condition, ".rds"))
      
  }
  
}
