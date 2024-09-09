
Exe.exp = function(model_name,
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
  
  system("rm -f dockerID *error.log *.log ExitStatusFile" )
  
  tag <<- Ablation
  Condition <<- therapy
  Exper <<- Exper
  r <<- replicates*cores
  model_name <<- model_name
  f_time <<- f_time
  time.step <<- time.step
  param_target <<- param_target
  distance_measure <<- distance_measure
  reference_data <<- reference_data
  multi_param <<- if (length(param_target) > 1) TRUE else FALSE
  
  transitions_fname = switch(tag,
                             "SuperAblated" = paste0(wd, "/net/SuperAblated.cpp"),
                             "Ablated" = paste0(wd, "/net/Ablated.cpp"),
                             "Unified" = paste0(wd, "/net/Unified.cpp"))
  
  model.generation(net_fname = paste0(wd, "/net/", net_fname, ".PNPRO"),
                   transitions_fname = transitions_fname,
                   fba_fname = paste0(wd, "/input/compiled_models/", fba_fname))
  
  system(paste0("mv ", net_fname, ".* ./net"))
  
  switch(Exper, 
         Model_Analysis = {
           parallel_processors = 1
           n_config = 1
           parameters_fname <<- "input/csv/ParametersListSinkHeme.csv"
         },
         Model_Sensitivity = {
           parallel_processors = cores
           n_config = parallel_processors*replicates
           parameters_fname <<- "input/csv/ParametersSensitivityAnalysis.csv"
         },
         FBA_Sensitivity = {
           parallel_processors = cores
           n_config = parallel_processors*replicates
           parameters_fname <<- "input/csv/ParametersSensitivityAnalysisBounds.csv"
         }
  )
  
  if( Exper == "Model_Analysis" ) {
    
    model.analysis(solver_fname = paste0(wd, "/net/", net_fname, ".solver"),
                   i_time = 0, f_time = f_time, s_time = time.step,
                   parallel_processors = parallel_processors,
                   n_config = n_config,
                   FVA = T,
                   fba_fname = paste0(wd, "/input/compiled_models/", fba_fname),
                   atol = atol, rtol = rtol, 
                   debug = debug,
                   parameters_fname = paste0(wd, "/", parameters_fname),
                   functions_fname = paste0(supp_function.dir, "/Functions.R"),
                   event_times = if (Condition != "NoDrug") event_times,
                   event_function = if (Condition != "NoDrug") "treat_generation")
    
    if ( FVA ) {
      
      model.sensitivity(folder_trace = paste0(wd, "/", net_fname, "_analysis/"),
                        FVA = T,
                        fva_gamma = 0.9,
                        debug = debugFVA,
                        parallel_processors = (parallel::detectCores() - 1),
                        fba_fname = paste0(wd, "/input/compiled_models/", fba_fname),
                        flux_fname = react,
                        functions_fname = paste0(supp_function.dir, "Functions.R"))
      
      system(paste("rm -r ./Results/FVA", tag, Condition, Exper, sep = "_"))
      system(paste("mv Model_sensitivity* ./Results/FVA", tag, Condition, Exper, sep = "_"))
      
      if (FVAcomparison) {
        
        source(paste0(supp_function.dir, "glpkVsSybil_4server.R"))
        
        # p = ((deltaVFAmin|deltaVFAmax)/(deltaFVA|deltaFBA))/(pFBA|pFVA)
        # 
        # ggsave(p, file = paste("FVAComp", tag, Condition, ".pdf", sep = ""), 
        #        width = 10, height = 10)
        # 
        # system(paste("cp ./FVAComp", tag, Condition, ".pdf ", "./Results/FVA", "_",
        #              tag, "_", Condition, "_", Exper, sep = ""))
        # 
        # system(paste("mv ./FVAComp", tag, Condition, ".pdf ", "./Results/FVA", "_",
        #              tag, "_", Condition, "_", Exper, sep = ""))
        
      }
      
    }
    
    system(paste("rm -r ./results/CDiff", tag, Condition, Exper, sep = "_"))
    system(paste0("mv ", net_fname, "_analysis* ./results/CDiff", "_", tag, "_", Condition, "_", Exper))
    
    source(paste0(wd, "/code/scripts/trace_fluxes_bounds.R"))
    source(paste0(wd, "/code/scripts/plotting_panorama.R"))
    
    writedata.csv(wd, tag, Condition, Exper, net_fname)
    PlottingPanorama(Exper, Condition, tag, pdf.name = "PanoramaPlot")
    
  }
  
  if( Exper == "Model_Sensitivity" ) {
    
    model.sensitivity(solver_fname = paste0(wd, "/net/", net_fname, ".solver"),
                      i_time = 0, 
                      f_time = f_time, 
                      s_time = time.step,
                      parallel_processors = parallel_processors,
                      n_config = n_config,
                      fba_fname = paste0(wd, "/input/compiled_models/", fba_fname),
                      atol = atol, 
                      rtol = rtol, 
                      debug = debug,
                      target_value = "CD",
                      distance_measure = distance_measure, 
                      reference_data = reference_data, 
                      parameters_fname = paste0(wd, "/", parameters_fname),
                      functions_fname = paste0(supp_function.dir, "/Functions.R"),
                      event_times = if (Condition != "NoDrug") event_times,
                      event_function = if (Condition != "NoDrug") "treat_generation")
    
    system(paste(
      "rm -r ./results/CDiff", tag, Condition, Exper, sep = "_"))
    system(paste0(
      "mv ", net_fname, "_", "analysis* ./results/CDiff", "_", tag, "_", Condition, "_", Exper))
    
    system(paste(
      "rm -r ./results/CDiff", tag, Condition, Exper, "sensitivity", sep = "_"))
    system(paste0(
      "mv ", net_fname, "_", "sensitivity* ./results/CDiff", "_", tag, "_", Condition, "_", Exper))
    
    source(paste0(wd, "/code/scripts/plotting_sensitivity.R"))
    
    PlottingSensitivityAnalysis(numb = r,
                                alph = 0.25,
                                multi_param = multi_param,
                                tag = tag,
                                Condition = Condition,
                                Exper = Exper,
                                net_fname = net_fname,
                                param_target = param_target,
                                colo1 <- c("black", "magenta", "gold"),
                                colo2 <- c("black", "#266867", "yellow"),
                                colo3 <- c("#ffd166", "#ee6c4d", "#293241"),
                                Na = 6.022e20,
                                c = 6.022e08,
                                fc = 1e-06)
    
    setwd(paste0(wd))
    
    ViolinPlotSensitivityAnalysis(Tempi = c(0, 16, 32, 48, 64, 72),
                                  multi_param = multi_param,
                                  tag = tag, 
                                  Condition = Condition, 
                                  Exper = Exper,
                                  net_fname = net_fname,
                                  per = 0.05,
                                  Na = 6.022e20,
                                  c = 6.022e08,
                                  fc = 1e-06)
    
  }
}
