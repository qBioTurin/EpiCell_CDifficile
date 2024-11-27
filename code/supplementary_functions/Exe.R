
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
  
  tag <<- Ablation
  Condition <<- therapy
  
  Exper <<- Exper
  r <<- replicates*cores
  
  model_name <<- model_name
  model_cat <<- model_cat
  
  f_time <<- f_time
  time.step <<- time.step
  
  param_target <<- param_target
  
  distance_measure <<- distance_measure
  reference_data <<- reference_data
  
  multi_param <<- if (length(param_target) > 1) TRUE else FALSE
  
  transitions_fname = switch(tag,
                             "SuperAblated" = paste0(wd, "/Net/SuperAblated.cpp"),
                             "Ablated" = paste0(wd, "/Net/Ablated.cpp"),
                             "Unified" = paste0(wd, "/Net/Unified.cpp"))
  
  model.generation(net_fname = paste0(wd, "/Net/", net_fname, ".PNPRO"),
                   transitions_fname = transitions_fname,
                   fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname))
  
  system(paste0("mv ", net_fname, ".* ./Net"))
  
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
    
    source(paste0(wd, supp_function.dir, "tracefluxesbounds.R"))
    source(paste0(wd, supp_function.dir, "PlottingAnalysisResults.R"))
    
    writedata.csv(wd, tag, Condition, Exper, net_fname)
    PlottingPanorama(Exper, Condition, tag, pdf.name = "PanoramaPlot")
    
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
    
    system(paste("rm -r ./results/CDiff", tag, Condition, Exper, sep = "_"))
    system(paste0("mv ", net_fname, "_", "analysis* ./results/CDiff", "_", tag, "_", Condition, "_", Exper))
    
    system(paste("rm -r ./results/CDiff", tag, Condition, Exper, "sensitivity", sep = "_"))
    system(paste0("mv ", net_fname, "_", "sensitivity* ./results/CDiff", "_", tag, "_", Condition, "_", Exper))
    
    net_fname <<- net_fname
    
    source(paste0(wd, supp_function.dir, "PlottingSensitivityAnalysis.R"))
    
    p = (pC/pBiom | pdrug/pI) / ((pcys/pval) | (pleu/ppro) | (pile/ptrp))
    Fl = pflux/pflux_special
    
    p3 = (pC3/pBiom3 | pdrug3/pI3) / ((pcys3/pval3) | (pleu3/ppro3) | (pile3/ptrp3))
    Fl3 = pflux3/pflux_special3
    
    pM = ((pCM + pdrugM) / (pIECsM + pBiomassM)) + pfluxM
    pM.median = ((pCM.median + pdrugM.median) / (pIECsM.median + pBiomassM.median)) + pfluxM.median
    
    ggsave(p, file = paste0(Exper, "_multiparam_places_", tag, "_", Condition, ".pdf"),  width = 25, height = 15)
    ggsave(p3, file = paste0(Exper, "_multiparam_places_3_", tag, "_", Condition, ".pdf"), width = 9, height = 10)
    
    ggsave(Fl, file = paste0(Exper, "_multiparam_fluxes_", tag, "_", Condition, ".pdf"), width = 15, height = 6)
    ggsave(Fl3, file = paste0(Exper, "_multiparam_fluxes_3_", tag, "_", Condition, ".pdf"), width = 15, height = 6)
    
    ggsave(ppheme, file = paste0(Exper, "_heme", tag, "_", Condition, ".pdf"), width = 15, height = 3.5)
    ggsave(pM, file = paste0("RankingPlot", tag, "_", Condition, ".pdf"), width = 10, height = 12)
    ggsave(pM.median, file = paste0("RankingPlot_median", tag, "_", Condition, ".pdf"), width = 10, height = 12)
    
    ggsave(pr, file = paste0("RankingBoxPlot", tag, "_", Condition, ".pdf"), width = 6, height = 6)
    
    setwd(paste0(wd))
    
    source(paste0(wd, supp_function.dir, "ViolinPlotSensitivityAnalysis.R"))
    
    pV1 = (V.trp | V.pro | V.val) / (V.ile | V.cys | V.leu)
    pV2 = ((V.Cdiff | V.Drug)/(V.Iecs | V.Biomass))
    
    pV3 = (V2.trp | V2.pro | V2.val) / (V2.ile | V2.cys | V2.leu)
    pV4 = ((V2.Cdiff | V2.Drug)/(V2.Iecs | V2.Biomass))
    
    ggsave(pV1, file = paste("ViolinPlot_AA_", tag, "_", Condition, ".pdf"), width = 20, height = 14)
    ggsave(pV2, file = paste("ViolinPlot_Places_", tag, "_" ,Condition, ".pdf"), width = 16, height = 12)
    ggsave(pV3, file = paste("ViolinPlot_AA_", tag, "_", Condition, "_Ver2.pdf"), width = 17, height = 12)
    ggsave(pV4, file = paste("ViolinPlot_Places_", tag, "_", Condition, "_Ver2.pdf"), width = 12, height = 10)
    
  }
  
}