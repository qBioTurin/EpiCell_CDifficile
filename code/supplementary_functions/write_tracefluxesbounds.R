
writedata.csv = function(wd, tag, Condition, Exper, net_fname) {
  
  # Set the directory where the results are stored
  results.dir = paste0(wd, "/results/")
  
  # Read in the trace file and select only the relevant columns
  subtrace = read.table(paste0(results.dir, "CDiff", "_", tag, "_", Condition, "_", 
                               Exper, "/", net_fname, "-analysis-1.trace"), header = T) %>%
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  # Read in the flux file and select only the relevant columns
  subflux = read.table(paste0(results.dir, "CDiff", "_", tag, "_", Condition, "_", 
                              Exper, "/", net_fname, "-analysis-1-0.flux"), header = T) %>% 
    dplyr::select(Time, !c(Obj_0, "X12DGR180ti_Lb":"sink_pheme_c_Ub")) %>% 
    tidyr::gather(key = "Reaction", value = "Flux", -Time)
  
  # Select only the relevant columns from the flux file
  subbounds = read.table(paste0(results.dir, "CDiff", "_", tag, "_", Condition, "_", 
                                Exper, "/", net_fname, "-analysis-1-0.flux"), header = T) %>% 
    dplyr::select(Time, "X12DGR180ti_Lb":"sink_pheme_c_Ub") %>% 
    tidyr::gather(key = "Reaction", value = "Flux", -Time)
  
  # Add the scenario and condition columns to the subtrace data frame
  subtrace$Scenario = tag
  subtrace$Condition = Condition
  
  # Add the scenario and condition columns to the subflux data frame
  subflux$Scenario = tag
  subflux$Condition = Condition
  
  # Add the scenario and condition columns to the subbounds data frame
  subbounds$Scenario = tag
  subbounds$Condition = Condition
  
  # Write the subtrace data frame to a CSV file
  write.csv(subtrace, file = paste0(results.dir, "CDiff_", tag, "_", Condition, "_", 
                                    Exper, "/SubTrace", tag, Condition, ".csv"))
  
  # Write the subflux data frame to a CSV file
  write.csv(subflux, file = paste0(results.dir, "CDiff_", tag, "_", Condition, "_", 
                                   Exper, "/SubFlux", tag, Condition, ".csv"))
  
  # Write the subbounds data frame to a CSV file
  write.csv(subbounds, file = paste0(results.dir, "CDiff_", tag, "_", Condition, "_", 
                                     Exper, "/SubBounds", tag, Condition, ".csv"))
  
}