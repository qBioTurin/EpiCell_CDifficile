
ModelAnalysisPlot=function(TracesPath, FluxPath, FluxVec) {
  
  library(dplyr)
  
  trace = read.table(TracesPath, header = T)
  
  subtrace = trace %>% select(Time, IECs, pheme_e, BiomassCD, CD,
                   pro_L_e, leu_L_e, ile_L_e, val_L_e, trp_L_e, cys_L_e, 
                   Drug, Damage) %>% 
    tidyr::gather(key = "Places", value = "Marking", -Time) %>% 
    ggplot() + ggplot2::geom_line(ggplot2::aes(x=Time, y=Marking, group = Places)) +
    ggplot2::facet_wrap(~Places, scales = "free_y") + ggplot2::scale_fill_brewer("Accent")
  
  subtrace
  
  flux = read.table(FluxPath, header = T, sep = ";")
  
  if(is.null(FluxVec)){
    
    # bla bla bla
    
    } else {
      
      subflux = flux %>% select(FluxVec, Time) %>%
        tidyr::gather(key = "Reaction", value = "Flux", -Time)
      
      ggplot(subflux, aes(Time, Reaction, height=Flux)) +
        ggridges::geom_density_ridges(stat = "binline", binwidth=1,
                                      draw_baseline = F)
    }
  
  subtrace
  
  }
