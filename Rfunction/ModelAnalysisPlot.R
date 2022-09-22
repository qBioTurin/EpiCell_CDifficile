
ModelAnalysisPlot=function(TracesPath, FluxPath, FluxVec) {
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  
  library(dplyr)
  library(ggplot2)
  
  trace = read.table(TracesPath, header = T)
  
  subtrace = trace %>% dplyr::select(Time, IECs, pheme_e, BiomassCD, CD,
                   pro_L_e, leu_L_e, ile_L_e, val_L_e, trp_L_e, cys_L_e, 
                   Drug, Damage) %>% tidyr::gather(key = "Places", value = "Marking", -Time)
  
  print(ggplot(subtrace, aes(x = Time, y = Marking)) + 
          geom_line(aes(color = Places)) + 
          facet_grid(Places ~ ., scales = "free_y") + 
          theme(legend.position = "none"), 
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  
  flux = read.csv(FluxPath, header=FALSE, sep=";")
  colnames(flux) = flux[1, ]
  flux = flux[-1, ]
  
  if(is.null(FluxVec)){
    
    FluxVec = colnames(flux)
    flux %>% dplyr::select(FluxVec, Time) %>% tidyr::gather(key = "Reaction", value = "Flux", -Time)
    
    subflux = flux %>% select(FluxVec, Time) %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
      ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x=Time, y=Flux, group = Reaction)) +
      ggplot2::facet_wrap(~Reaction, scales = "free_y") + ggplot2::scale_fill_brewer("Accent")
    
    print(subflux, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    
    } else {
      
      subflux = flux %>% dplyr::select(FluxVec, Time) %>%
        tidyr::gather(key = "Reaction", value = "Flux", -Time)
      
      subflux$Time = as.numeric(as.character(subflux$Time))
      subflux$Flux = as.numeric(as.character(subflux$Flux))
      
      library(ggplot2)

      
      print(ggplot(subflux, aes(x = Time, y = Flux)) + geom_line(aes(color = Reaction)) +
              facet_grid(Reaction ~ ., scales = "free_y") + 
              theme(legend.position = "none"), 
            vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
      
    }
  
  }
