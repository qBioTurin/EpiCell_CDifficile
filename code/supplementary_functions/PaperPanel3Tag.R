
colo1 <- c("black", "magenta", "gold")
colo2 <- c("black", "#266867", "yellow")
colo3 <- c("#ffd166", "#ee6c4d", "#293241")
coloTag = data.frame(ablated = "#EC6D67", unified = "#5351A2",superablated = "#293241")

net_name = c("ablated","unified","superablated")

fluxeInterest = c("EX_biomass_e", "sink_pheme_c", "EX_cys_L_e" , "EX_trp_L_e", "EX_val_L_e",
                  "EX_ile_L_e", "EX_leu_L_e", "EX_pro_L_e")
param_target = c("IECsDeath", "Death4Treat", "Detox")

alph = 0.25

Tempi = c(0, 16, 32, 48, 64, 72)

subtrace = do.call("rbind",lapply(net_name, function(j) {
  folder_traces = paste0(j,"_analysis/")
  load(paste(folder_traces,j, "-analysis.RData", sep = ""))
  id_param = c()
  for (p in 1:length(param_target)) {
    for(i in 1:length(config)) {
      if(param_target[p] == config[[i]][[1]][[1]]) {
        id_param[p] = i
      }
    }
  }
  
  
  failed.tr <<- c(); f.tr <<- 1
  
  list.Tr = list.files(pattern = "\\.trace$",path = folder_traces)
  
  subtrace = lapply(list.Tr, function(i,jj,folder){
    
    TracesPath = i
    
    n = as.numeric(gsub(pattern = paste("(", jj, "-analysis-)|(.trace)", sep = ""), 
                        replacement = "", x = i))
    
    trace = read.table(paste0(folder,"/",TracesPath), header = T)
    
    if(!all(sapply(trace, is.numeric))) {
      print(i)
      failed.tr[f.tr] <<- n
      f.tr <<- f.tr + 1
      return()
    }
    
    subtrace = trace %>% 
      dplyr::select(Time, IECs, Drug, Damage, pheme_e, 
                    pheme_c, BiomassCD, CD, pro_L_e, leu_L_e, 
                    ile_L_e, val_L_e, trp_L_e, cys_L_e) %>% 
      tidyr::gather(key = "Places", value = "Marking", -Time)
    
    subtrace = cbind(subtrace, config = rep(as.character(i), length(subtrace$Time)))
    
    for (p in 1:length(param_target)) {
      subtrace = cbind(subtrace, p = rep(config[[id_param[p]]][[n]][[3]], length(subtrace$Time)))
    }
    
    return(subtrace)
    
  }, jj = j,folder = folder_traces)
  
  subtrace = do.call(rbind, subtrace)
  names(subtrace) = c("Time", "Places", "Marking", "config", param_target)
  
  subtrace$Scenario = j
  return(subtrace)
} ) )


subflux = do.call("rbind",lapply(net_name, function(j){
  if(j == "superablated"){
    return(NULL)
  }else{
    folder_traces = paste0(j,"_analysis/")
    load(paste(folder_traces,j, "-analysis.RData", sep = ""))
    id_param = c()
    for (p in 1:length(param_target)) {
      for(i in 1:length(config)) {
        if(param_target[p] == config[[i]][[1]][[1]]) {
          id_param[p] = i
        }
      }
    }
    
    list.Fl = list.files(pattern = "\\.flux$",path = folder_traces)
    
    subtrace = lapply((list.Fl), function(i){
      
      TracesPath = paste0(folder_traces,i)
      
      n = as.numeric(gsub(pattern = paste("(", j, "-analysis-)|(-0.flux)", sep = ""),
                          replacement = "", x = i))
      
      flux = read.table(TracesPath, header = F)
      
      ReactionsNames = as.matrix(flux[1, -c(1, 2)])
      ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
      ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
      
      flux = read.table(TracesPath, header = T)
      
      colnames(flux) = c("Time", "Obj_0", ReactionsNames)
      
      subflux = flux[, c("Time", ReactionsNames)] %>%
        tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
        filter(Reaction %in% fluxeInterest )

      subflux$config = as.character(n)
      
      for (p in 1:length(param_target)) {
        subflux = cbind(subflux,
                        p = rep(config[[id_param[p]]][[n]][[3]], length(subflux$Time)))
      }
      
      return(subflux)
    })
    subtrace = do.call(rbind, subtrace)
    names(subtrace) = c("Time", "Reactions", "Fluxes", "config", param_target)
    subtrace$Scenario = j
    return(subtrace)
  }
}))

AA = filter(subtrace, Time %in% Tempi)

AA = filter(AA, AA$Places %in% 
              c("trp_L_e", "pro_L_e", "val_L_e", 
                "ile_L_e", "cys_L_e", "leu_L_e"))

AA$Time = as.character(AA$Time)

#################################################################################

for (j in net_name) {
  if(j == "superablated")
    assign(paste("EX_biomass_e", param_target[3], j, sep = "."),plot_spacer())
  else
  assign(paste("EX_biomass_e", param_target[3], j, sep = "."),
         filter(subflux[which(subflux$Reactions == "EX_biomass_e"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Fluxes, color = Detox)) +
           labs(x = "Time (h)",
                y = "EX_biomass(e) \n (mmol/gDW*h)",
                subtitle =  "EX_biomass(e) Flux",
                title = paste(j, " Experiments", sep = ""),
                colour = latex2exp::TeX(paste("Param: ", param_target[3], sep = ""))) +
           theme(plot.title = element_text(size = 14, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1b2838"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold"),
                 legend.position = "none") +
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo3,
                                  oob = scales::squish))
  
  assign(paste("BiomassCD", param_target[3], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "BiomassCD"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking, color = Detox)) + 
           labs(x = "Time (h)", 
                y = "Biomass CD",
                subtitle = "Biomass",
                title = paste(j, " Experiments", sep = ""),
                colour = paste("Param: ", param_target[3], sep = "")) + 
           theme(plot.title = element_text(size = 14, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1b2838"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo3,
                                  oob = scales::squish))
  assign(paste("CD", param_target[3], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "CD"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking, color = Detox)) + 
           labs(x = "Time (h)", 
                y = "C. difficile (cell)",
                subtitle = "CD",
                title = paste(j, " Experiments", sep = ""),
                colour = paste("Param: ", param_target[3], sep = "")) + 
           theme(plot.title = element_text(size = 14, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1b2838"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold"),
                 legend.position = "none") + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo3,
                                  oob = scales::squish))
  
  assign(paste("IECs", param_target[3], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "IECs"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking, color = Detox)) + 
           labs(x = "Time (h)", 
                y = "IECs (cell)",
                subtitle = "IECs",
                title = paste(j, " Experiments", sep = ""),
                colour = paste("Param: ", param_target[3], sep = "")) + 
           theme(plot.title = element_text(size = 14, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1b2838"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo3,
                                  oob = scales::squish))
  
  assign(paste("trp_L_e", param_target[1], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "trp_L_e"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking/(6.022e17/6.022e08), color = IECsDeath)) + 
           labs(x = "Time (h)", 
                y = latex2exp::TeX("\\textbf{Tryptophan ($\\mu$mol/mL)}"),
                subtitle = "trp_L_e",
                title = paste(j, " Experiments", sep = ""),
                colour = paste("Param: ", param_target[1], sep = "")) + 
           theme(plot.title = element_text(size = 14, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1b2838"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold"),
                 legend.position = "none") + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo1,
                                  oob = scales::squish))
  
  assign(paste("ile_L_e", param_target[1], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "ile_L_e"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking/(6.022e17/6.022e08), color = IECsDeath)) + 
           labs(x = "Time (h)", 
                y = latex2exp::TeX("\\textbf{Isoleucine ($\\mu$mol/mL)}"),
                subtitle = "ile_L_e",
                title = paste(j, " Experiments", sep = ""),
                colour = paste("Param: ", param_target[1], sep = "")) + 
           theme(plot.title = element_text(size = 14, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1b2838"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold"),
                 legend.position = "none") + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo1,
                                  oob = scales::squish))
  
  assign(paste("leu_L_e", param_target[1], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "leu_L_e"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking/(6.022e17/6.022e08), color = IECsDeath)) + 
           labs(x = "Time (h)", 
                y = latex2exp::TeX("\\textbf{Leucine ($\\mu$mol/mL)}"),
                subtitle = "leu_L_e",
                title = paste(j, " Experiments", sep = ""),
                colour = paste("Param: ", param_target[1], sep = "")) + 
           theme(plot.title = element_text(size = 14, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1b2838"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo1,
                                  oob = scales::squish))
  
  assign(paste("AA", j, sep = "."),
         filter(AA, Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking/(6.022e17/6.022e08))) +
           labs(x = "Time (h)", 
                y = latex2exp::TeX("\\textbf{L-amino acids ($\\mu$mol/mL)}"),
                subtitle = "L-amino acids",
                title = paste(j, " Experiments", sep = "")) +
           theme(text = element_text(size = 8, color = "black"),
                 panel.background = element_blank(),
                 panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
                 legend.position = "none",
                 plot.title = element_text(size = 12, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot", 
                 axis.text = element_text(size = 6, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) +
           geom_jitter(aes(col = Time), alpha = 0.35, size = 0.1) +
           geom_violin(scale = "width", alpha = 0.3, width = 0.5) +
           geom_boxplot(alpha = 0) +
           facet_wrap(~Places, nrow = 2) +
           ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")))
  
}

#################################################################################

ggsave((EX_biomass_e.Detox.superablated|BiomassCD.Detox.superablated|CD.Detox.superablated|IECs.Detox.superablated)/
         (EX_biomass_e.Detox.ablated|BiomassCD.Detox.ablated|CD.Detox.ablated|IECs.Detox.ablated) /
         (EX_biomass_e.Detox.unified|BiomassCD.Detox.unified|CD.Detox.unified|IECs.Detox.unified),
       file = paste("./Results/images4paper/PaperPanelA.pdf", sep = ""), width = 15, height = 5)


ggsave((trp_L_e.IECsDeath.superablated|ile_L_e.IECsDeath.superablated|leu_L_e.IECsDeath.superablated)/
         (trp_L_e.IECsDeath.ablated|ile_L_e.IECsDeath.ablated|leu_L_e.IECsDeath.ablated) / 
         (trp_L_e.IECsDeath.unified|ile_L_e.IECsDeath.unified|leu_L_e.IECsDeath.unified),
       file = paste("./Results/images4paper/PaperPanelB.pdf", sep = ""), width = 15, height = 6)


ggsave((AA.superablated|AA.ablated|AA.unified), 
       file = paste("./Results/images4paper/PaperPanelC.pdf", sep = ""), 
       width = 12, height = 4)

