
plotting_ablation_sensitivity = function(Exper,
                                         Condition,
                                         tag,
                                         param_target,
                                         colo1,
                                         colo2,
                                         colo3) {
  
  for (j in tag) {
    
    setwd(paste0(wd, "/results/CDiff", "_", j, "_", Condition, "_", Exper))
    load(file = paste0("subflux_" , j, "_", Condition, ".RData"))
    subflux = cbind(subflux, Scenario = rep(j, length(subflux$Time)))
    assign(paste0("subflux", j, Condition), subflux)
    
  }
  
  subflux = rbind(subfluxAblatedTherapy, subfluxUnifiedTherapy)
  
  for (j in tag) {
    
    assign(paste("pflux", param_target[1], j, sep = "."), 
           filter(subflux[which(subflux$Places %in%
                                  c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                                    "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")), ], Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking, color = IECsDeath)) +
             labs(x = "Time (h)", 
                  y = "Output fluxes (mmol/gDW*h)",
                  title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[1], sep = ""),
                  subtitle = paste("Scenario: " , j, sep = "")) + 
             theme(plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot",
                   axis.text = element_text(size = 12, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) + 
             geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
             scale_colour_gradientn(name = "",
                                    colors = colo1,
                                    oob = scales::squish))
    
    assign(paste("pflux", param_target[2], j, sep = "."), 
           filter(subflux[which(subflux$Places %in%
                                  c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                                    "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")), ], Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking, color = Death4Treat)) +
             labs(x = "Time (h)", 
                  y = "Output fluxes (mmol/gDW*h)",
                  title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[2], sep = ""),
                  subtitle = paste("Scenario: " , j, sep = "")) + 
             theme(plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot",
                   axis.text = element_text(size = 12, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) + 
             geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
             scale_colour_gradientn(name = "",
                                    colors = colo2,
                                    oob = scales::squish))
    
    assign(paste("pflux", param_target[3], j, sep = "."), 
           filter(subflux[which(subflux$Places %in%
                                  c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                                    "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")), ], Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking, color = Detox)) +
             labs(x = "Time (h)", 
                  y = "Output fluxes (mmol/gDW*h)",
                  title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[3], sep = ""),
                  subtitle = paste("Scenario: " , j, sep = "")) + 
             theme(plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot",
                   axis.text = element_text(size = 12, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) + 
             geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
             scale_colour_gradientn(name = "",
                                    colors = colo3,
                                    oob = scales::squish))
    
    assign(paste("pfluxspecial", param_target[1], j, sep = "."),
           filter(subflux[which(subflux$Places %in% 
                                  c("EX_biomass_e", "sink_pheme_c")), ], Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking, color = IECsDeath)) +
             labs(x = "Time (h)", 
                  y = "Output fluxes (mmol/gDW*h)",
                  title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[1], sep = ""),
                  subtitle = paste("Scenario: " , j, sep = "")) + 
             theme(plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot",
                   axis.text = element_text(size = 12, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) + 
             geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
             scale_colour_gradientn(name = "",
                                    colors = colo1,
                                    oob = scales::squish))
    
    assign(paste("pfluxspecial", param_target[2], j, sep = "."),
           filter(subflux[which(subflux$Places %in% 
                                  c("EX_biomass_e", "sink_pheme_c")), ], Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking, color = Death4Treat)) +
             labs(x = "Time (h)", 
                  y = "Output fluxes (mmol/gDW*h)",
                  title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[2], sep = ""),
                  subtitle = paste("Scenario: " , j, sep = "")) + 
             theme(plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot",
                   axis.text = element_text(size = 12, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) +
             geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
             scale_colour_gradientn(name = "",
                                    colors = colo2,
                                    oob = scales::squish))
    
    assign(paste("pfluxspecial", param_target[3], j, sep = "."),
           filter(subflux[which(subflux$Places %in% 
                                  c("EX_biomass_e", "sink_pheme_c")), ], Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking, color = Detox)) +
             labs(x = "Time (h)", 
                  y = "Output fluxes (mmol/gDW*h)",
                  title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[3], sep = ""),
                  subtitle = paste("Scenario: " , j, sep = "")) + 
             theme(plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot",
                   axis.text = element_text(size = 12, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) +
             geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
             scale_colour_gradientn(name = "",
                                    colors = colo3,
                                    oob = scales::squish))
    
  }
  
  ggsave(((pflux.Death4Treat.Ablated | pflux.IECsDeath.Ablated | pflux.Detox.Ablated) /
            (pflux.Death4Treat.Unified | pflux.IECsDeath.Unified | pflux.Detox.Unified)) / 
           ((pfluxspecial.Death4Treat.Ablated | pfluxspecial.IECsDeath.Ablated | pfluxspecial.Detox.Ablated) /
              (pfluxspecial.Death4Treat.Unified | pfluxspecial.IECsDeath.Unified | pfluxspecial.Detox.Unified)) , 
         file = paste0(wd, "/AblationAnalysis_Fluxes.pdf"), width = 28, height = 20)
  
  system(paste0("cp ", wd, "/AblationAnalysis_Fluxes.pdf ", wd, "/results/images4paper"))
  system(paste0("mv ", wd, "/AblationAnalysis_Fluxes.pdf ", wd, "/results/images4paper"))
  
}
