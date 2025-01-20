 
# Exper = "Model_Analysis"
# colConNoDrug = "darkred"
# colConTherapy = "darkcyan"
# coltag = "darkgray"
# tag = c("Ablated", "Unified")
# Condition = c("NoDrug", "Therapy")
# wd = wd

PlottingFigS4 = function(Exper,
                         colConNoDrug,
                         colConTherapy,
                         coltag, 
                         tag, 
                         Condition,
                         wd) {
  
  scientific_10 <- function(x) {
    text = gsub("e", " %*% 10^", scales::scientific_format()(x))
    text[1] = "0"
    parse(text=text)
  }
  
  # Create scale_x based on variable
  y_scale = scale_y_continuous(label = scientific_10)
  
  results.dir <- paste0(wd, "/results/")
  
  for (i in Condition) {
    for (j in tag) {
      assign(paste0("subtrace", j, i), readr::read_csv(paste0(results.dir, "CDiff_", j, "_", i, "_", Exper, "/SubTrace", j, i, ".csv"))[-1])
      assign(paste0("subflux", j, i), readr::read_csv(paste0(results.dir, "CDiff_", j, "_", i, "_", Exper, "/SubFlux", j, i, ".csv"))[-1])
    }
  }
  
  subtrace <- do.call(rbind, mget(ls(pattern = "subtrace")))
  subflux <- do.call(rbind, mget(ls(pattern = "subflux")))
  
  reactionsvec <- c("EX_biomass_e", "sink_pheme_c",
                    "EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e", 
                    "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")
  
  pl = ggplot(subtrace, aes(Time, Marking, colour = Condition, linetype = Scenario)) + 
    theme_bw() + geom_line(linewidth = 0.6) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("Unified" = "solid", "Ablated" = "dashed"), name = "Scenario") +
    labs(x = "Time (h)", y = "Marking") +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    facet_wrap(~Places, ncol = 3, scales = "free") + y_scale
  
  # Plot for fluxes
  pl_fluxes = ggplot(subflux %>% filter(Reaction %in% reactionsvec), 
                     aes(Time, Flux, colour = Condition, linetype = Scenario)) + 
    theme_bw() + 
    geom_line(linewidth = 0.6) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("Unified" = "solid", "Ablated" = "dashed"), name = "Scenario") +
    labs(x = "Time (h)", y = "Flux (mmol/gDW*h)",
         title = "Metabolic Flux Dynamics") +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = c(8, 16, 24, 32, 40, 48), 
               linetype = "dotted", color = "darkblue", lwd = 0.3) +
    facet_wrap(~Reaction, ncol = 3, scales = "free") + 
    y_scale
  
  return(list(pl, pl_fluxes))
}