
PlottingFigS4 = function(Exper,
                         colConNoDrug,
                         colConTherapy,
                         coltag, 
                         tag, 
                         Condition,
                         wd) {
  
  scientific_10 <- function(x) {parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x)))}
  results.dir <- paste0(wd, "/results/")
  
  for (i in Condition) {
    for (j in tag) {
      assign(paste0("subtrace", j, i), readr::read_csv(paste0(results.dir, "CDiff_", j, "_", i, "_", Exper, "/SubTrace", j, i, ".csv"))[-1])
      assign(paste0("subflux", j, i), readr::read_csv(paste0(results.dir, "CDiff_", j, "_", i, "_", Exper, "/SubFlux", j, i, ".csv"))[-1])
      assign(paste0("subbounds", j, i), readr::read_csv(paste0(results.dir, "CDiff_", j, "_", i, "_", Exper, "/SubBounds", j, i, ".csv"))[-1])
    }
  }
  
  subtrace <- do.call(rbind, mget(ls(pattern = "subtrace")))
  subflux <- do.call(rbind, mget(ls(pattern = "subflux")))
  subbounds <- do.call(rbind, mget(ls(pattern = "subbounds")))
  
  aa.placesvec <- c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  aa.reactionsvec <- c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e", "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")
  
  pl = ggplot(subtrace, 
              aes(Time, Marking, colour = Condition, linetype = Scenario)) + 
    theme_bw() + geom_line(linewidth = 0.6) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "Marking",
         title = paste("Model dynamics ", "(", "Unified VS Ablated", " | ", "NoDrug VS Therapy", ")", sep = "")
    ) +
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
    facet_wrap(~Places, ncol = 3,scales = "free")
  
  return(pl)
}