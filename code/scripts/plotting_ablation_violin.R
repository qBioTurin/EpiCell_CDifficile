
plotting_ablation_violin = function(Exper,
                                    Condition,
                                    tag,
                                    Tempi,
                                    wd) {
 
  setwd(wd)
  
  scientific_10 <- function(x) {parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))}
  
  for (j in tag) {
    
    setwd(paste0(wd, "/results/CDiff", "_", j, "_", Condition, "_", Exper))
    load(file = paste0("subtrace_" , j, "_", Condition, ".RData"))
    subtrace = cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
    assign(paste0("subtrace", j, Condition), subtrace)
    
  }
  
  subtrace = rbind(subtraceUnifiedTherapy, subtraceAblatedTherapy)
  subtrace = filter(subtrace, Time %in% Tempi)
  
  Cdiff = filter(subtrace, subtrace$Places == "CD")
  Cdiff$Time = as.character(Cdiff$Time)
  Biomass = filter(subtrace, subtrace$Places == "BiomassCD")
  Biomass$Time = as.character(Biomass$Time)
  IECs = filter(subtrace, subtrace$Places == "IECs")
  IECs$Time = as.character(IECs$Time)
  
  aa_vec = c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  
  AA = filter(subtrace, subtrace$Places %in% aa_vec)
  
  AA$Time = as.character(AA$Time)
  
  Na = 6.022e20
  c = 6.022e08
  
  for (j in tag) {
    
    assign(paste("C", j, sep = "."),
           filter(Cdiff, Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking)) +
             labs(x = "Time (h)", 
                  y = "CD Marking (cell)",
                  title = "Distribution of configurations | CD Place",
                  subtitle = paste("Scenario: " , j, sep = "")) +
             theme(text = element_text(size = 8, color = "black"),
                   panel.background = element_blank(),
                   panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
                   legend.position = "none",
                   plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot", 
                   axis.text = element_text(size = 10, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) +
             geom_jitter(aes(col = Time), alpha = 0.35) +
             geom_violin(scale = "width", alpha = 0.3, width = 1) +
             geom_boxplot(alpha = 0) +
             scale_y_continuous(label = scientific_10) +
             ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
             scale_x_discrete(name = "", label = rep(paste(paste(Tempi, "\n", sep = ""),
                                                           "(n = ", length(unique(Cdiff$config)), ")",
                                                           sep = ""), length(Tempi))))
    
    assign(paste("B", j, sep = "."),
           filter(Biomass, Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking)) +
             labs(x = "Time (h)", 
                  y = "Biomass Marking (pg)",
                  title = "Distribution of configurations | Biomass Place",
                  subtitle = paste("Scenario: " , j, sep = "")) +
             theme(text = element_text(size = 8, color = "black"),
                   panel.background = element_blank(),
                   panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
                   legend.position = "none",
                   plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot", 
                   axis.text = element_text(size = 10, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) +
             geom_jitter(aes(col = Time), alpha = 0.35) +
             geom_violin(scale = "width", alpha = 0.3, width = 1) +
             geom_boxplot(alpha = 0) +
             scale_y_continuous(label = scientific_10) +
             ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
             scale_x_discrete(name = "", label = rep(paste(paste(Tempi, "\n", sep = ""),
                                                           "(n = ", length(unique(Biomass$config)), ")",
                                                           sep = ""), length(Tempi))))
    
    assign(paste("I", j, sep = "."),
           filter(IECs, Scenario == j) %>%
             ggplot(aes(x = Time, y = Marking)) +
             labs(x = "Time (h)", 
                  y = "IECs Marking (cell)",
                  title = "Distribution of configurations | IECs Place",
                  subtitle = paste("Scenario: " , j, sep = "")) +
             theme(text = element_text(size = 8, color = "black"),
                   panel.background = element_blank(),
                   panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
                   legend.position = "none",
                   plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot", 
                   axis.text = element_text(size = 10, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) +
             geom_jitter(aes(col = Time), alpha = 0.35) +
             geom_violin(scale = "width", alpha = 0.3, width = 1) +
             geom_boxplot(alpha = 0) +
             scale_y_continuous(label = scientific_10) +
             ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
             scale_x_discrete(name = "", label = rep(paste(paste(Tempi, "\n", sep = ""),
                                                           "(n = ", length(unique(IECs$config)), ")",
                                                           sep = ""), length(Tempi))))
    assign(paste("AA", j, sep = "."),
           filter(AA, Scenario == j) %>%
             ggplot(aes(x = Time, y = (Marking/((Na/c)*1000)))) +
             labs(x = "Time (h)", 
                  y = latex2exp::TeX("\\textbf{aminoacids ($\\mu$mol)}"),
                  title = latex2exp::TeX("\\textbf{aminoacids ($\\mu$mol)}"),
                  subtitle = paste0("Scenario: " , j)) +
             theme(text = element_text(size = 8, color = "black"),
                   panel.background = element_blank(),
                   panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
                   legend.position = "none",
                   plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
                   plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
                   plot.title.position = "plot", 
                   axis.text = element_text(size = 5, color = "black"),
                   axis.title = element_text(size = 12, face = "bold")) +
             geom_jitter(aes(col = Time), alpha = 0.35, size = 0.1) +
             geom_violin(scale = "width", alpha = 0.3, width = 1) +
             geom_boxplot(alpha = 0) +
             facet_wrap(~Places, nrow = 2) +
             ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
             scale_x_discrete(
               name = "", 
               label = rep(paste0(
                 paste0(Tempi, "\n"), 
                 "(n = ", length(unique(AA$config)), ")"), length(Tempi))))
    
  }
 
  ggsave((C.Ablated/B.Ablated/I.Ablated/AA.Ablated) | (C.Unified/B.Unified/I.Unified/AA.Unified), 
         file = paste0(wd, "/ViolinAblationVSUnified.pdf"), width = 14, height = 16)
  
  system(paste0("cp ", wd, "/ViolinAblationVSUnified.pdf ", wd, "/results/images4paper"))
  system(paste0("mv ", wd, "/ViolinAblationVSUnified.pdf ", wd, "/results/images4paper"))
  
}
