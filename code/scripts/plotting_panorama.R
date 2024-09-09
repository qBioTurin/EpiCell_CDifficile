
PlottingPanorama <- function(Exper, Condition, tag, pdf.name) {
  
  # Set the directory where the results are stored
  results.dir = paste0(wd, "/results/")
  # Set scintific annotation function
  scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))}
  # Define color based on condition
  colConFirst = ifelse(Condition == "NoDrug", "darkred", "darkcyan")
  
  # Read in data
  subtrace = readr::read_csv(
    paste0(results.dir, "CDiff_", tag, "_", Condition, "_", 
           Exper, "/SubTrace", tag, Condition, ".csv"))[-1]
  subflux = readr::read_csv(
    paste0(results.dir, "CDiff_", tag, "_", Condition, "_", 
           Exper, "/SubFlux", tag, Condition, ".csv"))[-1]
  subbounds = readr::read_csv(
    paste0(results.dir, "CDiff_", tag, "_", Condition, "_", 
           Exper, "/SubBounds", tag, Condition, ".csv"))[-1]
  
  # Define vectors for amino acids
  aa.placesvec = c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  aa.reactionsvec = paste0("EX_", aa.placesvec)
  aa.boundsvec_Lb = paste0(aa.reactionsvec, "_Lb")
  aa.boundsvec_Ub = paste0(aa.reactionsvec, "_Ub")
  
  # Filter data by place
  Cdiff = dplyr::filter(subtrace, subtrace$Places == "CD")
  Biomass = dplyr::filter(subtrace, subtrace$Places == "BiomassCD")
  Drug = dplyr::filter(subtrace, subtrace$Places == "Drug")
  IECs = dplyr::filter(subtrace, subtrace$Places == "IECs")
  pheme_e = dplyr::filter(subtrace, subtrace$Places == "pheme_e")
  pheme_c = dplyr::filter(subtrace, subtrace$Places == "pheme_c")
  Toxin = dplyr::filter(subtrace, subtrace$Places == "Toxin")
  AA = dplyr::filter(subtrace, subtrace$Places %in% aa.placesvec)
  
  # Filter data by flux
  FluxBiom = dplyr::filter(subflux, subflux$Reaction == "EX_biomass_e")
  FluxPheme = dplyr::filter(subflux, subflux$Reaction == "sink_pheme_c")
  FluxAA = dplyr::filter(subflux, subflux$Reaction %in% aa.reactionsvec)
  
  # Filter data by bounds
  BoundsBiom = dplyr::filter(subbounds, subbounds$Reaction %in% c("EX_biomass_e_Ub", "EX_biomass_e_Lb"))
  BoundsSinkHeme = dplyr::filter(subbounds, subbounds$Reaction %in% c("sink_pheme_c_Ub", "sink_pheme_c_Lb"))
  BoundsAALb = dplyr::filter(subbounds, subbounds$Reaction %in% aa.boundsvec_Lb)
  BoundsAAUb = dplyr::filter(subbounds, subbounds$Reaction %in% aa.boundsvec_Ub)
  
  Na = 6.022e20
  c = 6.022e08
  pack = 1*(Na*(1/c))
  fc = 1e-06
    
  # Places
  C = ggplot(Cdiff, aes(Time, Marking)) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "C. difficile (cells)",
         title = paste("Model dynamics ", "(", tag, " | ", Condition, ")", sep = ""),
         subtitle = paste("Place: ", Cdiff[1, 2], " (cells)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70)
  
  B = ggplot(Biomass, aes(Time, Marking)) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "C. difficile average biomass (pg)",
         title = paste("Model dynamics ", "(", tag, " | ", Condition, ")", sep = ""),
         subtitle = paste("Place: ", Biomass[1, 2], " (pg)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 7, 1.5, label = "biomass max", vjust = 0.5, color = "black", size = 3) +
    geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 7, 0.6, label = "biomass mean", vjust = 0.5, color = "black", size = 3) +
    geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 7, 0.15, label = "biomass min", vjust = 0.5, color = "black", size = 3)
  
  D = ggplot(Drug, aes(Time, Marking*fc)) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "Drug concentration (µmol)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Place: ", Drug[1, 2], " (µmol)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_hline(yintercept = (5800*1e-06), linetype = "dashed", color = "darkgrey", lwd = 0.5) +
    annotate("text", 3, (6500*1e-06), label = "MIC value", vjust = 0.5, color = "black", size = 3)
  
  I = ggplot(IECs, aes(Time, Marking)) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "IECs (cell)",
         title = paste("Model dynamics ", "(", tag, " | ", Condition, ")", sep = ""),
         subtitle = paste("Place: ", IECs[1, 2], " (cell)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    scale_y_continuous(label = scientific_10)
  
  Hc = ggplot(pheme_c, aes(Time, (Marking*fc))) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "Intracell. heme concentration (µmol)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Place: ", pheme_c[1, 2], " (µmol)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70)
  
  He = ggplot(pheme_e, aes(Time, (Marking*fc))) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "Extr. heme concentration (µmol)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Place: ", pheme_e[1, 2], " (µmol)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_hline(yintercept = (1e+06*1e-06), linetype = "dashed", color = "darkred", lwd = 0.5) +
    annotate("text", 12, (0.95e+06*1e-06), label = "Hemorrhagic baseline", vjust = 0.5, color = "black", size = 3)
  
  A = ggplot(AA, aes(Time, (Marking/pack)*1000, color = Places)) + theme_bw() +
    geom_line(size = 0.8) +
    labs(x = "Time (h)", y = "Aminoacids concentration (µmol)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Place: aminoacids", " (C molecule)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    facet_wrap(~Places, nrow = 2)
  
  # Fluxes
  fB = ggplot(FluxBiom, aes(Time, Flux)) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "EX_biomass_e FBA flux (mmol/gDW*h)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Reaction: ", FluxBiom[1, 2], " (pg)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70)
  
  fH = ggplot(FluxPheme, aes(Time, Flux)) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "sink_pheme_c FBA flux (mmol/gDW*h)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Reaction: ", FluxPheme[1, 2], " (pmol)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70)
  
  fA = ggplot(FluxAA, aes(Time, Flux, color = Reaction)) + theme_bw() +
    geom_line(size = 0.8) +
    labs(x = "Time (h)", y = "EX_Aminoacid_e FBA flux (mmol/gDW*h)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Reaction: EX_Aminoacid_e", " (C molecule)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    facet_wrap(~Reaction, nrow = 2)
  
  # bounds
  bB = ggplot(BoundsBiom, aes(Time, Flux, linetype = Reaction)) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "EX_biomass_e Lb & Ub (mmol/gDW*h)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Double bounded flux (Lb =< v =< Ub)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70)
  
  LbA = ggplot(BoundsAALb, aes(Time, Flux, color = Reaction)) + theme_bw() +
    geom_line(size = 0.8) +
    labs(x = "Time (h)", y = "EX_Aminoacid_e Lb (mmol/gDW*h)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Double bounded flux (Lb =< v =< Ub)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    facet_wrap(~Reaction, nrow = 2)
  
  UbA = ggplot(BoundsAAUb, aes(Time, Flux, color = Reaction)) + theme_bw() +
    geom_line(size = 0.8) +
    labs(x = "Time (h)", y = "EX_Aminoacid_e Ub (mmol/gDW*h)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Double bounded flux (Lb =< v =< Ub)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    facet_wrap(~Reaction, nrow = 2)
  
  # bounds
  bH = ggplot(BoundsSinkHeme, aes(Time, Flux, linetype = Reaction)) + theme_bw() +
    geom_line(size = 0.8) + geom_line(color = colConFirst) +
    labs(x = "Time (h)", y = "sink_pheme_c Lb & Ub (mmol/gDW*h)",
         title = paste0("Model dynamics ", "(", tag, " | ", Condition, ")"),
         subtitle = paste0("Double bounded flux (Lb =< v =< Ub)")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.70) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.70)
  
  p = (((C|B|D)/(I|Hc|He|A))/(fB|fH|fA)/(bB|bH|LbA|UbA))
  
  ggsave(plot = p, filename = paste0(pdf.name, "_", tag, "_", Condition, ".pdf"), width = 24, height = 16)
  
  system(paste0("cp ./", pdf.name, "_", tag, "_", 
                Condition, ".pdf ", results.dir, "CDiff_", 
                tag, "_", Condition, "_", Exper))
  
  system(paste0("mv ./", pdf.name, "_", tag, "_", 
                Condition, ".pdf ", results.dir, "CDiff_", tag,
                "_", Condition, "_", Exper))
  
}
