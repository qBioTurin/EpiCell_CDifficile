
PlottingAblation = function(Exper,
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
  
  Cdiff <- dplyr::filter(subtrace, Places == "CD")
  Biomass <- dplyr::filter(subtrace, Places == "BiomassCD")
  Drug <- dplyr::filter(subtrace, Places == "Drug")
  IECs <- dplyr::filter(subtrace, Places == "IECs")
  pheme_e <- dplyr::filter(subtrace, Places == "pheme_e")
  pheme_c <- dplyr::filter(subtrace, Places == "pheme_c")
  AA <- dplyr::filter(subtrace, Places %in% aa.placesvec)
  
  FluxBiom <- dplyr::filter(subflux, Reaction == "EX_biomass_e")
  FluxPheme <- dplyr::filter(subflux, Reaction == "sink_pheme_c")
  FluxAA <- dplyr::filter(subflux, Reaction %in% aa.reactionsvec)
  
  BoundsBiom <- dplyr::filter(subbounds, Reaction %in% c("EX_biomass_e_Ub", "EX_biomass_e_Lb"))
  BoundsAALb <- dplyr::filter(subbounds, Reaction %in% paste0(aa.reactionsvec, "_Lb"))
  BoundsAAUb <- dplyr::filter(subbounds, Reaction %in% paste0(aa.reactionsvec, "_Ub"))
  
  # Plotting
  for (Con in Condition) {
    assign(paste("C", Con, sep = "."), 
           ggplot(dplyr::filter(Cdiff, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
             geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "C. difficile (cells)",
                  title = paste0("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")"),
                  subtitle = paste0("Place: ", Cdiff[1, 2], " (cells)")) +
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
             scale_y_continuous(label = scientific_10))
    
    assign(paste("B", Con, sep = "."), 
           ggplot(dplyr::filter(Biomass, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
             geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "C. difficile average biomass (pg)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Place: ", Biomass[1, 2], " (pg)", sep = "")) +
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
             geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.3) +
             annotate("text", 7, 1.5, label = "biomass max", vjust = 0.5, color = "black", size = 3) +
             geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.3) +
             annotate("text", 7, 0.6, label = "biomass mean", vjust = 0.5, color = "black", size = 3) +
             geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.3) +
             annotate("text", 7, 0.15, label = "biomass min", vjust = 0.5, color = "black", size = 3))
    
    assign(paste("D", Con, sep = "."), 
           ggplot(dplyr::filter(Drug, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
             geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "C. difficile # cells",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Place: ", Drug[1, 2], " (pmol)", sep = "")) +
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
             geom_hline(yintercept = 5800, linetype = "dashed", color = "darkgrey", lwd = 0.3) +
             annotate("text", 3, 6500, label = "MIC value", vjust = 0.5, color = "black", size = 3))
    
    # assign(paste("Txn", Con, sep = "."), 
    #        ggplot(dplyr::filter(Toxin, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
    #          geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
    #            colConNoDrug} else {colConTherapy})) +
    #          labs(x = "Time (h)", y = "Toxin concentration (pmol)",
    #               title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
    #               subtitle = paste("Place: ", Toxin[1, 2], " (pmol)", sep = "")) +
    #          theme(axis.text = element_text(size = 14), 
    #                axis.title = element_text(size = 14),
    #                plot.title = element_text(size = 16),
    #                plot.subtitle = element_text(size = 14)) + 
    #          geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    #          geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    #          geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    #          geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    #          geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    #          geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    #          geom_hline(yintercept = 5800, linetype = "dashed", color = "darkgrey", lwd = 0.3))
    
    assign(paste("I", Con, sep = "."), 
           ggplot(dplyr::filter(IECs, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
             geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "IECs # cells (cell)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Place: ", IECs[1, 2], " (cell)", sep = "")) +
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
             scale_y_continuous(label = scientific_10))
    
    assign(paste("He", Con, sep = "."), 
           ggplot(dplyr::filter(pheme_e, Condition == Con), aes(Time, Marking, colour = Scenario)) + theme_bw() +
             geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "Extr heme concentration (pmol)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Place: ", pheme_e[1, 2], " (pmol)", sep = "")) +
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
             geom_hline(yintercept = 1e+06, linetype = "dashed", color = "darkred", lwd = 0.3) +
             annotate("text", 10, 0.95e+06, label = "Hemorrhagic baseline", vjust = 0.5, color = "black", size = 3) +
             scale_y_continuous(label = scientific_10))
    
    assign(paste("A", Con, sep = "."), 
           ggplot(dplyr::filter(AA, Condition == Con), aes(Time, Marking, colour = Scenario)) +
             theme_bw() + geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "Aminoacids concentration (pmol)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Place: aminoacids", " (molecule)", sep = "")) +
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
             scale_y_continuous(label = scientific_10) +
             facet_wrap(~Places, nrow = 2))
    
    assign(paste("fB", Con, sep = "."), 
           ggplot(dplyr::filter(FluxBiom, Condition == Con), aes(Time, Flux, colour = Scenario)) + theme_bw() +
             geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "EX_biomass(e) FBA flux (mmol/gDW*h)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Reaction: ", FluxBiom[1, 2], " (pg)", sep = "")) +
             theme(axis.text = element_text(size = 14), 
                   axis.title=element_text(size = 14),
                   plot.title = element_text(size = 16),
                   plot.subtitle = element_text(size = 14)) + 
             geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3))
    
    assign(paste("fH", Con, sep = "."), 
           ggplot(dplyr::filter(FluxPheme, Condition == Con), aes(Time, Flux, colour = Scenario)) + theme_bw() +
             geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "sink_pheme(c) FBA flux (mmol/gDW*h)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Reaction: ", FluxPheme[1, 2], " (pmol)", sep = "")) +
             theme(axis.text = element_text(size = 14), 
                   axis.title=element_text(size = 14),
                   plot.title = element_text(size = 16),
                   plot.subtitle = element_text(size = 14)) + 
             geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3))
    
    assign(paste("fA", Con, sep = "."), 
           ggplot(dplyr::filter(FluxAA, Condition == Con), aes(Time, Flux, colour = Scenario)) + 
             theme_bw() + geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "EX_Aminoacid(e) FBA flux (mmol/gDW*h)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Reaction: EX_Aminoacid_e", " (molecule)", sep = "")) +
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
             scale_y_continuous(label = scientific_10) +
             facet_wrap(~Reaction, nrow = 2))
    
    ################################################################################
    
    assign(paste("bB", Con, sep = "."), 
           ggplot(dplyr::filter(BoundsBiom, Condition == Con), aes(Time, Flux, colour = Scenario, linetype = Reaction)) + 
             theme_bw() + geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "EX_biomass(e) Lb & Ub (mmol/gDW*h)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
             theme(axis.text = element_text(size = 14), 
                   axis.title=element_text(size = 14),
                   plot.title = element_text(size = 16),
                   plot.subtitle = element_text(size = 14)) + 
             geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
             geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3))
    
    assign(paste("LbA", Con, sep = "."), 
           ggplot(dplyr::filter(BoundsAALb, Condition == Con), aes(Time, Flux, colour = Scenario)) + 
             theme_bw() + geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "EX_Aminocid(e) Lb (mmol/gDW*h)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
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
             facet_wrap(~Reaction, nrow = 2))
    
    assign(paste("UbA", Con, sep = "."), 
           ggplot(dplyr::filter(BoundsAAUb, Condition == Con), aes(Time, Flux, colour = Scenario)) + 
             theme_bw() + geom_line(size = 0.8) + scale_color_manual(values = c(coltag, if(Con == "NoDrug") {
               colConNoDrug} else {colConTherapy})) +
             labs(x = "Time (h)", y = "EX_Aminoacid(e) Ub (mmol/gDW*h)",
                  title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", Con, ")", sep = ""),
                  subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
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
             facet_wrap(~Reaction, nrow = 2))
    
  }
  
  # Plotting all
  C.All = ggplot(Cdiff, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "C. difficile # cells",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", Cdiff[1, 2], " (cell)", sep = "")) +
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
    scale_y_continuous(label = scientific_10)
  
  B.All = ggplot(Biomass, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "C. difficile average biomass (pg)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", Biomass[1, 2], " (pg)", sep = "")) +
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
    geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.3) +
    annotate("text", 7, 1.5, label = "biomass max", vjust = 0.5, color = "black", size = 3) +
    geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.3) +
    annotate("text", 7, 0.6, label = "biomass mean", vjust = 0.5, color = "black", size = 3) +
    geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.3) +
    annotate("text", 7, 0.15, label = "biomass min", vjust = 0.5, color = "black", size = 3)
  
  D.All = ggplot(Drug, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "Drug concentration (pmol)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", Drug[1, 2], " (pmol)", sep = "")) +
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
    geom_hline(yintercept = 5800, linetype = "dashed", color = "darkgrey", lwd = 0.3) +
    annotate("text", 3, 6500, label = "MIC value", vjust = 0.5, color = "black", size = 3)
  
  # T.All = ggplot(Toxin, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
  #   geom_line(size = 0.8) + 
  #   scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  #   scale_linetype_manual(values = c("dashed", "solid")) +
  #   labs(x = "Time (h)", y = "Toxin concentration (pmol)",
  #        title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
  #        subtitle = paste("Place: ", Toxin[1, 2], " (pmol)", sep = "")) +
  #   theme(axis.text = element_text(size = 14), 
  #         axis.title=element_text(size = 14),
  #         plot.title = element_text(size = 16),
  #         plot.subtitle = element_text(size = 14)) + 
  #   geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3)
  
  I.All = ggplot(IECs, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "IECs # cells (cell)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", IECs[1, 2], " (cell)", sep = "")) +
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
    scale_y_continuous(label = scientific_10)
  
  He.All = ggplot(pheme_e, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "Extr heme concentration (pmol)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", pheme_e[1, 2], " (pmol)", sep = "")) +
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
    geom_hline(yintercept = 1e+06, linetype = "dashed", color = "darkred", lwd = 0.3) +
    annotate("text", 10, 0.95e+06, label = "Hemorrhagic baseline", vjust = 0.5, color = "black", size = 3) +
    scale_y_continuous(label = scientific_10)
  
  A.All = ggplot(AA, aes(Time, Marking, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "Aminoacids concentration (pmol)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: aminoacids", " (molecule)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    scale_y_continuous(label = scientific_10) +
    facet_wrap(~Places, nrow = 2)
  
  fB.All = ggplot(FluxBiom, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "EX_biomass(e) FBA flux (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Reaction: ", FluxBiom[1, 2], " (pg)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3)
  
  fH.All = ggplot(FluxPheme, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "sink_pheme(c) FBA flux (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Reaction: ", FluxPheme[1, 2], " (pmol)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3)
  
  fA.All = ggplot(FluxAA, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "EX_Aminoacid(e) FBA flux (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Reaction: EX_Aminoacid_e", " (molecule)", sep = "")) +
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
    facet_wrap(~Reaction, nrow = 2)
  
  bB.All = ggplot(BoundsBiom, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "EX_biomass(e) Bounds (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
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
    facet_wrap(~Reaction, nrow = 1)
  
  LbA.All = ggplot(BoundsAALb, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "EX_Aminoacid(e) Lb (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
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
    facet_wrap(~Reaction, nrow = 2)
  
  UbA.All = ggplot(BoundsAAUb, aes(Time, Flux, colour = Condition, linetype = Scenario)) + theme_bw() +
    geom_line(size = 0.8) + scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Time (h)", y = "EX_Aminoacid(e) Ub (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
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
    facet_wrap(~Reaction, nrow = 2)
  
  ######
  
  C.Unified = ggplot(dplyr::filter(Cdiff, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
    theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "C. difficile # cells",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", Cdiff[1, 2], " (cells)", sep = "")) +
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
    scale_y_continuous(label = scientific_10)
  
  B.Unified = ggplot(dplyr::filter(Biomass, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
    theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "C. difficile average biomass (pg)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", Biomass[1, 2], " (pg)", sep = "")) +
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
    geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.3) +
    annotate("text", 7, 1.5, label = "biomass max", vjust = 0.5, color = "black", size = 3) +
    geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.3) +
    annotate("text", 7, 0.6, label = "biomass mean", vjust = 0.5, color = "black", size = 3) +
    geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.3) +
    annotate("text", 7, 0.15, label = "biomass min", vjust = 0.5, color = "black", size = 3)
  
  D.Unified = ggplot(dplyr::filter(Drug, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
    theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "Drug concentration (pmol)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", Drug[1, 2], " (pmol)", sep = "")) +
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
    geom_hline(yintercept = 5800, linetype = "dashed", color = "darkgrey", lwd = 0.3) +
    annotate("text", 3, 6500, label = "MIC value", vjust = 0.5, color = "black", size = 3)
  
  # T.Unified = ggplot(dplyr::filter(Toxin, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
  #   theme_bw() +
  #   geom_line(size = 0.8) + 
  #   scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
  #   labs(x = "Time (h)", y = "Toxin concentration (pmol)",
  #        title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
  #        subtitle = paste("Place: ", Toxin[1, 2], " (pmol)", sep = "")) +
  #   theme(axis.text = element_text(size = 14), 
  #         axis.title=element_text(size = 14),
  #         plot.title = element_text(size = 16),
  #         plot.subtitle = element_text(size = 14)) + 
  #   geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
  #   geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3)
  
  I.Unified = ggplot(dplyr::filter(IECs, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
    theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "IECs # cells (cell)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", IECs[1, 2], " (cell)", sep = "")) +
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
    scale_y_continuous(label = scientific_10)
  
  He.Unified = ggplot(dplyr::filter(pheme_e, Scenario == "Unified"), aes(Time, Marking, colour = Condition)) + 
    theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "Extr heme concentration (pmol)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: ", pheme_e[1, 2], " (pmol)", sep = "")) +
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
    geom_hline(yintercept = 1e+06, linetype = "dashed", color = "darkred", lwd = 0.3) +
    annotate("text", 10, 0.95e+06, label = "Hemorrhagic baseline", vjust = 0.5, color = "black", size = 3) +
    scale_y_continuous(label = scientific_10)
  
  A.Unified = ggplot(dplyr::filter(AA, Scenario == "Unified"), 
                     aes(Time, (Marking/(6.022e20/6.022e08))*1e+09, colour = Condition)) + 
    theme_bw() + geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "Aminoacids concentration (pmol)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Place: aminoacids", " (molecule)", sep = "")) +
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
    scale_y_continuous(label = scientific_10) +
    facet_wrap(~Places, nrow = 2)
  
  ######
  
  fB.Unified = ggplot(dplyr::filter(FluxBiom, Scenario == "Unified"), aes(Time, Flux, colour = Condition)) + 
    theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "EX_biomass(e) FBA flux (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Reaction: ", FluxBiom[1, 2], " (pg)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3)
  
  fH.Unified = ggplot(dplyr::filter(FluxPheme, Scenario == "Unified"), aes(Time, Flux, colour = Condition)) + 
    theme_bw() +
    geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "sink_pheme(c) FBA flux (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Reaction: ", FluxPheme[1, 2], " (pmol)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3)
  
  fA.Unified = ggplot(dplyr::filter(FluxAA, Scenario == "Unified"), 
                      aes(Time, Flux, colour = Condition)) + 
    theme_bw() + geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "EX_Aminoacid(e) FBA flux (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Reaction: EX_Aminoacid_e", " (molecule)", sep = "")) +
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
    facet_wrap(~Reaction, nrow = 2)
  
  ######
  
  bB.Unified = ggplot(dplyr::filter(BoundsBiom, Scenario == "Unified"), 
                      aes(Time, Flux, colour = Condition, linetype = Reaction)) + 
    theme_bw() + geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "EX_biomass(e) Lb & Ub (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
    theme(axis.text = element_text(size = 14), 
          axis.title=element_text(size = 14),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14)) + 
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.3) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.3)
  
  LbA.Unified = ggplot(dplyr::filter(BoundsAALb, Scenario == "Unified"), 
                       aes(Time, Flux, colour = Condition)) + theme_bw() +
    geom_line(size = 0.8) + scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "EX_Aminoacid(e) Lb (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
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
    facet_wrap(~Reaction, nrow = 2)
  
  UbA.Unified = ggplot(dplyr::filter(BoundsAAUb, Scenario == "Unified"), 
                       aes(Time, Flux, colour = Condition)) + 
    theme_bw() + geom_line(size = 0.8) + 
    scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
    labs(x = "Time (h)", y = "EX_Aminoacid(e) Ub (mmol/gDW*h)",
         title = paste("Model dynamics ", "(", "Unified", " | ", "NoDrug VS Therapy", ")", sep = ""),
         subtitle = paste("Double bounded flux (Lb ≤ v ≤ Ub)", sep = "")) +
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
    facet_wrap(~Reaction, nrow = 2)
  
  # Define a list of input variables and output file names
  plots <- list(
    list(C.NoDrug, C.Therapy, C.Unified, C.All, "AblationAnalysis_CD.pdf"),
    list(B.NoDrug, B.Therapy, B.Unified, B.All, "AblationAnalysis_BiomassCD.pdf"),
    list(D.NoDrug, D.Therapy, D.Unified, D.All, "AblationAnalysis_Drug.pdf"),
    list(I.NoDrug, I.Therapy, I.Unified, I.All, "AblationAnalysis_IECs.pdf"),
    list(He.NoDrug, He.Therapy, He.Unified, He.All, "AblationAnalysis_pheme_e.pdf"),
    list(A.NoDrug, A.Therapy, A.Unified, A.All, "AblationAnalysis_AA.pdf"),
    list(fB.NoDrug, fB.Therapy, fB.Unified, fB.All, "AblationAnalysis_EX_biomass_e.pdf"),
    list(fH.NoDrug, fH.Therapy, fH.Unified, fH.All, "AblationAnalysis_sink_pheme_c.pdf"),
    list(fA.NoDrug, fA.Therapy, fA.Unified, fA.All, "AblationAnalysis_EX_aminocids_e.pdf"),
    list(bB.NoDrug, bB.Therapy, bB.Unified, bB.All, "AblationAnalysis_EX_biomass_Lb_Ub_e.pdf"),
    list(LbA.NoDrug, LbA.Therapy, LbA.Unified, LbA.All, "AblationAnalysis_EX_aminocids_e_Lb.pdf"),
    list(UbA.NoDrug, UbA.Therapy, UbA.Unified, UbA.All, "AblationAnalysis_EX_aminocids_e_Ub.pdf")
  )
  
  # Loop over the list of input variables and output file names
  for (p in plots) {
    ggsave((p[[1]] | p[[2]]) / (p[[3]] | p[[4]]),
           file = paste(p[[5]], sep = ""),
           width = 14, height = 8)
  }
  
  
  system(paste0("rm -r ./results/ScenarioComparison"))
  system(paste0("mkdir ./results/ScenarioComparison"))
  
  system(paste0("cp AblationAnalysis_*.pdf ", "./results/ScenarioComparison/"))
  system(paste0("mv AblationAnalysis_*.pdf ", "./results/ScenarioComparison/"))
  
}
