
plotting_panel_paper = function(Exper,
                                Condition,
                                tag,
                                param_target,
                                Tempi,
                                colo1,
                                colo2,
                                colo3,
                                wd) {
  
  Na = 6.022e20
  c = 6.022e08
  pack = 1*(Na*(1/c))
  fc = 1e-06
  
  coloTag = data.frame(SuperAblated = "darkred", Ablated = "#EC6D67", Unified = "#5351A2")
  nameTag = data.frame(SuperAblated = "Ablated", Ablated = "Partially Ablated", Unified = "Unified")
  
  aa_places <- c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  places <- c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", aa_places)
  units = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))
  
  subtrace <- do.call(rbind, lapply(tag, function(j) {
    setwd(file.path(wd, paste0("results/CDiff", "_", j, "_", Condition, "_", Exper)))
    load(file = paste0("subtrace_" , j, "_", Condition, ".RData"))
    cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
  }))
  
  setwd(wd)
  
  subflux <- do.call(rbind, lapply(tag, function(j) {
    setwd(file.path(wd, paste0("results/CDiff", "_", j, "_", Condition, "_", Exper)))
    load(file = paste0("subflux_" , j, "_", Condition, ".RData"))
    cbind(subflux, Scenario = rep(j, length(subflux$Time)))
  }))
  
  AA <- subtrace[subtrace$Time %in% Tempi & subtrace$Places %in% aa_places, ]
  AA$Time <- as.character(AA$Time)
  
  I <- subtrace[subtrace$Places == "IECs" & subtrace$Time %in% Tempi, ]
  I$Time <- as.character(I$Time)
  
  create_plot <- function(place, j, 
                          y_label, 
                          s = 12, st = 12, at = 12) {
    subtrace_filtered <- filter(subtrace[which(subtrace$Places == place), ], Scenario == j)
    plot <- ggplot(subtrace_filtered, 
                   aes(x = Time, y = if (place == "IECs") Marking 
                       else ((Marking/pack)*1000), color = IECsDeath)) + 
      labs(x = "Time (h)", y = y_label, 
           title = paste0(nameTag[, j], " Experiments"), 
           colour = paste0("Param: ", param_target)) + 
      theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
            plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
            plot.title.position = "plot",
            axis.text = element_text(size = at, color = "black"),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "none") + 
      geom_line(aes(group = config)) +
      scale_colour_gradientn(colors = colo1, oob = scales::squish)
    return(plot)
  }
  
  IECs_plots <- list()
  trp_L_e_plots <- list()
  leu_L_e_plots <- list()
  
  for (j in tag) {
    IECs_plots[[j]] <- create_plot(place = "IECs", j = j, y_label = paste0("IECs" , " ", "(cell)"))
    trp_L_e_plots[[j]] <- create_plot(place = "trp_L_e", j = j, y_label = paste0("Tryptophan" , " ", "(µmol)"))
    leu_L_e_plots[[j]] <- create_plot(place = "leu_L_e", j = j, y_label = paste0("Leucine" , " ", "(µmol)"))
  }
  
  create_plot <- function(place, j) {
    data <- if (place == "IECs") filter(I, Scenario == j) 
    else filter(AA, AA$Places == place, Scenario == j)
    
    y_label <- if (place == "IECs") paste0("IECs" , " ", "(cell)") 
    else if (place == "trp_L_e") paste0("Tryptophan" , " ", "(µmol)")
    else paste0("Leucine" , " ", "(µmol)")
    
    plot_name <- paste(place, "violin", param_target, j, sep = "_")
    
    plot <- data %>% 
      ggplot(aes(y = as.numeric(Time), 
                 x = if (place == "IECs") Marking 
                 else ((Marking/pack)*1000), group = Time)) +
      labs(y = "Time (h)", 
           x = y_label,
           title = paste0(nameTag[, j], " Experiments")) +
      theme(text = element_text(size = 8, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
            legend.position = "none",
            plot.title = element_text(size = 12, face = "bold", color = coloTag[, j]),
            plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 12, face = "bold")) +
      geom_jitter(aes(col = IECsDeath), alpha = 0.55, size = 1.5) +
      geom_violin(scale = "width", alpha = 0.5, width = 0.25) +
      geom_boxplot(alpha = 0) + 
      scale_colour_gradientn(colors = colo1, oob = scales::squish)
    
    if (place == "IECs") {
      # Store x-axis limits when place == "IECs"
      x_limits <<- ggplot_build(plot)$layout$panel_params[[1]]$x.range
      plot <- plot + scale_x_continuous(breaks = scales::pretty_breaks(n = 3), 
                                        labels = function(x) format(x, scientific = FALSE))
    } else if (!is.null(x_limits)) {
      # Use stored x-axis limits for all other plots
      plot <- plot + coord_cartesian(xlim = x_limits)
    }
    
    return(plot)
  }
  
  plot_list <- list()
  
  for (j in tag) {
    for (place in c("IECs", "leu_L_e", "trp_L_e")) {
      plot_name <- paste(place, "violin", param_target, j, sep = "_")
      plot_list[[plot_name]] <- create_plot(place, j)
    }
  }
  
  ######
  
  lineCon = data.frame(Therapy = "solid", NoDrug = "dashed")
  
  subtrace = list()
  subflux = list()
  
  for (tag in c("Ablated", "Unified")) {
    for (condition in c("NoDrug", "Therapy")) {
      file_prefix = paste0(wd, "/results/CDiff_", tag, "_", condition, "_", "Model_Analysis", "/")
      subtrace[[paste0(tag, condition)]] = readr::read_csv(paste0(file_prefix, "SubTrace", tag, condition, ".csv"))[-1]
      subflux[[paste0(tag, condition)]] = readr::read_csv(paste0(file_prefix, "SubFlux", tag, condition, ".csv"))[-1]
    }
  }
  
  subtrace = do.call(rbind, subtrace)
  subflux = do.call(rbind, subflux)
  
  Places.all = ggplot(subtrace, aes(Time, Marking, colour = Scenario, linetype = Condition)) +
    theme_bw() + geom_line(linewidth = 0.5) +
    scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
    scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
    theme(axis.text = element_text(size = 12), 
          axis.title = element_text(size = 12, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none", 
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          axis.ticks = element_line(size = 0.5)) +
    geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.25) +
    facet_wrap(~ Places, scales="free", ncol = 3)
  
  plot_place <- function(place, unit, subtrace, Sce, lineCon, coloTag) {
    ggplot(filter(subtrace, subtrace$Places == place & subtrace$Scenario == Sce), 
           aes(Time, 
               if (place %in% aa_places) ((Marking/pack)*1000)
               else if (place %in% c("Drug", "pheme_e", "pheme_c")) Marking*fc
               else Marking, 
               colour = Scenario, linetype = Condition)) + 
      theme_bw() + geom_line(size = 0.8) + 
      scale_color_manual(values = coloTag[, Sce]) +
      scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
      labs(x = "Time (h)", 
           y = paste0(place, " ", unit)) +
      theme(axis.text = element_text(size = 12), 
            axis.title = element_text(size = 12, face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none", 
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
            axis.ticks = element_line(size = 1.75)) + 
      geom_vline(xintercept = seq(8, 48, by = 8), linetype = "dotted", color = "darkblue", lwd = 0.75)
  }
  
  miniplot_list <- list()
  
  for (i in 1:length(places)) {
    miniplot_list[[places[i]]] <- plot_place(places[i], units[i], subtrace, "Unified", lineCon, coloTag)
  }

  plot_reaction <- function(subflux, subtrace, Sce, coloTag, lineCon) {
    # Define the reactions
    reactions = c("EX_biomass_e", "EX_trp_L_e", "EX_leu_L_e")
    
    # Initialize an empty list to store the plots
    plot_list <- list()
    
    # Loop over each reaction
    for (reaction in reactions) {
      # Generate the plot for the current reaction
      plot <- ggplot(filter(subflux, subflux$Reaction == reaction & subflux$Scenario == Sce),
                     aes(Time, Flux, colour = Scenario, linetype = Condition)) + 
        theme_bw() + geom_line(size = 0.8) + 
        scale_color_manual(values = coloTag[, "Unified"]) +
        scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
        labs(x = "Time (h)",
             y = paste0(reaction, " (mmol/gDW*h)")) +
        theme(axis.text = element_text(size = 12), 
              axis.title = element_text(size = 12, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none", 
              panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
              axis.ticks = element_line(size = 1.75)) + 
        geom_vline(xintercept = seq(8, 48, by = 8), linetype = "dotted", color = "darkblue", lwd = 0.75)
      
      # Add the plot to the list
      plot_list[[reaction]] <- plot
    }
    
    # Return the list of plots
    return(plot_list)
  }
  
  plot_list_reaction <- plot_reaction(subflux, subtrace, "Unified", coloTag, lineCon)
  
  #####
  
  mini_pl = c(miniplot_list, plot_list_reaction)
  
  for(i in 1:length(mini_pl) ) {
    ggsave(mini_pl[[i]], 
           file = paste0(wd, "/", names(mini_pl)[i], ".pdf"), 
           width = 6, height = 3)
    system(paste0("cp ", 
                  wd, "/", names(mini_pl)[i], ".pdf", " ", 
                  wd, "/results/images4paper"))
    system(paste0("mv ", 
                  wd, "/", names(mini_pl)[i], ".pdf", " ", 
                  wd, "/results/images4paper"))
  }
  
  ggsave((plot_list[["IECs_violin_IECsDeath_SuperAblated"]] |
            leu_L_e_plots[["SuperAblated"]] |
            trp_L_e_plots[["SuperAblated"]]) / 
           (plot_list[["IECs_violin_IECsDeath_Ablated"]] |
              leu_L_e_plots[["Ablated"]] |
              trp_L_e_plots[["Ablated"]]) / 
           (plot_list[["IECs_violin_IECsDeath_Unified"]] |
              leu_L_e_plots[["Unified"]] |
              trp_L_e_plots[["Unified"]]), 
         file = paste0(wd, "/Figure3A.pdf"), width = 8, height = 7)
  
  system(paste0("cp ", wd, "/Figure3A.pdf ", wd, "/results/images4paper"))
  system(paste0("mv ", wd, "/Figure3A.pdf ", wd, "/results/images4paper"))
  
  ggsave(Places.all, file = paste0(wd, "/Allplaces.pdf"), width = 12, height = 12)
  system(paste0("cp ", wd, "/Allplaces.pdf ", wd, "/results/images4paper"))
  system(paste0("mv ", wd, "/Allplaces.pdf ", wd, "/results/images4paper"))
  
}
