

plotting_ablation = function(Exper,
                             colConNoDrug,
                             colConTherapy,
                             coltag, 
                             tag, 
                             Condition,
                             wd) {
  
  results.dir <- paste0(wd, "/results/")
  
  Na = 6.022e20
  c = 6.022e08
  pack = 1*(Na*(1/c))
  fc = 1e-06
  
  # Define the file types and places
  file_types <- c("SubTrace", "SubFlux", "SubBounds")
  
  places <- c("CD", "IECs", "BiomassCD", 
              "Drug", "pheme_e", "pheme_c", 
              "trp_L_e", "pro_L_e", "val_L_e", 
              "ile_L_e", "cys_L_e", "leu_L_e")
  
  reactions <- c("EX_biomass_e", "sink_pheme_c", 
                 "EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e", 
                 "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")
  
  bounds_Lb <- paste0(reactions, "_Lb")
  bounds_Ub <- paste0(reactions, "_Ub")
  
  units_places = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))
  units_reactions = c(rep("(mmol/gDW*h)", length(reactions)))
  units_bounds_Lb = c(rep("(mmol/gDW*h)", length(reactions)))
  units_bounds_Ub = c(rep("(mmol/gDW*h)", length(reactions)))
  
  # Read the files and assign them to variables
  for (i in Condition) {
    for (j in tag) {
      for (file_type in file_types) {
        assign(paste0("sub", tolower(file_type), j, i), 
               readr::read_csv(
                 paste0(
                   results.dir, "CDiff_", j, "_", i, "_", 
                   Exper, "/", file_type, j, i, ".csv"))[-1])
      }
    }
  }
  
  # Combine the data frames
  subtrace <- do.call(rbind, mget(ls(pattern = "subtrace")))
  subflux <- do.call(rbind, mget(ls(pattern = "subflux")))
  subbounds <- do.call(rbind, mget(ls(pattern = "subbounds")))
  
  # Initialize lists to store the data frames
  places_list <- list()
  reactions_list <- list()
  bounds_Lb_list <- list()
  bounds_Ub_list <- list()
  
  # Filter the data frames
  for (place in places) {
    places_list[[place]] <- dplyr::filter(subtrace, Places == place)
  }
  
  for (reaction in reactions) {
    reactions_list[[reaction]] <- dplyr::filter(subflux, Reaction == reaction)
    bounds_Lb_list[[paste0(reaction, "_Lb")]] <- dplyr::filter(subbounds, Reaction == paste0(reaction, "_Lb"))
    bounds_Ub_list[[paste0(reaction, "_Ub")]] <- dplyr::filter(subbounds, Reaction == paste0(reaction, "_Ub"))
    }
  
  
  # Define the function to plot data frames
  plot_data_places <- function(df, 
                               place, 
                               unit, 
                               Scenario, 
                               coltag, 
                               Con, 
                               colConNoDrug, 
                               colConTherapy, 
                               fc, pack) {
    
    if ( place %in% c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e") ) {
      df$Marking = (df$Marking/pack)*1000
    } else if ( place %in% c("Drug", "pheme_e", "pheme_c") ) {
      df$Marking = df$Marking*fc
    }
    
    ggplot(dplyr::filter(df, Condition == Con), 
           aes(Time, Marking, colour = Scenario)) + 
      theme_bw() + geom_line(size = 0.8) + 
      scale_color_manual(
        values = c(coltag, if(Con == "NoDrug") colConNoDrug 
                   else colConTherapy)) +
      labs(x = "Time (h)", y = paste0(place, " ", unit),
           title = paste0("(", "Ablated VS Unified", " | ", Con, ")"),
           subtitle = paste0("Place: ", place, " ", unit)) +
      theme(axis.text = element_text(size = 14), 
            axis.title=element_text(size = 14),
            plot.title = element_text(size = 16),
            plot.subtitle = element_text(size = 14)) + 
      geom_vline(xintercept = seq(8, 48, by = 8), 
                 linetype = "dotted", color = "darkblue", lwd = 0.3)
  }
  
  # Define the function to plot data frames
  plot_data_reactions <- function(df, 
                                  reaction, 
                                  unit, 
                                  Scenario, 
                                  coltag, 
                                  Con, 
                                  colConNoDrug, 
                                  colConTherapy) {
    
    ggplot(dplyr::filter(df, Condition == Con), 
           aes(Time, Flux, colour = Scenario)) + 
      theme_bw() + geom_line(size = 0.8) + 
      scale_color_manual(
        values = c(coltag, if(Con == "NoDrug") colConNoDrug 
                   else colConTherapy)) +
      labs(x = "Time (h)", y = paste0(reaction, " ", unit),
           title = paste0("(", "Ablated VS Unified", " | ", Con, ")"),
           subtitle = paste0("Reaction: ", reactions, " ", unit)) +
      theme(axis.text = element_text(size = 14), 
            axis.title=element_text(size = 14),
            plot.title = element_text(size = 16),
            plot.subtitle = element_text(size = 14)) + 
      geom_vline(xintercept = seq(8, 48, by = 8), 
                 linetype = "dotted", color = "darkblue", lwd = 0.3)
  }
  
  # Define the function to plot data frames
  plot_data_bounds_Lb <- function(df, 
                                  bound_Lb, 
                                  unit, 
                                  Scenario, 
                                  coltag, 
                                  Con, 
                                  colConNoDrug, 
                                  colConTherapy) {
    
    ggplot(dplyr::filter(df, Condition == Con), 
           aes(Time, Flux, colour = Scenario)) + 
      theme_bw() + geom_line(size = 0.8) + 
      scale_color_manual(
        values = c(coltag, if(Con == "NoDrug") colConNoDrug 
                   else colConTherapy)) +
      labs(x = "Time (h)", y = paste0(bound_Lb, " ", unit),
           title = paste0("(", "Ablated VS Unified", " | ", Con, ")"),
           subtitle = paste0("Bound_Lb: ", bounds_Lb, " ", unit)) +
      theme(axis.text = element_text(size = 14), 
            axis.title=element_text(size = 14),
            plot.title = element_text(size = 16),
            plot.subtitle = element_text(size = 14)) + 
      geom_vline(xintercept = seq(8, 48, by = 8), 
                 linetype = "dotted", color = "darkblue", lwd = 0.3)
  }
  
  # Define the function to plot data frames
  plot_data_bounds_Ub <- function(df, 
                                  bound_Ub, 
                                  unit, 
                                  Scenario, 
                                  coltag, 
                                  Con, 
                                  colConNoDrug, 
                                  colConTherapy) {
    
    ggplot(dplyr::filter(df, Condition == Con), 
           aes(Time, Flux, colour = Scenario)) + 
      theme_bw() + geom_line(size = 0.8) + 
      scale_color_manual(
        values = c(coltag, if(Con == "NoDrug") colConNoDrug 
                   else colConTherapy)) +
      labs(x = "Time (h)", y = paste0(bound_Ub, " ", unit),
           title = paste0("(", "Ablated VS Unified", " | ", Con, ")"),
           subtitle = paste0("Bound_Ub: ", bounds_Ub, " ", unit)) +
      theme(axis.text = element_text(size = 14), 
            axis.title=element_text(size = 14),
            plot.title = element_text(size = 16),
            plot.subtitle = element_text(size = 14)) + 
      geom_vline(xintercept = seq(8, 48, by = 8), 
                 linetype = "dotted", color = "darkblue", lwd = 0.3)
  }
  
  # Apply the plot function
  apply_plot_place <- function(Con) {
    plots_list <- list()
    for (place in places) {
      plots_list[[place]] = 
        plot_data_places(df = places_list[[place]], 
                         place = place,
                         unit = units_places[places == place], 
                         Scenario = Scenario, 
                         coltag = coltag, 
                         Con = Con, 
                         colConNoDrug = colConNoDrug, 
                         colConTherapy = colConTherapy, 
                         fc = fc, pack = pack)
    }
    return(plots_list)
  }
  
  # Apply the plot function
  apply_plot_reaction <- function(Con) {
    plots_list <- list()
    for (reaction in reactions) {
      plots_list[[reaction]] = 
        plot_data_reactions(df = reactions_list[[reaction]], 
                            reaction = reaction,
                            unit = units_reactions[reactions == reaction], 
                            Scenario = Scenario, 
                            coltag = coltag, 
                            Con = Con, 
                            colConNoDrug = colConNoDrug, 
                            colConTherapy = colConTherapy)
    }
    return(plots_list)
  }
  
  # Apply the plot function
  apply_plot_bound_Lb <- function(Con) {
    plots_list <- list()
    for (bound_Lb in bounds_Lb) {
      plots_list[[bound_Lb]] = 
        plot_data_bounds_Lb(df = bounds_Lb_list[[bound_Lb]], 
                            bound_Lb = bound_Lb,
                            unit = units_bounds_Lb[bounds_Lb == bound_Lb], 
                            Scenario = Scenario, 
                            coltag = coltag, 
                            Con = Con,
                            colConNoDrug = colConNoDrug, 
                            colConTherapy = colConTherapy)
    }
    return(plots_list)
  }
  
  # Apply the plot function
  apply_plot_bound_Ub <- function(Con) {
    plots_list <- list()
    for (bound_Ub in bounds_Ub) {
      plots_list[[bound_Ub]] = 
        plot_data_bounds_Ub(df = bounds_Ub_list[[bound_Ub]], 
                            bound_Ub = bound_Ub,
                            unit = units_bounds_Ub[bounds_Ub == bound_Ub], 
                            Scenario = Scenario, 
                            coltag = coltag, 
                            Con = Con,
                            colConNoDrug = colConNoDrug, 
                            colConTherapy = colConTherapy)
    }
    return(plots_list)
  }
  
  plots_place_therapy <- apply_plot_place("Therapy")
  plots_place_no_drug <- apply_plot_place("NoDrug")
  plots_reaction_therapy <- apply_plot_reaction("Therapy")
  plots_reaction_no_drug <- apply_plot_reaction("NoDrug")
  plots_Lb_therapy <- apply_plot_bound_Lb("Therapy")
  plots_Lb_no_drug <- apply_plot_bound_Lb("NoDrug")
  plots_Ub_therapy <- apply_plot_bound_Ub("Therapy")
  plots_Ub_no_drug <- apply_plot_bound_Ub("NoDrug")
  
  # Plotting all
  plot_all_place = list()
  
  for ( name in 1:length(places_list) ) {
    
    df = dplyr::filter(places_list[[name]])
    
    if (names(places_list)[name] %in% c("CD", "IECs", "BiomassCD")) {
      df$Marking = df$Marking
    } else if (names(places_list)[name] %in% c("Drug", "pheme_e", "pheme_c")) {
      df$Marking = (df$Marking*fc)
    } else {
      df$Marking = ((df$Marking/pack)*1000)
    }
    
    plot_all_place[[names(places_list)[name]]] = 
      ggplot(df, aes(Time, Marking, colour = Condition, linetype = Scenario)) + 
      theme_bw() + geom_line(size = 0.8) + 
      scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
      scale_linetype_manual(values = c("dashed", "solid")) +
      labs(x = "Time (h)", y = paste0(names(places_list)[name], " ", units_places[name]),
           title = paste0("(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")"),
           subtitle = paste0("Place: ", names(places_list)[name], " ", units_places[name])) +
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
  }
  
  # Plotting all
  plot_all_reaction = list()
  
  for ( name in 1:length(reactions_list) ) {
    plot_all_reaction[[names(reactions_list)[name]]] = 
      ggplot(dplyr::filter(reactions_list[[name]]), 
             aes(Time, Flux, colour = Condition, linetype = Scenario)) + 
      theme_bw() + geom_line(size = 0.8) + 
      scale_color_manual(values = c(colConNoDrug, colConTherapy)) +
      scale_linetype_manual(values = c("dashed", "solid")) +
      labs(x = "Time (h)", y = paste0(names(reactions_list)[name], " ", units_reactions[name]),
           title = paste0("(", "Ablated VS Unified", " | ", "NoDrug VS Therapy", ")"),
           subtitle = paste0("Reaction: ", names(reactions_list)[name], " ", units_reactions[name])) +
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
  }
  
  place_1 = c("CD", "IECs", "BiomassCD", "Drug")
  place_2 = c("pheme_e", "pheme_c")
  place_3 = c("trp_L_e", "pro_L_e", "val_L_e", 
               "ile_L_e", "cys_L_e", "leu_L_e")
  
  reaction_1 = c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e", 
                 "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")
  
  reaction_2 = c("EX_biomass_e", "sink_pheme_c")
  
  bounds_Lb_1 = paste0(reaction_1, "_Lb")
  bounds_Lb_2 = paste0(reaction_2, "_Lb")
  bounds_Ub_1 = paste0(reaction_1, "_Ub")
  bounds_Ub_2 = paste0(reaction_2, "_Ub")
  
  # Organize your plots
  place_1_plots_therapy <- wrap_plots(plots_place_therapy[place_1], ncol = 2, nrow = 2)
  place_2_plots_therapy <- wrap_plots(plots_place_therapy[place_2], ncol = 2, nrow = 1)
  place_3_plots_therapy <- wrap_plots(plots_place_therapy[place_3], ncol = 3, nrow = 2)
  
  place_1_plots_no_drug <- wrap_plots(plots_place_no_drug[place_1], ncol = 2, nrow = 2)
  place_2_plots_no_drug <- wrap_plots(plots_place_no_drug[place_2], ncol = 2, nrow = 1)
  place_3_plots_no_drug <- wrap_plots(plots_place_no_drug[place_3], ncol = 3, nrow = 2)
  
  place_1_plots_therapy <- wrap_plots(plots_place_therapy[place_1], ncol = 2, nrow = 2)
  place_2_plots_therapy <- wrap_plots(plots_place_therapy[place_2], ncol = 2, nrow = 1)
  place_3_plots_therapy <- wrap_plots(plots_place_therapy[place_3], ncol = 3, nrow = 2)
  
  place_1_plots_no_drug <- wrap_plots(plots_place_no_drug[place_1], ncol = 2, nrow = 2)
  place_2_plots_no_drug <- wrap_plots(plots_place_no_drug[place_2], ncol = 2, nrow = 1)
  place_3_plots_no_drug <- wrap_plots(plots_place_no_drug[place_3], ncol = 3, nrow = 2)
  
  reaction_1_plots_therapy <- wrap_plots(plots_reaction_therapy[reaction_1], ncol = 3, nrow = 2)
  reaction_2_plots_therapy <- wrap_plots(plots_reaction_therapy[reaction_2], ncol = 2, nrow = 1)
  
  reaction_1_plots_no_drug <- wrap_plots(plots_reaction_no_drug[reaction_1], ncol = 3, nrow = 2)
  reaction_2_plots_no_drug <- wrap_plots(plots_reaction_no_drug[reaction_2], ncol = 2, nrow = 1)
  
  bounds_Lb_1_plots_therapy = wrap_plots(plots_Lb_therapy[bounds_Lb_1], ncol = 3, nrow = 2)
  bounds_Lb_2_plots_therapy = wrap_plots(plots_Lb_therapy[bounds_Lb_2], ncol = 2, nrow = 1)
  bounds_Ub_1_plots_therapy = wrap_plots(plots_Ub_therapy[bounds_Ub_1], ncol = 3, nrow = 2)
  bounds_Ub_2_plots_therapy = wrap_plots(plots_Ub_therapy[bounds_Ub_2], ncol = 2, nrow = 1)
  
  bounds_Lb_1_plots_no_drug = wrap_plots(plots_Lb_no_drug[bounds_Lb_1], ncol = 3, nrow = 2)
  bounds_Lb_2_plots_no_drug = wrap_plots(plots_Lb_no_drug[bounds_Lb_2], ncol = 2, nrow = 1)
  bounds_Ub_1_plots_no_drug = wrap_plots(plots_Ub_no_drug[bounds_Ub_1], ncol = 3, nrow = 2)
  bounds_Ub_2_plots_no_drug = wrap_plots(plots_Ub_no_drug[bounds_Ub_2], ncol = 2, nrow = 1)
  
  place_1_plots_all <- wrap_plots(plot_all_place[place_1], ncol = 2, nrow = 2)
  place_2_plots_all <- wrap_plots(plot_all_place[place_2], ncol = 2, nrow = 1)
  place_3_plots_all <- wrap_plots(plot_all_place[place_3], ncol = 3, nrow = 2)
  
  reaction_1_plots_all <- wrap_plots(plot_all_reaction[reaction_1], ncol = 3, nrow = 2)
  reaction_2_plots_all <- wrap_plots(plot_all_reaction[reaction_2], ncol = 2, nrow = 1)

  # Combine all plots
  p_therapy = ((place_1_plots_therapy | place_2_plots_therapy) / place_3_plots_therapy)
  p_no_drug = ((place_1_plots_no_drug | place_2_plots_no_drug) / place_3_plots_no_drug)
  p_all = ((place_1_plots_all | place_2_plots_all) / place_3_plots_all)
  
  r_therapy = (reaction_1_plots_therapy | reaction_2_plots_therapy)
  r_no_drug = (reaction_1_plots_no_drug | reaction_2_plots_no_drug)
  r_all = (reaction_1_plots_all | reaction_2_plots_all)
  
  Lb_therapy = (bounds_Lb_1_plots_therapy | bounds_Lb_2_plots_therapy)
  Lb_no_drug = (bounds_Lb_1_plots_no_drug | bounds_Lb_2_plots_no_drug)
  
  Ub_therapy = (bounds_Ub_1_plots_therapy | bounds_Ub_2_plots_therapy)
  Ub_no_drug = (bounds_Ub_1_plots_no_drug | bounds_Ub_2_plots_no_drug)
  
  # Combine all plots
  ggsave( ( p_all ) ,
          file = paste0(wd, "/places_all.pdf"), 
          width = 18, height = 12)
  
  ggsave( ( r_all ) ,
          file = paste0(wd, "/reactions_all.pdf"), 
          width = 25, height = 7)
  
  ggsave( p_therapy ,
          file = paste0(wd, "/places_therapy.pdf"), 
          width = 16, height = 10)
  
  ggsave( p_no_drug ,
          file = paste0(wd, "/places_no_drug.pdf"), 
          width = 16, height = 10)
  
  ggsave( r_therapy ,
          file = paste0(wd, "/reactions_therapy.pdf"), 
          width = 26, height = 6)
  
  ggsave( r_no_drug ,
          file = paste0(wd, "/reactions_no_drug.pdf"), 
          width = 26, height = 6)
  
  ggsave( (Lb_therapy / Lb_no_drug) ,
          file = paste0(wd, "/bounds_Lb.pdf"), 
          width = 26, height = 12)
  
  ggsave( (Ub_therapy / Ub_no_drug) ,
          file = paste0(wd, "/bounds_Ub_.pdf"), 
          width = 26, height = 12)
  
  system(paste0("rm -r ./results/ScenarioComparison"))
  system(paste0("mkdir ./results/ScenarioComparison"))
  
  system(paste0("cp *.pdf ", "./results/ScenarioComparison/"))
  system(paste0("mv *.pdf ", "./results/ScenarioComparison/"))
  
}
