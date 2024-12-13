library(ggnewscale)

# Exper = "Model_Sensitivity"
# Condition = "Therapy"
# tag = c("Ablated", "Unified")
# param_target = "IECsDeath"
# Tempi = c(0, 12, 24, 36, 48, 60)
# wd = wd
# subtrace = subtrace
# colConfigSets = colors_new_confParams

plotting_Fig2_paper = function(Exper,
                                Condition,
                                tag,
                                param_target,
                                Tempi,
                                wd,
                               subtrace,
                               colConfigSets) {
  
  Na = 6.022e20
  c = 6.022e08
  pack = 1*(Na*(1/c))
  fc = 1e-06
  
  coloTag = data.frame(Ablated = "darkred", ParAblated = "#EC6D67", Unified = "#5351A2")
  nameTag = data.frame(Ablated = "Ablated", ParAblated = "Partially Ablated", Unified = "Unified")
  
  aa_places <- c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  places <- c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", aa_places)
  units = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))
  
  trace <- do.call(rbind, lapply(tag, function(j) {
    trace = readRDS(file = paste0(wd, paste0("/results/CDiff", "_", j, "_", Condition, "_", Exper),"/subtrace_" , j, Condition, ".RDs"))
    cbind(trace, Scenario = rep(j, length(trace$Time)))
  }))
  
  # subflux <- do.call(rbind, lapply(tag, function(j) {
  #   subflux = readRDS(file = paste0(wd, paste0("/results/CDiff", "_", j, "_", Condition, "_", Exper),"/subflux_" , j, Condition, ".RDs"))
  #   cbind(subflux, Scenario = rep(j, length(subflux$Time)))
  # }))
  # 
  
  lineCon = data.frame(Therapy = "solid", NoDrug = "dashed")
  
  # Sample data for demonstration
  df_for_df <- data.frame(
    variable = c("CD", "IECs", "trp_L_e", "leu_L_e", "pheme_c", "sink_pheme_c", "pheme_e"),
    unit = c(1, 1, (1 / (1 * (6.022e+20 * (1 / 6.022e+08)))) * 1e+03, 
             (1 / (1 * (6.022e+20 * (1 / 6.022e+08)))) * 1e+03, 1e-03, 1, 1e-03),
    y_label = c("C. difficile (cell)", 
                "IECs (cell)", 
                "Tryptophan (mmol/mL)", 
                "Leucine (mmol/mL)", 
                "Heme (µmol)",
                "heme flux (mmol/gDW*h)",
                "Heme (µmol)"
    ),
    subtitle = c("C. difficile cells", 
                 "Intestinal epithelial cells", 
                 "Extracellular Tryptophan concentration", 
                 "Extracellular Leucine concentration", 
                 "Heme intracellular concentration",
                 "Heme intracellular metabolism",
                 "Heme extracellular concentration"),
    scale_x_min = c(0, 3e+05, 0, 0, 0, -0.0015, NA),
    scale_x_max = c(8e+08, 4.75e+05, 0.02, 0.15, 16, 0, NA),
    scale_y_min = c(0, NA, 0, 0, 0, -0.0013, NA),
    scale_y_max = c(1.1e+09, NA, 0.025, 0.2, 50, 0, NA),
    colo_l = c("#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff"),
    colo_m = c("#d09565ff", "#ff80e3ff", "#d09565ff", "#d09565ff", "#d09565ff", "#89ca66b6",  "#d09565ff"),
    colo_h = c("gold", "gold", "#ff8978ff", "#ff8978ff", "gold", "gold", "#ff8978ff"),
    trace = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  )
  
  # Create plotting function
  # variable = "pheme_c"
  # j = "Unified"
  generate_plot <- function(variable, j) {
    
    row <- df_for_df[df_for_df$variable == variable, ]
    
    if (row$trace) {
      df <- trace %>%
        dplyr::filter(Places == variable) %>%
        dplyr::select(-c("Death4Treat", "Detox")) %>%
        dplyr::mutate(Time = as.numeric(Time))
    } else {
      df <- subflux %>%
        dplyr::filter(Places == variable) %>%
        dplyr::select(-c("Death4Treat", "Detox")) %>%
        dplyr::mutate(Time = as.numeric(Time))
    }
    
    p1 <- dplyr::filter(df, Scenario == j) %>%
      dplyr::filter(Time %in% c(0, 4, 8, 12, 16, 20)) %>%
      dplyr::mutate(Time = as.character(Time)) %>%
      ggplot(aes(x = Marking*row$unit, y = as.numeric(Time), group = Time)) +
      theme_minimal() +
      labs(y = "Time (h)", 
           x = row$y_label,
           #subtitle = row$subtitle,
           #title = paste(nameTag[, j], " Experiments", sep = "")
           ) +
      theme(plot.title = element_text(size = 10, face = "bold", color = coloTag[, j]),
            plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 9, color = "black"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.key.size = unit(0.4, "cm"),
            legend.position = "bottom") + 
      geom_jitter(aes(col = IECsDeath), alpha = 0.7, size = 1) +
      geom_violin(scale = "width", alpha = 0.5, width = 2.5) +
      scale_x_continuous(limits = c(row$scale_x_min, row$scale_x_max)) +
      scale_colour_gradientn(colors = c(row$colo_l, row$colo_m, row$colo_h), oob = scales::squish)+
      ggnewscale::new_scale_color()+
      geom_point(data = subtrace %>%
                   dplyr::filter(Time %in% c(0, 4, 8, 12, 16, 20),
                                 Places == variable, tag == j,Places == variable, new_eps_value == "1e-6"),
                 aes(y = Time,x = Marking*row$unit, group = Time,col = ConfParams ),size = 3) +
      scale_color_manual(values = colConfigSets[unique(subtrace$ConfParams)])
      
    p2 <- dplyr::filter(df, Scenario == j) %>%
      ggplot(aes(x = Time, y = Marking*row$unit, color = IECsDeath)) + 
      theme_minimal() +
      labs(x = "Time (h)", 
           y = row$y_label,
           #subtitle = row$subtitle,
           #title = paste(nameTag[, j], " Experiments", sep = "")
           ) + 
      theme(plot.title = element_text(size = 10, face = "bold", color = coloTag[, j]),
            plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 9, color = "black"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.key.size = unit(0.4, "cm"),
            legend.position = "bottom") + 
      geom_line(aes(group = config), alpha = 0.65) +
      scale_colour_gradientn(colors = c(row$colo_l, row$colo_m, row$colo_h), oob = scales::squish)
      #scale_y_continuous(limits = c(row$scale_y_min, row$scale_y_max))
    
    #combined_plot <- (p1 + p2 + plot_layout(ncol = 2))
    return(list(p1,p2))
  }
  
  pl2C=(
    (
    generate_plot("IECs", "Ablated")[[1]]+labs(title = "Ablated") +
      generate_plot("pheme_c", "Ablated")[[1]]
  )|
    (
      generate_plot("IECs", "Unified")[[1]]+labs(title = "Unified") +
        generate_plot("pheme_c", "Unified")[[1]]
    ) 
  )+ plot_layout(guides = "collect")&theme(legend.position = "bottom")
  
  
  ####### fig2B
  # Define your col_settings
  col_settings = c("darkred", "#5351A2")
  
  # Define the plotting function
  generate_plots <- function(plot_data,trace,subtrace) {
  
    variable = plot_data$variable
    unit = plot_data$unit
    y_label = plot_data$y_label
    subtitle = plot_data$subtitle
    scale_y_min = plot_data$scale_y_min
    scale_y_max = plot_data$scale_y_max
    
    df <- trace %>%
      dplyr::filter(Places == variable) %>%
      dplyr::select(-c("Death4Treat", "Detox")) %>%
      dplyr::mutate(Time = as.numeric(Time)) %>%
      dplyr::filter(Scenario != "ParAblated")
    dfsub = subtrace%>%
      dplyr::filter(Places == variable) %>%
      dplyr::mutate(Time = as.numeric(Time)) %>%
      dplyr::filter(tag != "ParAblated", new_eps_value == "1e-6")
    
    p1 <- ggplot(df) +
      geom_line( aes(x = Time, y = Marking * unit, color = Scenario, group = interaction(config, Scenario)), alpha = 0.08) +
      labs(y = y_label, x = "Time (h)")+#, subtitle = subtitle, title = "Comparison of Experimental Settings") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot",
            axis.text = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, face = "bold")) +
      scale_colour_manual(values = col_settings) 
      #scale_y_continuous(limits = c(scale_y_min, scale_y_max))
    
    
    p2 <- ggplot(df) +
      geom_line( aes(x = Time, y = Marking * unit, group = config), alpha = 0.005) +
      geom_line(data = dfsub, aes(x = Time, y = Marking * unit,group = interaction(config, tag), color = ConfParams, linetype = tag),linewidth = 1) +
      scale_color_manual(values = colConfigSets[unique(dfsub$ConfParams)])+
      #stat_summary(aes(group = Scenario), fun = median, geom = "line", linetype = "solid", size = 0.5) +
      labs(y = y_label, x = "Time (h)",color = "Scenario",linetype = "Paramters")+ #, subtitle = subtitle, title = "Comparison of Experimental Settings") + 
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, face = "bold")) 
    
    p3 <- ggplot(df, aes(x = Time, y = Marking * unit, color = Scenario, linetype = Scenario)) +
      geom_line(aes(group = interaction(config, Scenario)), alpha = 0.0055) + 
      stat_summary(aes(group = Scenario), fun = median, geom = "line", linetype = "solid", size = 0.5) +
      labs(y = y_label, x = "Time (h)", subtitle = subtitle, title = "Comparison of Experimental Settings") + 
      theme_minimal() +
      theme(legend.position = "right",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, face = "bold")) +
      scale_colour_manual(values = col_settings) 
      #scale_y_continuous(limits = c(scale_y_min, scale_y_max))
    
    return(list(p1 = p1, p2 = p2, p3 = p3))
  }
  
  # Define the variables you want to plot
  variables_to_plot <- c("CD", "IECs", "pheme_c")
  
  # Initialize lists to store the plots
  p1_list <- list()
  p2_list <- list()
  p3_list <- list()
  
  # Generate the plots for each variable and store them in the lists
  for (var in variables_to_plot) {
    plot_data <- df_for_df %>%
      dplyr::filter(variable == var)
    
    plots <- generate_plots(plot_data,trace,subtrace)
    
    p1_list[[var]] <- plots$p1
    p2_list[[var]] <- plots$p2
    p3_list[[var]] <- plots$p3
  }
  
  # Combine the plots into a single plot using patchwork
  pl2B = (p1_list$CD + p2_list$CD +
                     p1_list$IECs + p2_list$IECs + 
                     p1_list$pheme_c + p2_list$pheme_c  + 
                     plot_layout(ncol = 2))
  
  
  return(list(pl2B=pl2B,pl2C=pl2C))
}
