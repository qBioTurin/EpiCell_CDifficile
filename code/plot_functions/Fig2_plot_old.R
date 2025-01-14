Exper = "Model_Sensitivity"
Condition = "Therapy"
tag = c("Ablated", "ParAblated", "Unified")
param_target = "IECsDeath"
Tempi = c(0, 12, 24, 36, 48, 60)
colo1 = c("black", "magenta", "gold")
colo2 = c("black", "#266867", "yellow")
colo3 = c("#ffd166", "#ee6c4d", "#293241")

plotting_Fig2_paper = function(Exper,
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
  
  coloTag = data.frame(Ablated = "darkred", ParAblated = "#EC6D67", Unified = "#5351A2")
  nameTag = data.frame(Ablated = "Ablated", ParAblated = "Partially Ablated", Unified = "Unified")
  
  aa_places <- c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  places <- c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", aa_places)
  units = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))
  
  subtrace <- do.call(rbind, lapply(tag, function(j) {
    subtrace = readRDS(file = paste0(wd, paste0("/results/CDiff", "_", j, "_", Condition, "_", Exper),"/subtrace_" , j, Condition, ".rds"))
    cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
  }))
  
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
    scale_x_min = c(0, 0, 0, 0, 0, -0.0015, NA),
    scale_x_max = c(8e+08, 4.75e+05, 0.02, 0.15, 50, 0, NA),
    scale_y_min = c(0, NA, 0, 0, 0, -0.0013, NA),
    scale_y_max = c(1.1e+09, NA, 0.15, 1, 150, 0, NA),
    colo_l = c("#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff"),
    colo_m = c("#d09565ff", "#ff80e3ff", "#d09565ff", "#d09565ff", "#d09565ff", "#89ca66b6",  "#d09565ff"),
    colo_h = c("gold", "gold", "#ff8978ff", "#ff8978ff", "gold", "gold", "#ff8978ff"),
    trace = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  )
  
  # Create plotting function
  generate_plot <- function(variable, j) {
    
    row <- df_for_df[df_for_df$variable == variable, ]
    
    if (row$trace) {
      df <- subtrace %>%
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
      labs(y = "Time (h)", x = row$y_label) +
      theme(plot.title = element_text(size = 10, face = "bold", color = coloTag[, j]),
            plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 9, color = "black"),
            axis.title = element_text(size = 12, face = "bold"),
            legend.key.size = unit(0.4, "cm"),
            legend.position = "bottom") + 
      geom_jitter(aes(col = IECsDeath), alpha = 0.7, size = 1) +
      geom_violin(scale = "width", alpha = 0.5, width = 2.5) +
      scale_colour_gradientn(colors = c(row$colo_l, row$colo_m, row$colo_h), oob = scales::squish)
    
    p2 <- dplyr::filter(df, Scenario == j) %>%
      ggplot(aes(x = Time, y = Marking*row$unit, color = IECsDeath)) + 
      theme_minimal() +
      labs(x = "Time (h)", y = row$y_label) + 
      theme(plot.title = element_text(size = 10, face = "bold", color = coloTag[, j]),
            plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 9, color = "black"),
            axis.title = element_text(size = 12, face = "bold"),
            legend.key.size = unit(0.4, "cm"),
            legend.position = "bottom") + 
      geom_line(aes(group = config), alpha = 0.65) +
      scale_colour_gradientn(colors = c(row$colo_l, row$colo_m, row$colo_h), oob = scales::squish)
    
    return(list(p1,p2))
  }
  
  pl2C=(
    (
    generate_plot("IECs", "Ablated")[[1]]+labs(title = "Ablated") +
    generate_plot("leu_L_e", "Ablated")[[2]]+
      generate_plot("trp_L_e", "Ablated")[[2]]+
      generate_plot("pheme_c", "Ablated")[[1]]  + plot_layout(ncol = 4)
  ) |
    (
      generate_plot("IECs", "ParAblated")[[1]]+labs(title = "Partially Ablated") +
        generate_plot("leu_L_e", "ParAblated")[[2]]+
        generate_plot("trp_L_e", "ParAblated")[[2]]+
        generate_plot("pheme_c", "ParAblated")[[1]]  + plot_layout(ncol = 4)
    )| 
    (
      generate_plot("IECs", "Unified")[[1]]+labs(title = "Unified") +
        generate_plot("leu_L_e", "Unified")[[2]]+
        generate_plot("trp_L_e", "Unified")[[2]]+
        generate_plot("pheme_c", "Unified")[[1]]  + plot_layout(ncol = 4)
    ) 
  )+ plot_layout(nrow = 3,guides = "collect")&theme(legend.position = "bottom")
  
  col_settings = c("darkred", "#5351A2")
  
  generate_plots <- function(variable, unit, y_label, subtitle, scale_y_min, scale_y_max) {
    
    df <- subtrace %>%
      dplyr::filter(Places == variable) %>%
      dplyr::select(-c("Death4Treat", "Detox")) %>%
      dplyr::mutate(Time = as.numeric(Time)) %>%
      dplyr::filter(Scenario != "ParAblated")
    
    p1 <- ggplot(df, aes(x = Time, y = Marking * unit, color = Scenario)) +
      geom_line(aes(group = interaction(config, Scenario)), alpha = 0.08) +
      labs(y = y_label, x = "Time (h)", subtitle = subtitle, title = "Comparison of Experimental Settings") +
      theme_minimal() +
      theme(legend.position = "right",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot",
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, face = "bold")) +
      scale_colour_manual(values = col_settings)
    
    p2 <- ggplot(df, aes(x = Time, y = Marking * unit, color = Scenario, linetype = Scenario)) +
      geom_line(aes(group = interaction(config, Scenario)), alpha = 0.0055) + 
      stat_summary(aes(group = Scenario), fun = median, geom = "line", linetype = "solid", size = 0.5) +
      labs(y = y_label, x = "Time (h)", subtitle = subtitle, title = "Comparison of Experimental Settings") + 
      theme_minimal() +
      theme(legend.position = "right",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, face = "bold")) +
      scale_colour_manual(values = col_settings)
    
    p3 <- ggplot(df, aes(x = Time, y = Marking * unit, group = Scenario, colour = Scenario)) +
      stat_summary(fun = median, geom = "line") +
      labs(y = y_label, x = "Time (h)", subtitle = subtitle, title = "Comparison of Experimental Settings") + 
      theme_minimal() +
      theme(legend.position = "right",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, face = "bold")) +
      scale_colour_manual(values = col_settings)
    
    pp3 <- ggplot_build(p3)
    
    x <- pp3$data[[1]]$x[pp3$data[[1]]$group == 1]
    y1 <- pp3$data[[1]]$y[pp3$data[[1]]$group == 1]
    y2 <- pp3$data[[1]]$y[pp3$data[[1]]$group == 2]
    
    df2 <- data.frame(x = x, ymin = pmin(y1, y2), ymax = pmax(y1, y2), side = (y1 < y2), ydiff = y2 - y1)
    
    df_median <- df %>%
      group_by(Time, Scenario) %>%
      summarize(median_value = median(Marking))
    
    g2 <- ggplot(df2) +
      stat_summary(data = df, aes(x = Time, y = Marking * unit, group = Scenario, colour = Scenario), 
                   fun = median, geom = "line", linetype = "solid", size = 0.5) + 
      geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = side), alpha = 0.1) +
      guides(alpha = FALSE, fill = FALSE) +
      theme_minimal() +
      labs(y = y_label, x = "Time (h)", subtitle = subtitle, title = "Comparison of Experimental Settings") + 
      theme(legend.position = "right",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, face = "bold")) +
      scale_colour_manual(values = col_settings) +
      scale_fill_manual(values = col_settings)
    
    g3 <- ggplot(df2) +
      geom_line(aes(x = x, y = abs(ydiff), colour = side)) +
      geom_area(aes(x = x, y = abs(ydiff), fill = side), alpha = 0.1) +
      guides(alpha = FALSE, fill = FALSE) +
      theme_minimal() +
      labs(y = "Absolute Difference", colour = "Difference", x = "Time (h)", 
           subtitle = subtitle, title = "Comparison of Experimental Settings") + 
      theme(legend.position = "right",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, face = "bold")) +
      scale_colour_manual(values = col_settings) +
      scale_fill_manual(values = col_settings)
    
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
    
    plots <- generate_plots(
      variable = plot_data$variable,
      unit = plot_data$unit,
      y_label = plot_data$y_label,
      subtitle = plot_data$subtitle,
      scale_y_min = plot_data$scale_y_min,
      scale_y_max = plot_data$scale_y_max
    )
    
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
