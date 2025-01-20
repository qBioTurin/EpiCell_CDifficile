library(ggnewscale)

# Exper = "Model_Sensitivity"
# Condition = "Therapy"
# tag = c("Ablated", "ParAblated", "Unified")
# param_target = "IECsDeath"
# Tempi = c(10, 20, 30, 40, 50, 60)
# wd = wd
# trajectories = trajectories
# colConfigSets = colors_new_confParams
# variables_to_plot = c("CD", "IECs", "pheme_c", "leu_L_e", "trp_L_e")

plotting_Fig2_paper = function(Exper,
                               Condition,
                               tag,
                               param_target,
                               Tempi,
                               wd,
                               trajectories,
                               colConfigSets,
                               variables_to_plot) {
  
  Na = 6.022e20
  c = 6.022e08
  pack = 1*(Na*(1/c))
  fc = 1e-06
  
  coloTag = data.frame(Ablated = "darkred", ParAblated = "#EC6D67", Unified = "#5351A2")
  nameTag = data.frame(Ablated = "Ablated", ParAblated = "Partially Ablated", Unified = "Unified")
  
  aa_places <- c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  places <- c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", aa_places)
  units = c(rep("(cell)", 2), "(pg)", rep("(mmol)", 9))
  
  subtrace <- do.call(rbind, lapply(tag, function(j) {
    subtrace = readRDS(file = paste0(wd, paste0("/results/CDiff", "_", j, "_", Condition, "_", Exper),"/subtrace_" , j, Condition, ".rds"))
    cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
  }))
  
  df_for_df <- data.frame(
    variable = c("CD", "IECs", "trp_L_e", "leu_L_e", "pheme_c", "sink_pheme_c", "pheme_e"),
    scale_x_min = c(NA, NA, NA, NA, NA, NA, NA),
    scale_x_max = c(1e09, NA, NA, NA, 700, NA, NA),
    scale_y_min = c(0, 0, NA, NA, 0, NA, 0),
    scale_y_max = c(1e09, NA, NA, NA, 700, NA, NA),
    unit = c(1, 1, (1 / (1 * (6.022e+20 * (1 / 6.022e+08)))) * 1e+03,
             (1 / (1 * (6.022e+20 * (1 / 6.022e+08)))) * 1e+03, 1e-03, 1, 1e-03),
    y_label = c("log(C. difficile) (cell)",
                "log(IECs) (cell)",
                "log(Tryptophan) (mmol/mL)",
                "log(Leucine) (mmol/mL)",
                "log(Heme) (µmol)",
                "log(heme flux) (mmol/gDW*h)",
                "log(Heme) (µmol)"),
    subtitle = c("C. difficile cells",
                 "Intestinal epithelial cells",
                 "Extracellular Tryptophan concentration",
                 "Extracellular Leucine concentration",
                 "Heme intracellular concentration",
                 "Heme intracellular metabolism",
                 "Heme extracellular concentration"),
    colo_l = c("#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff", "#373332ff"),
    colo_m = c("#89ca66b6", "#ff80e3ff", "#d09565ff", "#d09565ff", "#d09565ff", "#89ca66b6",  "#d09565ff"),
    colo_h = c("gold", "gold", "#ff8978ff", "#ff8978ff", "gold", "gold", "#ff8978ff"),
    trace = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  )
  
  variable = "CD"
  j = "Ablated"
  
  generate_plot <- function(variable, j, coloTag) {
    
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
    
    scientific_10 <- function(x) {
      text = gsub("e", " %*% 10^", scales::scientific_format()(x))
      text[1] = "0"
      parse(text=text)
    }
    
    # Create scale_x based on variable
    x_scale <- if(variable == "CD") {
      scale_x_continuous(limits = c(row$scale_x_min, row$scale_x_max), 
                         label=scientific_10, breaks = (seq(0, row$scale_x_max, length.out = 3)))
    } else {
      scale_x_continuous(limits = c(row$scale_x_min, row$scale_x_max))
    }
    
    p1 <- dplyr::filter(df, Scenario == j) %>%
      dplyr::filter(Time %in% Tempi) %>%
      dplyr::mutate(Time = as.character(Time)) %>%
      ggplot(aes(x = Marking*row$unit, y = as.numeric(Time), group = Time)) +
      theme_minimal() +
      labs(y = "Time (h)", x = gsub("log\\(([^)]+)\\)\\s*\\(([^)]+)\\)", "\\1 (\\2)", row$y_label)) +
      theme(plot.title = element_text(size = 15, face = "bold", color = coloTag[, j]),
            plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.key.size = unit(0.4, "cm"),
            legend.position = "bottom") + 
      geom_jitter(aes(col = IECsDeath), alpha = 0.7, size = 1) +
      geom_violin(scale = "width", alpha = 0.5, width = 2.5) +
      scale_colour_gradientn(colors = c(row$colo_l, row$colo_m, row$colo_h), oob = scales::squish) +
      ggnewscale::new_scale_color() +
      geom_point(data = trajectories %>% 
                   dplyr::filter(Time %in% Tempi,
                                 Places == variable, 
                                 tag == j, 
                                 new_eps_value == "1e-6"),
                 aes(y = Time, x = Marking*row$unit, group = Time, col = ConfParams), size = 3) +
      scale_color_manual(values = colConfigSets[unique(trajectories$ConfParams)]) +
      x_scale
    
    p2 <- dplyr::filter(df, Scenario == j) %>%
      ggplot() + 
      geom_line(aes(x = Time, y = Marking*row$unit, color = IECsDeath, group = config), alpha = 0.1) +
      scale_colour_gradientn(colors = c(row$colo_l, row$colo_m, row$colo_h), oob = scales::squish) +
      new_scale_color() +
      geom_line(data = trajectories %>% 
                  dplyr::filter(Places == variable, 
                                tag == j,
                                new_eps_value == "1e-6") %>%
                  distinct(Time, ConfParams, tag, .keep_all = TRUE),
                aes(x = Time, y = Marking*row$unit,
                    group = interaction(ConfParams, tag),
                    color = ConfParams,
                    linetype = tag),
                linewidth = 1) +
      scale_color_manual(values = colConfigSets, name = "Sets") +
      scale_linetype_manual(values = c("Unified" = "solid", "Ablated" = "dashed"), name = "Scenario") +
      theme_minimal() +
      labs(x = "Time (h)", 
           y = gsub("log\\(([^)]+)\\)\\s*\\(([^)]+)\\)", "\\1 (\\2)", row$y_label),
           color = "Configuration", 
           linetype = "Scenario") + 
      theme(plot.title = element_text(size = 15, face = "bold", color = coloTag[, j]),
            plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.key.size = unit(0.4, "cm"),
            legend.position = "bottom")
    
    return(list(p1, p2))
  }
  
  pl2C = (
    ((generate_plot("CD", "Ablated", coloTag)[[1]] + labs(title = "Ablated")) |
       generate_plot("leu_L_e", "Ablated", coloTag)[[2]] | 
       generate_plot("trp_L_e", "Ablated", coloTag)[[2]] |
       generate_plot("pheme_c", "Ablated", coloTag)[[1]]
     ) /
      ((generate_plot("CD", "ParAblated", coloTag)[[1]] + labs(title = "ParAblated")) |
         generate_plot("leu_L_e", "ParAblated", coloTag)[[2]] | 
         generate_plot("trp_L_e", "ParAblated", coloTag)[[2]] |
         generate_plot("pheme_c", "ParAblated", coloTag)[[1]]
      ) /
      ((generate_plot("CD", "Unified", coloTag)[[1]] + labs(title = "Unified")) |
         generate_plot("leu_L_e", "Unified", coloTag)[[2]] | 
         generate_plot("trp_L_e", "Unified", coloTag)[[2]] |
         generate_plot("pheme_c", "Unified", coloTag)[[1]])
    ) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  # Define your col_settings
  col_settings = c(coloTag["Ablated"][1, ], coloTag["Unified"][1, ])
  
  # plot_data <- df_for_df %>% dplyr::filter(variable == "CD")
  
  generate_plots <- function(plot_data, trajectories, subtrace) {
    
    variable = plot_data$variable
    unit = plot_data$unit 
    y_label = plot_data$y_label
    # subtitle = plot_data$subtitle
    scale_y_min = plot_data$scale_y_min
    scale_y_max = plot_data$scale_y_max
    
    scientific_10 <- function(x) {
      text = gsub("e", " %*% 10^", scales::scientific_format()(x))
      text[1] = "0"
      parse(text=text)
    }
    
    # Create scale_x based on variable
    y_scale <- if(variable == "CD") {
      scale_y_continuous(label = scientific_10, breaks = seq(scale_y_min, scale_y_max, length.out = 6))
    } else {
      NULL
    }
    
    df <- subtrace %>%
      dplyr::filter(Places == variable) %>%
      dplyr::select(-c("Death4Treat", "Detox")) %>%
      dplyr::mutate(Time = as.numeric(Time)) %>%
      dplyr::filter(Scenario != "ParAblated")
    
    dfsub = trajectories %>%
      dplyr::filter(Places == variable) %>%
      dplyr::mutate(Time = as.numeric(Time)) %>%
      dplyr::filter(tag != "ParAblated", new_eps_value == "1e-6")
    
    transform_y <- function(marking, unit) {
      if(variable == "CD") {
        return(marking * unit)
      } else {
        return(log(abs(marking * unit) + 1))
      }
    }
    
    # Modify y_label based on variable
    if(variable == "CD") {
      y_label = gsub("log\\(([^)]+)\\)\\s*\\(([^)]+)\\)", "\\1 (\\2)", y_label)
    }
    
    # Base plot p1
    p1 <- ggplot(df) +
      geom_line(
        aes(x = Time, y = transform_y(Marking, unit), 
            color = Scenario, group = interaction(config, Scenario)), alpha = 0.05) +
      labs(y = y_label, x = "Time (h)") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot",
            axis.text = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, face = "bold")) +
      scale_colour_manual(values = col_settings) + y_scale
    
    if(variable == "IECs") {
      p1 <- p1 + geom_vline(xintercept = 20, linetype = "dotted", color = "darkgrey")
    }
    
    # Base plot p2
    p2 <- ggplot(df) +
      geom_line(aes(x = Time, y = transform_y(Marking, unit), 
                    group = config, color = Scenario), alpha = 0.005) +
      scale_color_manual(values = c("Ablated" = "darkred", "Unified" = "#5351A2"), name = "Scenario") +
      new_scale_color() +
      geom_line(data = dfsub |> distinct(Time, ConfParams, tag, .keep_all = TRUE), 
                aes(x = Time, y = transform_y(Marking, unit),
                    group = interaction(ConfParams, tag), 
                    color = ConfParams, linetype = tag), linewidth = 0.75) +
      scale_color_manual(values = colConfigSets, name = "Sets") +
      scale_linetype_manual(values = c("Unified" = "solid", "Ablated" = "dashed"), name = "Scenario") +
      labs(y = y_label, x = "Time (h)", color = "Configuration", linetype = "Scenario") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, face = "bold")) + y_scale
    
    # Base plot p3
    p3 <- ggplot(df, aes(x = Time, y = transform_y(Marking, unit), 
                         color = Scenario, linetype = Scenario)) +
      geom_line(aes(group = interaction(config, Scenario)), alpha = 0.01) + 
      stat_summary(aes(group = Scenario), fun = median, geom = "line", linetype = "solid", linewidth = 0.75) +
      labs(y = y_label, x = "Time (h)", 
           # subtitle = subtitle, 
           # title = "Median trajectories in each setting"
           ) + 
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold", color = "black"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#2a475e"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, face = "bold")) +
      scale_colour_manual(values = col_settings) + y_scale
    
    if(variable == "IECs") {
      p2 <- p2 + geom_vline(xintercept = 20, linetype = "solid", color = "black", size = 0.35)
    }
    
    return(list(p1 = p1, p2 = p2, p3 = p3))
  }
  
  # Initialize lists to store the plots
  p1_list <- list()
  p2_list <- list()
  p3_list <- list()
  
  # Generate the plots for each variable and store them in the lists
  for (var in variables_to_plot) {
    plot_data <- df_for_df %>%
      dplyr::filter(variable == var)
    
    plots <- generate_plots(plot_data, trajectories, subtrace)
    
    p1_list[[var]] <- plots$p1
    p2_list[[var]] <- plots$p2
    p3_list[[var]] <- plots$p3
  }
  
  # Combine the plots into a single plot using patchwork
  pl2B = (p3_list$CD + p2_list$CD +
            p3_list$IECs + p2_list$IECs + 
            p3_list$pheme_c + p2_list$pheme_c + 
            plot_layout(ncol = 2))
  
  return(list(pl2B = pl2B, pl2C = pl2C))
}
