
library(dplyr)
library(ggplot2)
library(patchwork)

# Define your col_settings
col_settings = c("darkred", "#5351A2")

# Create your df_for_df data frame
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
  scale_y_min = c(0, 0, 0, 0, 0, 0, 0),
  scale_y_max = c(1e+09, 5e+05, 0.017, 0.14, 48, NA, NA)
)

# Define the plotting function
generate_plots <- function(variable, unit, y_label, subtitle, scale_y_min, scale_y_max) {
  
  df <- subtrace %>%
    dplyr::filter(Places == variable) %>%
    dplyr::select(-c("Death4Treat", "Detox")) %>%
    dplyr::mutate(Time = as.numeric(Time)) %>%
    dplyr::filter(Scenario != "Ablated")
  
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
    scale_colour_manual(values = col_settings) +
    scale_y_continuous(limits = c(scale_y_min, scale_y_max))
  
  
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
    scale_colour_manual(values = col_settings) +
    scale_y_continuous(limits = c(scale_y_min, scale_y_max))
  
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
    scale_colour_manual(values = col_settings) +
    scale_y_continuous(limits = c(scale_y_min, scale_y_max))
  
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
    scale_fill_manual(values = col_settings) +
    scale_y_continuous(limits = c(scale_y_min, scale_y_max))
  
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
    scale_fill_manual(values = col_settings) +
    scale_y_continuous(limits = c(scale_y_min, scale_y_max))
  
  # combined_plot <- (p1 + p2 + p3 + g2 + g3 + plot_layout(ncol = 5))
  # 
  # ggsave(combined_plot, 
  #        file = paste0(wd, "/fig_3/", variable, "_diff.pdf"), 
  #        width = 22.5, height = 3, limitsize = FALSE)
  
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
combined_plot = (p1_list$CD + p2_list$CD +
                   p1_list$IECs + p2_list$IECs + 
                   p1_list$pheme_c + p2_list$pheme_c  + 
                   plot_layout(ncol = 2))

ggsave(combined_plot, 
       file = paste0(wd, "/fig_3/fig_3_diff.pdf"), 
       width = 9, height = 9, limitsize = FALSE)

