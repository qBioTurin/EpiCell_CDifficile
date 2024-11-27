
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(grid)
library(patchwork)

# setting working directory
wd = "/home/raucello/EpiCell_CDifficile"
supp_function.dir = "/code/supplementary_functions/"

####

Exper = "Model_Sensitivity"
Condition = "Therapy"
tag = c("SuperAblated", "Ablated", "Unified")
param_target = c("IECsDeath", "Death4Treat", "Detox")

injection_times = c(0, 5, 10, 15, 20, 25)

coloTag = data.frame(SuperAblated = "darkred", Ablated = "#EC6D67", Unified = "#5351A2")
nameTag = data.frame(SuperAblated = "Ablated", Ablated = "Partially Ablated", Unified = "Unified")

# for (j in tag) {
#   
#   setwd(paste0(wd, "/results/CDiff", "_", j, "_", Condition, "_", Exper))
#   
#   load(file = paste0("subtrace_" , j, Condition, ".RData"))
#   subtrace = cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
#   
#   saveRDS(subtrace, file = paste0("subtrace_unfolded_", j, "_", Condition, ".rds"))
# }
# 
# for (j in tag) {
#   
#   setwd(paste0(wd, "/results/CDiff", "_", j, "_", Condition, "_", Exper))
#   
#   load(file = paste0("subflux_" , j, Condition, ".RData"))
#   subflux = cbind(subflux, Scenario = rep(j, length(subflux$Time)))
#   
#   saveRDS(subflux, file = paste0("subflux_unfolded_", j, "_", Condition, ".rds"))
# }

subtrace = rbind(
  readRDS(file = paste0(wd, "/results/CDiff", "_", "SuperAblated", "_", "Therapy", "_", Exper, "/subtrace_unfolded_", "SuperAblated", "_", Condition, ".rds")),
  readRDS(file = paste0(wd, "/results/CDiff", "_", "Ablated", "_", "Therapy", "_", Exper, "/subtrace_unfolded_", "Ablated", "_", Condition, ".rds")), 
  readRDS(file = paste0(wd, "/results/CDiff", "_", "Unified", "_", "Therapy", "_", Exper, "/subtrace_unfolded_", "Unified", "_", Condition, ".rds")))

subflux = rbind(
  readRDS(file = paste0(wd, "/results/CDiff", "_", "SuperAblated", "_", "Therapy", "_", Exper, "/subflux_unfolded_", "SuperAblated", "_", Condition, ".rds")),
  readRDS(file = paste0(wd, "/results/CDiff", "_", "Ablated", "_", "Therapy", "_", Exper, "/subflux_unfolded_", "Ablated", "_", Condition, ".rds")), 
  readRDS(file = paste0(wd, "/results/CDiff", "_", "Unified", "_", "Therapy", "_", Exper, "/subflux_unfolded_", "Unified", "_", Condition, ".rds")))

####

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
    labs(y = "Time (h)", 
         x = row$y_label,
         subtitle = row$subtitle,
         title = paste(nameTag[, j], " Experiments", sep = "")) +
    theme(plot.title = element_text(size = 10, face = "bold", color = coloTag[, j]),
          plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
          plot.title.position = "plot", 
          axis.text = element_text(size = 9, color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          legend.key.size = unit(0.4, "cm")) + 
    geom_jitter(aes(col = IECsDeath), alpha = 0.7, size = 1) +
    geom_violin(scale = "width", alpha = 0.5, width = 2.5) +
    scale_x_continuous(limits = c(row$scale_x_min, row$scale_x_max)) +
    scale_colour_gradientn(colors = c(row$colo_l, row$colo_m, row$colo_h), oob = scales::squish)
  
  p2 <- dplyr::filter(df, Scenario == j) %>%
    ggplot(aes(x = Time, y = Marking*row$unit, color = IECsDeath)) + 
    theme_minimal() +
    labs(x = "Time (h)", 
         y = row$y_label,
         subtitle = row$subtitle,
         title = paste(nameTag[, j], " Experiments", sep = "")) + 
    theme(plot.title = element_text(size = 10, face = "bold", color = coloTag[, j]),
          plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
          plot.title.position = "plot", 
          axis.text = element_text(size = 9, color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          legend.key.size = unit(0.4, "cm")) + 
    geom_line(aes(group = config), alpha = 0.65) +
    scale_colour_gradientn(colors = c(row$colo_l, row$colo_m, row$colo_h), oob = scales::squish) +
    scale_y_continuous(limits = c(row$scale_y_min, row$scale_y_max))
  
  combined_plot <- (p1 + p2 + plot_layout(ncol = 2))
  
  ggsave(combined_plot, file = paste0(wd, "/fig_3/", variable, "_", j, ".pdf"), width = 8, height = 3)
}

# Loop over variables and tags to generate plots
for (variable in df_for_df$variable) {
  for (j in tag) {
    generate_plot(variable, j)
  }
}

# Our research introduces a novel modeling paradigm integrating Ordinary Differential Equations (ODEs) and Flux Balance Analysis (FBA) into a unified meta-formalism. This approach combines the precision of mechanism-based models with the comprehensive scope of constraint-based ones. Our study encompasses hybrid model development, computational implementation, solution techniques, and parameter estimation.
# 
# The framework was utilized to analyze microbial growth variations concerning key metabolite concentrations. A unique feature of this paradigm is its ability to continuously compute system evolution, accounting for all metabolic reactions and environmental changes. Specifically, we investigated the interaction between *Clostridioides difficile* (C. difficile), intestinal epithelial cells (IECs), and gut metabolites during antibiotic treatment.
# 
# We aim to model a biological scenario where the composition of the gut growth medium significantly influences C. difficile toxin production. These toxins damage the epithelial cell barrier and cause IEC apoptosis. When nutrients are scarce, the bacteria produce toxins that release host-derived nutrients through cytotoxic activity. Upon IEC death, nutrients such as amino acids are released into the gut environment, becoming available for the bacteria. We also model the action of the antibiotic metronidazole, which induces oxidative stress and pathogen cell death. Without the drug, successful C. difficile colonization leads to increased bacterial cells, IEC death, gut lining damage, and increased heme concentration due to red blood cell lysis. 
# 
# Antibiotic injections reduce bacterial numbers at specific time intervals (e.g., 8, 16, 24 hours). Despite this reduction, dynamically computing the FBA throughout the simulation reveals mechanisms of antibiotic resistance. This mechanism, involving heme release from gut barrier leakage, promotes chronic C. difficile infection even when antibiotics are effective.
# 
# Three experimental settings were proposed to assess the model:
# 1. **Ablated experiment**: FBA computed only at the start.
# 2. **Partial-ablated experiment**: FBA computed during antibiotic injections.
# 3. **Unified experiment**: FBA continuously solved with ODEs.
# 
# Comparing these settings shows different system trajectories. In all settings, IEC cell numbers change linearly with the "IECsDeath" rate, reflecting toxic gene regulation due to low nutrients. The unified model shows a significant decrease in cell numbers due to increased cell death compared to the ablated settings. Tryptophan and leucine dynamics, linked to the bacterium's cytotoxicity ("IECsDeath"), differ across settings. Ablation experiments suggest infection clearance reduces IEC damage, while unified experiments reveal C. difficile entering nutrient-dependent stationary phases during treatment. Only the unified setting uncovers a heme-mediated antibiotic-resistant strain.
# 
# The aim is to visualize the results effectively. For example, considering the dynamics of Leucine, we can use a violin plot to show the distribution of dynamics based on the variation of the parameter "IECsDeath" at selected time points. Additionally, the same distribution can be visualized using a line plot. Below are some R plots illustrating these visualization attempts, where "j" represents the experimental settings:
# 
# The question is which could be another way to visualize these results, I wish a plot that merges all three different experimental settings onto a single graph to highlight the differences

# Define custom color gradients for each scenario
color_gradients <- list(
  "SuperAblated" = c("darkgrey", "#ffa498ff", "#7c0f00ff"),
  "Ablated" = c("darkgrey", "#ffdedaff", "#ffa498ff"),
  "Unified" = c("darkgrey", "lightblue", "#5351A2")
)

# Function to create a plot for a specific scenario
create_plot <- function(variable, j) {
  
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
  
  df %>%
    dplyr::filter(Time %in% c(0, 4, 8, 12, 16, 20, 24), Scenario == j) %>%
    dplyr::mutate(Time = as.character(Time)) %>%
    ggplot(aes(x = ordered(Time, levels = unique(Time)), 
               y = Marking * row$unit)) +
    labs(y = row$y_label, 
         x = "Time (h)",
         subtitle = row$subtitle,
         title = paste0(nameTag[, j], " Experiments")) +
    theme_minimal() +
    theme(
      text = element_text(size = 10, color = "black"),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "#2a475e"),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    ) +
    scale_y_continuous(limits = c(row$scale_x_min, row$scale_x_max)) +
    geom_jitter(aes(color = IECsDeath), alpha = 0.15, size = 0.4) +
    scale_color_gradientn(colors = if (j == "Unified") {
      color_gradients[["Unified"]]
    } else if (j == "Ablated") {
      color_gradients[["Ablated"]]
    } else {
      color_gradients[["SuperAblated"]]
    }) +
    geom_boxplot(alpha = 0.25, width = 0.7)
}

for (variable in df_for_df$variable){
  # Create plots for each scenario
  plot_super_ablated <- create_plot(variable, "SuperAblated")
  plot_ablated <- create_plot(variable, "Ablated")
  plot_unified <- create_plot(variable, "Unified")
  
  # Combine the plots using patchwork
  combined_plot <- plot_super_ablated + plot_ablated + plot_unified + plot_layout(ncol = 3)
  
  ggsave(combined_plot, file = paste0(wd, "/fig_3/", variable, "_all_scenario_fig_3.pdf"), 
         width = 10, height = 2.5, limitsize = FALSE)
}
