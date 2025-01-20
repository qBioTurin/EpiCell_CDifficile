
#' Configuration Analysis for Close-to-Median Trajectories
#' 
#' @return List containing:
#'   - plot: ggplot visualization of key metrics
#'   - Config: Data frame of configuration parameters
#'
configuration.closeMedian <- function(tag, 
                                      Condition, 
                                      subtrace = NULL, 
                                      numberOfSets = 2,
                                      places2plot = c("CD", "IECs", "BiomassCD", "pheme_c")) {
  
  # Constants
  CONSTANTS <- list(
    AVOGADRO = 6.022e20,
    CONVERSION = 6.022e08,
    PACKING = 1,
    FACTOR = 1e-06
  )
  
  EXPERIMENT <- "Model_Sensitivity"
  
  MEASUREMENTS <- list(
    amino_acids = c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e"),
    units = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))
  )
  
  places <- c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", 
              MEASUREMENTS$amino_acids)
  
  # Load data if not provided
  if (is.null(subtrace)) {
    subtrace <- load_subtrace_data(tag, Condition, EXPERIMENT, wd)
  }
  
  # Process data into groups
  grouped_data <- create_grouped_data(subtrace,numberOfSets)
  
  # Transform data for visualization
  processed_data <- transform_data_for_plotting(grouped_data)
  
  # Create visualization
  plot_output <- create_visualization(processed_data, places2plot)
  
  # Generate configuration summary
  config_summary <- summarize_configurations(processed_data, subtrace)
  
  return(list(plot = plot_output, Config = config_summary))
}

#' Helper function to load subtrace data
load_subtrace_data <- function(tag, Condition, Experiment, wd) {
  do.call(rbind, lapply(tag, function(j) {
    file_path <- file.path(wd, "results", 
                           paste0("CDiff_", j, "_", Condition, "_", Experiment),
                           paste0("subtrace_", j, Condition, ".rds"))
    subtrace <- readRDS(file_path)
    cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
  }))
}

#' Create grouped data based on IECs measurements
create_grouped_data <- function(subtrace,numberOfSets) {
  subtrace %>%
    filter(Places == "IECs", round(Time,digits = 1) == 20) %>%
    group_by(Scenario) %>%
    mutate(group = ntile(Marking, numberOfSets)) %>%
    group_by(Scenario, group) %>%
    mutate(median_value = median(Marking)) %>%
    slice_min(abs(Marking - median_value), with_ties = FALSE)
}

#' Transform data for visualization
transform_data_for_plotting <- function(grouped_data) {
  do.call(rbind, lapply(1:nrow(grouped_data), function(j) {
    config_file <- grouped_data$config[j]
    df <- readr::read_delim(config_file) %>%
      mutate(
        Config = config_file,
        ConfigNumb = gsub(".*-([0-9]+)\\.trace", "\\1", basename(config_file)),
        Scenario = grouped_data$Scenario[j],
        ConfParams = grouped_data$group[j]
      )
  })) %>%
    tidyr::gather(-ConfParams, -ConfigNumb, -Config, -Time, -Scenario,
                  key = "Places", value = "Marking")
}

#' Create visualization plot
create_visualization <- function(processed_data, places2plot) {
  processed_data %>%
    filter(Places %in% places2plot) %>%
    ggplot() +
    geom_line(aes(x = Time, y = Marking, 
                  group = ConfParams, 
                  linetype = Scenario, 
                  col = as.factor(ConfParams))) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_wrap(~Places, scales = "free")
}

#' Summarize configuration parameters
summarize_configurations <- function(processed_data, subtrace) {
  config_data <- processed_data %>%
    filter(Places %in% c("CD", "IECs", "BiomassCD")) %>%
    select(Scenario, Config, ConfParams) %>%
    distinct() %>%
    mutate(ConfParams = paste0("Set ", ConfParams))
  
  params_data <- subtrace %>%
    filter(config %in% unique(config_data$Config)) %>%
    select(Detox, Death4Treat, IECsDeath, Scenario, config) %>%
    distinct()
  
  merge(params_data, config_data, 
        by.x = c("Scenario", "config"),
        by.y = c("Scenario", "Config")) %>%
    select(-config)
}
