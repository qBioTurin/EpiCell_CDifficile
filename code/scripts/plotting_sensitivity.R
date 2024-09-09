
PlottingSensitivityAnalysis = function(numb, 
                                       alph, 
                                       multi_param,
                                       tag, 
                                       Condition, 
                                       Exper,
                                       net_fname,
                                       param_target,
                                       colo1,
                                       colo2,
                                       colo3,
                                       Na,
                                       c,
                                       fc) {
 
  library(doParallel)
  library(foreach)
  
  failed.tr = c()
  f.tr = 1
  
  if (multi_param) {
    setwd(paste0(wd, "/results/CDiff_", tag, "_", Condition, "_", Exper))
  } else {
    setwd(paste0(wd, "/results/CDiff_", tag, "_", Condition, "_", Exper))
  }
  
  load(paste0(net_fname, "-analysis.RData"))
  
  # Function to find the index of a target parameter in a configuration list
  find_param_index = function(param_target, config, multi_param = FALSE) {
    # Initialize the id_param variable
    id_param <- ifelse(multi_param, vector(mode = "list", length(param_target)), NA)
    # Loop over the configuration list
    for (i in seq_along(config)) {
      # Check if the target parameter(s) match the current configuration parameter
      if (multi_param) {
        id_param[sapply(param_target, function(p) p == config[[i]][[1]][[1]])] <- i
      } else if (param_target == config[[i]][[1]][[1]]) {
        id_param <- i
        break
      }
    }
    
    return(id_param)
  }
  
  id_param = unlist(find_param_index(param_target, config, multi_param))
  
  if (multi_param){
    registerDoParallel(cores = numb)
    
    subtrace = foreach(
      i = sample(list.files(pattern = "\\.trace$"), numb), .combine = rbind) %dopar% {
        TracesPath <- i
        n <- as.numeric(gsub(
          pattern = paste0("(", net_fname, "-analysis-)|(.trace)"), replacement = "", x = i))
        trace = read.table(TracesPath, header = TRUE)
        
        if ( !all(sapply(trace, is.numeric)) ) {
          print(i)
          failed.tr[f.tr] <<- n
          f.tr <<- f.tr + 1
          
          return()
        }
        
        subtrace <- trace %>% 
          dplyr::select(Time, IECs, Drug, CD, Damage, pheme_e, 
                        pheme_c, BiomassCD, pro_L_e, leu_L_e,
                        ile_L_e, val_L_e, trp_L_e, cys_L_e) %>%
          tidyr::gather(key = "Places", value = "Marking", -Time)
        
        subtrace = cbind(subtrace, config = rep(as.character(i), length(subtrace$Time)))
        
        for (p in 1:length(param_target)) {
          subtrace <- cbind(subtrace, p = rep(config[[id_param[p]]][[n]][[3]], length(subtrace$Time)))
        }
        
        return(subtrace)
      }
    
    stopImplicitCluster()
    
    registerDoParallel(cores = numb)
    
    subflux <- foreach(i = sample(list.files(pattern = "\\.flux$"), numb), .combine = rbind) %dopar% {
      
      TracesPath <- i
      n <- as.numeric(gsub(pattern = paste0("(", net_fname, "-analysis-)|(-0.flux)"), replacement = "", x = i))
      
      # Get the names of the columns from the first few lines of the file
      first_lines <- read.table(TracesPath, nrows = 1)
      df <- data.frame(strings = unlist(first_lines))
      colnames_all <- df$strings
      
      ReactionsNames = c("EX_pro_L_e", "EX_leu_L_e", "EX_ile_L_e", 
                         "EX_val_L_e", "EX_trp_L_e", "EX_cys_L_e", 
                         "sink_pheme_c", "EX_biomass_e")
      
      # Determine which columns to read
      cols_to_read <- c("Time", ReactionsNames)
      colClasses <- ifelse(colnames_all %in% cols_to_read, NA, "NULL")
      
      # Read in only the necessary columns
      flux <- read.table(TracesPath, header = TRUE, colClasses = colClasses)
      
      if (TracesPath %in% paste0(net_fname, "-analysis-", failed.tr, "-0.flux")) {
        print(TracesPath)
        return(NULL)
      }
      
      subflux <- flux %>%
        tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
        dplyr::filter(!grepl("_Ub|_Lb", Reaction)) %>%
        cbind(config = rep(as.character(i), length(.$Time)))
      
      for (p in 1:length(param_target)) {
        subflux <- cbind(subflux, p = rep(config[[id_param[p]]][[n]][[3]], length(subflux$Time)))
      }
      
      return(subflux)
    }
    
    stopImplicitCluster()
    
  } else {
    
    # Define a function to process trace files
    process_trace <- function(i) {
      n <- as.numeric(gsub(pattern = paste0("(", net_name, "-analysis-)|(.trace)"), replacement = "", x = i))
      trace <- read.table(i, header = TRUE)
      subtrace <- trace %>%
        dplyr::select(Time, IECs, Drug, Damage, pheme_e, pheme_c, BiomassCD, CD, 
                      pro_L_e, leu_L_e, ile_L_e, val_L_e, trp_L_e, cys_L_e) %>%
        tidyr::gather(key = "Places", value = "Marking", -Time) %>%
        cbind(config = rep(as.character(i), length(trace$Time)), 
              param_sen = rep(config[[id_param]][[n]][[3]], length(trace$Time)))
      return(subtrace)
    }
    
    # Define a function to process flux files
    process_flux <- function(i) {
      n <- as.numeric(gsub(pattern = paste0("(", net_name, "-analysis-)|(-0.flux)"), replacement = "", x = i))
      flux <- read.table(i, header = FALSE)
      ReactionsNames <- gsub("\\(|\\)", "", gsub("\\(", "_", flux[1, -c(1, 2)]))
      colnames(flux) <- c("Time", "Obj_0", ReactionsNames)
      subflux <- flux[, c("Time", ReactionsNames)] %>%
        tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
        dplyr::filter(!grepl("_Ub|_Lb", Reaction)) %>%
        cbind(config = rep(as.character(i), length(flux$Time)), param_sen = rep(config[[id_param]][[n]][[3]], length(flux$Time)))
      return(subflux)
    }
    
    # Apply the functions to the files
    subtrace <- lapply(sample(list.files(pattern = "\\.trace$"), numb), process_trace)
    subflux <- lapply(sample(list.files(pattern = "\\.flux$"), numb), process_flux)
    
  }
  
  if(multi_param) {
    names(subtrace) <- 
      names(subflux) <- 
      c("Time", "Places", "Marking", "config", param_target)
  }
  
  save(subtrace, file = paste0("subtrace_" , tag, "_", Condition, ".RData"))
  save(subflux, file = paste0("subflux_" , tag, "_", Condition, ".RData"))
  
  if ( multi_param ) {
    
    EX_amino = c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                 "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")
    other_react = c("EX_biomass_e", "sink_pheme_c")
    
    amino = c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
    places = c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", amino)
    units = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))
    df_units = data.frame(places = places, units = units)
    
    # Function to generate plots
    generate_plots <- function(df, 
                               place, 
                               units, 
                               param_target, 
                               colo1, 
                               colo2, 
                               colo3) {
      
      p1 = ggplot(df, aes(x = Time, y = Marking, color = IECsDeath)) + 
        xlab("Time (hour)") + ylab(paste0(place, " ", units)) +
        geom_line(aes(group = config)) + 
        ggtitle(paste0(place, " | ", param_target[1])) +
        scale_colour_gradientn(name = "", colors = colo1, oob = scales::squish)
      
      p2 = ggplot(df, aes(x = Time, y = Marking, color = Death4Treat)) + 
        xlab("Time (hour)") + ylab(paste0(place, " ", units)) +
        geom_line(aes(group = config)) + 
        ggtitle(paste0(place, " | ", param_target[2])) +
        scale_colour_gradientn(name = "", colors = colo2, oob = scales::squish)
      
      p3 = ggplot(df, aes(x = Time, y = Marking, color = Detox)) + 
        xlab("Time (hour)") + ylab(paste0(place, " ", units)) +
        geom_line(aes(group = config)) + 
        ggtitle(paste0(place, " | ", param_target[3])) +
        scale_colour_gradientn(name = "", colors = colo3, oob = scales::squish)
      
      return(list(p1, p2, p3))
    }
    
    # function combing the plots
    combine_plots <- function(plots_list) {
      Reduce(`|`, plots_list)
    }
    
    # in the model aminoacids amounts are measured as (c_molecule) i.e. molecule packet
    pack = 1*(Na*(1/c))
    
    # Generate and combine plots for each dataframe
    for ( i in 1:length(places) ) {
      if ( places[i] %in% places[4:6] ) {
        df = filter(subtrace, subtrace$Places == places[i])
        df$Marking = df$Marking*fc # from pmol to micromol
      } else if( places[i] %in% amino ) {
        df = filter(subtrace, subtrace$Places == places[i])
        df$Marking = ((df$Marking)/pack)*1000 # from c_molecule to micromol
      } else {
        df = filter(subtrace, subtrace$Places == places[i])
      }
      
      plots = generate_plots(df, 
                             places[i], 
                             units[i], 
                             param_target, 
                             colo1, 
                             colo2, 
                             colo3)
      
      assign(paste0("p", places[i]), combine_plots(plots))
    }
    
    generate_flux_plot <- function(data, 
                                   places_filter, 
                                   color_var, 
                                   param_target, 
                                   colo) {
      
      ggplot(data[which(data$Places %in% places_filter), ],
             aes(x = Time, y = Marking, color = get(color_var))) +
        xlab("Time (h)") + 
        ylab("Output fluxes (mmol/gDW*h)") +
        geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
        ggtitle(paste0("Fluxes (mmol/gDW*h) | ", param_target)) +
        scale_colour_gradientn(name = "", colors = colo, oob = scales::squish)
    }
    
    pflux1 = generate_flux_plot(subflux, EX_amino, "IECsDeath", param_target[1], colo1)
    pflux2 = generate_flux_plot(subflux, EX_amino, "Death4Treat", param_target[2], colo2)
    pflux3 = generate_flux_plot(subflux, EX_amino, "Detox", param_target[3], colo3)
    
    pflux_special1 = generate_flux_plot(subflux, other_react, "IECsDeath", param_target[1], colo1)
    pflux_special2 = generate_flux_plot(subflux, other_react, "IECsDeath", param_target[2], colo2)
    pflux_special3 = generate_flux_plot(subflux, other_react, "IECsDeath", param_target[3], colo3)
  }
  
  # Set working directory
  setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper, "/", net_fname, "_", "sensitivity"))
  # Load ranking data
  load(paste0("ranking_", net_fname, "-sensitivity.RData"))
  # Clean up rank$id column
  rank$id <- stringr::str_remove(stringr::str_remove(rank$id, paste0(net_fname, "-analysis-")), ".trace")
  rank$tr <- numeric(length(rank$id))
  
  # Calculate median value and categorize rank$measure as "low" or "high"
  value <- median(rank$measure)
  rank <- mutate(rank, tr = case_when(measure < value ~ "low", measure >= value ~ "high"))
  names(rank)[2] <- "config"
  
  # Create plot for ranking measure distribution
  pr = ggplot(rank, aes(x = "ranking exper", y = measure)) +
    geom_boxplot(width = 0.4, fill = "white") +
    geom_jitter(aes(color = tr, shape = tr), width = 0.25, size = 2.5) +
    scale_color_manual(values = c("#00AFBB", "#E7B800")) +
    labs(subtitle = sprintf(
      "low; measure < %.1f | high; measure >= %.1f", value, value), x = NULL) +
    ggtitle("Ranking measure distribution")
  
  # Clean up rank$config and subtrace$config columns
  rank$config <- stringr::str_remove(
    stringr::str_remove(rank$config, paste0(net_fname, "-analysis-")), ".trace")
  subtrace$config <- stringr::str_remove(
    stringr::str_remove(subtrace$config, paste0(net_fname, "-analysis-")), ".trace")
  subtrace$measure <- numeric(length(subtrace$config))
  # Merge subtrace and subflux with rank to add measure information
  subtrace <- merge(subtrace, rank[, c(1, 2)], by = "config")
  subflux$config <- stringr::str_remove(
    stringr::str_remove(subflux$config, paste0(net_fname, "-analysis-")), "-0.flux")
  subflux$measure <- numeric(length(subflux$config))
  subflux <- merge(subflux, rank[, c(1, 2)], by = "config")
  # Get the ID with the minimum measure value and write it to ConfigID.txt file
  id <- rank$config[which(rank$measure == min(rank$measure))]
  writeLines(id, "ConfigID.txt")
  
  # Extract trace and flux paths based on the ID
  TracesPath <- paste0(net_fname, "-analysis-", id, ".trace")
  FluxPath <- paste0(net_fname, "-analysis-", id, "-0.flux")
  
  # Extract numerical values from ID
  n <- as.numeric(gsub(
    pattern = paste0("(", net_fname, "-analysis-)|(.trace)"), replacement = "", x = id))
  f <- as.numeric(gsub(
    pattern = paste0("(", net_fname, "-analysis-)|(-0.flux)"), replacement = "", x = id))
  
  # Set working directory
  setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper))
  
  # Read trace and flux data
  traceSelected <- read.table(TracesPath, header = T)
  fluxSelected <- read.table(FluxPath, header = T)
  # Reshape traceSelected and fluxSelected data
  traceSelected <- tidyr::gather(traceSelected, key = "Places", value = "Marking", -Time)
  fluxSelected <- tidyr::gather(fluxSelected, key = "Places", value = "Marking", -Time)
  
  places = c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", amino)
  reactions = c("EX_biomass_e", "sink_pheme_c",
                "EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")
  
  # Extract specific data frames for different Places
  Cdiff <- subtrace[subtrace$Places == "CD", ]
  Drug <- subtrace[subtrace$Places == "Drug", ]
  Iecs <- subtrace[subtrace$Places == "IECs", ]
  Biomass <- filter(subtrace, Places == "BiomassCD")
  
  pheme_e <- filter(subtrace, Places == "pheme_e")
  pheme_c <- filter(subtrace, Places == "pheme_c")
  
  trp <- filter(subtrace, Places == "trp_L_e")
  pro <- filter(subtrace, Places == "pro_L_e")
  val <- filter(subtrace, Places == "val_L_e")
  ile <- filter(subtrace, Places == "ile_L_e")
  cys <- filter(subtrace, Places == "cys_L_e")
  leu <- filter(subtrace, Places == "leu_L_e")
  
  flux = filter(subflux, subflux$Places %in% reactions)
  
  create_plot <- function(data, 
                          x_col, 
                          y_col, 
                          color_col, 
                          places_col, 
                          id_col, 
                          title_prefix) {
    
    ggplot(data, aes_string(x = x_col, y = y_col, color = color_col)) +
      xlab("Time (hour)") + ylab(y_col) +
      geom_line(aes_string(group = "config"), alpha = alph) + 
      geom_line(data = filter(
        traceSelected, traceSelected[[places_col]] == id_col), col = "red") + 
      ggtitle(paste0(title_prefix, " | config selected: ", id_col)) + 
      scale_colour_gradientn(name = "", colors = c("#00AFBB", "#E7B800"), 
                             oob = scales::squish)
  }
  
  pCM <- create_plot(Cdiff, "Time", "Marking", "measure.y", "Places", "CD", "CD Place")
  pdrugM <- create_plot(Drug, "Time", "Marking", "measure.y", "Places", "Drug", "Drug Place")
  pIECsM <- create_plot(Iecs, "Time", "Marking", "measure.y", "Places", "IECs", "IECs")
  pBiomassM <- create_plot(Biomass, "Time", "Marking", "measure.y", "Places", "BiomassCD", "Biomass")
  pfluxM <- ggplot(flux, aes(x = Time, y = Marking, color = measure.y)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle("Fluxes (mmol/gDW*h)") +
    scale_colour_gradientn(name = "", colors = c("#00AFBB", "#E7B800"), oob = scales::squish)
  
  # Calculate median
  median <- median(Cdiff$measure.y)
  
  # CD Place
  pCM.median <- ggplot(Cdiff, aes(x = Time, y = Marking, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
    xlab("Time (hour)") + ylab("C. difficile (cell)") +
    geom_line(aes(group = config), alpha = alph) + 
    geom_line(data = filter(traceSelected, traceSelected$Places == "CD"), col = "red") + 
    ggtitle("CD Place") +
    scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800"))
  
  # Drug Place
  pdrugM.median <- ggplot(Drug, aes(x = Time, y = Marking*fc, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
    xlab("Time (hour)") + ylab(paste0("Drug", " ", "(µmol)")) +
    geom_line(aes(group = config), alpha = alph) + 
    geom_line(data = filter(traceSelected, traceSelected$Places == "Drug"), col = "red") + 
    ggtitle("Drug Place") +
    scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800"))
  
  # IECs Place
  pIECsM.median <- ggplot(Iecs, aes(x = Time, y = Marking, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
    xlab("Time (hour)") + ylab("IECs (cell)") +
    geom_line(aes(group = config), alpha = alph) + 
    geom_line(data = filter(traceSelected, traceSelected$Places == "IECs"), col = "red") + 
    ggtitle("IECs Place") +
    scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800"))
  
  # Biomass Place
  pBiomassM.median <- ggplot(Biomass, aes(x = Time, y = Marking, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
    xlab("Time (hour)") + ylab("BiomassCD (pg)") +
    geom_line(aes(group = config), alpha = alph) + 
    geom_line(data = filter(traceSelected, traceSelected$Places == "BiomassCD"), col = "red") + 
    ggtitle("Biomass Place") +
    scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800"))
  
  # Sensitivity Analysis: Fluxes
  pfluxM.median <- ggplot(flux, aes(x = Time, y = Marking, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config), alpha = alph) +
    facet_wrap(~Places, scale = "free") + 
    ggtitle("Fluxes (mmol/gDW*h)") + 
    scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800")) + 
    facet_wrap(~Places, scale = "free")
  
  p = (pCD/pBiomassCD | pDrug/pIECs) / ((pcys_L_e / pval_L_e) | (pleu_L_e/ppro_L_e) | (pile_L_e/ptrp_L_e))
  p_heme = (ppheme_e | ppheme_c) 
  Fl = (pflux1 | pflux2 | pflux3) / (pflux_special1 | pflux_special2 | pflux_special3)
  pM = ((pCM + pdrugM) / (pIECsM + pBiomassM)) + pfluxM
  pM.median = ((pCM.median + pdrugM.median) / (pIECsM.median + pBiomassM.median)) + pfluxM.median
  
  setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper, "/", net_fname, "_", "sensitivity"))
  
  ggsave(p, file = paste0(
    Exper, "_multiparam_places_", tag, "_", Condition, ".pdf"), width = 32, height = 7)
  ggsave(Fl, file = paste0(
    Exper, "_multiparam_fluxes_", tag, "_", Condition, ".pdf"), width = 15, height = 6)
  ggsave(ppheme_e, file = paste0(
    Exper, "_heme", tag, "_", Condition, ".pdf"), width = 12, height = 2)
  ggsave(pM, file = paste0(
    "RankingPlot", tag, "_", Condition, ".pdf"), width = 9, height = 10)
  ggsave(pM.median, file = paste0(
    "RankingPlot_median", tag, "_", Condition, ".pdf"), width = 9, height = 10)
  ggsave(pr, file = paste0(
    "RankingBoxPlot", tag, "_", Condition, ".pdf"), width = 6, height = 4)
  
}

ViolinPlotSensitivityAnalysis = function(Tempi,
                                         multi_param,
                                         tag, 
                                         Condition, 
                                         Exper,
                                         net_fname,
                                         per,
                                         Na,
                                         c,
                                         fc) {
  
  if (multi_param) {
    setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper))
  } else {
    setwd(paste90(wd, "/results/CDiff_", tag, "_", Condition, "_", Exper))
  }
  
  load(paste0(net_fname, "-analysis.RData"))
  
  if (multi_param) {
    id_param = sapply(param_target, function(x) which(sapply(config, function(y) y[[1]][[1]] == x)))
  } else {
    id_param = which(sapply(config, function(x) x[[1]][[1]] == param_target))
  }
  
  setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper))
  
  load(file = paste0("subtrace_" , tag, "_", Condition, ".RData"))
  load(file = paste0("subflux_" , tag, "_", Condition, ".RData"))
  
  subtrace = dplyr::filter(subtrace, Time %in% Tempi)
  subtrace$config = stringr::str_remove(subtrace$config, paste0(net_fname, "-analysis-"))
  subtrace$config = stringr::str_remove(subtrace$config, ".trace")
  
  # Create a list of your dataframe names
  amino = c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  df_names = c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", amino)
  marking_units = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))
  
  # Create a list of dataframes
  df_list = lapply(df_names, function(x) dplyr::filter(subtrace, subtrace$Places == x))
  
  # Assign dataframes to variables without using list2env
  for(i in seq_along(df_list)) {
    assign(df_names[i], df_list[[i]])
  }
  
  # Define the function
  process_dataframe = function(df, df_name, pack = 1*(Na*(1/c))) {
    # Handle special cases
    if (df_name %in% df_names[3:5]) {
      df$Marking = df$Marking*fc # from pmol to micromol
    } else if (df_name %in% amino) {
      df$Marking = ((df$Marking)/pack)*1000 # from c_molecule to micromol
    } else {
      df$Marking = as.double(df$Marking)
    }
    
    df$Time = as.character(df$Time)
    
    # Calculate the number of dots to display
    num_dots <- round(nrow(df) * per)
    # Sample a subset of the data for jitter
    subset_df <- df[sample(nrow(df), num_dots), ]
    # Calculate the mean for each time
    means <- aggregate(Marking ~ Time, df, mean)
    # Return the processed dataframe and means
    list(df = df, subset_df = subset_df, means = means)
  }
  
  # Use lapply to apply the function to each dataframe
  processed_list <- lapply(df_names, function(df_name) {
    df <- get(df_name)
    process_dataframe(df, df_name)
  })
  
  violin_list = list()
  
  for ( i in 1:length(df_names) ) {
    
    violin_list[[i]] = processed_list[[i]]$df %>%
      ggplot(aes(x = Time, y = Marking)) +
      labs(x = "Time (h)", 
           y = paste0("' ", df_names[i], " '", " Marking ", marking_units[i]),
           title = paste0("Distribution of configurations | ","' ", df_names[i], " '", " Place"),
           subtitle = paste("Total observations: ", 
                            length(unique(processed_list[[i]]$df$config)))) +
      theme(text = element_text(size = 8, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
            legend.position = "none",
            plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
            plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
            plot.title.position = "plot", 
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 12, face = "bold")) +
      # Use subset data for jitter
      geom_jitter(data = processed_list[[i]]$subset_df, aes(col = Time), alpha = 0.35) +
      # Use all data for violin
      geom_violin(scale = "width", alpha = 0.3, width = 1) +
      geom_boxplot(alpha = 0) +  # Use all data for boxplot
      geom_point(data = processed_list[[i]]$means, aes(y = Marking), 
                 color = "darkred", size = 4) +  # Add mean points
      ggrepel::geom_label_repel(data = processed_list[[i]]$means, 
                                aes(label = paste("mean =", round(Marking, 2)), 
                                    fill = "white"), size = 3.5,
                                color = 'black', box.padding = unit(0.45, "lines"), 
                                segment.color = 'grey50') + 
      ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
      scale_x_discrete(name = "", label = rep(paste0(
        paste0(Tempi, "\n"),"(n = ", length(unique(processed_list[[i]]$df$config)), ")"), 
        length(Tempi)))
    
  }
  
  pV1 = (violin_list[[7]] | violin_list[[8]] | violin_list[[9]]) / 
    (violin_list[[10]] | violin_list[[11]] | violin_list[[12]])
  pV2 = ((violin_list[[2]] | violin_list[[3]])/(violin_list[[1]] | violin_list[[4]]))
  pV3 = (violin_list[[5]] | violin_list[[6]])
  
  setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper, "/", net_fname, "_", "sensitivity"))
  
  ggsave(pV1, file = paste0("ViolinPlot_AA_", tag, "_", Condition, ".pdf"), width = 15, height = 8)
  ggsave(pV2, file = paste0("ViolinPlot_Places_", tag, "_", Condition, ".pdf"), width = 9, height = 9)
  ggsave(pV3, file = paste0("ViolinPlot_Heam_", tag, "_", Condition, ".pdf"), width = 12, height = 5)
  
}
