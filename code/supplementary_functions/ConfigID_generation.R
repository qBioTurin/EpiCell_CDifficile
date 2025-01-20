
csvFileUpdating.minRankConfigID = function(wd, net_fname) {
  
  tag = "Unified"
  Exper = "Model_Sensitivity"
  Condition = "Therapy"
  
  # Load ranking data
  load(paste0(wd,"/results/CDiff", "_", tag, "_", Condition, "_", Exper, "/", net_fname, "_", "sensitivity/ranking_", net_fname, "-sensitivity.RData"))
  
  # Clean up rank$id column
  rank$id <- stringr::str_remove(stringr::str_remove(rank$id, paste0(net_fname, "-analysis-")), ".trace")
  rank$tr <- numeric(length(rank$id))
  
  # Calculate median value and categorize rank$measure as "low" or "high"
  value <- median(rank$measure)
  rank <- mutate(rank, tr = case_when(measure < value ~ "low", measure >= value ~ "high"))
  names(rank)[2] <- "config"
  
  # Create plot for ranking measure distribution
  pr <- ggplot(rank, aes(x = "ranking exper", y = measure)) +
    geom_boxplot(width = 0.4, fill = "white") +
    geom_jitter(aes(color = tr, shape = tr), width = 0.25, size = 2.5) +
    scale_color_manual(values = c("#00AFBB", "#E7B800")) +
    labs(subtitle = sprintf("low; measure < %.1f | high; measure >= %.1f", value, value), x = NULL) +
    ggtitle("Ranking measure distribution")
  
  # Clean up rank$config and subtrace$config columns
  rank$config <- stringr::str_remove(stringr::str_remove(rank$config, paste0(net_fname, "-analysis-")), ".trace")
  
  # Get the ID with the minimum measure value and write it to ConfigID.txt file
  id <- rank$config[which(rank$measure == min(rank$measure))]

  print(paste0("#### Configuration ID: ", id))

  param_target = c("Detox", "Death4Treat", "IECsDeath")
  
  load(paste0(wd,"/results/CDiff", "_", tag, "_", Condition, "_", Exper, "/", net_fname, "-analysis.RData"))

  id_param <- sapply(param_target, function(x) {
    which(sapply(config, function(y) y[[1]][[1]] == x))
  })
  
  for (p in seq_along(param_target)) {
    assign(param_target[p], config[[id_param[p]]][[as.numeric(id)]][[3]])
  }
  
  print(paste0("#### Modifying csv file:\n ", wd, "/input/csv/ParametersListSinkHeme.csv"))
  
  csv <- read.csv(file = paste0(wd, "/input/csv/ParametersListSinkHeme.csv"), header = F, quote = "")
  
  par = c(paste0("c; IECsDeath; ", IECsDeath),
          paste0("g; Death4Treat; ", Death4Treat),
          paste0("c; Detox; ", Detox))
  
  for(p in 1:length(id_param)) {csv[id_param[p], ] = par[p]}
  
  write.table(csv, file = paste0(wd, "/input/csv/ParametersListSinkHeme.csv"), col.names = FALSE, row.names = F, quote = F)
  write.table(csv, file = paste0(wd, "/input/csv/ParametersListSinkHeme_eps.csv"), col.names = FALSE, row.names = F, quote = F)
  print(paste0("#### File Updated ######"))
  
  return(list(plot = pr, configID = id))
}
