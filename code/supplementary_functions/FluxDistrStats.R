
distr_stats = function(result.dir,
                       fbasol.file,
                       file.standard_dev,
                       thrsd,
                       linethrsd) {
  
  fbasol = readRDS(fbasol.file)
  
  # Most variable fluxes
  fbasol = as.data.frame(fbasol)
  
  # Define function to calculate statistics for each Reaction
  get_stats = function(r) {
    
    d = dplyr::filter(fbasol, Reaction == r)
    Flux = tidyr::gather(d, key, value, -c(Reaction))$value
    Flux = as.numeric(Flux, digits = 12)
    sd = sd(Flux)
    E = mean(Flux)
    median = median(Flux)
    variance = var(Flux)
    c(r, as.double(sd), E, median, variance)
    
  }
  
  # use multicore, set to the number of our cores
  registerDoParallel((detectCores()))
  
  data = foreach(r = unique(fbasol$Reaction), 
                 .combine = rbind) %dopar% {
                   get_stats(r)
                 }
  
  colnames(data) = c("react.id", "sd", "E", "Median", "Variance")
  row.names(data) = NULL
  data = as.data.frame(data)
  data$sd = as.double(data$sd); data$E = as.double(data$E)
  data$Median = as.double(data$Median); data$Variance = as.double(data$Variance) 
  
  save(data, file = file.standard_dev)
  
  names(all.react)[1] = names(data)[1]
  
  subdata = merge(x = data, y = all.react, by = "react.id", all.x = TRUE)
  
  mode = unique(as.double(subdata$E))[which.max(
    tabulate(match(as.double(subdata$E), unique(as.double(subdata$E)))))]
  
  # Plotting fluxes variance among configurations and differences compare to reference
  title.size = 8
  subtitle.size = 6
  text.size = 8 
  axistext.size = 6
  axistitle.size = 6
  legendsize = 0.25
  point.size = 0.5
  
  subdata.sd = dplyr::filter(subdata, as.double(subdata$sd) > thrsd)
  
  psdAll = ggplot(subdata, aes(x = as.double(sd), fill = "")) + 
    labs(x = "sd (mmol/gDW*h)", y = "Counts",
         title = paste0("SD Distr. | ", "number of reaction = ",
                        length(subdata$react.id)),
         subtitle = paste0("dotted line : sd = ", linethrsd, "\n", 
                           paste0("config = ", (dim(param_values)[1] + 1)))) +
    theme(text = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold")) + 
    geom_vline(xintercept = linethrsd, linetype = "dotted", color = "darkred", lwd = 0.5) +
    geom_histogram(aes(y = ..count..), alpha = 0.5, 
                   position = "identity", binwidth = 9,
                   col = "#651e3e", fill = "#651e3e")
  
  psd = ggplot(subdata.sd, aes(x = as.double(sd), fill = "")) + 
    labs(x = "sd (mmol/gDW*h)", y = "Counts",
         title = paste0("SD Distr. | ", "number of reaction = ",
                        length(subdata.sd$react.id)),
         subtitle = paste0("sd < ", thrsd, " was filtered out", "\n", 
                           "dotted line : sd = ", linethrsd, "\n" , 
                           paste0("config = ", (dim(param_values)[1] + 1)))) +
    theme(text = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold")) + 
    geom_vline(xintercept = linethrsd, linetype = "dotted", color = "darkred", lwd = 0.5) +
    geom_histogram(aes(y = ..count..), alpha = 0.5, 
                   position = "identity", binwidth = 9,
                   col = "#ceb5b7", fill = "#ceb5b7")
  
  psd_wrapAll = ggplot(subdata, aes(as.double(sd))) +
    geom_histogram(aes(y = (..count..), fill = ReactionDistrict)) +
    facet_wrap("ReactionDistrict", scales = "free") +
    labs(x = "sd (mmol/gDW*h)", y = "Counts",
         title = paste("SD Distr. | ", "number of reaction = ", 
                       length(subdata$react.id), sep = ""),
         subtitle = paste0("Number of ractions (core) = ", 
                           sum(subdata$ReactionDistrict == "core"), " \n",
                           "Number of ractions (boundary) = ", 
                           sum(subdata$ReactionDistrict == "boundary"), " \n", 
                           paste0("config = ", (dim(param_values)[1] + 1)))) +
    theme(text = element_text(size = text.size, color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          legend.key.size = unit(legendsize, 'cm'),
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold"))
  
  psd_wrap = ggplot(subdata.sd, aes(as.double(sd))) +
    geom_histogram(aes(y = (..count..), fill = ReactionDistrict)) +
    facet_wrap("ReactionDistrict", scales = "free") +
    labs(x = "sd (mmol/gDW*h)", y = "Counts",
         title = paste0("SD Distr. | ", "number of reaction = ", 
                        length(subdata.sd$react.id)),
         subtitle = paste0("Number of ractions (core) = ", 
                           sum(subdata.sd$ReactionDistrict == "core"), " \n",
                           "Number of ractions (boundary) = ", 
                           sum(subdata.sd$ReactionDistrict == "boundary"), " \n", 
                           paste0("config = ", (dim(param_values)[1] + 1)))) +
    theme(text = element_text(size = text.size, color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          legend.key.size = unit(legendsize, 'cm'),
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold"))
  
  ggsave((psdAll|psd_wrapAll)/(psd|psd_wrap), 
         file = paste0(result.dir, "/sd.png"), width = 7, height = 4)
  
  ggsave((psdAll|psd_wrapAll)/(psd|psd_wrap), 
         file = paste0(result.dir, "/sd.pdf"), width = 7, height = 4)
  
  subdata.sd_boundary = dplyr::filter(subdata.sd, subdata.sd$ReactionDistrict == "boundary")
  
  ssd = dplyr::filter(subflux.boun, subflux.boun$Reaction %in% subdata.sd_boundary$react.id)
  
  react.subset = ggplot(ssd, aes(x = Flux, fill = ordered(Reaction, levels = unique(Reaction)))) +
    geom_density(aes(y = ..count..), alpha = 0.5, position = "identity", col = c, fill = c) +
    ggplot2::labs(x = "Flux (mmol/gDW*h)", y = "Count",
                  title = paste0("Reactions with sd > ", thrsd),
                  subtitle = paste0("Number of reactions: ", length(unique(ssd$Reaction)), " \n", 
                                    "Number of configurations: ", length(unique(ssd$config)))) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Condition block")) +
    theme(text = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
          plot.title = element_text(size = 14, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold")) +
    facet_wrap(~ Reaction, scales = "free")
  
  dplyr::filter(subdata, react.id %in% targetrnx)[, c(1, 2)]
  
  ggsave(react.subset, file = paste0(result.dir, "/react_subsetsd.png"), 
         width = sqrt(length(unique(ssd$Reaction)))*1.75, height = sqrt(length(unique(ssd$Reaction))*1.75))
  
  ggsave(react.subset, file = paste0(result.dir, "/react_subsetsd.pdf"),
         width = sqrt(length(unique(ssd$Reaction)))*1.75, height = sqrt(length(unique(ssd$Reaction))*1.75))
  
}