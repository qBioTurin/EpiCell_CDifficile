library(dplyr)
library(ggplot2)
library(patchwork)
library(epimod)
library(grid)
library(gtable)
library(stringr)
library("ggplotify")

changecolorsFacet = function(pt,colors){
  # Convert the plot to a gtable object
  g <- ggplotGrob(pt)
  
  # Locate the strip grobs for the x-axis (top strips)
  strip_idx <- which(grepl("strip-t", g$layout$name))
  
  # Customize the strip backgrounds
  for (i in seq_along(strip_idx)) {
    # Get the corresponding facet label
    facet_label <- g$grobs[[strip_idx[i]]]$children
    
    strip_grob <- g$grobs[[strip_idx[i]]]
    label_grob <- strip_grob$grobs[[1]]$children[[which(grepl("text", names(strip_grob$grobs[[1]]$children)))]]
    facet_label <- label_grob$children[[1]]$label
    
    # Set the fill color based on the mapping
    fill_color <- colors[facet_label]
    
    # Replace the background rect
    g$grobs[[strip_idx[i]]]$grobs[[1]]$children[[1]]$gp$fill <- fill_color
  }
  
  return(as.ggplot(g))
  
}
ModelAnalysisPlot=function(TracesPath, FluxPath, FluxVec,new_eps_value,config,tag) {
  
  trace = read.table(TracesPath, header = T)
  
  subtrace = trace %>% tidyr::gather(key = "Places", value = "Marking", -Time)
  
  
  flux = read.csv(FluxPath, header=FALSE, sep=" ")
  colnames(flux) = flux[1, ]
  flux = flux[-1, ]
  
  if(is.null(FluxVec)){
    
    FluxVec = colnames(flux)
    flux %>% dplyr::select(FluxVec, Time) %>% tidyr::gather(key = "Reaction", value = "Flux", -Time)
    
    subflux = flux %>% select(FluxVec, Time) %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
      ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x=Time, y=Flux, group = Reaction)) +
      ggplot2::facet_wrap(~Reaction, scales = "free_y") + ggplot2::scale_fill_brewer("Accent")
    
    print(subflux, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    
  } else {
    
    subflux = flux[,c(FluxVec, "Time")] %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time)
    
    subflux$Time = as.numeric(as.character(subflux$Time))
    subflux$Flux = as.numeric(as.character(subflux$Flux))
  }
  
  subflux$new_eps_value = subtrace$new_eps_value = new_eps_value
  subflux$config = subtrace$config = config
  subflux$tag = subtrace$tag = tag
  
  return(list(subflux = subflux,subtrace = subtrace))
}
update_eps_value <- function(file_path, new_eps_value) {
  # Read the file into a vector of lines
  lines <- readLines(file_path)
  
  # Find the line containing "eps;"
  eps_line_index <- grep("^g;\\s*eps;", lines)
  
  # Check if the eps line exists
  if (length(eps_line_index) == 0) {
    stop("No line with 'eps' found in the file.")
  }
  
  # Extract the specific line
  eps_line <- lines[eps_line_index]
  
  # Modify the last value after the last semicolon
  modified_eps_line <- sub("(g;\\s*eps;.*?).*", paste0("\\1 ", new_eps_value), eps_line)
  
  # Update the line in the lines vector
  lines[eps_line_index] <- modified_eps_line
  
  # Optionally, write the modified lines back to the file
  writeLines(lines, file_path)
  
  # Return the modified content as a vector of lines
  return(lines)
}
# Function to set working directory to the closest folder with a specific name
set_closest_directory <- function(target_name) {
  wd = getwd()
  if( basename(wd) != target_name){
    # Get all directories relative to the current working directory
    all_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
    
    # Filter directories with the specific name
    matching_dirs <- all_dirs[basename(all_dirs) == target_name]
    
    # Check if there are matches
    if (length(matching_dirs) == 0) {
      stop("No folder with the name '", target_name, "' found in the current working directory or its subdirectories.")
    }
    
    # Find the directory with the shortest path (closest to the current working directory)
    closest_dir <- matching_dirs[which.min(nchar(matching_dirs))]
    
    # Set the working directory
    setwd(closest_dir)
  }
  # Return the new working directory
  return(getwd())
}
wd <- set_closest_directory("EpiCell_CDifficile")
setwd(wd)

system("rm -f dockerID *error.log *.log ExitStatusFile" )

if(!dir.exists("./results")) dir.create("./results")
if(!dir.exists("./results/TimingEval")) dir.create("./results/TimingEval")

######## Variables definition:
Condition <- "Therapy"
model.type = "Clostridium"
model.name = "CD196HemeSink"
net_fname = "EpitCellDifficileHemeSink"
fba_fname = paste0(model.name, ".txt")
atol = 1e-06
rtol = 1e-06
time.step = 1
f_time = 17
# f_time = 72
event_times = c(8, 16)
# event_times = c(8, 16, 24, 32, 40, 48)
distance_measure = "ReferenceM"
parameters_fname <- "input/csv/ParametersListSinkHeme_eps.csv"
supp_function.dir = "/code/supplementary_functions/"

######################


# #### parameters closest to the median trances in FIG2B ####
# source("~/EpiCell_CDifficile/code/supplementary_functions/MedianParamsConfiguration.R")
# 
# resParams = configuration.closeMedian(Condition = "Therapy",
#                                       tag = c( "Unified"))
# 
# resParams$plot
# saveRDS(resParams,file = "paramsProve.RDs") ## Da caricare se non si ha la sensitivity
resParams <- readRDS("~/git/EpiCell_CDifficile/paramsProve.RDs")
############

# epstimes = c("1e-6","1e-4","1e-2")
# numbConfig = 1:dim(resParams$Config)[1]
# tag = c("Ablated", "Unified")
# paramsgrid = rbind(expand.grid("1e-6",numbConfig,"Ablated"), expand.grid(epstimes,numbConfig,"Unified"))

tag = c("Ablated", "ParAblated", "Unified")
paramsgrid = expand.grid("1e-6",1,tag)

MultipleAnalysis = lapply(seq_along(paramsgrid[,1]),
                          function(c,paramsConfig,paramsgrid){
                            new_eps_value = paramsgrid[c,"Var1"]
                            nconfig = paramsgrid[c,"Var2"]
                            tag = paramsgrid[c,"Var3"]
                            updated_content = update_eps_value(parameters_fname, new_eps_value = new_eps_value)
                            
                            # Print updated content
                            cat(updated_content, sep = "\n")
                            
                            ######## Model creation
                            transitions_fname = paste0(wd, "/Net/",tag,".cpp")
                            
                            model.generation(net_fname = paste0(wd, "/Net/", net_fname, ".PNPRO"),
                                             transitions_fname = transitions_fname,
                                             fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname))
                            
                            system(paste0("mv ", net_fname, ".* ./Net"))
                            ######################
                            
                            execution_start <- Sys.time()
                            model.analysis(solver_fname = paste0(wd, "/Net/", net_fname, ".solver"),
                                           i_time = 0, f_time = f_time, s_time = time.step,
                                           fba_fname = paste0(wd, "/input/CompiledModels/", fba_fname),
                                           atol = atol, rtol = rtol, debug = T,
                                           ini_v = as.numeric(paramsConfig[nconfig,c("Detox","Death4Treat","IECsDeath")]),
                                           parameters_fname = paste0(wd, "/", parameters_fname),
                                           functions_fname = paste0(wd, supp_function.dir, "Functions.R"),
                                           event_times = if (Condition != "NoDrug") event_times,
                                           event_function = if (Condition != "NoDrug") "treat_generation")
                            execution_end <- Sys.time()
                            
                            execution_time = as.numeric(difftime(execution_end, execution_start, units = "secs"))
                            
                            # Read the file
                            file_path <- list.files(path = wd, pattern = "\\.log$", full.names = TRUE)
                            
                            file.rename(file_path, paste0(tag, ".cout"))

                            # # lines <- readLines(file_path)
                            # # file.remove(file_path)
                            # # Filter lines containing "Total memory used" or "Total time required"
                            # extracted_lines <- grep("Total memory used|Total time required", lines, value = TRUE)
                            # 
                            # # Extract time and memory separately
                            # times <- as.numeric(sub(".*Total time required: ([0-9]+)s.*", "\\1", extracted_lines[grepl("time required", extracted_lines)]))
                            # memory <- as.numeric(sub(".*Total memory used: ([0-9]+)KB.*", "\\1", extracted_lines[grepl("memory used", extracted_lines)]))
                            # 
                            # # Create a data frame
                            # result_df <- data.frame(Time_s = sum(times), Memory_KB = sum(memory),
                            #                         GlobalExecution_time = execution_time, eps = new_eps_value,
                            #                         params = paramsConfig[nconfig,"ConfParams"], Scenario = tag)
                            
                            resFolder = paste0("./results/TimingEval/CDiff", tag, Condition,
                                               new_eps_value, gsub(pattern = " ",replacement = "",x = paramsConfig[nconfig,"ConfParams"]), collapse  = "_")
                            
                            if(dir.exists(resFolder)) system(paste("rm -r ",resFolder))
                            system(paste0("mv ", net_fname, "_analysis* ",resFolder))
                            
                            traces = ModelAnalysisPlot(paste0(resFolder,"/EpitCellDifficileHemeSink-analysis-1.trace"),
                                                       paste0(resFolder,"/EpitCellDifficileHemeSink-analysis-1-0.flux"), 
                                                       FluxVec = c(
                                                         "EX_biomass_e", "sink_pheme_c", "EX_cys_L_e", 
                                                         "EX_trp_L_e", "EX_ile_L_e", "EX_pro_L_e", "EX_leu_L_e", "EX_val_L_e"),
                                                       new_eps_value,
                                                       paramsConfig[nconfig,"ConfParams"],
                                                       tag
                            )
                            return(list(traces = traces))
                          },paramsConfig = resParams$Config,paramsgrid)


# timing = do.call(rbind, lapply(MultipleAnalysis,"[[",2) )
flux = do.call(rbind, lapply(lapply(MultipleAnalysis,"[[",1),"[[",1) )
subtrace = do.call(rbind, lapply(lapply(MultipleAnalysis,"[[",1),"[[",2) )
# subtrace$new_eps_value = factor(subtrace$new_eps_value,levels = epstimes)
# flux$new_eps_value = factor(flux$new_eps_value,levels = epstimes)
#subtrace$config = factor(subtrace$config, levels = numbConfig)
#flux$config = factor(flux$config, levels = numbConfig)
subtrace$ConfParams = subtrace$config

ggplot( subtrace%>%filter(Time < 15) )+
  geom_line(aes(x = Time, y = Marking, linetype = tag, col = tag ))+
  facet_wrap(Places~new_eps_value,scales = "free")+
  theme_bw()+
  theme(
    plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
    plot.title.position = "plot", 
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    strip.text.x = element_text(size = 10, face = "bold", colour = "white"),
    strip.text.y = element_text(size = 10, face = "bold", colour = "black"),
    strip.background.y = element_rect( fill = "white"),
    legend.position = "top")+
  labs(x = "Time (h)",y = "Quantity", linetype = "Scenario") +
  geom_vline(xintercept = 8)

saveRDS(list(trace = subtrace,flux = flux, timing = timing),file = "tracesFluxes.RDs")
tracesFluxes <- readRDS(paste0(wd,"/tracesFluxes.RDs"))
trace= tracesFluxes$trace
flux= tracesFluxes$flux

# Define color mappings for facets

colors_new_confParams <- c("#3B9AB2", "#89ca66b6", "#972D15")
names(colors_new_confParams) <- unique(trace$ConfParams)

colors_new_eps_value <- rev(grey.colors(5))#c("#F21A00", "#3B9AB2" , "#D8A499", "#EBCC2A" ,"#C6CDF7")
names(colors_new_eps_value) <- unique(trace$new_eps_value)

source("~/EpiCell_CDifficile/code/plot_functions/Fig2_plot.R")
Fig2 = plotting_Fig2_paper(Exper = "Model_Sensitivity",
                    Condition = "Therapy",
                    tag = c("Ablated", "Unified"),
                    param_target = "IECsDeath",
                    Tempi = c(0, 12, 24, 36, 48, 60),
                    wd = wd,
                    subtrace = trace,
                    colConfigSets = colors_new_confParams)

ggsave(plot = Fig2$pl2B,filename = "Figures/Fig2B.pdf",width = 8,height = 10)
ggsave(plot = Fig2$pl2C,filename = "Figures/Fig2C.pdf",width = 15,height = 4.5)



pt = ggplot(trace%>%filter(Places %in% c("CD","IECs","BiomassCD")))+
  geom_line(aes(x = Time, y = Marking, linetype = tag, col = ConfParams ))+
  facet_grid(Places~new_eps_value,scales = "free")+
  scale_color_manual(values = colors_new_confParams[unique(trace$ConfParams)],guide="none")+
  theme_bw()+
  theme(
        plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 15, face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 10, face = "bold", colour = "white"),
        strip.text.y = element_text(size = 10, face = "bold", colour = "black"),
        strip.background.y = element_rect( fill = "white"),
        legend.position = "top")+
  labs(x = "Time (h)",y = "Quantity", linetype = "Scenario")

g = changecolorsFacet(pt,colors_new_eps_value)
g

# Calculate percentage difference

df_diff <- trace %>%
  filter(Places %in% c("CD","IECs","BiomassCD"),tag!="ParAblated") %>%
  mutate(new_eps_value = as.character(new_eps_value) ) %>% 
  mutate(new_eps_value = str_replace(new_eps_value, "1e-6", "baseline")) %>%
  group_by(Time,tag,config) %>%
  tidyr::spread(key = new_eps_value, value = Marking) %>%
  tidyr::gather(-Time,-Places,-baseline,-config, -tag,-ConfParams, key = "new_eps_value", value = "Marking") %>%
  mutate(
    perc_diff = 100 * (Marking - baseline) / baseline        # % Difference
  )


pt_perc = ggplot(df_diff %>% filter(tag == "Unified"))+
  geom_bar(aes(x = as.factor(Time), y = perc_diff, group = tag, fill = new_eps_value),
           stat = "identity", position = "dodge")+
  facet_grid(Places~ConfParams,scales = "free")+
  theme_bw()+
  scale_fill_manual(values = colors_new_eps_value[unique(df_diff$new_eps_value)])+
  theme(
        plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 15, face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 10, face = "bold", colour = "white"),
        strip.text.y = element_text(size = 10, face = "bold", colour = "black"),
        strip.background.y = element_rect( fill = "white"),
        legend.position = "none")+
  labs(x = "Time (h)",y = "% difference with 1e-6")+ 
  scale_x_discrete(breaks = seq(0,max(df_diff$Time),20))


g2= changecolorsFacet(pt_perc,colors_new_confParams)

Fig2D = (g|g2) & theme(legend.position = "bottom")
ggsave(plot = Fig2D,filename = "Figures/Fig2D.pdf",width = 20,height = 5)


#################

pf = ggplot(flux)+
  geom_line(aes(x = Time, y = Flux, col = new_eps_value,linetype = tag))+
  facet_grid(config~Reaction,scales = "free")+
  theme_bw()

# Scale Memory_KB to match the range of Time_s
range_time <- range(c(timing$Time_s,timing$GlobalExecution_time))
range_memory <- range(timing$Memory_KB)
scale_factor <- diff(range_time) / diff(range_memory)
timing$Memory_KB_scaled <- (timing$Memory_KB - min(range_memory)) * scale_factor + min(range_time)

timing$ID = 1
timing$eps = factor(timing$eps,levels = epstimes)
# Plot
p <- ggplot(timing, aes(x = eps ,group = ID) ) +
  geom_line(aes(y = Time_s, color = "Time (s)"), linewidth = 1) +  # Time_s line
  geom_line(aes(y = GlobalExecution_time, color = "Global Time (s)"), linewidth = 1, linetype = "dashed") +  # Time_s line
  geom_line(aes(y = Memory_KB_scaled, color = "Memory (KB)"), linewidth = 1) +  # Scaled Memory_KB line
  scale_y_continuous(
    name = "Time (s)",
    sec.axis = sec_axis(~ (. - min(range_time)) / scale_factor + min(range_memory), 
                        name = "Memory (KB)")
  ) +
  scale_color_manual(values = c("Time (s)" = "blue", "Memory (KB)" = "red","Global Time (s)" = "blue")) +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "blue"),  # Secondary axis styling
    axis.text.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "red"),  # Secondary axis styling
    axis.text.y.right = element_text(color = "red"),
    legend.title = element_blank()
  ) +
  labs(x = "Epsilon (eps)", color = "Variable")

(pf|pt)/p+ plot_layout(heights  = c(2, 1),guides = 'collect')


