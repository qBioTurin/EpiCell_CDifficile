
# libraries required
library(dplyr)
library(tidyr)
library(epiR)
library(tibble)
library(ggplot2)
library(readr)
library(sensitivity)
library(pander)
library(RColorBrewer)
library(patchwork)
library(fbar)
library(foreach)
library(doParallel)
library(fbar)
library(stringr)

# setting working directory
wd = getwd()
setwd(wd)

## functions and methods ##

## FBA processing methods
source(paste0(wd, "/code/processing_functions/ExtractEx.R"))
source(paste0(wd, "/code/processing_functions/FBAmodel_Metadata.R"))
source(paste0(wd, "/code/processing_functions/createReactionString.R"))
source(paste0(wd, "/code/processing_functions/modelorg2tsv.R"))
source(paste0(wd, "/code/processing_functions/read_excel_all_sheets.R"))
source(paste0(wd, "/code/processing_functions/Eval_ReactionExpression.R"))
source(paste0(wd, "/code/processing_functions/findReactEq.R"))

## SA processing methods
source(paste0(wd, "/SA_FBA/code/functions/class_generation.R"))
source(paste0(wd, "/SA_FBA/code/functions/FBAgreatmodeClass.R"))
source(paste0(wd, "/SA_FBA/code/functions/readMat.R"))
source(paste0(wd, "/SA_FBA/code/functions/ParallelFBA_sens.R"))

# computing flux distributions stats
source(paste0(wd, "/code/supplementary_functions/FluxDistrStats.R"))

# setting GEMM tags
model.type = "Clostridium"
model.name = "CD196HemeSink"

# setting results directory
result_dir = paste0(wd, "/input/", model.type, "/", model.name, "/global_sens_analysis")

# setting analysis FBA files names
mat_file = paste0(wd, "/input/", model.type, "/", model.name, "/", model.name, ".mat")
model_file = paste0(wd, "/input/", model.type, "/", model.name, "/", model.name, ".RData")
all_react_file = paste0(wd, "/input/", model.type, "/", model.name, "/all_react.rds")

# generate S3 R object FBA_model
model = FBA4Greatmod.generation(fba_mat = mat_file)

files = c("geneAssociation.rds", 
          "geni.rds", "officialName.rds", 
          "subsystem.rds", "genesFromGeneAss.rds",
          "met_KEGGID.rds", "rxn_KEGGID.rds")

for(f in files) {
  system(paste0("cd ", wd, " && ", "mv ", f, " ", wd, "/input/", model.type, "/", model.name, "/", f))
}

# saving RData FBA model
save(model, file = model_file)

# generate FBA model metadata and saving FBA model reaction data
bigg.path = paste0(wd, "/input/", model.type, "/", model.name, "/BIGGdata_", model.name, "_react.tsv")

metadataEQ(model = model, 
           model.name = model.name,
           prefix = paste("BIGGdata", model.name, sep = "_"),
           suffix = "tsv",
           extMetFlag = "b", 
           fielddelim = "\t", 
           entrydelim = ", ", 
           makeClosedNetwork = FALSE, 
           onlyReactionList = FALSE, 
           minimalSet = TRUE, 
           fpath = paste0(wd, "/input/", model.type, "/", model.name))

all_react = FBAmodel.metadata(model = model, 
                              model.name = model.name, 
                              wd = wd, 
                              bigg.path = bigg.path)

saveRDS(all_react, all_react_file)

# saving biochemical reaction equations
BIGGdata = read.delim2(bigg.path)

saveRDS(BIGGdata$equation, file = paste0(wd, "/input/", model.type, "/", model.name, "/equation.rds"))

# We will investigate the sensitivty of the reaction objective ("biomass205") 
# flux solution by varying the D parameters within the bounds (-10.0, 0.0) 
# defined by the problem shown below.
# 
# I’m going to use saltelli to generate a collection of parameter values to use for FBA. 
# Then sobol to analyze the results and generate the sensitivity analysis. 
# To use sample.saltelli, you have to provide the problem definition in a specific format.

# moving files
system(paste0("cd ", wd, "/input/", model.type, "/", model.name, " && ", 
              "cp", " geneAssociation.rds ", " all_react.rds ", " officialName.rds ", 
              " subsystem.rds ", " equation.rds ", model.name, ".RData", " ", wd, "/SA_FBA/input", " && ", 
              "mv ", wd, "/SA_FBA/input/", model.name, ".RData", " ", wd, "/SA_FBA/input/FBAmodel.RData"))

# script calling
source(paste0(wd, "/SA_FBA/code/SA_FBA_R.R"))

# moving files
system(paste0("cd ", wd, "/SA_FBA", " && ", 
              "rm -r ", wd, "/input/", model.type, "/", model.name, "/global_sens_analysis", " && ",
              "cp -R ./res ", wd, "/input/", model.type, "/", model.name, "/global_sens_analysis"))

# freeing up memory
gc()

param_values = readRDS(paste0(result_dir, "/param_values.rds"))
df_ind = readRDS(paste0(result_dir, "/df_ind.rds"))

df_ind_filtered = dplyr::filter(df_ind, sol == "Y" & value > 0.01)
df_ind_filtered$param = as.integer(df_ind_filtered$param)

# model annotation
all_react = readRDS(all_react_file)
br = all_react[[2]]
all_react = all_react[[1]]

br$param = seq(br$react.id)

p = left_join(df_ind_filtered, br[, c(2, 6)], by = "param")
p = p[order(p$value, decreasing = TRUE), ]

p_f = ggplot(dplyr::filter(p, order == "first"),
             aes(x = reorder(react.id, -value), y = value, fill = react.id)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.1, position = position_dodge(0.5)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(p$react.id)))) +
  labs(x = "", y = "First-order sensitivity index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, family = "Helvetica"),
        axis.text.y = element_text(face = "plain", size = 10, family = "Helvetica"),
        axis.title.x = element_text(face = "plain", size = 10, family = "Helvetica"),
        axis.title.y = element_text(face = "bold", size = 9, family = "Helvetica"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(size = 0.35)) 

p_t = ggplot(dplyr::filter(p, order == "total"), 
             aes(x = reorder(react.id, -value), y = value, fill = react.id)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.1, position = position_dodge(0.5)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(p$react.id)))) +
  labs(x = "", y = "Total-order sensitivity index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, family = "Helvetica"),
        axis.text.y = element_text(face = "plain", size = 10, family = "Helvetica"),
        axis.title.x = element_text(face = "plain", size = 10, family = "Helvetica"),
        axis.title.y = element_text(face = "bold", size = 9, family = "Helvetica"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(size = 0.35))

ggsave((p_f / p_t), file = paste0(wd, "/SA_C_diff.pdf"), width = 7, height = 6)
system(paste0("cp ", wd, "/SA_C_diff.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/SA_C_diff.pdf ", wd, "/results/images4paper"))

png(paste0(result_dir, "/key_r_indices.png"), height = 25*length(p$value), width = 650)
gridExtra::grid.arrange(gridExtra::tableGrob(p))
dev.off()

# plotting indices

# New facet label names for supp variable
supp.labs = c("First-order indices", "Total-order indices")
names(supp.labs) = c("first", "total")

a = ggplot(dplyr::filter(df_ind, sol == "Y"), 
           aes(x = "", y = value, fill = order, color = order)) + 
  geom_boxplot(alpha = 0.2, size = 0.25) +
  coord_flip() +
  labs(title = paste0("Sensitivity indices evaluation | ", "config = " , nrow(param_values)), 
       subtitle = paste0("number of parameters : ", length(unique(df_ind$param))),
       x = "Output Y (biomass)", 
       y = "") +
  theme(legend.position="none",
        axis.ticks = element_line(size = 0.25),
        strip.background = element_blank(),
        line =  element_line(size = 0.75, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        axis.title = element_text(size = 6, face = "bold"),
        strip.text.x = element_text(size = 6, face = "bold", color = "#2A475E", hjust = 0.5),
        text = element_text(size = 8,  color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
        plot.title = element_text(size = 8, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 7, face = "bold", color = "#1b2838"),
        plot.title.position = "plot") +
  facet_wrap(~order, labeller = labeller(order = supp.labs)) +
  scale_fill_manual(values=c('blue', 'purple')) +
  scale_color_manual(values=c('black', 'black'))

b = ggplot(dplyr::filter(df_ind, sol == "Y"), 
           aes(x = value, fill = order, color = order)) + 
  geom_histogram(binwidth = 0.05, alpha = 0.2, size = 0.25) +
  labs(x = "Index Value", 
       y = "Count") +
  theme(legend.position="none",
        axis.ticks = element_line(size = 0.25),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        line =  element_line(size = 0.75, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        axis.title = element_text(size = 6, face = "bold"),
        text = element_text(size = 8,  color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
        plot.title = element_text(size = 8, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 6, face = "bold", color = "#1b2838"),
        plot.title.position = "plot") +
  facet_wrap(~order) + 
  scale_fill_manual(values=c('blue', 'purple')) +
  scale_color_manual(values=c('#1b2838', '#1b2838'))

ggsave(a/b, file = paste0(result_dir, "/indices.png"), width = 3, height = 2.5)
ggsave(a/b, file = paste0(result_dir, "/indices.pdf"), width = 3, height = 2.5)

mat_list = split(param_values,
                 rep(1:(nrow(param_values) %/% chunks + 1),
                     each = chunks, length.out = nrow(param_values)))

for(i in 1:length(mat_list)) {

  bounds_set = readRDS(paste0(result_dir, "/bounds_set", "_", i, ".rds"))

  bounds_set_lb = bounds_set[, 1:(ncol(bounds_set)/2)]
  bounds_set_ub = bounds_set[, ((ncol(bounds_set)/2) + 1):ncol(bounds_set)]
  colnames(bounds_set_lb) = NULL
  colnames(bounds_set_ub) = NULL

  fbasol = readRDS(paste0(result_dir, "/fbasol", "_", i, ".rds"))
  colnames(fbasol) = NULL

  if(i == 1) {

    bounds_set_lb_t = bounds_set_lb
    bounds_set_ub_t = bounds_set_ub
    fbasol_t = fbasol

  } else {

    bounds_set_lb_t = cbind(bounds_set_lb_t, bounds_set_lb)
    bounds_set_ub_t = cbind(bounds_set_ub_t, bounds_set_ub)
    fbasol_t = cbind(fbasol_t, fbasol[, -1])

  }
}

colnames(bounds_set_lb_t) = paste("lb", 1:ncol(bounds_set_lb_t), sep = ".")
colnames(bounds_set_ub_t) = paste("ub", 1:ncol(bounds_set_ub_t), sep = ".")
saveRDS(cbind(bounds_set_lb_t, bounds_set_ub_t), bounds_set_file)

colnames(fbasol_t) = c("Reaction", paste("config", 2:ncol(fbasol_t), sep = "."))
saveRDS(fbasol_t, paste0(result_dir, "/fbasol.rds"))

# list_files = list.files(result_dir)
# 
# file.remove(paste0(result_dir, "/", list_files[grepl("bounds_set", list_files)]))
# v_fbasol = paste0(result_dir, "/", list_files[grepl("fbasol", list_files)])
# file.remove(v_fbasol[-length(v_fbasol)])

if( identical(which(all_react$React_ID == "sink_pheme_c"), integer(0)) ) {
  
  targetrnx = c("EX_trp_L_e", "EX_pro_L_e", 
                "EX_val_L_e", "EX_ile_L_e", 
                "EX_cys_L_e", "EX_leu_L_e",
                "EX_lys_L_e",  "EX_biomass_e")
  
} else {
  
  targetrnx = c("EX_trp_L_e", "EX_pro_L_e", 
                "EX_val_L_e", "EX_ile_L_e", 
                "EX_cys_L_e", "EX_leu_L_e", 
                "EX_biomass_e", "sink_pheme_c")
  
}

### Uncertainty of flux boundaries on a FBA model affects its objective

fbasol = readRDS(paste0(result_dir, "/fbasol.rds"))
fbasol_t = readRDS(paste0(result_dir, "/fbasol_template.rds"))

subflux = data.frame(Reaction = rep(fbasol[, 1], (dim(fbasol)[2] - 1)),
                     Flux = gather(as.data.frame(fbasol), key, 
                                   value, -c(Reaction))$value,
                     config = str_remove(gather(as.data.frame(fbasol), key, 
                                                value, -c(Reaction))$key, "config."))

subflux = cbind(subflux, 
                data.frame(lapply(param_values[, which(br$react.id %in% targetrnx)], 
                                  rep, each = dim(fbasol)[1])))

subflux$Flux = as.double(subflux$Flux)

df_plot = data.frame(ass.param = paste0("P.", which(br$react.id %in% targetrnx)),
                     react.index = which(br$react.id %in% targetrnx),
                     Reaction = br$react.id[which(br$react.id %in% targetrnx)],
                     col.data = seq(colnames(subflux[4:dim(subflux)[2]])) + 3,
                     col = brewer.pal(length(targetrnx),"Spectral"))

subflux_t = cbind(data.frame(fbasol_t$abbreviation,
                             fbasol_t$flux,
                             rep(0, dim(fbasol)[1])),
                  data.frame(lapply(fbasol_t$lowbnd[df_plot$react.index],
                                    rep, each = dim(fbasol)[1])))

colnames(subflux_t) = colnames(subflux)

subflux = rbind(subflux, subflux_t)

l = "#cc002d"; m = "#fffdbf"; h = "#112764"
  c = "#653592"; pn = "darkred"
    
  title.size = 5; subtitle.size = 5; text.size = 5; axistext.size = 4; 
  axistitle.size = 5; legendsize = 0.25; point.size = 0.5
  
  PN.plots.violin = list()
  PN.plots.count = list()
  
  names(df_ind)[3] = "react.index"
  df_ind$react.index = as.integer(df_ind$react.index)
  
  # PN-related reactions indices
  p = left_join(
    dplyr::filter(df_ind, sol == "Y" & 
                    parameter %in% df_plot$ass.param)[, c(2, 3, 5)],
    df_plot, by = "react.index")[, c(1, 3, 4, 5)]
  
  p = p[order(p$value, decreasing = TRUE), ]
  
  png(paste0(result_dir, "/PN_r_indices.png"), height = 25*length(p$value), width = 350)
  gridExtra::grid.arrange(gridExtra::tableGrob(p))
  dev.off()
  
  for (i in 1:length(df_plot$Reaction)) {
    
    tmp = subflux %>% dplyr::filter(Reaction == df_plot$Reaction[i])
    tmp = dplyr::filter(tmp, config != "0")
    
    PN.plots.violin[[i]] = ggplot(dplyr::filter(tmp, config != "0"), 
                                  aes(x = Reaction, y = Flux)) +
      labs(x = df_plot$Reaction[i], y = "Flux (mmol/gDW*h)",
           title = paste0(df_plot$Reaction[i]),
           subtitle = paste0("number of reactions varied: ", 
                             length(param_values)),
           colour = names(tmp)[i + 3]) +
      theme(text = element_text(size = text.size, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
            plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
            plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
            plot.title.position = "plot",
            legend.key.size = unit(legendsize, 'cm'),
            axis.text = element_text(size = axistext.size, color = "black"),
            axis.title = element_text(size = axistitle.size, face = "bold")) +
      scale_color_gradient2(
        midpoint = mean(dplyr::filter(tmp, config != "0")[, names(tmp)[i + 3]]),
        low = l, 
        mid = m, 
        high = h) +
      geom_boxplot(alpha = 0.5) +
      geom_violin(scale = "width", alpha = 0.3, width = 1) +
      scale_x_discrete(name = "", 
                       label = rep( paste0(paste0(tmp$Reaction, "\n"),
                                           "(n = ", (dim(param_values)[1]), ")"),
                                    length(tmp$Reaction)))
    
    PN.plots.count[[i]] = ggplot(tmp, 
                                 aes(x = Flux, fill = ordered(Reaction, levels = unique(Reaction)))) + 
      geom_histogram(aes(y = ..count..), alpha = 0.5, bins = 25, fill = df_plot$col[i]) +
      labs(x = "Flux (mmol/gDW*h)", y = "Count",
           subtitle = paste0("number of reactions varied: ", length(param_values)),
           title = paste0("Fluxes distribution | ", df_plot$Reaction[i])) +
      guides(fill = guide_legend(title = "Condition block")) +
      theme(text = element_text(size = 8, color = "black"),
            panel.background = element_blank(),
            legend.position = "none",
            panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
            plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
            plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
            plot.title.position = "plot", 
            axis.text = element_text(size = axistext.size, color = "black"),
            axis.title = element_text(size = axistitle.size, face = "bold"))
    
  }
  
  names(PN.plots.violin) = df_plot$Reaction
  names(PN.plots.count) = df_plot$Reaction
  
  ggsave(((PN.plots.violin[[1]] | PN.plots.violin[[2]] | 
             PN.plots.violin[[3]] | PN.plots.violin[[4]]) / 
            (PN.plots.violin[[5]] | PN.plots.violin[[6]] | 
               PN.plots.violin[[7]] | PN.plots.violin[[8]])),
         file = paste0(result_dir, "/FBA_Sens1.png"), 
         width = 5, height = 3.5)
  
  ggsave((PN.plots.count[[1]] | PN.plots.count[[2]] | 
            PN.plots.count[[3]] | PN.plots.count[[4]]) / 
           (PN.plots.count[[5]] | PN.plots.count[[6]] | 
              PN.plots.count[[7]] | PN.plots.count[[8]]),
         file = paste0(result_dir, "/FBA_Sens2.png"), 
         width = 6, height = 3)
  
  ggsave(((PN.plots.violin[[1]] | PN.plots.violin[[2]] | 
             PN.plots.violin[[3]] | PN.plots.violin[[4]]) / 
            (PN.plots.violin[[5]] | PN.plots.violin[[6]] | 
               PN.plots.violin[[7]] | PN.plots.violin[[8]])),
         file = paste0(result_dir, "/FBA_Sens1.pdf"), 
         width = 5, height = 3)
  
  ggsave((PN.plots.count[[1]] | PN.plots.count[[2]] | 
            PN.plots.count[[3]] | PN.plots.count[[4]]) / 
           (PN.plots.count[[5]] | PN.plots.count[[6]] | 
              PN.plots.count[[7]] | PN.plots.count[[8]]),
         file = paste0(result_dir, "/FBA_Sens2.pdf"), 
         width = 6, height = 3.5)
  
  subflux = dplyr::filter(subflux, config != "0")[, c(1, 2, 3)]
  
  subflux.boun = 
    dplyr::filter(subflux, subflux$Reaction %in% 
                    dplyr::filter(all_react, ReactionDistrict == "boundary")$React_ID)
  
  dall_wrap = ggplot(subflux.boun, 
                     aes(x = Flux, fill = ordered(Reaction, levels = unique(Reaction)))) + 
    geom_histogram(aes(y = ..count..), alpha = 0.25, 
                   position = "identity", bin = 30,
                   binwidth = 0.05, col = c, fill = c) +
    labs(x = "Flux (mmol/gDW*h)", y = "Count") +
    guides(fill = guide_legend(title = "Condition block")) +
    theme(text = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold")) +
    facet_wrap(~ Reaction, scales = "free")
  
  ggsave(dall_wrap, 
         file = paste(result_dir, "/FBA_Sens_allreact.pdf", sep = ""), 
         width = 15, height = 15, limitsize = FALSE)
  
  gc()
  
  distr_stats(result_dir,
              paste0(result_dir, "/fbasol.rds"),
              file.standard_dev = paste0(result_dir, "/data.RData"),
              thrsd = 10,
              linethrsd = 100)
  