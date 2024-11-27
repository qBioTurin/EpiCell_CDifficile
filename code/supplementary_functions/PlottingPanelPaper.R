
Exper = "Model_Sensitivity"
Condition = "Therapy"
tag = c("SuperAblated", "Ablated", "Unified")
param_target = c("IECsDeath", "Death4Treat", "Detox")

Tempi = c(0, 12, 24, 36, 48, 60)

coloTag = data.frame(SuperAblated = "darkred", Ablated = "#EC6D67", Unified = "#5351A2")
nameTag = data.frame(SuperAblated = "Ablated", Ablated = "Partially Ablated", Unified = "Unified")

colo1 <- c("black", "magenta", "gold")
colo2 <- c("black", "#266867", "yellow")
colo3 <- c("#ffd166", "#ee6c4d", "#293241")

for (j in tag) {
  
  setwd(paste0(wd, "/results/CDiff", "_", j, "_", Condition, "_", Exper))
  
  load(file = paste0("subtrace_" , j, Condition, ".RData"))
  subtrace = cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
  assign(paste0("subtrace", j, Condition), subtrace)
  
  load(file = paste0("subflux_" , j, Condition, ".RData"))
  subflux = cbind(subflux, Scenario = rep(j, length(subflux$Time)))
  assign(paste0("subflux", j, Condition), subflux)
  
}

subtrace = rbind(subtraceSuperAblatedTherapy,
                 subtraceAblatedTherapy, 
                 subtraceUnifiedTherapy)

subflux = rbind(subfluxAblatedTherapy, subfluxUnifiedTherapy)

AA = filter(subtrace, Time %in% Tempi)

AA = filter(AA, AA$Places %in% 
              c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e"))

AA$Time = as.character(AA$Time)

I = filter(subtrace, subtrace$Places == "IECs")
I = filter(I, Time %in% c(0, 4, 8, 12, 16, 20))
I$Time = as.character(I$Time)

s = 12; st = 12; at = 12

##### 

# setwd("~/EpiCell_CDifficile")
# multi_param = T
# tag = "SuperAblated"
# Condition = "Therapy"
# Exper = "Model_Sensitivity"
# param_target = c("IECsDeath", "Death4Treat", "Detox")
# net_name = "superablated"
# numb = 1000
# 
# if (multi_param) {
#   setwd(paste("./Results/CDiff", tag, Condition, Exper, sep = "_"))
# } else {
#   setwd(paste("./Results/CDiff_", tag, "_", 
#               Condition, "_", Exper, sep = ""))
# }
# 
# load(paste("./", net_name, "-analysis.RData", sep = ""))
# 
# if (multi_param) {
#   id_param = c()
#   for (p in 1:length(param_target)) {
#     for(i in 1:length(config)) {
#       if(param_target[p] == config[[i]][[1]][[1]]) {
#         id_param[p] = i
#       }
#     }
#   }
# } else {
#   for(i in 1:length(config)) {
#     if(param_target == config[[i]][[1]][[1]]) {
#       id_param = i
#     }
#   }
# }
# 
# failed.tr <<- c(); f.tr <<- 1
# list.Tr = list.files(pattern = "\\.trace$")
# 
# subtrace = lapply(sample(list.Tr, numb), function(i){
#   TracesPath = i
#   n = as.numeric(gsub(pattern = paste("(", net_name, "-analysis-)|(.trace)", sep = ""), 
#                       replacement = "", x = i))
#   trace = read.table(TracesPath, header = T)
#   if(!all(sapply(trace, is.numeric))) {
#     print(i)
#     failed.tr[f.tr] <<- n
#     f.tr <<- f.tr + 1
#     return()
#   }
#   subtrace = trace %>% 
#     dplyr::select(Time, IECs, Drug, Damage, pheme_e, 
#                   pheme_c, BiomassCD, CD, pro_L_e, leu_L_e, 
#                   ile_L_e, val_L_e, trp_L_e, cys_L_e) %>% 
#     tidyr::gather(key = "Places", value = "Marking", -Time)
#   subtrace = cbind(subtrace, config = rep(as.character(i), length(subtrace$Time)))
#   
#   for (p in 1:length(param_target)) {
#     subtrace = cbind(subtrace, p = rep(config[[id_param[p]]][[n]][[3]], length(subtrace$Time)))
#   }
#   
#   return(subtrace)
#   
# })
# subtrace = do.call(rbind, subtrace)
# 
# if(multi_param) {
#   names(subtrace) = c("Time", "Places", "Marking", "config", param_target)
# }
# 
# setwd(paste("~/EpiCell_CDifficile/Results/CDiff", tag, Condition, Exper, sep = "_"))
# save(subtrace, file = paste("subtrace_" , tag, Condition, ".RData", sep = ""))

j = "Ablated"

assign(paste("EX_biomass_e", param_target[3], j, sep = "."),
       filter(subflux[which(subflux$Places == "EX_biomass_e"), ], Scenario == j) %>%
         ggplot(aes(x = Time, y = Marking, color = Detox)) +
         labs(x = "Time (h)", 
              y = "EX_biomass(e) \n (mmol/gDW*h)",
              subtitle =  "EX_biomass(e) Flux",
              title = paste(nameTag[, j], " Experiments", sep = ""),
              colour = latex2exp::TeX(paste("Param: ", param_target[3], sep = ""))) + 
         theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
               plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
               plot.title.position = "plot",
               axis.text = element_text(size = 12, color = "black"),
               axis.title = element_text(size = 12, face = "bold"),
               legend.position = "none") + 
         geom_line(aes(group = config)) +
         scale_colour_gradientn(colors = colo3,
                                oob = scales::squish))

j = "Unified"

assign(paste("EX_biomass_e", param_target[3], j, sep = "."),
       filter(subflux[which(subflux$Places == "EX_biomass_e"), ], Scenario == j) %>%
         ggplot(aes(x = Time, y = Marking, color = Detox)) +
         labs(x = "Time (h)", 
              y = "EX_biomass(e) \n (mmol/gDW*h)",
              subtitle =  "EX_biomass(e) Flux",
              title = paste(nameTag[, j], " Experiments", sep = ""),
              colour = latex2exp::TeX(paste("Param: ", param_target[3], sep = ""))) + 
         theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
               plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
               plot.title.position = "plot",
               axis.text = element_text(size = 12, color = "black"),
               axis.title = element_text(size = 12, face = "bold")) + 
         geom_line(aes(group = config)) +
         scale_colour_gradientn(colors = colo3,
                                oob = scales::squish))

######

j = "SuperAblated"

assign(paste("CD", param_target[3], j, sep = "."), 
       filter(subtrace[which(subtrace$Places == "CD"), ], Scenario == j) %>%
         ggplot(aes(x = Time, y = Marking, color = Detox)) + 
         labs(x = "Time (h)", 
              y = "C. difficile (cell)",
              subtitle = "CD",
              title = paste(nameTag[, j], " Experiments", sep = ""),
              colour = paste("Param: ", param_target[3], sep = "")) + 
         theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
               plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
               plot.title.position = "plot",
               axis.text = element_text(size = at, color = "black"),
               axis.title = element_text(size = 12, face = "bold"),
               legend.position = "none") + 
         geom_line(aes(group = config)) +
         scale_colour_gradientn(colors = colo3,
                                oob = scales::squish))

j = "Ablated"

assign(paste("CD", param_target[3], j, sep = "."), 
       filter(subtrace[which(subtrace$Places == "CD"), ], Scenario == j) %>%
         ggplot(aes(x = Time, y = Marking, color = Detox)) + 
         labs(x = "Time (h)", 
              y = "C. difficile (cell)",
              subtitle = "CD",
              title = paste(nameTag[, j], " Experiments", sep = ""),
              colour = paste("Param: ", param_target[3], sep = "")) + 
         theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
               plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
               plot.title.position = "plot",
               axis.text = element_text(size = at, color = "black"),
               axis.title = element_text(size = 12, face = "bold"),
               legend.position = "none") + 
         geom_line(aes(group = config)) +
         scale_colour_gradientn(colors = colo3,
                                oob = scales::squish))

j = "Unified"

assign(paste("CD", param_target[3], j, sep = "."), 
       filter(subtrace[which(subtrace$Places == "CD"), ], Scenario == j) %>%
         ggplot(aes(x = Time, y = Marking, color = Detox)) + 
         labs(x = "Time (h)", 
              y = "C. difficile (cell)",
              subtitle = "CD",
              title = paste(nameTag[, j], " Experiments", sep = ""),
              colour = paste("Param: ", param_target[3], sep = "")) + 
         theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
               plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
               plot.title.position = "plot",
               axis.text = element_text(size = at, color = "black"),
               axis.title = element_text(size = 12, face = "bold")) + 
         geom_line(aes(group = config)) +
         scale_colour_gradientn(colors = colo3,
                                oob = scales::squish))

for (j in tag) {
  
  assign(paste("IECs", param_target[3], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "IECs"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking, color = Detox)) + 
           labs(x = "Time (h)", 
                y = "IECs (cell)",
                subtitle = "IECs",
                title = paste(nameTag[, j], " Experiments", sep = ""),
                colour = paste("Param: ", param_target[3], sep = "")) + 
           theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = at, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo3,
                                  oob = scales::squish))
  
  assign(paste("trp_L_e", param_target[1], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "trp_L_e"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking/(6.022e17/6.022e08), color = IECsDeath)) + 
           labs(x = "Time (h)", 
                y = latex2exp::TeX("\\textbf{Tryptophan ($\\mu$mol/mL)}"),
                subtitle = "trp_L_e",
                title = paste(nameTag[, j], " Experiments", sep = ""),
                colour = paste("Param: ", param_target[1], sep = "")) + 
           theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = at, color = "black"),
                 axis.title = element_text(size = 12, face = "bold"),
                 legend.position = "none") + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo1,
                                  oob = scales::squish))
  
  assign(paste("leu_L_e", param_target[1], j, sep = "."), 
         filter(subtrace[which(subtrace$Places == "leu_L_e"), ], Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking/(6.022e17/6.022e08), color = IECsDeath)) + 
           labs(x = "Time (h)", 
                y = latex2exp::TeX("\\textbf{Leucine ($\\mu$mol/mL)}"),
                subtitle = "leu_L_e",
                title = paste(nameTag[, j], " Experiments", sep = ""),
                colour = paste("Param: ", param_target[1], sep = "")) + 
           theme(plot.title = element_text(size = s, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = st, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot",
                 legend.position = "none",
                 axis.text = element_text(size = at, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) + 
           geom_line(aes(group = config)) +
           scale_colour_gradientn(colors = colo1,
                                  oob = scales::squish))
  
  assign(paste("leuViolin", param_target[1], j, sep = "."),
         filter(AA, AA$Places == "leu_L_e", Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking/(6.022e17/6.022e08))) +
           labs(x = "Time (h)", 
                y = latex2exp::TeX("\\textbf{Leucine ($\\mu$mol/mL)}"),
                subtitle = "leu_L_e",
                title = paste(nameTag[, j], " Experiments", sep = "")) +
           theme(text = element_text(size = 8, color = "black"),
                 panel.background = element_blank(),
                 panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
                 legend.position = "none",
                 plot.title = element_text(size = 12, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot", 
                 axis.text = element_text(size = at, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) +
           geom_jitter(aes(col = IECsDeath), alpha = 0.35, size = 0.1) +
           geom_violin(scale = "width", alpha = 0.3, width = 0.5) +
           geom_boxplot(alpha = 0) +
           scale_colour_gradientn(colors = colo1, oob = scales::squish))
  
  assign(paste("trpViolin", param_target[1], j, sep = "."),
         filter(AA, AA$Places == "trp_L_e", Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking/(6.022e17/6.022e08))) +
           labs(x = "Time (h)", 
                y = latex2exp::TeX("\\textbf{Tryptophan ($\\mu$mol/mL)}"),
                subtitle = "trp_L_e",
                title = paste(nameTag[, j], " Experiments", sep = "")) +
           theme(text = element_text(size = 8, color = "black"),
                 panel.background = element_blank(),
                 panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
                 legend.position = "none",
                 plot.title = element_text(size = 12, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot", 
                 axis.text = element_text(size = at, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) +
           geom_jitter(aes(col = IECsDeath), alpha = 0.35, size = 0.1) +
           geom_violin(scale = "width", alpha = 0.3, width = 0.5) +
           geom_boxplot(alpha = 0) +
           scale_colour_gradientn(colors = colo1, oob = scales::squish))
  
  assign(paste("IECsViolin", param_target[1], j, sep = "."),
         filter(I, Scenario == j) %>%
           ggplot(aes(x = as.numeric(Time), y = Marking, group = Time)) +
           labs(x = "Time (h)", 
                y = "IECs (cell)",
                subtitle = "Intestinal Epithelial Cells",
                title = paste(nameTag[, j], " Experiments", sep = "")) +
           theme(text = element_text(size = 8, color = "black"),
                 panel.background = element_blank(),
                 panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
                 legend.position = "none",
                 plot.title = element_text(size = 12, face = "bold", color = coloTag[, j]),
                 plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot", 
                 axis.text = element_text(size = at, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) +
           geom_jitter(aes(col = IECsDeath), alpha = 0.55, size = 1.5) +
           geom_violin(scale = "width", alpha = 0.5, width = 0.25) +
           geom_boxplot(alpha = 0) + scale_y_continuous(limits = c(300000, 480000)) +
           scale_colour_gradientn(colors = colo1, oob = scales::squish))
  # ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")))
  
}

#################################################################################

Exper = "Model_Analysis"
tag = c("Ablated", "Unified")
con = c("NoDrug", "Therapy")

lineCon = data.frame(Therapy = "solid", NoDrug = "dashed")

for (i in con) {
  for (j in tag) {
    
    subtrace = readr::read_csv(paste0(wd, "/results/CDiff_", j, "_", i, "_", Exper, "/SubTrace", j, i, ".csv"))[-1]
    subflux = readr::read_csv(paste0(wd, "/results/CDiff_", j, "_", i, "_", Exper, "/SubFlux", j, i, ".csv"))[-1]
    
    assign(paste0("subtrace", j, i), subtrace)
    assign(paste0("subflux", j, i), subflux)
    
  }
}

subtrace = rbind(subtraceAblatedNoDrug, subtraceAblatedTherapy,
                 subtraceUnifiedNoDrug, subtraceUnifiedTherapy)

subflux = rbind(subfluxAblatedNoDrug, subfluxAblatedTherapy,
                subfluxUnifiedNoDrug, subfluxUnifiedTherapy)

Places.all = ggplot(subtrace, aes(Time, Marking, colour = Scenario, linetype = Condition)) +
  theme_bw() + geom_line(linewidth = 0.5) +
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 0.5)) +
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.25) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.25) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.25) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.25) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.25) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.25) +
  facet_wrap(~ Places, scales="free", ncol = 3)

CD.all = ggplot(filter(subtrace, subtrace$Places == "CD"), 
                aes(Time, Marking, colour = Scenario, linetype = Condition)) + 
  theme_bw() + geom_line(size = 0.8) + 
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  labs(x = "Time (h)", 
       y = "C. difficile (cell)") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 1.75)) + 
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  scale_y_continuous(label = scientific_10)

pheme_e.all = ggplot(filter(subtrace, subtrace$Places == "pheme_e"), 
                aes(Time, Marking*1e-06, colour = Scenario, linetype = Condition)) + 
  theme_bw() + geom_line(size = 0.8) + 
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  labs(x = "Time (h)", 
       y = latex2exp::TeX("\\textbf{Extr. heam ($\\mu$mol/mL)}")) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 1.75)) + 
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_hline(yintercept = 1e+06*1e-06, linetype = "dashed", color = "darkred", lwd = 0.3) +
  annotate("text", 13, 0.95e+06*1e-06, label = "Hemorrhagic baseline", 
           vjust = 0.5, color = "black", size = 3.5)

pheme_c.all = ggplot(filter(subtrace, subtrace$Places == "pheme_c"), 
                     aes(Time, Marking*1e-06, colour = Scenario, linetype = Condition)) + 
  theme_bw() + geom_line(size = 0.8) + 
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  labs(x = "Time (h)",
       y = latex2exp::TeX("\\textbf{Intra. heam ($\\mu$mol/mL)}")) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 1.75)) + 
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.75)

EX_biomass.all = ggplot(filter(subflux[which(subflux$Reaction == "EX_biomass_e"), ]), 
                     aes(Time, Flux, colour = Scenario, linetype = Condition)) + 
  theme_bw() + geom_line(size = 0.8) + 
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  labs(x = "Time (h)",
       y = "EX_biomass(e) (mmol/gDW*h)") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 1.75)) + 
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.75)

EX_trp_L_e.all = ggplot(filter(subflux[which(subflux$Reaction == "EX_trp_L_e"), ]), 
                        aes(Time, Flux, colour = Scenario, linetype = Condition)) + 
  theme_bw() + geom_line(size = 0.8) + 
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  labs(x = "Time (h)",
       y = "EX_trp_L(e) (mmol/gDW*h)") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 1.75)) + 
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  scale_y_continuous(label = scientific_10)


drug.all = ggplot(filter(subtrace, subtrace$Places == "Drug"), 
                     aes(Time, Marking, colour = Scenario, linetype = Condition)) + 
  theme_bw() + geom_line(size = 0.8) + 
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  labs(x = "Time (h)",
       y = "Drug (pmol/mL)") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 1.75)) + 
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_hline(yintercept = 5800, linetype = "dashed", color = "darkgrey", lwd = 0.3) +
  annotate("text", 62, 6500, label = "MIC value", vjust = 0.5, color = "black", size = 3.5) +
  scale_y_continuous(label = scientific_10)

IECs.all = ggplot(filter(subtrace, subtrace$Places == "IECs"), 
                  aes(Time, Marking, colour = Scenario, linetype = Condition)) + 
  theme_bw() + geom_line(size = 0.8) + 
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  labs(x = "Time (h)",
       y = "IECs (cell)") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 1.75)) + 
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  scale_y_continuous(label = scientific_10)

biomass.all = ggplot(filter(subtrace, subtrace$Places == "BiomassCD"), 
                  aes(Time, Marking, colour = Scenario, linetype = Condition)) + 
  theme_bw() + geom_line(size = 0.8) + 
  scale_color_manual(values = c(coloTag[, "Ablated"], coloTag[, "Unified"])) +
  scale_linetype_manual(values = c(lineCon[, "NoDrug"], lineCon[, "Therapy"])) +
  labs(x = "Time (h)",
       y = "C. difficile average biomass (pg)") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(size = 1.75)) + 
  geom_vline(xintercept = 8, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 16, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 32, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 40, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_vline(xintercept = 48, linetype = "dotted", color = "darkblue", lwd = 0.75) +
  geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.3) +
  annotate("text", 62, 1.7, label = "biomass max", vjust = 0.5, color = "black", size = 3) +
  geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.3) +
  annotate("text", 62, 0.6, label = "biomass mean", vjust = 0.5, color = "black", size = 3) +
  geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.3) +
  annotate("text", 62, 0.18, label = "biomass min", vjust = 0.5, color = "black", size = 3)

##### 

ggsave((plot_spacer()|EX_biomass_e.Detox.Ablated|EX_biomass_e.Detox.Unified) /
         (CD.Detox.SuperAblated|CD.Detox.Ablated|CD.Detox.Unified),
  file = paste0(wd, "/FigSupp2.pdf"), width = 14, height = 5)

ggsave(((trpViolin.IECsDeath.SuperAblated|trp_L_e.IECsDeath.SuperAblated)/
          (trpViolin.IECsDeath.Ablated|trp_L_e.IECsDeath.Ablated)/
          (trpViolin.IECsDeath.Unified|trp_L_e.IECsDeath.Unified)) |
         ((leuViolin.IECsDeath.SuperAblated|leu_L_e.IECsDeath.SuperAblated)/
            (leuViolin.IECsDeath.Ablated|leu_L_e.IECsDeath.Ablated)/
            (leuViolin.IECsDeath.Unified|leu_L_e.IECsDeath.Unified)),
       file = paste0(wd, "/Figure3A.pdf"), width = 12, height = 8)

ggsave((IECsViolin.IECsDeath.SuperAblated|
          IECsViolin.IECsDeath.Ablated|
          IECsViolin.IECsDeath.Unified),
       file = paste0(wd, "/Figure3B.pdf"), width = 12, height = 6)

system(paste0("cp ", wd, "/FigSupp2.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/FigSupp2.pdf ", wd, "/results/images4paper"))

system(paste0("cp ", wd, "/Figure3A.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/Figure3A.pdf ", wd, "/results/images4paper"))

system(paste0("cp ", wd, "/Figure3B.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/Figure3B.pdf ", wd, "/results/images4paper"))

###### 

height = 3

ggsave(CD.all, file = paste0(wd, "/PaperCD.pdf"), width = 6, height = height)
system(paste0("cp ", wd, "/PaperCD.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/PaperCD.pdf ", wd, "/results/images4paper"))

ggsave(pheme_e.all, file = paste0(wd, "/PaperHeam.pdf"), width = 6, height = height)
system(paste0("cp ", wd, "/PaperHeam.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/PaperHeam.pdf ", wd, "/results/images4paper"))

ggsave(EX_biomass.all, file = paste0(wd, "/PaperEXBiom.pdf"), width = 6, height = height)
system(paste0("cp ", wd, "/PaperEXBiom.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/PaperEXBiom.pdf ", wd, "/results/images4paper"))

ggsave(drug.all, file = paste0(wd, "/PaperDrug.pdf"), width = 6, height = height)
system(paste0("cp ", wd, "/PaperDrug.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/PaperDrug.pdf ", wd, "/results/images4paper"))

ggsave(IECs.all, file = paste0(wd, "/PaperIECs.pdf"), width = 6, height = height)
system(paste0("cp ", wd, "/PaperIECs.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/PaperIECs.pdf ", wd, "/results/images4paper"))

ggsave(EX_trp_L_e.all, file = paste0(wd, "/PaperEX_trp.pdf"), width = 6, height = height)
system(paste0("cp ", wd, "/PaperEX_trp.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/PaperEX_trp.pdf ", wd, "/results/images4paper"))

ggsave(biomass.all, file = paste0(wd, "/PaperBiomass.pdf"), width = 6, height = height)
system(paste0("cp ", wd, "/PaperBiomass.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/PaperBiomass.pdf ", wd, "/results/images4paper"))

ggsave(Places.all, file = paste0(wd, "/Allplaces.pdf"), width = 12, height = 12)
system(paste0("cp ", wd, "/Allplaces.pdf ", wd, "/results/images4paper"))
system(paste0("mv ", wd, "/Allplaces.pdf ", wd, "/results/images4paper"))
