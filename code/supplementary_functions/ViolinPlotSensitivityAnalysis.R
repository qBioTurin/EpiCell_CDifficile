
# for debug
#
# setwd(wd)
# replicates = 5
# r = replicates*cores
# multi_param = T
# tag = "Unified"
# Condition = "Therapy"
# Exper = "Model_Sensitivity"
# param_target = c("IECsDeath", "Death4Treat", "Detox")
# multi_param = T
# net_fname = "EpitCellDifficileHemeSink"

Tempi = c(0, 16, 32, 48, 64, 72)

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

load(file = paste0("subtrace_" , tag, Condition, ".RData"))
load(file = paste0("subflux_" , tag, Condition, ".RData"))

subtrace = dplyr::filter(subtrace, Time %in% Tempi)

subtrace$config = stringr::str_remove(subtrace$config, paste0(net_fname, "-analysis-"))
subtrace$config = stringr::str_remove(subtrace$config, ".trace")

Cdiff = subtrace[which(subtrace$Places == "CD"), ]
Cdiff$Marking = as.double(Cdiff$Marking)

Drug = subtrace[which(subtrace$Places == "Drug"), ]
Drug$Marking = as.double(Drug$Marking)*1e-09

Iecs = subtrace[which(subtrace$Places == "IECs"), ]
Iecs$Marking = as.double(Iecs$Marking)

Biomass = subtrace[which(subtrace$Places == "BiomassCD"), ]
Biomass$Marking = as.double(Biomass$Marking)

pheme_e = dplyr::filter(subtrace, subtrace$Places == "pheme_e")
pheme_e$Marking = as.double(pheme_e$Marking)
pheme_c = dplyr::filter(subtrace, subtrace$Places == "pheme_c")
pheme_c$Marking = as.double(pheme_c$Marking)

trp = dplyr::filter(subtrace, subtrace$Places == "trp_L_e")
trp$Marking = as.double(trp$Marking)/((6.022e20/6.022e08)*1000)
pro = dplyr::filter(subtrace, subtrace$Places == "pro_L_e")
pro$Marking = as.double(pro$Marking)/((6.022e20/6.022e08)*1000)
val = dplyr::filter(subtrace, subtrace$Places == "val_L_e")
val$Marking = as.double(val$Marking)/((6.022e20/6.022e08)*1000)
ile = dplyr::filter(subtrace, subtrace$Places == "ile_L_e")
ile$Marking = as.double(ile$Marking)/((6.022e20/6.022e08)*1000)
cys = dplyr::filter(subtrace, subtrace$Places == "cys_L_e")
cys$Marking = as.double(cys$Marking)/((6.022e20/6.022e08)*1000)
leu = dplyr::filter(subtrace, subtrace$Places == "leu_L_e")
leu$Marking = as.double(leu$Marking)/((6.022e20/6.022e08)*1000)

setwd(paste0(wd, "/results/CDiff_", tag, "_", Condition, "_", Exper))

##### plotting ranking

V.Cdiff = ggstatsplot::ggbetweenstats(data = Cdiff, x = Time, y = Marking) +
  labs(x = "Time (h)", 
       y = "CD Marking (cell)",
       title = "Distribution of configurations | CD Place") +
  theme(text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
    plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
    plot.title.position = "plot", 
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
  )

Cdiff$Time = as.character(Cdiff$Time)

V2.Cdiff = Cdiff %>%
  ggplot(aes(x = Time, y = Marking)) +
  labs(x = "Time (h)", 
       y = "CD Marking (cell)",
       title = "Distribution of configurations | CD Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.Drug = ggstatsplot::ggbetweenstats(data = Drug, x = Time, y = Marking) +
  labs(x = "Time (h)", 
       y = "Drug Marking (mmol)",
       title = "Distribution of configurations | Drug Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

Drug$Time = as.character(Drug$Time)

V2.Drug = Drug  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", 
       y = "Drug Marking (mmol)",
       title = "Distribution of configurations | Drug Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.Iecs = ggstatsplot::ggbetweenstats(data = Iecs, x = Time, y = Marking) +
  labs(x = "Time (h)", 
       y = "IECs Marking (cell)",
       title = "Distribution of configurations | IECs Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

Iecs$Time = as.character(Iecs$Time)

V2.Iecs = Iecs  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", 
       y = "IECs Marking (cell)",
       title = "Distribution of configurations | IECs Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.Biomass = ggstatsplot::ggbetweenstats(data = Biomass, x = Time, y = Marking) +
  labs(x = "Time (h)", 
       y = "Biomass Marking (pg)",
       title = "Distribution of configurations | Biomass Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

Biomass$Time = as.character(Biomass$Time)

V2.Biomass = Biomass  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", 
       y = "Biomass Marking (pg)",
       title = "Distribution of configurations | Biomass Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.pheme_c = ggstatsplot::ggbetweenstats(data = pheme_c, x = Time, y = Marking) +
  labs(x = "Time (h)", 
       y = "intracell. heme Marking (pmol)",
       title = "Distribution of configurations | pheme_c Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

pheme_c$Time = as.character(pheme_c$Time)

V2.pheme_c = pheme_c  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", 
       y = "intracell. heme Marking (pmol)",
       title = "Distribution of configurations | pheme_c Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.pheme_e = ggstatsplot::ggbetweenstats(data = pheme_e, x = Time, y = Marking) +
  labs(x = "Time (h)", y = "extracell. heme Marking (pmol)",
       title = "Distribution of configurations | pheme_e Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

pheme_e$Time = as.character(pheme_e$Time)

V2.pheme_e = pheme_e  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", 
       y = "extracellular. heme Marking (pmol)",
       title = "Distribution of configurations | pheme_e Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

######

V.trp = ggstatsplot::ggbetweenstats(data = trp, x = Time, y = Marking) +
  labs(x = "Time (h)", y = "Trp_e Marking (mmol)",
       title = "Distribution of configurations | Trp Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

trp$Time = as.character(trp$Time)

V2.trp = trp  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", 
       y = "Trp_e Marking (mmol)",
       title = "Distribution of configurations | Trp Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.pro = ggstatsplot::ggbetweenstats(data = pro, x = Time, y = Marking) +
  labs(x = "Time (h)", y = "Pro_e Marking (mmol)",
       title = "Distribution of configurations | Pro Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

pro$Time = as.character(pro$Time)

V2.pro = pro  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", y = "Pro_e Marking (mmol)",
       title = "Distribution of configurations | Pro Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.val = ggstatsplot::ggbetweenstats(data = val, x = Time, y = Marking) +
  labs(x = "Time (h)", y = "Val_e Marking (mmol)",
       title = "Distribution of configurations | Val Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

val$Time = as.character(val$Time)

V2.val = val  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", y = "Val_e Marking (mmol)",
       title = "Distribution of configurations | Val Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.ile = ggstatsplot::ggbetweenstats(data = val, x = Time, y = Marking) +
  labs(x = "Time (h)", y = "Ile_e Marking (mmol)",
       title = "Distribution of configurations | Ile Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

ile$Time = as.character(ile$Time)

V2.ile = ile  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", y = "Ile_e Marking (mmol)",
       title = "Distribution of configurations | Ile Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.cys = ggstatsplot::ggbetweenstats(data = cys, x = Time, y = Marking) +
  labs(x = "Time (h)", y = "Cys_e Marking (mmol)",
       title = "Distribution of configurations | Cys Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

cys$Time = as.character(cys$Time)

V2.cys = cys  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", y = "Cys_e Marking (mmol)",
       title = "Distribution of configurations | Cys Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))

V.leu = ggstatsplot::ggbetweenstats(data = leu, x = Time, y = Marking) +
  labs(x = "Time (h)", y = "Leu_e Marking (mmol)",
       title = "Distribution of configurations | Leu Place") +
  theme(text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12)
  )

leu$Time = as.character(leu$Time)

V2.leu = leu  %>%
  ggplot(aes(x= Time, y = Marking)) +
  labs(x = "Time (h)", y = "Cys_e Marking (mmol)",
       title = "Distribution of configurations | Cys Place") +
  theme(text = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
        plot.subtitle = element_text(size = 9, face = "bold", color = "#1b2838"),
        plot.title.position = "plot", 
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  geom_jitter(aes(col = Time), alpha = 0.35) +
  geom_violin(scale = "width", alpha = 0.3, width = 1) +
  geom_boxplot(alpha = 0) +
  ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
  scale_x_discrete(name = "", 
                   label = rep(paste(paste(Tempi, "\n", sep = ""),
                                     "(n = ", length(unique(Cdiff$config)), ")", 
                                     sep = ""), length(Tempi)))
