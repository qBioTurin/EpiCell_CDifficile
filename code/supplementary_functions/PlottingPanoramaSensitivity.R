
setwd("~/EpiCell_CDifficile")

tag = "Unified"
Condition = "Therapy"
Exper = "Model_Sensitivity"

net_name = "EpitCellDifficileHemeSink"

param_target = c("IECsDeath", "Death4Treat", "Detox")

numb = 1000
multi_param = T
alph = 0.25

Tempi = c(0, 16, 32, 48, 64, 72)

colo1 <- c("black", "magenta", "gold")
colo2 <- c("black", "#266867", "yellow")
colo3 <- c("#ffd166", "#ee6c4d", "#293241")

if (multi_param) {
  setwd(paste("./Results/CDiff", tag, Condition, Exper, sep = "_"))
} else {
  setwd(paste("./Results/CDiff_", tag, "_", 
              Condition, "_", Exper, sep = ""))
}

load(paste(net_name, "-analysis.RData", sep = ""))

if (multi_param) {
  id_param = c()
  for (p in 1:length(param_target)) {
    for(i in 1:length(config)) {
      if(param_target[p] == config[[i]][[1]][[1]]) {
        id_param[p] = i
      }
    }
  }
} else {
  for(i in 1:length(config)) {
    if(param_target == config[[i]][[1]][[1]]) {
      id_param = i
    }
  }
}

load(paste("~/EpiCell_CDifficile/Results/CDiff_", tag, 
           "_Therapy_Model_Sensitivity/subtrace_", tag, 
           "Therapy.RData", sep = ""))

Cdiff = subtrace[which(subtrace$Places == "CD"), ]
Drug = subtrace[which(subtrace$Places == "Drug"), ]
Iecs = subtrace[which(subtrace$Places == "IECs"), ]
Biomass = filter(subtrace, subtrace$Places == "BiomassCD")

pheme_e = filter(subtrace, subtrace$Places == "pheme_e")
pheme_c = filter(subtrace, subtrace$Places == "pheme_c")

trp = filter(subtrace, subtrace$Places == "trp_L_e")
pro = filter(subtrace, subtrace$Places == "pro_L_e")
val = filter(subtrace, subtrace$Places == "val_L_e")
ile = filter(subtrace, subtrace$Places == "ile_L_e")
cys = filter(subtrace, subtrace$Places == "cys_L_e")
leu = filter(subtrace, subtrace$Places == "leu_L_e")

subtrace = filter(subtrace, Time %in% Tempi)
AA = filter(subtrace, subtrace$Places %in% 
              c("trp_L_e", "pro_L_e", "val_L_e", 
                "ile_L_e", "cys_L_e", "leu_L_e"))

AA$Time = as.character(AA$Time)

if (multi_param) {
  
  pC1 = ggplot(Cdiff, aes(x = Time, y = Marking, color = IECsDeath)) + 
    xlab("Time (hour)") + ylab("C. difficile (cell)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("CD Place | target param: ", param_target[1], sep = "")) +
    scale_colour_gradientn(name = "", 
                           colors = colo1,
                           oob = scales::squish)
  
  pC2 = ggplot(Cdiff, aes(x = Time, y = Marking, color = Death4Treat)) + 
    xlab("Time (hour)") + ylab("C. difficile (cell)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("CD Place | target param: ", param_target[2], sep = "")) +
    scale_colour_gradientn(name = "", 
                           colors = colo2,
                           oob = scales::squish)
  
  pC3 = ggplot(Cdiff, aes(x = Time, y = Marking, color = Detox)) + 
    xlab("Time (hour)") + ylab("C. difficile (cell)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("CD Place | target param: ", param_target[3], sep = "")) +
    scale_colour_gradientn(name = "", 
                           colors = colo3,
                           oob = scales::squish)
  
  ppheme1 = ggplot(pheme_c, aes(x = Time, y = Marking, color = IECsDeath)) + 
    xlab("Time (hour)") + ylab("pheme_c (pmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("pheme_c Place | param: ", param_target[1], sep = "")) +
    scale_colour_gradientn(name = "", 
                           colors = colo1,
                           oob = scales::squish)
  
  ppheme2 = ggplot(pheme_c, aes(x = Time, y = Marking, color = Death4Treat)) + 
    xlab("Time (hour)") + ylab("pheme_c (pmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("pheme_c Place | param: ", param_target[2], sep = "")) +
    scale_colour_gradientn(name = "", 
                           colors = colo2,
                           oob = scales::squish)
  
  ppheme3 = ggplot(pheme_c, aes(x = Time, y = Marking, color = Detox)) + 
    xlab("Time (hour)") + ylab("pheme_c (pmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("pheme_c Place | param: ", param_target[3], sep = "")) +
    scale_colour_gradientn(name = "", 
                           colors = colo3,
                           oob = scales::squish)
  
  
  ppheme = (ppheme1 + ppheme2 + ppheme3)
  
  pC = (pC1 | pC2 | pC3)
  
  pBiom1 = ggplot(Biomass, aes(x = Time, y = Marking, color = IECsDeath)) + 
    geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 1.5, label = "biomass_max", vjust = 0.5, color = "black", size = 2.5) +
    geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 0.6, label = "biomass_mean", vjust = 0.5, color = "black", size = 2.5) +
    geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 0.02, label = "biomass_min", vjust = 0.5, color = "black", size = 2.5) + 
    xlab("Time (hour)") + ylab("Biomass (pg)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("Biomass Place | target param: ", param_target[1], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pBiom2 = ggplot(Biomass, aes(x = Time, y = Marking, color = Death4Treat)) + 
    geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 1.5, label = "biomass_max", vjust = 0.5, color = "black", size = 2.5) +
    geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 0.6, label = "biomass_mean", vjust = 0.5, color = "black", size = 2.5) +
    geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 0.02, label = "biomass_min", vjust = 0.5, color = "black", size = 2.5) + 
    xlab("Time (hour)") + ylab("Biomass (pg)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("Biomass Place | target param: ", param_target[2], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pBiom3 = ggplot(Biomass, aes(x = Time, y = Marking, color = Detox)) + 
    geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 1.5, label = "biomass_max", vjust = 0.5, color = "black", size = 2.5) +
    geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 0.6, label = "biomass_mean", vjust = 0.5, color = "black", size = 2.5) +
    geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 0.02, label = "biomass_min", vjust = 0.5, color = "black", size = 2.5) + 
    xlab("Time (hour)") + ylab("Biomass (pg)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("Biomass Place | target param: ", param_target[3], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pBiom = (pBiom1 | pBiom2 | pBiom3)

  aa = filter(subtrace, subtrace$Places %in% 
                c("trp_L_e", "pro_L_e", "val_L_e",
                  "ile_L_e", "cys_L_e", "leu_L_e"))
  
  paa1 = ggplot(aa, aes(x = Time, y = Marking, color = IECsDeath)) +
    xlab("Time (h)") + ylab("AAs (mmol) (mmol)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle("Sensitivity Analysis: AAs (mmol)") +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  paa2 = ggplot(aa, aes(x = Time, y = Marking, color = Death4Treat)) +
    xlab("Time (h)") + ylab("AAs (mmol) (mmol)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle("Sensitivity Analysis: AAs (mmol)") +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  paa3 = ggplot(aa, aes(x = Time, y = Marking, color = Detox)) +
    xlab("Time (h)") + ylab("AAs (mmol) (mmol)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle("Sensitivity Analysis: AAs (mmol)") +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  paa = (paa1 | paa2 | paa3)
  
  pI1 = ggplot(Iecs, aes(x = Time, y = Marking, color = IECsDeath)) + 
    xlab("Time (hour)") + ylab("IECs (cell)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("IECs Place | target param: ", param_target[1], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pI2 = ggplot(Iecs, aes(x = Time, y = Marking, color = Death4Treat)) + 
    xlab("Time (hour)") + ylab("IECs (cell)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("IECs Place | target param: ", param_target[2], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pI3 = ggplot(Iecs, aes(x = Time, y = Marking, color = Detox)) + 
    xlab("Time (hour)") + ylab("IECs (cell)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("IECs Place | target param: ", param_target[3], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pI = (pI1 | pI2 | pI3)
  
  pdrug1 = ggplot(Drug, aes(x = Time, y = Marking*1e-09, color = IECsDeath)) +
    geom_hline(yintercept = 5.8e-06, linetype = "dashed", color = "darkgrey", lwd = 0.5) +
    annotate("text", 4, 6.3e-06, label = "MIC value", vjust = 0.5, color = "black", size = 2) + 
    xlab("Time (hour)") + ylab("Drug (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("Drug Place | target param: ", param_target[1], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pdrug2 = ggplot(Drug, aes(x = Time, y = Marking*1e-09, color = Death4Treat)) +
    geom_hline(yintercept = 5.8e-06, linetype = "dashed", color = "darkgrey", lwd = 0.5) +
    annotate("text", 4, 6.3e-06, label = "MIC value", vjust = 0.5, color = "black", size = 2) + 
    xlab("Time (hour)") + ylab("Drug (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("Drug Place | target param: ", param_target[2], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pdrug3 = ggplot(Drug, aes(x = Time, y = Marking*1e-09, color = Detox)) +
    geom_hline(yintercept = 5.8e-06, linetype = "dashed", color = "darkgrey", lwd = 0.5) +
    annotate("text", 4, 6.3e-06, label = "MIC value", vjust = 0.5, color = "black", size = 2) + 
    xlab("Time (hour)") + ylab("Drug (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("Drug Place | target param: ", param_target[3], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pdrug = (pdrug1 | pdrug2 | pdrug3)
  
  ptrp1 = ggplot(trp, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Trp_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("trp_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  ptrp2 = ggplot(trp, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Trp_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("trp_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  ptrp3 = ggplot(trp, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Trp_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("trp_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  ptrp = (ptrp1 | ptrp2 | ptrp3)
  
  ppro1 = ggplot(pro, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Pro_L (mmol)") +
    geom_line(aes(group = config)) + ggtitle("pro_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  ppro2 = ggplot(pro, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Pro_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("pro_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  ppro3 = ggplot(pro, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Pro_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("pro_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  ppro = (ppro1 | ppro2 | ppro3)
  
  pval1 = ggplot(val, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Val_L (mmol)") +
    geom_line(aes(group = config)) + ggtitle("val_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pval2 = ggplot(val, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Val_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("val_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  
  pval3 = ggplot(val, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Val_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("val_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pval = (pval1 | pval2 | pval3)
  
  pile1 = ggplot(ile, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Ile_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("ile_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pile2 = ggplot(ile, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Ile_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("ile_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pile3 = ggplot(ile, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Ile_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("leu_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pile = (pile1 | pile2 | pile3)
  
  pleu1 = ggplot(leu, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Leu_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("leu_L_e Place") +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pleu2 = ggplot(leu, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Leu_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("leu_L_e Place") +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pleu3 = ggplot(leu, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Leu_L (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("leu_L_e Place") +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pleu = (pleu1 | pleu2 | pleu3)
  
  pcys1 = ggplot(cys, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Cys_L (mmol)") + 
    geom_line(aes(group = config)) + 
    ggtitle("cys_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pcys2 = ggplot(cys, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Cys_L (mmol)") + 
    geom_line(aes(group = config)) + 
    ggtitle("cys_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pcys3 = ggplot(cys, aes(x = Time, (Marking/(6.022e20/6.022e08)), 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Cys_L (mmol)") + 
    geom_line(aes(group = config)) + 
    ggtitle("cys_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  
  pcys = (pcys1 | pcys2 | pcys3)
  
}

#################################################################################

pAA = ggplot(AA, aes(x = Time, y = Marking/(6.022e20/6.022e08))) +
         labs(x = "Time (h)", 
              y = "Aminoacids (mmol)",
              subtitle = "Aminoacids (mmol)",
              title = paste(tag, "Experiments ", sep = "")) +
         theme(text = element_text(size = 8, color = "black"),
               panel.background = element_blank(),
               panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
               legend.position = "none",
               plot.title = element_text(size = 12, face = "bold", color = "#2a475e"),
               plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold", color = "#1b2838"),
               plot.title.position = "plot", 
               axis.text = element_text(size = 5, color = "black"),
               axis.title = element_text(size = 12, face = "bold")) +
         geom_jitter(aes(col = Time), alpha = 0.35, size = 0.1) +
         geom_violin(scale = "width", alpha = 0.3, width = 0.3) +
         geom_boxplot(alpha = 0) +
         facet_wrap(~Places, nrow = 2) +
         ggplot2::scale_color_manual(values = brewer.pal(length(Tempi), "Set2")) +
         scale_x_discrete(name = "", label = rep(paste(paste(Tempi, "\n", sep = ""),
                                                       "(n = ", length(unique(AA$config)), ")",
                                                       sep = ""), length(Tempi)))

#################################################################################

setwd("~/EpiCell_CDifficile")

p = (pC/pBiom | pdrug/pI) / ((pcys/pval) + (pleu/ppro) | (pile/ptrp))

ggsave(p, file = paste("MultilinePlot_", tag, "_", Exper, "_", Condition, ".pdf", sep = ""), 
       width = 32, height = 18)

system(paste("cp ", paste("MultilinePlot_", tag, "_", Exper, "_", Condition, ".pdf ", sep = ""),
"./Results/images4paper/", sep = ""))
system(paste("mv ", paste("MultilinePlot_", tag, "_", Exper, "_", Condition, ".pdf ", sep = ""),
             "./Results/images4paper/", sep = ""))


ggsave(pAA, file = paste("AA.", tag, ".pdf", sep = ""), width = 7, height = 4)

system(paste("cp ",  paste("AA.", tag, ".pdf ", sep = ""), "./Results/images4paper/", sep = ""))
system(paste("mv ",  paste("AA.", tag, ".pdf ", sep = ""), "./Results/images4paper/", sep = ""))
