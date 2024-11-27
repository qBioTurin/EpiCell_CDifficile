
library(doParallel)
library(foreach)

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

alph = 0.25
numb = r

failed.tr <<- c()
f.tr <<- 1

if (multi_param) {
  setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper))
} else {
  setwd(paste90(wd, "/results/CDiff_", tag, "_", Condition, "_", Exper))
}

load(paste0(net_fname, "-analysis.RData"))

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

if (multi_param){
  
  registerDoParallel(cores = detectCores())
  subtrace <- foreach(i = sample(list.files(pattern = "\\.trace$"), numb), .combine = rbind) %dopar% {
    
    TracesPath <- i
    
    n <- as.numeric(gsub(pattern = paste0("(", net_fname, "-analysis-)|(.trace)"), replacement = "", x = i))
    
    trace <- read.table(TracesPath, header = TRUE)
    
    if (!all(sapply(trace, is.numeric))) {
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
    
    subtrace <- cbind(subtrace, config = rep(as.character(i), length(subtrace$Time)))
    
    for (p in 1:length(param_target)) {
      subtrace <- cbind(subtrace, p = rep(config[[id_param[p]]][[n]][[3]], length(subtrace$Time)))
    }
    
    return(subtrace)
  }
  stopImplicitCluster()
  
  registerDoParallel(cores = detectCores())
  subflux <- foreach(i = sample(list.files(pattern = "\\.flux$"), numb), .combine = rbind) %dopar% {
    TracesPath <- i
    n <- as.numeric(gsub(pattern = paste("(", net_fname, "-analysis-)|(-0.flux)", sep = ""), replacement = "", x = i))
    flux <- read.table(TracesPath, header = FALSE)
    ReactionsNames <- gsub("\\(|\\)", "", gsub("\\(", "_", flux[1, -c(1, 2)]))
    flux <- read.table(TracesPath, header = TRUE)
    if (TracesPath %in% paste0(net_fname, "-analysis-", failed.tr, "-0.flux")) {
      print(TracesPath)
      return(NULL)
    }
    colnames(flux) <- c("Time", "Obj_0", ReactionsNames)
    subflux <- flux[, c("Time", ReactionsNames)] %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time)
    obj <- flux[, c("Time", "Obj_0")]
    subbounds <- dplyr::filter(subflux, grepl("_Ub|_Lb", Reaction))
    obj_bounds <- dplyr::filter(subbounds, grepl("biomass205", Reaction))
    subflux <- dplyr::filter(subflux, !grepl("_Ub|_Lb", Reaction))
    subflux <- cbind(subflux, config = rep(as.character(i), length(subflux$Time)))
    for (p in 1:length(param_target)) {
      subflux <- cbind(subflux, p = rep(config[[id_param[p]]][[n]][[3]], length(subflux$Time)))
    }
    return(subflux)
  }
  stopImplicitCluster()
  
} else {
  
  subtrace <- lapply(sample(list.files(pattern = "\\.trace$"), numb), function(i) {
    TracesPath <- i
    n <- as.numeric(gsub(pattern = paste0("(", net_name, "-analysis-)|(.trace)"), replacement = "", x = i))
    trace <- read.table(TracesPath, header = TRUE)
    subtrace <- trace %>%
      dplyr::select(Time, IECs, Drug, Damage, 
                    pheme_e, pheme_c, BiomassCD, CD,
                    pro_L_e, leu_L_e, ile_L_e, val_L_e, trp_L_e, cys_L_e) %>%
      tidyr::gather(key = "Places", value = "Marking", -Time)
    subtrace <- cbind(subtrace, 
                      config = rep(as.character(i), length(subtrace$Time)), 
                      param_sen = rep(config[[id_param]][[n]][[3]], length(subtrace$Time)))
    
    return(subtrace)
  })
  
  subflux <- lapply(sample(list.files(pattern = "\\.flux$"), numb), function(i) {
    TracesPath <- i
    n <- as.numeric(gsub(pattern = paste("(", net_name, "-analysis-)|(-0.flux)", sep = ""), replacement = "", x = i))
    flux <- read.table(TracesPath, header = FALSE)
    ReactionsNames <- gsub("\\(|\\)", "", gsub("\\(", "_", flux[1, -c(1, 2)]))
    colnames(flux) <- c("Time", "Obj_0", ReactionsNames)
    subflux <- flux[, c("Time", ReactionsNames)] %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time)
    obj <- flux[, c("Time", "Obj_0")]
    subbounds <- dplyr::filter(subflux, grepl("_Ub|_Lb", Reaction))
    obj_bounds <- dplyr::filter(subbounds, grepl("biomass205", Reaction))
    subflux <- dplyr::filter(subflux, !grepl("_Ub|_Lb", Reaction))
    subflux <- cbind(subflux, 
                     config = rep(as.character(i), length(subflux$Time)), 
                     param_sen = rep(config[[id_param]][[n]][[3]], length(subflux$Time)))
    
    return(subflux)
  })
  
}

if(multi_param) {
  names(subtrace) = c("Time", "Places", "Marking", "config", param_target)
}

if(multi_param) {
  names(subflux) = c("Time", "Places", "Marking", "config", param_target)
}

setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper))

save(subtrace, file = paste0("subtrace_" , tag, Condition, ".RData"))
save(subflux, file = paste0("subflux_" , tag, Condition, ".RData"))

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

if (multi_param == F) {
  
  pC = ggplot(Cdiff, aes(x = Time, y = Marking, color = param_sen)) + 
    xlab("Time (hour)") + ylab("C. difficile (cell)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("CD Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "", colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)), oob = scales::squish)
  
  pBiom = ggplot(Biomass, aes(x = Time, y = Marking, color = param_sen)) + 
    geom_hline(yintercept = 1.57872, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 1.5, label = "biomass_max", vjust = 0.5, color = "black", size = 2.5) +
    geom_hline(yintercept = 0.496328, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 0.6, label = "biomass_mean", vjust = 0.5, color = "black", size = 2.5) +
    geom_hline(yintercept = 0.06738716, linetype = "dotted", color = "darkred", lwd = 0.5) +
    annotate("text", 5, 0.02, label = "biomass_min", vjust = 0.5, color = "black", size = 2.5) + 
    xlab("Time (hour)") + ylab("Biomass (pg)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("Biomass Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  pI = ggplot(Iecs, aes(x = Time, y = Marking, color = param_sen)) + 
    xlab("Time (hour)") + ylab("IECs (cell)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("IECs Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  pdrug = ggplot(Drug, aes(x = Time, y = Marking*1e-09, color = param_sen)) +
    geom_hline(yintercept = 5.8e-06, linetype = "dashed", color = "darkgrey", lwd = 0.5) +
    annotate("text", 4, 6.3e-06, label = "MIC value", vjust = 0.5, color = "black", size = 2) + 
    xlab("Time (hour)") + ylab("Drug (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("Drug Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  ptrp = ggplot(trp, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, color = param_sen)) +
    xlab("Time (hour)") + ylab("Trp_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("trp Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  ppro = ggplot(pro, aes(x = Time, (Marking/(6.022e20/6.022e08)), color = param_sen)) +
    xlab("Time (hour)") + ylab("Pro_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("pro_L_e Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  pval = ggplot(val, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, color = param_sen)) +
    xlab("Time (hour)") + ylab("Val_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("val_L_e Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  pile = ggplot(ile, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, color = param_sen)) +
    xlab("Time (hour)") + ylab("Ile_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("ile_L_e Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  pleu = ggplot(leu, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, color = param_sen)) +
    xlab("Time (hour)") + ylab("Leu_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("leu_L_e Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  pcys = ggplot(cys, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, color = param_sen)) +
    xlab("Time (hour)") + ylab("Cys_L (μmol)") + 
    ggtitle(paste("cys_L_e Place | target param: ", 
                  param_target, sep = "")) +
    geom_line(aes(group = config)) + 
    ggtitle("cys_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen), max(subtrace$param_sen)),
                           oob = scales::squish)
  
  pphemec = ggplot(pheme_c, aes(x = Time, y = Marking*1e-09, color = param_sen)) + 
    xlab("Time (hour)") + ylab("pheme_c (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("pheme_c Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  ppheme = ggplot(pheme_e, aes(x = Time, y = Marking*1e-09, color = param_sen)) + 
    xlab("Time (hour)") + ylab("pheme_e (mmol)") +
    geom_line(aes(group = config)) + 
    ggtitle(paste("pheme_e Place | target param: ", param_target, sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = c("black", "magenta", "gold"),
                           limits = c(min(subtrace$param_sen),
                                      max(subtrace$param_sen)),
                           oob = scales::squish)
  
  pflux = ggplot(subflux[which(subflux$Reaction %in% 
                                 c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                                   "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")), ],
                 aes(x = Time, y = Flux, color = param_sen)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") + 
    geom_line(aes(group = config)) + facet_wrap(~Reaction, scale = "free") + 
    ggtitle("Sensitivity Analysis: Fluxes (mmol/gDW*h)") +
    scale_colour_gradientn(name = "", colors = c("black", "magenta", "gold"),
                           limits = c(min(subflux$param_sen), max(subflux$param_sen)),
                           oob = scales::squish)
  
  pflux_special = ggplot(subflux[which(subflux$Reaction %in% 
                                         c("EX_biomass_e", "sink_pheme_c")), ],
                         aes(x = Time, y = Flux, color = param_sen)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config)) + facet_wrap(~Reaction, scale = "free") + 
    ggtitle("Sensitivity Analysis: Fluxes (mmol/gDW*h)") +
    scale_colour_gradientn(name = "", colors = c("black", "magenta", "gold"),
                           limits = c(min(subflux$param_sen), max(subflux$param_sen)),
                           oob = scales::squish)
  
}

if (multi_param) {
  
  colo1 <- c("black", "magenta", "gold")
  colo2 <- c("black", "#266867", "yellow")
  colo3 <- c("#ffd166", "#ee6c4d", "#293241")
  
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
  
  pC = (pC1 | pC2)
  
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
  
  pBiom = (pBiom1 | pBiom2)
  
  pflux1 = ggplot(subflux[which(subflux$Places %in% 
                                  c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                                    "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")), ],
                  aes(x = Time, y = Marking, color = IECsDeath)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
    ggtitle(paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[1], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pflux2 = ggplot(subflux[which(subflux$Places %in% 
                                  c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                                    "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")), ],
                  aes(x = Time, y = Marking, color = Death4Treat)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle(paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[2], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pflux3 = ggplot(subflux[which(subflux$Places %in% 
                                  c("EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                                    "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e")), ],
                  aes(x = Time, y = Marking, color = Detox)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle(paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[3], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pflux = (pflux1 | pflux2)
  
  pflux_special1 = ggplot(subflux[which(subflux$Places %in% 
                                          c("EX_biomass_e", "sink_pheme_c")), ],
                          aes(x = Time, y = Marking, color = IECsDeath)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle(paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[1], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pflux_special2 = ggplot(subflux[which(subflux$Places %in%
                                          c("EX_biomass_e", "sink_pheme_c")), ],
                          aes(x = Time, y = Marking, color = Death4Treat)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle(paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[2], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pflux_special3 = ggplot(subflux[which(subflux$Places %in%
                                          c("EX_biomass_e", "sink_pheme_c")), ],
                          aes(x = Time, y = Marking, color = Detox)) +
    xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle(paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[3], sep = "")) +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pflux_special = (pflux_special1 | pflux_special2)
  
  aa = filter(subtrace, subtrace$Places %in% 
                c("trp_L_e", "pro_L_e", "val_L_e",
                  "ile_L_e", "cys_L_e", "leu_L_e"))
  
  paa1 = ggplot(aa, aes(x = Time, y = Marking, color = IECsDeath)) +
    xlab("Time (h)") + ylab("AAs (mmol)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle("Sensitivity Analysis: AAs (mmol)") +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  paa2 = ggplot(aa, aes(x = Time, y = Marking, color = Death4Treat)) +
    xlab("Time (h)") + ylab("AAs (mmol)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle("Sensitivity Analysis: AAs (mmol)") +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  paa3 = ggplot(aa, aes(x = Time, y = Marking, color = Detox)) +
    xlab("Time (h)") + ylab("AAs (mmol)") +
    geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
    ggtitle("Sensitivity Analysis: AAs (mmol)") +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  paa = (paa1 | paa2)
  
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
  
  pI = (pI1 | pI2)
  
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
  
  pdrug = (pdrug1 | pdrug2)
  
  ptrp1 = ggplot(trp, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Trp_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("trp_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  ptrp2 = ggplot(trp, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Trp_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("trp_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  ptrp3 = ggplot(trp, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Trp_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("trp_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  ptrp = (ptrp1 | ptrp2)
  
  ppro1 = ggplot(pro, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Pro_L (μmol)") +
    geom_line(aes(group = config)) + ggtitle("pro_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  ppro2 = ggplot(pro, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Pro_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("pro_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  ppro3 = ggplot(pro, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Pro_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("pro_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  ppro = (ppro1 | ppro2)
  
  pval1 = ggplot(val, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Val_L (μmol)") +
    geom_line(aes(group = config)) + ggtitle("val_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pval2 = ggplot(val, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Val_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("val_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  
  pval3 = ggplot(val, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Val_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("val_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pval = (pval1 | pval2)
  
  pile1 = ggplot(ile, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Ile_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("ile_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pile2 = ggplot(ile, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Ile_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("ile_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pile3 = ggplot(ile, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Ile_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("leu_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pile = (pile1 | pile2)
  
  pleu1 = ggplot(leu, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Leu_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("leu_L_e Place") +
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pleu2 = ggplot(leu, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Leu_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("leu_L_e Place") +
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pleu3 = ggplot(leu, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Leu_L (μmol)") +
    geom_line(aes(group = config)) + 
    ggtitle("leu_L_e Place") +
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  pleu = (pleu1 | pleu2)
  
  pcys1 = ggplot(cys, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = IECsDeath)) +
    xlab("Time (hour)") + ylab("Cys_L (μmol)") + 
    geom_line(aes(group = config)) + 
    ggtitle("cys_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo1,
                           oob = scales::squish)
  
  pcys2 = ggplot(cys, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Death4Treat)) +
    xlab("Time (hour)") + ylab("Cys_L (μmol)") + 
    geom_line(aes(group = config)) + 
    ggtitle("cys_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo2,
                           oob = scales::squish)
  
  pcys3 = ggplot(cys, aes(x = Time, (Marking/(6.022e20/6.022e08))*1000, 
                          color = Detox)) +
    xlab("Time (hour)") + ylab("Cys_L (μmol)") + 
    geom_line(aes(group = config)) + 
    ggtitle("cys_L_e Place") + 
    scale_colour_gradientn(name = "",
                           colors = colo3,
                           oob = scales::squish)
  
  
  pcys = (pcys1 | pcys2)
  
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
pr <- ggplot(rank, aes(x = "ranking exper", y = measure)) +
  geom_boxplot(width = 0.4, fill = "white") +
  geom_jitter(aes(color = tr, shape = tr), width = 0.25, size = 2.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  labs(subtitle = sprintf("low; measure < %.1f | high; measure >= %.1f", value, value), x = NULL) +
  ggtitle("Ranking measure distribution")

# Clean up rank$config and subtrace$config columns
rank$config <- stringr::str_remove(stringr::str_remove(rank$config, paste0(net_fname, "-analysis-")), ".trace")
subtrace$config <- stringr::str_remove(stringr::str_remove(subtrace$config, paste0(net_fname, "-analysis-")), ".trace")
subtrace$measure <- numeric(length(subtrace$config))

# Merge subtrace and subflux with rank to add measure information
subtrace <- merge(subtrace, rank[, c(1, 2)], by = "config")
subflux$config <- stringr::str_remove(stringr::str_remove(subflux$config, paste0(net_fname, "-analysis-")), "-0.flux")
subflux$measure <- numeric(length(subflux$config))
subflux <- merge(subflux, rank[, c(1, 2)], by = "config")

# Get the ID with the minimum measure value and write it to ConfigID.txt file
id <- rank$config[which(rank$measure == min(rank$measure))]
writeLines(id, "ConfigID.txt")

# Extract trace and flux paths based on the ID
TracesPath <- paste0(net_fname, "-analysis-", id, ".trace")
FluxPath <- paste0(net_fname, "-analysis-", id, "-0.flux")

# Extract numerical values from ID
n <- as.numeric(gsub(pattern = paste0("(", net_fname, "-analysis-)|(.trace)"), replacement = "", x = id))
f <- as.numeric(gsub(pattern = paste0("(", net_fname, "-analysis-)|(-0.flux)"), replacement = "", x = id))

# Set working directory
setwd(paste0(wd, "/results/CDiff", "_", tag, "_", Condition, "_", Exper))

# Read trace and flux data
traceSelected <- read.table(TracesPath, header = T)
fluxSelected <- read.table(FluxPath, header = T)

# Reshape traceSelected and fluxSelected data
traceSelected <- tidyr::gather(traceSelected, key = "Places", value = "Marking", -Time)
fluxSelected <- tidyr::gather(fluxSelected, key = "Places", value = "Marking", -Time)

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

flux = filter(subflux, subflux$Places %in% 
                c("EX_biomass_e", "sink_pheme_c",
                  "EX_trp_L_e", "EX_pro_L_e", "EX_val_L_e",
                  "EX_ile_L_e", "EX_cys_L_e", "EX_leu_L_e"))

pCM = ggplot(Cdiff, aes(x = Time, y = Marking, color = measure.y)) +
  xlab("Time (hour)") + ylab("C. difficile (cell)") +
  geom_line(aes(group = config), alpha = alph) + 
  geom_line(data = filter(traceSelected, traceSelected$Places == "CD"), col = "red") + 
  ggtitle(paste("CD Place | config selected: ", id, sep = "")) + 
  scale_colour_gradientn(name = "", colors = c("#00AFBB", "#E7B800"), oob = scales::squish)

pdrugM = ggplot(Drug, aes(x = Time, y = Marking*1e-09, color = measure.y)) +
  geom_hline(yintercept = 5.8e-06, linetype = "dashed", color = "darkgrey", lwd = 0.5) +
  annotate("text", 4, 6.3e-06, label = "MIC value", vjust = 0.5, color = "black", size = 2) +
  xlab("Time (hour)") + ylab("Drug (mmol)") +
  geom_line(aes(group = config), alpha = alph) + 
  geom_line(data = filter(traceSelected, traceSelected$Places == "Drug"), col = "red") +
  ggtitle(paste("Drug Place | config selected: ", id, sep = "")) + 
  scale_colour_gradientn(name = "", colors = c("#00AFBB", "#E7B800"), oob = scales::squish)

pIECsM = ggplot(Iecs, aes(x = Time, y = Marking, color = measure.y)) +
  xlab("Time (hour)") + ylab("IECs (cell)") +
  geom_line(aes(group = config), alpha = alph) + 
  geom_line(data = filter(traceSelected, traceSelected$Places == "IECs"), col = "red") +
  ggtitle(paste("IECs | config selected: ", id, sep = "")) + 
  scale_colour_gradientn(name = "", colors = c("#00AFBB", "#E7B800"), oob = scales::squish)

pBiomassM = ggplot(Biomass, aes(x = Time, y = Marking, color = measure.y)) +
  xlab("Time (hour)") + ylab("Biomass (pg)") +
  geom_line(aes(group = config), alpha = alph) + 
  geom_line(data = filter(traceSelected, traceSelected$Places == "BiomassCD"), col = "red") + 
  ggtitle(paste("Biomass | config selected: ", id, sep = "")) + 
  scale_colour_gradientn(name = "", colors = c("#00AFBB", "#E7B800"), oob = scales::squish)

pfluxM = ggplot(flux, aes(x = Time, y = Marking, color = measure.y)) +
  xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
  geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") + 
  ggtitle("Sensitivity Analysis: Fluxes (mmol/gDW*h)") +
  scale_colour_gradientn(name = "", colors = c("#00AFBB", "#E7B800"), oob = scales::squish)

#####

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
pdrugM.median <- ggplot(Drug, aes(x = Time, y = Marking, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
  xlab("Time (hour)") + ylab("C. difficile (cell)") +
  geom_line(aes(group = config), alpha = alph) + 
  geom_line(data = filter(traceSelected, traceSelected$Places == "Drug"), col = "red") + 
  ggtitle("Drug Place") +
  scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800"))

# IECs Place
pIECsM.median <- ggplot(Iecs, aes(x = Time, y = Marking, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
  xlab("Time (hour)") + ylab("C. difficile (cell)") +
  geom_line(aes(group = config), alpha = alph) + 
  geom_line(data = filter(traceSelected, traceSelected$Places == "IECs"), col = "red") + 
  ggtitle("IECs Place") +
  scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800"))

# Biomass Place
pBiomassM.median <- ggplot(Biomass, aes(x = Time, y = Marking, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
  xlab("Time (hour)") + ylab("C. difficile (cell)") +
  geom_line(aes(group = config), alpha = alph) + 
  geom_line(data = filter(traceSelected, traceSelected$Places == "BiomassCD"), col = "red") + 
  ggtitle("Biomass Place") +
  scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800"))

# Sensitivity Analysis: Fluxes
pfluxM.median <- ggplot(flux, aes(x = Time, y = Marking, color = cut(measure.y, c(min(measure.y)-1, median, max(measure.y))))) +
  xlab("Time (h)") + ylab("Output fluxes (mmol/gDW*h)") +
  geom_line(aes(group = config), alpha = alph) +
  facet_wrap(~Places, scale = "free") + 
  ggtitle("Sensitivity Analysis: Fluxes (mmol/gDW*h)") + 
  scale_color_manual(labels = c("low", "high"), "Groups", values = c("#00AFBB", "#E7B800")) + 
  facet_wrap(~Places, scale = "free")
