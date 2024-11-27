
setwd("~/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink")
load(file = "./CD196HemeSinkDiet.RData")

BIGGdata_CD196HemeSink_react <- 
  read.delim2("~/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/BIGGdata_CD196HemeSink_react.tsv")

rxnNames <- readxl::read_excel("./rxnNames.xlsx", col_names = FALSE)
CD196HemeSink_BiGGdata = cbind(rxnNames, BIGGdata_CD196HemeSink_react)
colnames(CD196HemeSink_BiGGdata)[1] <- "rxnNames"

ReactionsNames = unlist(model.mat@react_id)
# reactant's nomenclature in accordance with VMH
ReagentsNames = unlist(model.mat@met_id)

S = as.matrix(model.mat@S)
ReagentsIndex = which(ReagentsNames %in% 
                        c("isocapr[c]", "isocapr[e]", 
                          "isoval[c]", "isoval[e]",
                          "leu_L[c]", "leu_L[e]"))

ReactionIndex = c()

for(i in 1:length(ReagentsIndex)) {
  print(which(S[ReagentsIndex[i], ] != 0))
  ReactionIndex = append(ReactionIndex, which(S[ReagentsIndex[i], ] != 0))
}

reactList = ReactionsNames[ReactionIndex]

reactList = gsub("\\(", replacement = "_", reactList)
reactList = gsub("\\)", replacement = "", reactList)

setwd("~/EpiCell_CDifficile")

Exper = "Model_Sensitivity"
Condition = "Therapy"
tag = c("Ablated", "Unified")
param_target = c("IECsDeath", "Death4Treat", "Detox")

Tempi = c(0, 12, 24, 36, 48, 60)

coloTag = data.frame(SuperAblated = "darkred", Ablated = "#EC6D67", Unified = "#5351A2")
nameTag = data.frame(SuperAblated = "Ablated", Ablated = "Partially Ablated", Unified = "Unified")

colo1 <- c("black", "magenta", "gold")
colo2 <- c("black", "#266867", "yellow")
colo3 <- c("#ffd166", "#ee6c4d", "#293241")

for (j in tag) {
  
  setwd(paste("~/EpiCell_CDifficile/Results/CDiff", j, Condition, Exper, sep = "_"))
  
  load(file = paste("subflux_" , j, Condition, ".RData", sep = ""))
  subflux = cbind(subflux, Scenario = rep(j, length(subflux$Time)))
  assign(paste("subflux", j, Condition, sep = ""), subflux)
  
}

subflux = rbind(subfluxAblatedTherapy, subfluxUnifiedTherapy)

Th = 0.05

for (j in tag) {
  
  assign(paste("pflux", param_target[1], j, sep = "."), 
         filter(subflux[which(subflux$Places %in% 
                                reactList & subflux$Marking < Th & subflux$Marking > -Th), ], 
                Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking, color = IECsDeath)) +
           labs(x = "Time (h)", 
                y = "Output fluxes (mmol/gDW*h)",
                title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[1], sep = ""),
                subtitle = paste("Scenario: " , j, sep = "")) + 
           theme(plot.subtitle = element_text(size = 12, face = "bold", color = coloTag[, j]),
                 plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) + 
           geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
           scale_colour_gradientn(name = "",
                                  colors = colo1,
                                  oob = scales::squish))
  
  assign(paste("pflux", param_target[2], j, sep = "."), 
         filter(subflux[which(subflux$Places %in% 
                                reactList & subflux$Marking < Th & subflux$Marking > -Th), ], 
                Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking, color = IECsDeath)) +
           labs(x = "Time (h)", 
                y = "Output fluxes (mmol/gDW*h)",
                title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[2], sep = ""),
                subtitle = paste("Scenario: " , j, sep = "")) + 
           theme(plot.subtitle = element_text(size = 12, face = "bold", color = coloTag[, j]),
                 plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) + 
           geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
           scale_colour_gradientn(name = "",
                                  colors = colo2,
                                  oob = scales::squish))
  
  assign(paste("pflux", param_target[3], j, sep = "."), 
         filter(subflux[which(subflux$Places %in% 
                                reactList & subflux$Marking < Th & subflux$Marking > -Th), ], 
                Scenario == j) %>%
           ggplot(aes(x = Time, y = Marking, color = IECsDeath)) +
           labs(x = "Time (h)", 
                y = "Output fluxes (mmol/gDW*h)",
                title =  paste("Sensitivity Analysis: Fluxes (mmol/gDW*h) | ", param_target[3], sep = ""),
                subtitle = paste("Scenario: " , j, sep = "")) + 
           theme(plot.subtitle = element_text(size = 12, face = "bold", color = coloTag[, j]),
                 plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#2a475e"),
                 plot.title.position = "plot",
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, face = "bold")) + 
           geom_line(aes(group = config)) + facet_wrap(~Places, scale = "free") +
           scale_colour_gradientn(name = "",
                                  colors = colo3,
                                  oob = scales::squish))
  
}

p = (pflux.Death4Treat.Unified / pflux.Death4Treat.Ablated) |
  (pflux.IECsDeath.Unified / pflux.IECsDeath.Ablated) |
  (pflux.Detox.Unified / pflux.Detox.Ablated)

setwd("~/EpiCell_CDifficile")
ggsave(p, file = paste("SelectedReact.pdf", sep = ""), width = 30, height = 14)
system(paste("cp SelectedReact.pdf ", "./Results/images4paper/", sep = ""))
system(paste("mv SelectedReact.pdf ", "./Results/images4paper/", sep = ""))
