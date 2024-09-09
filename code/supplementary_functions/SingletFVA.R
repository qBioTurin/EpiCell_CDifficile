
# Differential reaction fluxes between NoDrug and Therapy models. 

#' To compare the fluxes between the NoDrug and Therapy models, 
#' the comparison of individual fluxes could yield inaccurate results given the
#' large alternative optimal solution space. 
#'  
#' We compared the flux distributions per reaction rather than the value provided 
#' by one solution. To obtain the flux probability distributions, we performed a
#' FVA of both the NoDrug and Therapy model.
#'   

load("./Results/FVA_Unified_NoDrug_Model_Analysis/metadataFVAUnifiedNoDrug.RData")
dfFVANoDrug = dfFVA
load("./Results/FVA_Unified_Therapy_Model_Analysis/metadataFVAUnifiedTherapy.RData")
dfFVATherapy = dfFVA

t = 48

tollFBAsol = 1e-10
tolldiff = 1e-03

path_model = "./Input/CDmodels/CD196HemeSink/CD196HemeSink.RData"
path_diet = "./Input/Diets/diets.RData"

################################################################################

for (co in conditions) {
  
  load(path_model)
  
  Vcd = ((pi*(0.7)^2)/4)*(((4*0.7)/6) + 9) # μm^3
  gDW_CDmax = (1/3)*Vcd*1.3 # [pg]
  
  Na = 6.02214154e+20
  c = 6.022e8
  
  load(path_diet)
  
  d = diets[["EU_average"]]
  
  ex <- sybil::findExchReact(model.mat)
  ex = ex@react_id
  
  value = d$`Flux Value`[which(d$Reaction %in% ex)]
  name = d$Reaction[which(d$Reaction %in% ex)]
  
  index = c()
  
  for (i in 1:length(name)) {index[i] = which(model.mat@react_id == name[i])}
  
  react_id_ex = which(model.mat@react_id %in% unlist(ex))
  lbNOTdiet = -10
  
  TracesPath = paste("./Results/CDiffFBA", co, "_analysis", 
                     "/EpitCellDifficileHemeSink-analysis-1.trace", sep = "")
  
  trace = read.table(TracesPath, header = F)
  ReactionsNames = trace[1, ]
  ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
  ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
  
  trace = read.table(TracesPath, header = T)
  colnames(trace) = ReactionsNames
  
  subtrace = trace %>% 
    dplyr::select(Time, pheme_c, BiomassCD, CD, 
                  pro_L_e, leu_L_e, ile_L_e, 
                  val_L_e, trp_L_e, cys_L_e) %>% 
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  Biom = subtrace$Marking[which(subtrace$Places == "BiomassCD" & subtrace$Time == t)]
  
  if ((gDW_CDmax - Biom) > 1e-03) {
    uppbndBiomass = (gDW_CDmax - Biom)
  } else {
    uppbndBiomass = 1e-03
  }
  
  lowbndBiomass = lbNOTdiet
  
  model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]] =
    rep(lbNOTdiet, length(model.mat@lowbnd[react_id_ex[which(model.mat@lowbnd[react_id_ex] != 0)]]))
  
  nBac = subtrace$Marking[which(subtrace$Places == "CD" & subtrace$Time == t)]
  
  model.mat@lowbnd[index] = (value*1e-04)/(nBac*Biom*1e-12) # (mmol/per community)
  
  Ub = 10
  aaUb = 1e-14
  
  model.mat@uppbnd[react_id_ex[which(model.mat@uppbnd[react_id_ex] != 0)]] =
    rep(Ub, length(model.mat@uppbnd[react_id_ex[which(model.mat@uppbnd[react_id_ex] != 0)]]))
  
  model.mat@uppbnd[which(model.mat@react_id %in% 
                           c("EX_pro_L(e)", "EX_leu_L(e)",
                             "EX_ile_L(e)", "EX_val_L(e)", 
                             "EX_trp_L(e)", "EX_cys_L(e)"))] = rep(aaUb, 6)
  
  Met = subtrace$Marking[which(subtrace$Places %in% 
                                 c("pro_L_e", "leu_L_e", "ile_L_e", 
                                   "val_L_e", "trp_L_e", "cys_L_e") & 
                                 subtrace$Time == t)]
  
  aaLb = model.mat@lowbnd[which(model.mat@react_id %in% 
                                  c("EX_pro_L(e)", "EX_leu_L(e)",
                                    "EX_ile_L(e)", "EX_val_L(e)", 
                                    "EX_trp_L(e)", "EX_cys_L(e)"))]
  
  aaLb = -((Met*c)/Na)/(nBac*Biom*1e-12)
  
  model.mat@lowbnd[which(model.mat@react_id %in% 
                           c("EX_pro_L(e)", "EX_leu_L(e)",
                             "EX_ile_L(e)", "EX_val_L(e)", 
                             "EX_trp_L(e)", "EX_cys_L(e)"))] = aaLb
  
  heme = subtrace$Marking[which(subtrace$Places == "pheme_c" & 
                                  subtrace$Time == t)]
  
  Lbheme = -(heme*1e-09)/(nBac*Biom*1e-12)
  
  Ubheme = 1e-04;
  
  model.mat@uppbnd[which(model.mat@react_id == "sink_pheme(c)")] = Ubheme
  model.mat@lowbnd[which(model.mat@react_id == "sink_pheme(c)")] = Lbheme
  
  model.mat@uppbnd[which(model.mat@react_id == "EX_biomass(e)")] = uppbndBiomass
  model.mat@lowbnd[which(model.mat@react_id == "EX_biomass(e)")] = lowbndBiomass
  
  assign(paste("model", co, sep = ""), model.mat)
  Irr = sybil::mod2irrev(model.mat)
  assign(paste("model", co, "Irr", sep = ""), Irr)
  fv = sybil::fluxVar(Irr)
  
  assign(paste("lpMaxSol", co, sep = ""), sybil::maxSol(fv, "lp_obj"))
  assign(paste("lpMinSol", co, sep = ""), sybil::minSol(fv, "lp_obj"))
  
}

df = data.frame(lpMaxSolTherapy, lpMaxSolNoDrug, 
                Maxdiff = abs((lpMaxSolTherapy - lpMaxSolNoDrug)/(lpMaxSolNoDrug)),
                lpMinSolTherapy, lpMinSolNoDrug, 
                Mindiff = abs((lpMinSolTherapy - lpMinSolNoDrug)/(lpMinSolNoDrug)),
                React_ID = fv@react@react_id, 
                ReactionType = c(rep("A", (Irr@react_num/2)), 
                                 rep("B", (Irr@react_num/2))),
                PNmodel = c(rep("Y", (Irr@react_num/2)), 
                            rep("N", (Irr@react_num/2))))

a = dplyr::filter(df, grepl("EX_", React_ID))
a$ReactionType = rep("Exchange", length(a$ReactionType))

b = dplyr::filter(df, grepl("DM_", React_ID))
b$ReactionType = rep("Demand/Sink", length(b$ReactionType))

c = dplyr::filter(df, grepl("sink_", React_ID))
c$ReactionType = rep("Demand/Sink", length(c$ReactionType))

d = dplyr::filter(df, grepl("tra", React_ID))
d$ReactionType = rep("Transcription", length(d$ReactionType))

e = dplyr::filter(df, grepl("rep", React_ID))
e$ReactionType = rep("Replication", length(e$ReactionType))

f = dplyr::filter(df, !grepl("biomass205|rep|tra|EX_|DM_|sink_", React_ID))
f$ReactionType = rep("Internals/Transporters", length(f$ReactionType))

h = dplyr::filter(df, grepl("biomass205", React_ID))
h$ReactionType = rep("Objective", length(h$ReactionType))

df = rbind(a, b, c, d, e, f, h)

# Model in irreversible format (forward/backward reactions separated):
#   
# forward reactions have to be marked with an _f in the id, backward with an _b. 
# Reverse reactions (i.e. reactions which were originally only allowed to carry negative 
#                    flux must have a _r tag at the end of the reaction id. 
# It is further assumed, that there is an exact 1:1 match between _f and _b reactions.

ReactVec = c("EX_biomass(e)_f", "sink_pheme(c)_b", 
             "EX_pro_L(e)_b", "EX_leu_L(e)_b", "EX_ile_L(e)_b", 
             "EX_val_L(e)_b", "EX_trp_L(e)_b", "EX_cys_L(e)_b")

pn1 = dplyr::filter(df, React_ID %in% ReactVec)
pn1$PNmodel = rep("Y", length(pn1$PNmodel))
pn2 = dplyr::filter(df, !(React_ID %in% ReactVec))
pn2$PNmodel = rep("N", length(pn2$PNmodel))

df[df$React_ID %in% ReactVec, ] = pn1
df[!(df$React_ID %in% ReactVec), ] = pn2

df$lpMaxSolTherapy[which(df$lpMaxSolTherapy < 0)] = 
  numeric(length(df$lpMaxSolTherapy[which(df$lpMaxSolTherapy < 0)]))

df$lpMaxSolNoDrug[which(df$lpMaxSolNoDrug < 0)] = 
  numeric(length(df$lpMaxSolNoDrug[which(df$lpMaxSolNoDrug < 0)]))

df$lpMaxSolTherapy[which(df$lpMaxSolTherapy < tollFBAsol)] = 
  numeric(length(df$lpMaxSolTherapy[which(df$lpMaxSolTherapy < tollFBAsol)]))

df$lpMaxSolNoDrug[which(df$lpMaxSolNoDrug < tollFBAsol)] = 
  numeric(length(df$lpMaxSolNoDrug[which(df$lpMaxSolNoDrug < tollFBAsol)]))

df = filter(df, df$Maxdiff > tolldiff)

df = filter(df, !(df$lpMaxSolTherapy == 0 & df$lpMaxSolNoDrug == 0))

df <- df[order(df$Maxdiff), ]

# ex <- sybil::findExchReact(Irr)
# ex = ex@react_id
# 
# df = filter(df, df$React_ID %in% ex)

setwd("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Results")

scientific_10 <- function(x) {parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))}

m = min(df$Maxdiff[df$React_ID %in% ReactVec])
data = df %>% filter((Maxdiff > 31.5) | (React_ID %in% ReactVec))

pdf(paste("FVAplot_t", t, ".pdf", sep = ""), width = 7, height = 5.5, bg = "transparent")

ggplot(aes(x = lpMaxSolTherapy, y = lpMaxSolNoDrug), data = df) + theme_bw() +
  geom_point(aes(color = ReactionType, shape = ReactionType), size = 2.3) +
  ggrepel::geom_label_repel(data = data, 
                            size = 2.75, family = 'Times', segment.size = 0.5,
                            box.padding = 0.4, nudge_y = 0.4, alpha = 0.75, 
                            aes(label = React_ID),
                            max.overlaps = 20, segment.color="darkgrey",
                            label.padding = unit(0.1, "lines"),
                            label.size = NA, seed = 1234, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, size = 0.2, linetype="solid", color="black") +
  ggtitle("Flux Variability Analysis (mmol/gDW*h)") +
  xlab("log scale FV CDmodel Therapy (mmol/gDW*h)") + 
  ylab("log scale FV CDmodel No Drug (mmol/gDW*h)") +
  labs(subtitle = paste("time = ", as.character(t), sep = "")) +
  theme(plot.title = element_text(color="black", size = 12, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.83, 0.25),
        legend.background = element_rect(size = 0.2, fill = "white", color = "black")) + 
  scale_x_log10(label = scientific_10) + scale_y_log10(label = scientific_10)

dev.off()

summary(df$Maxdiff)
# pdf("test.pdf", width = 7, height = 5.5, bg = "transparent")
# m<-mean(df$Maxdiff)
# std<-sqrt(var(df$Maxdiff))
# hist(df$Maxdiff,prob=T, main="Max difference")
# curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE)
# dev.off()
