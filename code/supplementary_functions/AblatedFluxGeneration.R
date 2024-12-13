pla=c("EX_biomass_e" = "BiomassCD", "EX_cys_L_e" = "cys_L_e","EX_trp_L_e" = "trp_L_e", "EX_val_L_e" = "val_L_e",
"EX_ile_L_e" = "ile_L_e",
"EX_leu_L_e" = "leu_L_e",
"EX_pro_L_e" = "pro_L_e",
"sink_pheme_c_in" = "pheme_c")


f = flux %>% filter(new_eps_value=="1e-6",Time == 0,tag =="Unified") %>% select(-config,-tag) %>% distinct()
f[f$Reaction == "EX_phe_L_e","Flux"] = -0.001172480768218138
f[f$Flux == -0.001172480768218138,"Reaction"] = "sink_pheme_c_in" 

merge(f%>% mutate(Places = pla[Reaction]),
trace %>% filter(new_eps_value=="1e-6",Time == 0,Places %in% pla) %>% select(-ConfParams,-config,-tag) %>% distinct()) %>% mutate(Division = Flux/Marking)


Epit<- readr::read_table("results/TimingEval/CDiffAblatedTherapy1e-6Set2/EpitCellDifficileHemeSink-analysis-1-0.flux")
Epit %>% select(EX_biomass_e,sink_pheme_c,EX_cys_L_e,EX_leu_L_e) %>% distinct()

EpitU<- readr::read_table("results/TimingEval/CDiffUnifiedTherapy1e-6Set2/EpitCellDifficileHemeSink-analysis-1-0.flux")
EpitU  %>% select(Time,EX_biomass_e,sink_pheme_c,EX_cys_L_e,EX_leu_L_e) %>% distinct() 

EpitP<- readr::read_table("results/TimingEval/CDiffParAblatedTherapy1e-6Set2/EpitCellDifficileHemeSink-analysis-1-0.flux")
EpitP  %>% select(Time,EX_biomass_e,sink_pheme_c,EX_cys_L_e,EX_leu_L_e) %>% distinct() 
