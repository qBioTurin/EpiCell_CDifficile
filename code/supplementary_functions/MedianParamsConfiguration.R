Condition = "Therapy"
tag = c("Ablated", "Unified")

configuration.closeMedian = function(tag,Condition,subtrace=NULL){
  Exper = "Model_Sensitivity"

  Na = 6.022e20
  c = 6.022e08
  pack = 1*(Na*(1/c))
  fc = 1e-06
  
  aa_places <- c("trp_L_e", "pro_L_e", "val_L_e", "ile_L_e", "cys_L_e", "leu_L_e")
  places <- c("CD", "IECs", "BiomassCD", "Drug", "pheme_e", "pheme_c", aa_places)
  units = c(rep("(cell)", 2), "(pg)", rep("(µmol)", 9))
  
  if(is.null(subtrace)){
  subtrace <- do.call(rbind, lapply(tag, function(j) {
    subtrace = readRDS(file = paste0(wd, paste0("/results/CDiff", "_", j, "_", Condition, "_", Exper),"/subtrace_" , j, Condition, ".RDs"))
    cbind(subtrace, Scenario = rep(j, length(subtrace$Time)))
  }))
  }
  # 3 discrete intervals
  
  dfgroups = subtrace %>%
    filter(Places == "IECs", Time == 20) %>% 
    group_by(Scenario) %>%
    mutate(group = ntile(Marking, 3))%>%
    group_by(Scenario,group) %>%
    mutate(median_value = median(Marking)) %>%
    slice_min(abs(Marking - median_value),with_ties = F) 
  
  df= do.call(rbind,
              lapply(1:dim(dfgroups)[1],function(j){
                jj = dfgroups$config[j]
                df = readr::read_delim(jj)
                df$Config=dfgroups$config[j]
                df$ConfigNumb = gsub(pattern = ".*-([0-9]+)\\.trace", "\\1",x = basename(dfgroups$config[j]))
                df$Scenario = dfgroups$Scenario[j]
                df$ConfParams = dfgroups$group[j]
                df
              })
  ) %>% tidyr::gather(-ConfParams,-ConfigNumb,-Config,-Time,-Scenario,key = "Places",value="Marking")
  
  
  pl = ggplot(df%>%filter(Places %in% c("CD","IECs","BiomassCD","pheme_e")))+
    geom_line(aes(x=Time,y=Marking, group = ConfParams, linetype=Scenario, col = as.factor(ConfParams)) )+
    theme_bw()+theme(legend.position = "bottom")+
    facet_wrap(~Places,scales = "free")
  
  
  Config = df%>%filter(Places %in% c("CD","IECs","BiomassCD")) %>% select(Scenario,Config,ConfParams) %>% distinct() %>%
    mutate(ConfParams = paste0("Set ", ConfParams))
  paramsConfig = subtrace %>% filter(config %in% unique(Config$Config) ) %>% select(Detox,Death4Treat,IECsDeath,Scenario,config) %>% distinct() 
  paramsConfig= merge(paramsConfig,Config, by.x = c("Scenario","config"), by.y = c("Scenario","Config")) %>%select(-config)
  groupsConfig = list(plot = pl, Config = paramsConfig)
  
  return(groupsConfig)
}








