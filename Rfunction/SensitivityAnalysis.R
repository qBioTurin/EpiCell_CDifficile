
SensitivityAnalysis = function(seed, combinations, 
                               algo, model, verbose, 
                               sink, save,
                               distr, mean, sd, 
                               minBound, maxBound) {
  
  # seed = integer to set seed
  # combinations = integer to set the number of bounds sets
  # algo = "fba"; one objective, biomass maximization
  # algo = "mtf", double objective; biomass maximization and minimun total flux minimization
  # model = path to genome scale metabolic reconstruction in .RData format
  # verbose = logical, T then print
  # sink = logical, T if model with sink_heme
  # distr = "rnorm" character vector to set distribution
  # distr = "runif" character vector to set distribution
  # mean = double for rnorm()
  # sd = double for rnorm()
  # maxBound = double for maximun bounds possible
  # minBound =  double for minimum bounds possible
  
  set.seed(seed)
  combinations = combinations
  
  scenario = as.character(1:combinations)
  
  if (algo == "fba") {
    source("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Rfunction/flux_balance_fba.R")
  } else {
    source("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Rfunction/flux_balance.R")
    }
  
  load(paste("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/",
         model, "/", model, ".RData", sep = ""))
  
  ex <- sybil::findExchReact(model.mat)
  ex = ex@react_id
  
  EX_biomass = c()
  EX_ile = c(); EX_pro = c(); EX_leu = c(); EX_val = c(); EX_cys = c(); EX_trp = c();  
  Sink_heme = c()
  hyd4 = c()
  
  FBA.allScenarios = list()
  BoundsSet = list()
  
  for (i in 1:combinations) {
    
    lowbndTemplate = model.mat@lowbnd[which(model.mat@react_id %in% ex)]
    uppbndTemplate = model.mat@uppbnd[which(model.mat@react_id %in% ex)]
    
    if(distr == "rnorm") {
      
      # test
      # plot(density(rnorm(length(lowbndTemplate), mean = -100, sd = pi*25)))
      # plot(density(rnorm(length(uppbndTemplate), mean = 100, sd = pi*25)))
      
      x = rnorm(length(lowbndTemplate), mean = -mean, sd = sd)
      x[which(x >= 0.0)] = rep(-minBound, length(x[which(x >= 0.0)]))
      x[which(x < -maxBound)] = rep(-maxBound, length(x[which(x < -maxBound)]))
      lowbnd = x; # plot(density(lowbnd))
      lowbnd[which(lowbndTemplate == 0.0)] = rep(0.0, length(which(lowbndTemplate == 0.0)))
      
      y = rnorm(length(lowbndTemplate), mean = mean, sd = sd)
      y[which(y > maxBound)] = rep(maxBound, length(y[which(y > maxBound)]))
      y[which(y <= 0.0)] = rep(minBound, length(y[which(y <= 0.0)]))
      uppbnd = y; # plot(density(uppbnd))
      
      } else {
        
        lowbnd = runif(length(lowbndTemplate), min = -maxBound, max = -minBound)
        uppbnd = runif(length(lowbndTemplate), min = minBound, max = maxBound)
        
        }
    
    df = data.frame(lowbnd, uppbnd)
    BoundsSet[[i]] = df
    
  }
  
  if (save == T) {
      save(BoundsSet, 
       file = paste("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/",
                    model, "/BoundsSet", algo, distr, ".RData", sep = ""))
  }
    
  for (i in 1:length(BoundsSet)) {
    
    load(paste("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/",
               model, "/", model, ".RData", sep = ""))
    
    if (verbose == T) {
      print(i)
    }
    
    model.mat@lowbnd[which(model.mat@react_id %in% ex)] = BoundsSet[[i]]$lowbnd
    model.mat@uppbnd[which(model.mat@react_id %in% ex)] = BoundsSet[[i]]$uppbnd
    
    sol = flux_balance(model = model.mat, typename = model.mat@mod_desc, diet = F)
    
    FBA.allScenarios[[i]] = list(sol = sol, lb = sybil::lowbnd(model.mat), ub = sybil::uppbnd(model.mat))
    
    if (algo == "mtf") {
      v = sybil::getFluxDist(sol[[1]])
      EX_biomass[i] = FBA.allScenarios[[i]][[1]][[1]]@alg_par[['wtobj']]
    } else {
      v = sol[[1]][["fluxes"]]
      EX_biomass[i] = v[which(model.mat@react_id == 'EX_biomass(e)')]
    }
    
    EX_pro[i] = v[which(model.mat@react_id == 'EX_pro_L(e)')]
    EX_leu[i] = v[which(model.mat@react_id == 'EX_leu_L(e)')]
    EX_ile[i] = v[which(model.mat@react_id == 'EX_ile_L(e)')]
    EX_val[i] = v[which(model.mat@react_id == 'EX_val_L(e)')]
    EX_cys[i] = v[which(model.mat@react_id == 'EX_cys_L(e)')]
    EX_trp[i] = v[which(model.mat@react_id == 'EX_trp_L(e)')]
    
    Sink_heme[i] = v[which(model.mat@react_id == 'sink_pheme(c)')]
    
    hyd4[i] = v[which(model.mat@react_id == 'HYD4')]

    load(paste("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/",
               model, "/", model, ".RData", sep = ""))
    
  }
  
  if (sink == T) {
    SensitivityData = data.frame(EX_biomass, EX_pro, EX_leu, EX_val, EX_ile, EX_cys, EX_trp,
                                 Sink_heme, hyd4, BoundsSet = scenario)
  } else {
    SensitivityData = data.frame(EX_biomass, EX_pro, EX_leu, EX_val, EX_ile, EX_cys, EX_trp,
                                 hyd4, BoundsSet = scenario)
    }

  SensitivityData = tidyr::gather(SensitivityData, key, Value, -BoundsSet)
  
  if (save == T) {
      save(SensitivityData, 
       file = paste("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/",
                    model, "/SensitivityData", algo, distr, ".RData", sep = ""))
  }
  
  return(SensitivityData)
}