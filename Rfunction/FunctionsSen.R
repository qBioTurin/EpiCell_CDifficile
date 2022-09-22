
init.gen <- function(n_file, file, d, zero, zeroH) {
  
  yini.names <- readRDS(n_file)
  
  # the determination of heme concentrations in different conditions range 
  # from 20 to >350 μmol/L as occurring, e.g., in hemolytic events.
  
  if (zeroH == T ){
    initHeme = 1e-12 # [pmol/mL]
  } else {
    initHeme = (0.01*1e-03*1e+09)/100 # [pmol/mL]
  }
  
  
  Vcd = ((pi*(0.5)^2)/4)*(((4*0.5)/6) + 5.5) # micrometer^3
  gDW_CDmean = (1/3)*Vcd*1.3 # picograms
  
  MOI = 100
  
  Na = 6.022e+20; c = 602200000
  
  load(file) # fluxes in mmol/h*person
  
  diet = diets[[d]]
  
  vec_iniAA = (-1*diets$EU_average[match(
    c("EX_pro_L(e)", "EX_leu_L(e)", 
      "EX_ile_L(e)", "EX_val_L(e)", 
      "EX_trp_L(e)", "EX_cys_L(e)"),
    diets$EU_average$Reaction), ]$`Flux Value`)*1e-04
  # [mmol/h*person]
  
  vec_iniAA = (vec_iniAA*Na)/c
  
  # y_ini = c(IECs, Damage, 
  #           BiomassCD, CD, 
  #           Drug, 
  #           pheme_e, pheme_c, 
  #           pro_L_e, leu_L_e, ile_L_e, val_L_e, trp_L_e, cys_L_e,
  #           pro_L_v, leu_L_v, ile_L_v, val_L_v, trp_L_v, cys_L_v)
  
  if (d == "no diet") {
    vec_iniAA = c(0.36, 0.57, 0.32, 0.58, 0.098, 0.25)*1e-06 # [mol/mL]
    vec_iniAA = (vec_iniAA*Na)/c
    # vec_MMAA = c(115.13, 131.17, 131.17, 117.151, 204.23, 121.16) # [g/mol]
  }
  
  if (zero == T) {
    vec_iniAA = vec_iniAA*1e-15
  }
  
  y_ini <- c(4.57e05, 0, 
             gDW_CDmean, 4.57e05*MOI, 
             0, initHeme, 0, vec_iniAA,
             0, 0, 0, 0, 0, 0)
  
  names(y_ini) <- c("IECs", "Damage", 
                    "BiomassCD", "CD", 
                    "Drug", 
                    "pheme_e", "pheme_c",
                    "pro_L_e", "leu_L_e", "ile_L_e", "val_L_e", "trp_L_e", "cys_L_e",
                    "pro_L_v", "leu_L_v", "ile_L_v", "val_L_v", "trp_L_v", "cys_L_v")
  
  y_ini = y_ini[yini.names]
  
  return(y_ini)
  
}

ComputeBiomassBac = function(diameter, len) {
  
  Vcd = ((pi*(diameter)^2)/4)*(((4*diameter)/6) + len) # μm^3
  mass = (1/3)*Vcd*1.3 # [pg]
  
  return(mass)
}

InflammationFunction = function(gamma) {
  
  W = 0.1 # villus width [mm]
  L = 1 # villus height [mm]
  EM = 13 # surface amplification due to microvilli [unit]
  
  VIarea = 2*pi*(W/2)*L*EM # [mm^2]
  
  Na = 6.022e20 # [molecule]
  
  Vrbc = 11000 # red blood cell velocity [mm/h]
  RHOrbc = 6e+06 # [cell/microliter]
  HeamRBC = 5e+09 # [molecule]
  Dv = 0.03 # [mm]
  Aw = 113 # [mm^2]
  
  VA = pi*((Dv/2))^2 # vassel section [mm^2]
	
  Inflammation = gamma*VA*(Vrbc*RHOrbc*((HeamRBC/Na)*1e+09))*(Aw/VIarea) # [pmol]
  
  return(Inflammation)
  
  }

EvalDiet = function(file, diet, ex, c, HemeDiet) {
  
  if (is.null(diet)) {
    
    Mpheme = 616.5*1e-12 # [g/pmol]
    
    # What percentage of haem iron is absorbed? (15% to 35%)
    
    # The median intake of Total, Heme, and Nonheme iron was:
    # 23.0, 0.75, and 22.0 [mg/day*person] among men, and:
    # 20.0, 0.63, and 19.1 [mg/day*person] among women
    
    if (HemeDiet == "poor") {
      vheme = 0.75*1e-03*1e-04 # [g/(day*person)]
      } else {
        vheme = 0.75*100*1e-03*1e-04 # [g/(day*person)]
      }
    
    vD = (1/24)*((vheme/Mpheme)) # [pmol/h]
    
    } else {
      
      Na = 6.022e20 # [molecule]
      c = 6.022e8 # [unit]
      
      # diet = readr::read_delim(paste0("./Input/Diets/vmh/", diet, ".tsv", sep=""), 
      #                         "\t", escape_double = FALSE, trim_ws = TRUE)
      
      load(file)
      
      diet = diets[[diet]]
      flux = -1*diet[["Flux Value"]][which(diet[["Reaction"]] %in% ex)]*1e-04 # [mmol/h*person]
      # reason behind "1e-04": flux of nutrient are given per person, then 1e-04
      # could be a good scaler indicative for mL orger of magnitude
      vD = (flux*(Na/c))
      
      }
  
  return(vD)
}

EvalTransport = function(yield) {

	IECs_t0 = 4.57e05
	vT = yield/IECs_t0
	
	return(vT)
}

# EvalTreat = function(dose){
#   
#   MM = 171.16*1e-12 # [g/pmol]
#   
#   # dose: Adults, 0.5 or 0.75 or 1 [g] 3 times a day for 5 to 10 days
#   
#   dose = dose*3 # [g/day*person]
#   
# 	ther = (1/24)*((dose/MM))*1e-03 #[pmol/mL*h^-1]
# 	
# 	return(ther)
# 
# }

# EvalEmax = function(p0, p24){
# 
#   # 1e+06*cfu/ml;
# 
#   # DOI: 10.1093/jac/39.suppl_2.23
# 
#   MM = 171.16*1e-12 # [g/pmol]
# 
#   MIC = 8e-06 # [g/mL]
#   drug0 = MIC/MM # [pmol/mL]
# 
#   Emax = (abs(10^p24 - 10^p0))/(24*drug0*10^p0)
# 
#   return(Emax)
# 
# }

EvalDetox = function(HsmA){
  
  # doi: 10.1289/ehp.856453.
  
  # k for reaction: PFe3 + RO2· -> [PFe3] + RO2-
  # it is: 
  # k = 2.5e+08*60 # [1/mmol*min]
  kHsmA = 2.5e+17*60*60 # [1/pmol*h]
  
  # HsmA [pmol] # HsmA = 1e-20
  
  # [pmol/h]
  detox = HsmA*kHsmA
  
  return(detox)
  
}

EvalEfflux = function(MDET){
  
  # MDET; Multi Drug Efflux Transporter
  
  # MDET = runif(1, min = 0.0, max = 1e-12)
  
  # https://doi.org/10.1101/2020.01.16.909507
  # k for reaction: drug_c -> drug_e
  # k is 0.0235 [1/mol*s] -> 1420 [1/mmol*min]
  
  kMDET = 1420*60*1e+09 # [1/pmol*h]
  
  # MDET [pmol]
  # is about 2.35e-15/1e+06 [mol] -> 2.35e-09 [pmol]
  
  # [pmol/h]
  efflux = MDET*kMDET
  
  return(efflux)
  
}

EvalVmax = function(V, gamma, k) {
  
  # link: https://doi.org/10.1016/j.cell.2005.06.025
  
  # V = 3.1 [pmol/min/μg protein] -> 
  V = V*60*1e+06 # [pmol/h/g protein] 
  
  # gamma = proportion of mass protein which is HEMEti = 1e-03 (ref?)
  # k = proportion of reversibility [unit]
  
  gDW_CDmean = (1/3)*(((pi*(0.5)^2)/4)*(((4*0.5)/6) + 5.5))*1.3e-12 # [g]
  
  vmax = k*(V*gDW_CDmean*gamma) # [pmol/gTransporter*h]
}
  
MapParameter = function(optim_v, indx) {
  return(optim_v[indx])
}

treat_generation = function(marking, time) {
  
  MM = 171.16*1e-12 # [g/pmol]
  
  #   MIC = 1e-06 # [g/mL]
  #   drug0 = (MIC/MM)*1e-09 # [mmol/mL] about 5.8e-06 [mmol/mL]
  
  # dose = 0.5 # [g/day] (prescribed dose)
  dose = 0.05 # dose for MIC
  
  dose = 8*(1/24)*((dose/MM))*1e-04 # [pmol/mL*8*h^-1]
  
  new_marking = marking[1:19]
  new_marking[10] = marking[10] + dose
  
  return(new_marking)
}
 
targetFN = function(output) {
  return(output[, "CD"])
}

##################################################
##### How to compute aa's arcs multiplicity ######
##################################################

# M = c(0.11513, 0.13117, 0.13117, 0.117151, 0.20423, 0.12116) # [g/mmol]
# S = c(0.41248, 0.54554, 0.28608, 0.352610, 0.13306, 0.046571) # [unit]
# c = 6.022e8 # [molecule]
# Na = 6.022e20 # [molecule]
# DW_IECs = 1e-9 # [g/cell]

# car = round(((DW_IECs*S)/M)*(Na/c), 0) # ([g/cell]*[unit]/[g/mmol])*([molecule]/[molecule])
