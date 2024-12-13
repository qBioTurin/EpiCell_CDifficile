
init.gen <- function(n_file, file, d, zero, zeroH) {
  yini.names = readRDS(n_file)
  
  # The determination of heme concentrations in different conditions range from 20 to > 350 μmol/L as occurring, 
  # e.g., in hemolytic events.
  
  if (zeroH == T) {
    initHeme = 1e-12 # (pmol/mL)
  } else {
    initHeme = 100 # (μM)
    
    # (https://doi.org/10.1038/s41419-020-02950-8)
    # quick conversion: 1 (μM) = 1.0e-06 (mmol/mL) 
    # initHeme * 1.0e-06  = (mmol/mL)
    # conversion: 1 (mmol/mL) = 1e+09 (pmol/mL) 
    # initHeme * 1.0e-06 * 1e+09 (pmol/mL)
    
    initHeme = initHeme * 1.0e-06 * 1e+09 # (pmol/mL)
    # Note: hemorragic heme = 1000 (μM) = 1000 * 1.0e-06 * 1e+09 (pmol/mL)
  }
  
  Vcd = ((pi*(0.5)^2)/4)*(((4*0.5)/6) + 5.5) # (μm^3)
  gDW_CDmean = round((1/3)*Vcd*1.3, 4) # (pg)
  
  MOI = 100
  
  Na = 6.02214154e+20 # (molecule/mmol)
  c = 6.022e+08
  
  load(file) # fluxes (mmol/h*person)
  
  diet = diets[[d]]
  
  Delta = 1e-04
  
  vec_iniAA = (diets$EU_average[match(
    c("EX_cys_L_e", "EX_leu_L_e", 
      "EX_ile_L_e", "EX_val_L_e", 
      "EX_trp_L_e", "EX_pro_L_e"),
    diets$EU_average$Reaction), ]$`Flux Value`)*Delta # (mmol)
  
  # In the model aminoacids amounts are measured as (C_molecules)
  #
  # Literature information is units of (mmol)
  # quick conversion: 1 (mmol) = 6.02214154+20 (molecule)
  # 6.02214154+20 (molecule) = Na (Avogadro number reduced)
  
  vec_iniAA = vec_iniAA*(Na*(1/c))
  
  # y_ini = c(IECs, Damage, 
  #           BiomassCD, CD, 
  #           Drug, 
  #           pheme_e, pheme_c, 
  #           pro_L_e, leu_L_e, ile_L_e, val_L_e, trp_L_e, cys_L_e,
  #           pro_L_v, leu_L_v, ile_L_v, val_L_v, trp_L_v, cys_L_v)
  
  if (d == "no diet") {
    vec_iniAA = c(0.25, 0.57, 0.32, 0.58, 0.098, 0.36) # (μmol/mL)
    
    # in the model AAs are measured as (C_molecules)
    # literature information as (μmol)
    # quick conversion: 1 μmol = 6.0221415e+17 (molecule)
    
    vec_iniAA = vec_iniAA*(6.0221415e+17*(1/c)) # (C_molecule)
    # x = mean(250005874, 570013393, 320007519, 580013628, 98002303, 360008459) (C_molecule)
    # y = mean(3.130003e+11, 1.953705e+12,1.205098e+12, 1.574386e+12, 1.897888e+11, 2.242356e+12) (C_molecule)
    # Delta = (x/y)
    # vec_MMAA = c(115.13, 131.17, 131.17, 117.151, 204.23, 121.16) # [g/mol]
  }
  
  if (zero == T) {
    vec_iniAA = vec_iniAA*1e-15 # (C_molecule)
  }
  
  y_ini <- c(4.57e05, 0, 
             gDW_CDmean, 4.57e05*MOI, 
             0, 
             initHeme, 1e+03, 
             # gut comp aminoacids
             vec_iniAA,
             # blood comp aminoacids
             0, 0, 0, 0, 0, 0)
  
  names(y_ini) <- c("IECs", "Damage", 
                    "BiomassCD", "CD", 
                    "Drug", 
                    "pheme_e", "pheme_c",
                    "cys_L_e", "leu_L_e", "ile_L_e", "val_L_e", "trp_L_e", "pro_L_e",
                    "cys_L_v", "leu_L_v", "ile_L_v", "val_L_v", "trp_L_v", "pro_L_v")
  
  y_ini = y_ini[yini.names]
  
  return(y_ini)
  
}

ComputeBiomassBac = function(diameter, len) {
  
  Vcd = ((pi*(diameter)^2)/4)*(((4*diameter)/6) + len) # (μm^3)
  mass = round((1/3)*Vcd*1.3, 4) # (pg)
  
  return(mass)
  
}

InflammationFunction = function(G) {
  
  W = 0.1 # villus width (mm)
  L = 1 # villus height (mm)
  EM = 13 # surface amplification due to microvilli
  
  VIarea = 2*pi*(W/2)*L*EM # (mm^2)
  
  Na = 6.02214154e+20 # (molecule/mmol)
  
  Vrbc = 11000 # red blood cell velocity (mm/h)
  RHOrbc = 6e+06 # [cell/μL]
  HeamRBC = 5e+09 # (molecule/cell)
  Dv = 0.03 # (mm)
  Aw = 113 # (mm^2)
  
  VA = pi*((Dv/2))^2 # vassel section (mm^2)
  
  Inflammation = G*VA*(Vrbc*RHOrbc*((HeamRBC/Na)*1e+09))*(Aw/VIarea) # (pmol)
  
  return(Inflammation)
  
}

EvalDiet = function(file, diet, ex, HemeDiet) {
  
  if (is.null(diet)) {
    
    # Mpheme = 616.5 (g/mol)
    Mpheme = 616.5*1e-12 # (g/pmol)
    
    # What percentage of haem iron is absorbed? (15% to 35%)
    
    # The median intake of Total, Heme, and Nonheme iron was:
    # 23.0, 0.75, and 22.0 (mg/day*person) among men, and:
    # 20.0, 0.63, and 19.1 (mg/day*person) among women
    
    if (HemeDiet == "median") {
      vheme = ((1/24)*1.5)*1e-03 # (g/day)
    } else {
      vheme = ((1/24)*750)*1e-03 # (g/day)
    }
    
    vD = (vheme/Mpheme) # (pmol)
    
  } else {
    
    Na = 6.02214154e+20 # (molecule/mmol)
    c = 6.022e8
    
    # diet = readr::read_delim(paste0("./Input/Diets/vmh/", diet, ".tsv", sep=""), 
    #                         "\t", escape_double = FALSE, trim_ws = TRUE)
    
    load(file)
    
    Delta = 1e-04
    
    diet = diets[[diet]]
    flux = -1*diet[["Flux Value"]][which(diet[["Reaction"]] %in% ex)]*Delta # (mmol)
    vD = flux*1e-09 # (pmol)
    
  }
  
  return(vD)
}

EvalTransport = function(yield) {
  
  IECs_t0 = 4.57e05
  vT = yield/IECs_t0
  
  return(vT)
}

EvalDetox = function(HsmA){
  
  HsmA = runif(n = 1, min = 0, max = 1e-24)
  
  # doi: 10.1289/ehp.856453.
  
  # k for reaction: PFe3 + RO2· -> [PFe3] + RO2-
  # it is: 
  # k = 2.5e+08*60 # [1/mmol*min]
  kHsmA = 2.5e+17*60*60 # [1/pmol*h]
  
  # HsmA [pmol] # HsmA = 1e-20 # [pmol/bac]
  
  # [pmol/h]
  detox = HsmA*kHsmA
  
  return(detox)
  
}

EvalEfflux = function(MDET, unit){
  
  # MDET = Multi Drug Efflux Transporter
  
  # https://doi.org/10.1101/2020.01.16.909507
  # k for reaction: drug_c -> drug_e
  V_max = 2.35e-11 # (1/mol*s)*1e+09*cell
  V_max = V_max/1e+09 # (1/cell*mol*s)
  # 2.35e-20 (1/cell*mol*s) -> 
  # 2.35e-20*1000*3600 (1/cell*mmol*min) -> 
  V_max = (V_max/1e+09)*1000*3600*60*1e+09 # (1/cell*pmol*h)
  
  ## ref (1) ##
  # https://doi.org/10.1101/2020.01.16.909507
  MDET.1 = 2.35e-12  # (mol/1e+09*cell) of transporter
  
  MDET.1 = MDET.1*1e+12 # (pmol/1e+09*cell) of transporter
  MDET.1 = MDET.1/1e+09 # (pmol/cell) of transporter
  
  ## ref (2) ##
  # https://www.nature.com/articles/srep21909
  # 1 (mmol) = 6.0221415e+20 (molecule)
  # MDET = 100 (molecule/cell) -> 1.660538863127e-19 (mmol/cell)
  # MDET = 1.660538863127e-19*1e+09 (mmol/1e+09*cell)
  MDET.2 = MDET # (mmol/1e+09*cell)
  
  MDET.2 = MDET.2*1e+12 # (pmol/1e+09*cell) of transporter
  MDET.2 = MDET.2/1e+09 # (pmol/cell) of transporter
  
  V_max = (V_max*MDET.2)/MDET.1 # (1/cell*pmol*h)
  
  if(unit == "molecule") {
    V_max = 602214150000*V_max # (1/cell*molecule*h)
  }
  
  efflux = V_max # (according to "unit")
  
  return(efflux)
  
}

EvalVmax = function(V, gamma, k) {
  
  # link: https://doi.org/10.1016/j.cell.2005.06.025
  
  # V = 3.1 [pmol/min/μg protein] -> 
  V = V*60*1e+06 # [pmol/h/g protein] 
  
  # gamma = proportion of mass protein which is HEMEti = 1e-03 (ref?)
  # k = proportion of reversibility [unit]
  
  Vcd = ((pi*(0.5)^2)/4)*(((4*0.5)/6) + 5.5) # (μm^3)
  gDW_CDmean = round((1/3)*Vcd*1.3, 4) # (pg)
  gDW_CDmean = gDW_CDmean*1e-12 # (pg)
  
  vmax = k*(V*gDW_CDmean*gamma) # (pmol/gTransporter*h)
  
  return(vmax)
  
}

MapParameter = function(optim_v, indx) {
  return(optim_v[indx])
}

treat_generation = function(marking, time, unit = "pmol") {
  
  MM = 171.16*1e-12 # (g/pmol)
  
  # ref (1) (DOI: 10.1371/journal.pone.0053757)
  MIC.conc = 1e-06 # (g/mL)
  MIC.conc = (MIC.conc/MM)*1e-09 # (mmol)
  MIC.conc = MIC.conc*1e+09 # (pmol)
  
  # dose = 0.0 # (g/day) - no dose -
  # dose = 0.5 # (g/day) - prescribed dose -
  dose = 0.05 # (g/day) - MIC dose -
  # dose = 0.01 # (g/day) - SubMIC dose -
  
  Delta = 1e-04 # (unit)
  
  dose = 8*(1/24)*((dose/MM))*Delta # (pmol/8*h)
  
  if (unit == "molecule") {
    dose = dose*602214150000 # (molecule/8*h)
    MIC.conc = MIC.conc*602214150000 # (molecule/8*h)
  }
  
  new_marking = marking[1:19]
  new_marking[10] = marking[10] + dose
  
  return(new_marking)
}

# targetFN = function(output) {
#   return(output[, "CD"])
# }

ReferenceM = function(reference, output) {
  mic.value = 6300; cd.value = 0.5e+08
  E1 = abs(output[output$Time == 36, "CD"] - cd.value)
  E2 = abs(output[output$Time == 72, "IECs"] - output[output$Time == 0, "IECs"])
  E3 = abs(output[output$Time == 72, "Drug"] - mic.value)
  return(sum(E1, E2, E3))
}

##################################################
##################################################

# M = c(0.11513, 0.13117, 0.13117, 0.117151, 0.20423, 0.12116) # [g/mmol]
# S = c(0.41248, 0.54554, 0.28608, 0.352610, 0.13306, 0.046571) # [unit]
# c = 6.022e8 # [molecule]
# Na = 6.022e20 # [molecule]
# DW_IECs = 1e-9 # [g/cell]

# car = round(((DW_IECs*S)/M)*(Na/c), 0) # ([g/cell]*[unit]/[g/mmol])*([molecule]/[molecule])
