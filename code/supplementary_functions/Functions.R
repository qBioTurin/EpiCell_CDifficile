
init.gen <- function(n_file, zero, zeroH) {
  
  # for debug:
  # 
  # n_file = paste0(wd, "/input/NAMES.RDS")
  # d = "EU_average"
  # zero = F
  # zeroH = F
  #
  # init = init.gen(n_file,
  #                 file,
  #                 d,
  #                 zero,
  #                 zeroH)
  
  yini.names = readRDS(n_file)
  
  # initial Intestinal Epithelial Cells number
  IECs = 457000 # (cell)
  MOI = 100
  # initial C. diff cells number
  nBac = MOI*457000 # (cell)
  # single cell C.diff mean observed biomass
  # gDW_CDmean = ComputeBiomassBac(diameter = 0.5, len = 5.5)
  gDW_CDmean = 0.05727 #(pg)
    
  if (zeroH == T) {
    initHeme_e = 1e-12 # (pmol/mL)
  } else {
    pheme_e.i = 100 # (μM)
    # first conversion: 1 (μM) = 1.0e-06 (mmol/mL)
    pheme_e.i = pheme_e.i*1.0e-06 # (mmol/mL)
    # second conversion: 1 (mmol/mL) = 1.0e+09 (pmol/mL)
    pheme_e.i = pheme_e.i*1e+09 # (pmol/mL)
    initHeme_e = pheme_e.i # (pmol/mL)
  }
  
  initHeme_c = 1e-03 # (pmol/cell)
  initHeme_c = initHeme_c*nBac # (pmol)
  
  # scaled avogadro number per mmol
  Na = 6.022e+20 # (molecule/mmol)
  # in the model aminoacids amounts are measured as (c_molecule) i.e. molecule packet
  c = 6.022e+08 # (molecule)
  # Literature information is units of (mmol)
  # quick conversion: 1 (mmol) = 6.02214154+20 (molecule)
  # 6.02214154+20 (molecule) = Na (Avogadro number reduced)
  # packaging molecule
  pack = 1*(Na*(1/c))
  
  # aa concentration from literature
  iniAA_bib = c(0.25, 0.57, 0.32, 0.58, 0.098, 0.36) # (μmol/mL)
  iniAA_bib = iniAA_bib*0.001  # (mmol/mL)
  iniAA_bib = iniAA_bib*pack  # (c_molecule)
  
  if (zero == T) {
    vec_iniAA = iniAA_bib*1e-25 # (c_molecule)
  }
  
  y_ini <- c(IECs, 0, 
             gDW_CDmean, nBac, 
             0, 
             initHeme_e, initHeme_c, 
             # gut comp aminoacids
             iniAA_bib,
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

EvalDiet = function(file, d, ex, HemeDiet) {
  
  Na = 6.022e+20 # (molecule/mmol)
  c = 6.022e+08 # (molecule)
  pack = 1*(Na*(1/c)) # (c_molecule)
  
  mm_heam = 616.5 # (g/mol)
  Mpheme = mm_heam*1e-12 # (g/pmol)
  
  Delta_aa = 5e-04
  Delta_heme = 1e-04
  
  if (d == "No_diet") {
    
    mean_men  = ((20.0 + 0.63 + 19.1)/3) # (mg/day*person)
    mean_women = (((23.0 + 0.75 + 22.0)/3)) # (mg/day*person)
    
    if (HemeDiet == "median") {
      mean_heme = (mean_men + mean_women)/2 # (mg/day*person)
      mean_heme = (mean_heme*1e-03)/24 # (g/h*person)
      vheme = mean_heme/mm_heam ## (mmol/h*person)
    } else {
      mean_heme = 750 # (mg/day*person)
      mean_heme = (mean_heme*1e-03)/24 # (g/h*person)
      vheme = mean_heme/mm_heam # (mmol/h*person)
    }
    
    vheme = vheme*Delta_heme # (mmol/h*ml)
    vD = vheme*1e+9 # (pmol/h*ml)
    
  } else {
    
    load(file)
    
    diet = diets[[d]]
    flux = diet[["Flux Value"]][which(diet[["Reaction"]] %in% ex)] # (mmol/h*person)
    
    flux = flux*Delta_aa # (mmol/h*ml)
    vD = flux*pack # (c_molecule/h*ml)
    
  }
  return(vD)
}

ComputeBiomassBac = function(diameter, len) {
  
  # The R function ComputeBiomassBac calculates the biomass (pg) of a bacterium given: 
  # 1) diameter (μm)
  # 2) length (μm)
  
  # The bacterium is assumed to be shaped like a cylinder with hemispherical ends
  Vcd = ((pi*(diameter)^2)/4)*(((4*diameter)/6) + len) # (μm^3)
  
  # where Vcd is volume of the bacillus (i.e. rod-shaped) C. difficile 
  # The formula used here is the volume of a cylinder (πr²h), 
  # where r is the radius (diameter/2), and h is the height. 
  # However, the height is calculated as the sum of the length of the cylinder (len) 
  # and the diameter of the two hemispheres at the ends of the cylinder
  
  # The density of a bacterial cell, when considering its wet mass, is approximately 1 g/cm³
  rho_bac = 1 # (pg/μm^3)
  # calucutating the mass of the bacterium
  mass = round((1/3)*Vcd*rho_bac, 4) # (pg)
  
  # The dry weight of a prototypical bacterial cell, 
  # is usually considered to be about 10-20% of its wet weight
  mass = mass*0.15 # (pg)
  
  return(mass)
  
}

InflammationFunction = function(G) {
  
  # fenestration parameter
  G = 0.06 # (unit)
  
  W = 0.1 # villus width (mm)
  L = 1 # villus height (mm)
  EM = 13 # surface amplification due to microvilli
  
  VIarea = 2*pi*(W/2)*L*EM # (mm^2)
  
  Vrbc = 11000 # red blood cell velocity (mm/h)
  RHOrbc = 6e+06 # (cell/μL)
  HeamRBC = 5e+09 # (molecule/cell)
  Dv = 0.03 # (mm)
  Aw = 113 # (mm^2)
  
  VA = pi*((Dv/2))^2 # vassel section (mm^2)
  Na = 6.022e+20  # (molecule/mmol)
  heam_m = ((HeamRBC/Na)*1e+09) # (mmol)
  # Amount of hemoglobin (in mmol) passing through the blood vessel per unit time.
  g = G*VA*(Vrbc*RHOrbc*heam_m)
  # proportion of the blood vessel wall that is in contact with the villi.
  A = (Aw/VIarea)
  # inflammation level
  Inflammation = g*A # (pmol/h)
  
  return(Inflammation)
  
}

EvalTransport = function(yield) {
  
  # initial Intestinal Epithelial Cells number
  IECs_t0 = 457000 # (cell)
  vT = yield/IECs_t0
  
  return(vT)
}

EvalDetox = function(k){
  
  k = runif(n = 1, min = 0, max = 1e-04) # ((pmol/mL)*h)
  
  Vcd_mean = 1.145372 # (μm^3)
  Vcd_mean_ml = Vcd_mean*1e-12 # (mL)
  
  # the enzyme concentration HsmA
  HsmA = 2e+04 # (pmol/L)
  HsmA = HsmA*1e-03 # (pmol/mL)
  HsmA = HsmA*Vcd_mean_ml # (pmol/cell)
  
  detox = k*HsmA # (1/cell*h)
  
  return(detox)
  
}

EvalEfflux = function(){
  
  # Define the number of molecules
  MDET = 2.35e-12 # (mol/1e+09*cell)
  MDET = MDET*1e-09 # (mol/cell)
  MDET = MDET*1e+12 # (pmol/cell)
  
  # https://doi.org/10.1101/2020.01.16.909507
  V_max = 2.35e-11 # (1/mol*s)*1e+09*cell
  V_max = V_max*1e-09 # (1/cell*mol*s)
  V_max = V_max*60*60 # (1/cell*mol*h)
  V_max = V_max*1e+12 # (1/cell*pmol*h)
  
  # Calculate efflux rate
  efflux = V_max * MDET 
  
  return(efflux)
  
}

EvalVmax = function(k) {
  
  gDW_CDmean = 0.05727 # (pg)
  
  Vmax = 3.1 # (pmol/min*μg_protein)
  Vmax = Vmax*60 # (pmol/h*μg_protein)
  Vmax = Vmax*1e-06 # (pmol/h*pg_protein)
  
  # Calculate gamma
  prop_of_prot_in_DW = 0.55 # (unit)
  # Convert the average protein content from percentage to mass
  protein_mass = gDW_CDmean*prop_of_prot_in_DW # (pg)
  
  MW_prot = 50 # (kDa)
  MW_prot = MW_prot*1e+03 # (Da)
  avg_protein_weight = MW_prot*1.66e-24*1e+12 # (pg)
  
  num_proteins = protein_mass/avg_protein_weight
  
  frac = 0.001 # (unit)
  num_transporter_proteins = num_proteins*frac
  
  gamma = num_transporter_proteins / num_proteins # (unit)
  
  gDW_protein = gDW_CDmean*gamma # (pg_of_protein/cell)
  
  vmax = k*Vmax*gDW_protein
  
  return(vmax)
  
}

MapParameter = function(optim_v, indx) {
  return(optim_v[indx])
}

treat_generation = function(marking, time) {
  
  mm_drug = 171.16 # (g/mol)
  mm_drug = mm_drug*1e-03 # (g/mmol)
  
  # ref (1) (DOI: 10.1371/journal.pone.0053757)
  MIC.conc = 1e-06 # (g/mL)
  MIC.conc = (MIC.conc/mm_drug) # (mmol)
  
  # dose = 0.0 # (g/day*person) - no dose -
  # dose = 0.5 # (g/day*person) - prescribed dose -
  dose = 0.05 # (g/day*person) - MIC dose -
  # dose = 0.01 # (g/day*person) - SubMIC dose -
  
  dose = dose/24 # (g/h*person)
  dose = dose/mm_drug # (mmol/h*person)
  
  Delta_drug = round(MIC.conc/dose, 4)
  
  dose = dose*Delta_drug # (mmol/h*mL)
  dose = 8*dose*1e+09 #  (pmol/8*h)
  
  MIC.conc = MIC.conc*1e+09 # (pmol)
  
  new_marking = marking[1:19]
  new_marking[10] = marking[10] + dose
  
  return(new_marking)
}

ReferenceM = function(reference, output) {
  mic.value = 5842.5 # (pmol)
  cd.value = 0.5e+08 # (cell)
  E1 = abs(output[output$Time == 36, "CD"] - cd.value)
  E2 = abs(output[output$Time == 72, "IECs"] - output[output$Time == 0, "IECs"])
  E3 = abs(output[output$Time == 72, "Drug"] - mic.value)
  return(sum(E1, E2, E3))
}

# targetFN = function(output) {
#   return(output[, "CD"])
# }
