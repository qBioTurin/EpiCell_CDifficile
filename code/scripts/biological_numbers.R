
####### Initial conditions

gut_V = 1 # (ml)

essential_aa = c("EX_pro_L_e", "EX_leu_L_e", "EX_ile_L_e", 
                 "EX_val_L_e", "EX_trp_L_e", "EX_cys_L_e")

# single cell C.diff minimum observed biomass
gDW_CDmin = ComputeBiomassBac(diameter = 0.3, len = 2) # (pg)
# single cell C.diff mean observed biomass
gDW_CDmean = ComputeBiomassBac(diameter = 0.5, len = 5.5) # (pg); diameter(μm), len (μm)
# # single cell C.diff maximun observed biomass
gDW_CDmax = ComputeBiomassBac(diameter = 0.7, len = 9) # (pg)

# initial Intestinal Epithelial Cells number
IECs = 457000 # (cell)
# Multiplicity of infection (MOI)
MOI = 100 # (bac)
# initial C. diff cells number
nBac = MOI*457000 # (cell)

# The bacterium is assumed to be shaped like a cylinder with hemispherical ends
diameter_mean = 0.5 # (μm)
len_mean = 5.5 # (μm)
Vcd = ((pi*(diameter_mean)^2)/4)*(((4*diameter_mean)/6) + len_mean) # (μm^3)

# the assumption of 100 μM ( from -> ( https://doi.org/10.1038/s41419-020-02950-8 ) ) 
# heme concentration in the human colon seems a reasonable starting point 
# for a model or simulation, but it’s important to note that actual heme concentrations 
# can vary widely depending on a number of factors
# (e.g. diet, gastrointestinal disorders, and individual differences in metabolism and absorption)
# [ from -> ( https://doi.org/10.1158/0008-5472.CAN-14-2554 ) ] 
pheme_e.i = 100 # (μM)
# first conversion: 1 (μM) = 1e-06 (mmol/mL)
pheme_e.i = pheme_e.i*1e-06 # (mmol/mL)
# second conversion: 1 (mmol/mL) = 1.0e+09 (pmol/mL)
pheme_e.i = pheme_e.i*1e+09 # (pmol/mL)

# Note: hemorragic heme = 1000 (μM) = 1000 * 1.0e-06 * 1e+09 (pmol/mL)
# The determination of heme concentrations in different conditions range from 20 to > 350 μmol/L as occurring

initHeme_c = 1e-03 # (pmol/cell)
initHeme_c = initHeme_c*nBac # (pmol)

# calculating the heme concentration per bacterium. 
# the negative sign indicates consumption flux, meaning that the heme can be uptaken by the bacteria.
pheme_c.lb = -(initHeme_c*1e-09)/(nBac*gDW_CDmean*1e-12) # (mmol/gDW*h)
# " pheme_c.lb "value is then set as the lower bound for the heme concentration in the FBA model. 
# This means that in the simulations run using this model, the heme concentration cannot go below this value.

# scaled avogadro number per mmol
Na = 6.022e+20 # (molecule/mmol)
# in the model aminoacids amounts are measured as (c_molecule) i.e. molecule packet
c = 6.022e+08 # (molecule)
# Literature information is units of (mmol)
# quick conversion: 1 (mmol) = 6.02214154+20 (molecule)
# 6.02214154+20 (molecule) = Na (Avogadro number reduced)
# packaging molecule
pack = 1*(Na*(1/c)) # (molecule)

# aa concentration from literature
iniAA_bib = c(0.25, 0.57, 0.32, 0.58, 0.098, 0.36) # (μmol/mL)
iniAA_bib = iniAA_bib*0.001  # (mmol/mL)
iniAA_bib = iniAA_bib*pack  # (c_molecule)

iniAA.lb = - (iniAA_bib/pack)/((nBac*(gDW_CDmean*1e-12)))

####### Dietary intake metabolites ( EvalDiet )

# scaled avogadro number per mmol
Na = 6.022e+20 # (molecule/mmol)
c = 6.022e+08 # (molecule)
pack = 1*(Na*(1/c)) # (molecule)

# Rescaling factor from (mmol/h*person) to (mmol/h*ml)

# aminoacid concentration from literature
iniAA_bib = c(0.25, 0.57, 0.32, 0.58, 0.098, 0.36) # (μmol/mL)
iniAA_bib = iniAA_bib*0.001  # (mmol/mL)

# heme concentration from literature
pheme_e.i = 100 # (μM)
# first conversion: 1 (μM) = 1e-06 (mmol/mL)
pheme_e.i = pheme_e.i*1e-06 # (mmol/mL)

P <<- -180

# Introducing the Delta conversion:
#
# Approach:
#
# 1) Data Retrieval:
# the dietary intake of 6 essential amino acids (vec_iniAA) from a nutritional database.
#
# 2) Literature-based Amino Acid Concentrations:
# the amino acid concentrations in the gut from the literature (iniAA_bib)
#
# 3) Assumption of Steady State Concentration:
# the assuming, at steady state, the concentration of amino acids in the gut is equal 
# to the ratio of the literature-based concentration to the intake rate. 
# This is a reasonable assumption under the assumption that the body is in a state 
# where the intake and removal of amino acids are balanced, resulting in a stable concentration.
#
# 4) Conversion Factor (Delta):
# introduction of a conversion factor (Delta) to convert the dietary intake rate 
# from mmol/h*person to mmol/h*mL:

# aminoacid flux from nutrition data
file = paste0(wd, "/input/diets/diets.RData")
load(file)

EU_average = diets[["EU_average"]]
fluxAA = EU_average[match(essential_aa, EU_average$Reaction), ]$`Flux Value` # (mmol/h*person)

Delta_aa = round(mean((iniAA_bib/fluxAA)), 4)

fluxAA = fluxAA*Delta_aa # (mmol/h*ml)
vAA = fluxAA*pack # (c_molecule/h*ml)

# heme flux from nutrition data

mm_heam = 616.5 # (g/mol)
mm_heam = mm_heam*0.001 # (g/mmol)

HemeDiet = "median"

if (HemeDiet == "median") {
  # The median intake of Total, Heme, and Nonheme iron was:
  # 23.0, 0.75, and 22.0 (mg/day*person) among men, and:
  # 20.0, 0.63, and 19.1 (mg/day*person) among women
  mean_men  = ((20.0 + 0.63 + 19.1)/3) # (mg/day*person)
  mean_women = (((23.0 + 0.75 + 22.0)/3)) # (mg/day*person)
  mean_heme = (mean_men + mean_women)/2 # (mg/day*person)
  mean_heme = (mean_heme*1e-03)/24 # (g/h*person)
  vheme = mean_heme/mm_heam ## (mmol/h*person)
} else {
  mean_heme = 750 # (mg/day*person)
  mean_heme = (mean_heme*1e-03)/24 # (g/h*person)
  vheme = mean_heme/mm_heam # (mmol/h*person)
}

Delta_heme = round(pheme_e.i/vheme, 4)
vheme = vheme*Delta_heme # (mmol/h*ml)
vheme = vheme*1e+9 # (pmol/h*ml)

# Unit of Delta: converting dietary intake from mmol/h*person to mmol/mL, which is a concentration. 
# This assumes that the rate of intake is evenly distributed across the 1 mL of gut volume.
# If this assumption holds, the conversion factor Delta could be a valid approach.

######## Compute bacterial biomass

# The R function ComputeBiomassBac calculates the biomass (pg) of a bacterium given: 
# 1) diameter (μm)
# 2) length (μm)

diameter_mean = 0.5 # (μm)
len_mean = 5.5 # (μm)
# The bacterium is assumed to be shaped like a cylinder with hemispherical ends
Vcd = ((pi*(diameter_mean)^2)/4)*(((4*diameter_mean)/6) + len_mean) # (μm^3)

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

######## Inflammation function

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

# " Inflammation " is then used to compute rate of inflammation
# Calculate the rate based on the percentage of damage
DAMAGEmax = 457000 # (cell)
DamagePlace <- seq(0, DAMAGEmax, by = 1000) # (cell)
PercDamage = DamagePlace / DAMAGEmax # (unit)
rate = PercDamage * Inflammation # (pmol/h)

inflam = ggplot(data.frame(DamagePlace, rate), aes(x=DamagePlace, y=rate)) +
  geom_line() +
  labs(x = "Damage Place (cell)", 
       y = "Rate of Inflammation (pmol/h)", 
       title="Inflammation Rate")

######## EvalTransport

# Given values
yield_aa = 0.99 # (unit)
yield_heme = 0.15 # (unit)

IECs_t0 = 457000 # (cell)

pheme_e.i = 100 # (μM)
pheme_e.i = pheme_e.i*1e-06 # (mmol/mL)
pheme_e.i = pheme_e.i*1e+09 # (pmol/mL)

iniAA_bib = c(0.25, 0.57, 0.32, 0.58, 0.098, 0.36) # (μmol/mL)
iniAA_bib = iniAA_bib*0.001  # (mmol/mL)
Na = 6.022e+20 # (molecule/mmol)
c = 6.022e+08 # (molecule)
pack = 1*(Na*(1/c)) # (molecule)
iniAA_bib = iniAA_bib*pack # (c_molecule)
aa_e.i = iniAA_bib[1] # (c_molecule)

vT_aa = yield_aa/IECs_t0 # (1/cell)
vT_heme = yield_heme/IECs_t0 # (1/cell)

IECs = seq(0, 457000, by = 1000) # (cell)

transport_aa = vT_aa*aa_e.i*IECs # (c_molecule/h)
transport_heme = vT_heme*pheme_e.i*IECs # (pmol/h)

# Plotting
t_aa = ggplot(data.frame(IECs, transport_aa), aes(x=IECs, y=transport_aa)) +
  geom_line() +
  labs(x = "IECs (cell)", 
       y = "Transport (c_molecule/h)", 
       title="Transport of Amino Acids")

t_heme = ggplot(data.frame(IECs, transport_heme), aes(x=IECs, y=transport_heme)) +
  geom_line() +
  labs(x = "IECs (cell)", 
       y = "Transport (pmol/h)", 
       title="Transport of heme")

trasport = t_aa + t_heme
  
######## EvalDetox

# This system repurposes heme to counteract antimicrobial oxidative stress responses
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7486240)
# This research provides insights into how C. difficile exploits the toxic heme within 
# the inflamed gut as a shield against antimicrobial compounds

# \text{HsmA} + \text{pheme_c} + \text{Drug} \xrightarrow[\text{}]{k_{\text{HsmA}}} \text{HsmA}

# This equation suggests that the HsmA protein interacts with heme and the drug, 
# and then returns to its original state, presumably after having altered 
# the heme and/or drug in some way to reduce their toxicity or reactivity. 
# This is consistent with the findings of the research paper, 
# which suggest that HsmA uses heme to protect C. difficile from oxidative stress.

# In enzymatic reactions, the rate constant can be described using the Michaelis-Menten kinetics. 
# The general reaction scheme of an enzyme-catalyzed reaction is as follows:

# E + S → k1 [ES] → k2 E + P

# In mass action kinetics, the rate of a reaction is proportional 
# to the product of the concentrations of the reactants. 
# Each reactant concentration is raised to the power of its stoichiometric coefficient.

# Define the parameters

# For the rate constant k, we can refer to the general knowledge of reaction rate constants. 
# For a unimolecular step, a unimolecular rate constant has an upper limit of k1 ≤ ~10^13 s^-1. 
# For a bimolecular step, a bimolecular rate constant has an upper limit of k2 ≤ ~10^10 M^-1 s^-1

# Define the function to calculate detox_rate
calculate_detox_rate <- function(nBac, k) {
  
  Vcd_mean = 1.145372 # (μm^3)
  Vcd_mean_ml = Vcd_mean*1e-12 # (mL)
  
  # the enzyme concentration HsmA
  HsmA = 2e+04 # (pmol/L)
  HsmA = HsmA*1e-03 # (pmol/mL)
  HsmA = HsmA*Vcd_mean_ml # (pmol/cell)
  # 1 pmol are 602214150000 molecules
  mol = HsmA*602214150000 # (molecule)
  
  # Define the mass action equation
  detox = k*HsmA # (1/cell*h)
  
  initHeme_c = 45700 # (pmol/mL)
  
  mm_drug = 171.16 # (g/mol)
  mm_drug = mm_drug*1e-03 # (g/mmol)
  
  MIC.conc = 5842.487 # (pmol/mL)
  MIC.conc = MIC.conc*1e-09 # (mmol/mL)
  
  dose = 0.05 # - MIC dose - (g/day*person)
  dose = dose/24 # (g/h*person)
  dose = dose/mm_drug # (mmol/h*person)
  
  Delta_drug = round(MIC.conc/dose, 4)
  
  dose = dose*Delta_drug # (mmol/h*mL)
  dose = 8*dose*1e+09 #  (pmol/8*h)
  
  detox_rate = detox*dose*initHeme_c*nBac
  return(detox_rate)
}

k_values <- runif(n = 24, min = 0, max = 1e-04)

# Create a data frame to store results
results <- data.frame()
nBac <- seq(0, 45700000, by = 1000) # (cell)

# Calculate detox_rate for each k and store in the data frame
for (k in k_values) {
  detox_rate <- calculate_detox_rate(nBac, k)
  results <- rbind(results, data.frame(nBac, detox_rate, k))
}

# Plot the results
detox_plot <- ggplot(results, aes(x=nBac, y=detox_rate, group = k, color=k)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Number of Bacteria (cell)", 
       y = "Detoxification Rate (pmol*pmol/h)",
       title="Detox Rate") +
  geom_hline(aes(yintercept=48687.39), color="blue", linetype="dashed")+
  scale_color_gradient(low = "blue", high = "red", name = "Detox")

######## EvalEfflux

# Multi Drug Efflux Transporter (MDET) 
# https://www.nature.com/articles/srep21909
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
# The efflux rate is calculated as the product of:
# - V_max: The maximum efflux velocity per transporter, 
# representing the maximum rate at which a single transporter can move the drug out of the cell.
# - MDET: The amount of active transporter per cell, 
# indicating the availability of transporters to perform the efflux process.

Drug = 48687.39 # (pmol)

nBac = seq(0, 45700000, by = 1000) # (cell)

efflux_rate = Drug*nBac*efflux # (pmol/h)

efflux_plot = ggplot(data.frame(nBac, efflux_rate), aes(x=nBac, y=efflux_rate)) +
  geom_line() +
  labs(x = "Number of Bacteria (cell)", 
       y = "Efflux Rate (pmol/h)", 
       title="Transport of drug")

####### EvalVmax

# link: https://doi.org/10.1016/j.cell.2005.06.025
# Subtraction of the endogenous uptake from the heme transport-induced uptake 
# revealed a saturable component with an apparent
# Vmaxapp of 3.1 pmol/min/μg protein and Kmapp of 125 μM heme

k = 0.9 # (unit)

gDW_CDmean = 0.05727 # (pg)

Vmax = 3.1 # (pmol/min*μg_protein)
Vmax = Vmax*60 # (pmol/h*μg_protein)
Vmax = Vmax*1e-06 # (pmol/h*pg_protein)

# Given that the average protein content of a bacterial cell is around 55% of the dry weight 
# and the average molecular weight of proteins in the cell is around 50 kDa,

# Calculate gamma
prop_of_prot_in_DW = 0.55 # (unit)
# Convert the average protein content from percentage to mass
protein_mass = gDW_CDmean*prop_of_prot_in_DW # (pg)

# Convert the average molecular weight of proteins from kDa to Da (1 kDa = 1000 Da) 
# and then to pg (1 Da = 1.66 x 10^-24 g):
MW_prot = 50 # (kDa)
MW_prot = MW_prot*1e+03 # (Da)
avg_protein_weight = MW_prot*1.66e-24*1e+12 # (pg)

num_proteins = protein_mass/avg_protein_weight

# Assume the number of transporter proteins is a small fraction of the total number of proteins
frac = 0.001 # (unit)
num_transporter_proteins = num_proteins*frac

# If we assume that the transporter protein has a similar weight to the average protein, 
# then the mass of the transporter protein in a unit of bacterial mass (gamma) 
# can be calculated as the ratio of the number of transporter proteins to 
# the total number of proteins in the cell.
gamma = num_transporter_proteins / num_proteins # (unit)

gDW_protein = gDW_CDmean*gamma # (pg_of_protein/cell)

rate_intake = k*Vmax*gDW_protein
rate_efflux = (1 - k)*Vmax*gDW_protein

nBac = seq(0, 45700000, by = 1000) # (cell)
pheme_e.i = 1e+05 # (pmol)
pheme_c.i = 45700 # (pmol)

# Mass Action Kinetics
intake = rate_intake*nBac*pheme_e.i
efflux = rate_efflux*nBac*pheme_c.i

# Create the ggplot objects
plot_intake <- ggplot(data.frame(nBac, intake), aes(x=nBac, y=intake)) +
  geom_line(color="blue") +
  labs(x="Number of Bacteria (cell)", 
       y="Heme Intake Rate (pmol/h)", 
       title="Heme Intake")

plot_efflux <- ggplot(data.frame(nBac, efflux), aes(x=nBac, y=efflux)) +
  geom_line(color="red") +
  labs(x="Number of Bacteria (cell)", 
       y="Heme Efflux Rate (pmol/h)", 
       title="Heme Efflux")

# Combine the plots side by side using patchwork
combined_plot <- plot_intake + plot_efflux

####### treat_generation

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

######### aminoacids releasing arcs cardinality

# the number of molecules of each amino acid per cell

# Define constants
Na = 6.022e+20 # Scaled Avogadro number per mmol (molecule/mmol)
c = 6.022e+08 # Packaging unit (molecule)
pack = 1*(Na*(1/c)) # Packaging molecule (molecule/c_molecule)

# Molar mass of each amino acid
vec_MMAA = c(115.13, 131.17, 131.17, 117.151, 204.23, 121.16) # (g/mol)
vec_MMAA = vec_MMAA*1e-03 #  (g/mmol)
# Stoichiometry
stoich = c(0.41248, 0.54554, 0.28608, 0.352610, 0.13306, 0.046571) # (unitless)
# Dry weight of IECs
DW_IECs = 1e-09 # (g/cell)
aa_prop = DW_IECs*stoich # (g/cell)

# Number of molecules of each amino acid per cell
car = round(((aa_prop)/vec_MMAA)*pack, 0) # (c_molecule/cell)

######## FBA rates conversion

Na = 6.022e+20 # (molecule/mmol)
c = 6.022e+08 # (molecule)
pack = 1*(Na*(1/c)) # (molecule)

MWbio = 1 # (g/mmol)

nBac = 45700000 # (cell)

gDW_CDmean = 0.05727 # (pg)
Biom = gDW_CDmean
Biom_g = Biom*1e-12 # (g)
total_Biom_g = nBac*Biom_g # (g)

flux_biom = 0.2 # flux (fill from FBA estimates) (mmol/gDW*h)
flux_biom = flux_biom*MWbio # (g/gDW*h)
flux_biom = flux_biom*1e+12 # (pg/gDW*h)

flux_aa = rep(-1000, 6) # flux (fill from FBA estimates) (mmol/gDW*h)
flux_heme = -100 # flux (fill from FBA estimates) (mmol/gDW*h)

r_aa = flux_aa*pack*total_Biom_g # (c_molecule/h)
r_heme = flux_heme*1e+09*total_Biom_g # (pmol/h)
r_biom = flux_biom*Biom*1e-12 # (pg/h)

# aa concentration from literature
iniAA_bib = c(0.25, 0.57, 0.32, 0.58, 0.098, 0.36) # (μmol/mL)
iniAA_bib = iniAA_bib*0.001  # (mmol/mL)
iniAA_bib = iniAA_bib*pack  # (c_molecule)

initHeme_c = 1e-03 # (pmol/cell)
initHeme_c = initHeme_c*nBac # (pmol)

######## FBA bounds conversion

# The lower bound (Lb) in Flux Balance Analysis (FBA) is not necessarily limited
# by the concentration of the metabolite in the environment. 
# Instead, it represents the rate at which a metabolite is being consumed or produced
# in the system.

# the high negative value for Lb indicates that the met is being consumed rapidly
# by the bacteria. This makes sense given the small biomass (t_Biom_g) and the relatively
# large amount of met.

# the unit for the lower bound in FBA is mmol/gDW/h, which is a rate, not a concentration.
# So, a high absolute value for the lower bound simply means that the metabolite
# is being consumed or produced at a high rate relative to the biomass.

calculate_Lb_aa <- function(Biom, Met, nBac) {
  
  Na = 6.022e+20 # (molecule/mmol)
  c = 6.022e+08 # (molecule)
  pack = 1*(Na*(1/c)) # (molecule)
  
  epsilon = 1e-24 # (pg)
  gDW_CDmin = 0.0077 # (pg)
  
  # Calculate Lb
  Lb = ifelse(nBac < 1, 
              -( (Met/pack)/((Biom*1e-12) + epsilon) ), 
              -( (Met/pack)/((nBac*(Biom*1e-12)) + epsilon) ) )
  
  Lb = ifelse(Biom <= gDW_CDmin,
              ( Lb * Biom/gDW_CDmin ),
              Lb )
  
  Lb = ifelse(Lb <= -30,
              - 30,
              Lb )
  
  return(Lb) # (mmol/gWD*h)
}

# nBac = 0.5 # (cell)
nBac = 45700000 # (cell)
gDW_CDmax = 0.18216 # (pg)

# Create a sequence of values for Biom and Met
Biom <- seq(0, gDW_CDmax, length.out = 12) # (pg)

# aa concentration from literature
iniAA_bib = c(0.25, 0.57, 0.32, 0.58, 0.098, 0.36) # (μmol/mL)
iniAA_bib = iniAA_bib*0.001  # (mmol/mL)
iniAA_bib = iniAA_bib*pack  # (c_molecule)

Met <- seq(0, mean(iniAA_bib), length.out = 1000) # (c_molecule)

# Create an empty list to store plots
plots <- list()
Lb_max = Lb_min = c()

for( i in 1:length(Biom) ) {
  Lb_max[i] <- max(sapply(Met, calculate_Lb_aa, Biom = Biom[i], nBac))
  Lb_min[i] <- min(sapply(Met, calculate_Lb_aa, Biom = Biom[i], nBac))
}

# Create an empty data frame to store all Lb values for different Biom
df_all <- data.frame()

# Loop over Biom values
for( i in 1:length(Biom) ) {
  # Calculate Lb for each Met value
  Lb <- sapply(Met, calculate_Lb_aa, Biom = Biom[i], nBac)
  
  # Create a temporary data frame for plotting
  df_temp <- data.frame(Biom = Biom[i], Met = Met, Lb = Lb)
  
  # Bind the temporary data frame to the main data frame
  df_all <- rbind(df_all, df_temp)
}

# Create a ggplot object for the multiline plot
p_aa <- ggplot(df_all, aes(x = Met, y = Lb, group = Biom, color = Biom)) +
  geom_line() +
  scale_color_gradient(low = "blue", high = "red", name = "Biom (pg)") +
  labs(x = "aa_e (c_molecule)", y = "Lb (mmol/gDW*h)", title = "EX_aa_L_e Lower Bound") +
  theme(legend.position = "bottom")

calculate_Lb_pheme <- function(Biom, Met, nBac) {
  
  epsilon = 1e-24 # (pg)
  gDW_CDmin = 0.0077 # (pg)
  
  # Calculate Lb
  Lb = ifelse(nBac < 1, 
              -( (Met*1e-09)/((Biom*1e-12) + epsilon) ), 
              -( (Met*1e-09)/((nBac*(Biom*1e-12)) + epsilon) ) )
  
  Lb = ifelse(Biom <= gDW_CDmin,
              ( Lb * Biom/gDW_CDmin ),
              Lb )
  
  Lb = ifelse(Lb <= -3,
              - 3,
              Lb )
  
  return(Lb) # (mmol/gWD*h)
}

# nBac = 0.5 # (cell)
nBac = 45700000 # (cell)
gDW_CDmax = 0.18216 # (pg)

# Create a sequence of values for Biom and Met
Biom <- seq(0, gDW_CDmax, length.out = 12) # (pg)

initHeme_c = 1e-03 # (pmol/cell)
initHeme_c = initHeme_c*nBac # (pmol)
Met = seq(0, initHeme_c, length.out = 1000) # (pmol)

# Create an empty list to store plots
plots <- list()
Lb_max = Lb_min = c()

for( i in 1:length(Biom) ) {
  Lb_max[i] <- max(sapply(Met, calculate_Lb_pheme, Biom = Biom[i], nBac))
  Lb_min[i] <- min(sapply(Met, calculate_Lb_pheme, Biom = Biom[i], nBac))
}

# Create an empty data frame to store all Lb values for different Biom
df_all <- data.frame()

# Loop over Biom values
for( i in 1:length(Biom) ) {
  # Calculate Lb for each Met value
  Lb <- sapply(Met, calculate_Lb_pheme, Biom = Biom[i], nBac)
  
  # Create a temporary data frame for plotting
  df_temp <- data.frame(Biom = Biom[i], Met = Met, Lb = Lb)
  
  # Bind the temporary data frame to the main data frame
  df_all <- rbind(df_all, df_temp)
}

# Create a ggplot object for the multiline plot
p_heme <- ggplot(df_all, aes(x = Met, y = Lb, group = Biom, color = Biom)) +
  geom_line() +
  scale_color_gradient(low = "blue", high = "red", name = "Biom (pg)") +
  labs(x = "pheme_c (pmol)", y = "Lb (mmol/gDW*h)", title = "sink_pheme_c Lower Bound") +
  theme(legend.position = "bottom")

######## FBA bounds conversion (biomass)

# ref -> https://doi.org/10.1016/j.ymben.2021.10.012
# In this effort we considered the molecular weight (MW) of the biomass metabolite 
# equal 1 (g dry weight (gDW) mmol-1) (https://doi.org/10.1038/nprot.2009.203), 
# which is needed to correctly map biomass flux to specific growth rate 
# (https://doi.org/10.1093/bioinformatics/btx453).

# The approach is setting the growth limit for the biomass production reaction 
# considering the maximum observed dry weight biomass (gDW_CDmax) and the 
# average dry weight biomass (Biom) (Biom could be different from gDW_CDmean) 
# of the bacteria population. The upper bound (Ub) for the biomass production reaction 
# is set based on these values



required <- 0.14
max <- 0.2
stiff <- 150
tB <- 1e-12

UbFunction <- function(Biom, nBac) {
  
  Ub <- 0.0
  
  if (nBac < 1) {
    Ub <- Biom * nBac
  } else {
    if ((gDW_CDmax - Biom) > tB) {
      Ub <- (gDW_CDmax - Biom)
    } else {
      Ub <- tB
    }
    
    Ub <- ifelse(Ub > max, max, Ub)
    
    if (Biom < gDW_CDmin) {
      Ub <- 0.0
    } else if (Biom == gDW_CDmin) {
      Ub <- max
    }
    
    Ub <- Ub / (1 + exp(-stiff * (Biom - gDW_CDmin)))
    Ub <- Ub * (1 - 1 / (1 + exp(-stiff * (Biom - required))))
  }
  
  return(Ub)
}

Biom_values <- seq(0, gDW_CDmax, length.out = 100)
nBac_values <- seq(0, nBac, length.out = 100)
grid <- expand.grid(Biom = Biom_values, nBac = nBac_values)
grid$Ub <- mapply(UbFunction, grid$Biom, grid$nBac)

# Convert the grid to a matrix for plotting
Ub_matrix <- matrix(grid$Ub, nrow = length(Biom_values), byrow = TRUE)

library(plotly)

# Assuming Ub_matrix is the matrix of Ub values obtained from applying the UbFunction to the grid of Biom and nBac values.

# Define the color scale
colorscale <- list(
  list(0, "blue"),
  list(0.5, "yellow"),
  list(1, "red")
)

# Create the 3D surface plot with the color scale
p <- plot_ly(x = ~Biom_values, y = ~nBac_values, z = ~Ub_matrix, type = "surface",
             surfacecolor = ~Ub_matrix, colorscale = colorscale) %>%
  layout(title = "3D Surface Plot of Ub Function",
         scene = list(
           xaxis = list(title = "Biom"),
           yaxis = list(title = "nBac"),
           zaxis = list(title = "Ub", colorscale = colorscale)
         ))

# Create the heatmap
h = ggplot(grid, aes(x = Biom, y = nBac, fill = Ub)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                       midpoint = median(grid$Ub), 
                       limit = c(min(grid$Ub), 
                                 max(grid$Ub)), space = "Lab", name="Ub") +
  theme_minimal() +
  labs(x = "Biom", y = "nBac", title = "Heatmap of Ub Function") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Generate a sequence of Biom values
Biom_values <- seq(0, 0.1821, length.out = 100)

# Choose a set of nBac values for which you want to plot Ub
nBac_levels <-  seq(0, nBac, by = 1000000)

# Create a data frame to store the results
results <- data.frame()

# Calculate Ub for each nBac level across all Biom values
for(nBac in nBac_levels) {
  Ub_values <- sapply(Biom_values, UbFunction, nBac = nBac)
  results <- rbind(results, data.frame(Biom = Biom_values, Ub = Ub_values, nBac = as.factor(nBac)))
}

# Plot the results
p2 = ggplot(results, aes(x = Biom, y = Ub, color = nBac, group = nBac)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Biom", y = "Ub", color = "nBac", title = "Ub vs Biom for different nBac levels") +
  theme(legend.position = "right")

p_bounds_1 = (p_aa + p_heme)
p_bounds_2 = (h / p2)

######## Starv transition

# (DOI: - 10.1038/s41598-019-55041-w)
# a growth rate of 0.20 (1/h) (stardard for C. diff)
# requires biomass production = 0.14 (mmol/gDW*h) 
# to make further comparison (https://doi.org/10.1371/journal.pcbi.1008137)
# a growth rate of 0.20 (1/h) (stardard for C. diff)
# requires glucose uptake = 30 (mmol/gDW*h)
# also, (https://www.nature.com/articles/s41598-019-55041-w)

gDW_CDmax = 0.18216 # (pg)
gDW_CDmin = 0.00777 # (pg)
gDW_CDmean = 0.05727 # (pg)

MWbio = 1
gDW_CDmean_g = gDW_CDmean*1e-12 # (g)

EX_B_starv = 0.14 # (mmol/gDW*h)
max_rate = 0.2 # (mmol/gDW*h)

EX_B_starv = EX_B_starv*MWbio # (g/gDW*h)
EX_B_starv = EX_B_starv*1e+12 # (pg/gDW*h)
EX_B_starv = EX_B_starv*gDW_CDmean_g

max_rate = max_rate*MWbio # (g/gDW*h)
max_rate = max_rate*1e+12 # (pg/gDW*h)
max_rate = max_rate*gDW_CDmean_g

logistic_rate <- function(Biom) {
  rate <- max_rate * ((Biom / gDW_CDmax)**3) / ((Biom / gDW_CDmax)**3 + max_rate)
  return(rate)
}

# Example usage:
Biom <- seq(0, gDW_CDmax, length.out = 1000)
rate <- logistic_rate(Biom)

p_starv = ggplot(data.frame(Biom, rate), aes(x=Biom, y=rate)) +
  geom_line() +
  geom_vline(aes(xintercept=gDW_CDmax), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=gDW_CDmean), color="green", linetype="dashed") +
  geom_vline(aes(xintercept=gDW_CDmin), color="red", linetype="dashed") +
  geom_hline(aes(yintercept=max_rate), color="red", linetype="dashed") +
  geom_hline(aes(yintercept=EX_B_starv), color="green", linetype="dashed") +
  labs(x = "Biomass (pg)", y = "Starvation Rate (pg/h)", title="Starvation Rate")

p_starv

####### Death4Treat

calculate_D4T_rate <- function(Death4Treat, DrugPlace) {
  
  MIC.conc <- 1e-06 # (g/mL)
  
  mm_drug <- 171.16 # (g/mol)
  mm_drug <- mm_drug*1e-03 # (g/mmol)
  
  dose <- 0.05 # (g/day*person)
  dose <- dose/24 # (g/h*person)
  dose <- dose/mm_drug # (mmol/h*person)
  
  MIC.conc <- (MIC.conc/mm_drug) # (mmol)
  
  Delta_drug <- round(MIC.conc/dose, 4)
  
  dose <- dose*Delta_drug # (mmol/h*mL)
  dose <- 8*dose*1e+09 # (pmol)
  
  MIC.conc <- MIC.conc*1e+09 # (pmol)
  
  nBac <- 45700000 # (cell)
  
  # the rate increases as the ratio of the 
  # drug concentration to its MIC increases. 
  # This could model the increased efficacy of the drug as its 
  # concentration surpasses the MIC.
  
  D4T_rate <- DrugPlace*nBac*Death4Treat*(DrugPlace/MIC.conc)
  return(D4T_rate)
}

Death4Treat_values <- runif(n = 24, min = 0, max = 1e-12)
DrugPlace <- seq(0, dose, length.out = 1000)

# Prepare data for ggplot
data_for_plot <- data.frame()
for(i in 1:length(Death4Treat_values)) {
  D4T_rates <- calculate_D4T_rate(Death4Treat_values[i], DrugPlace)
  data_for_plot <- rbind(data_for_plot, data.frame(DrugPlace, D4T_rate = D4T_rates, Death4Treat = Death4Treat_values[i]))
}

# Plotting with ggplot
D4T_plot = ggplot(
  data_for_plot, 
  aes(x = DrugPlace, 
      y = D4T_rate, 
      group = Death4Treat, 
      color = Death4Treat)) +
  geom_line() +
  labs(x = "DrugPlace", y = "Death4Treat Rate (cell*pmol/h)", title = "Death4Treat Rate") + 
  theme_minimal() +
  geom_vline(aes(xintercept=5842.487), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=48687.39), color="red", linetype="dashed") +
  scale_color_gradient(low = "blue", high = "red", name = "Death4Treat")

####### Duplication

# Define the function to calculate dup_rate
calculate_dup_rate <- function(Biom, nBac) {
  
  gDW_CDmax = 0.18216 # (pg)
  gDW_CDmin = 0.00777 # (pg)
  gDW_CDmean = 0.05727 # (pg)
  
  # The maximum number of cells in 1 mL can vary greatly depending on the type of cell. 
  # For example, a saturated bacterial culture can contain about 10^9 cells per mL
  # (http://book.bionumbers.org/what-is-the-concentration-of-bacterial-cells-in-a-saturated-culture/)
  nBacMax = 1e+09 # (cell/mL)
  
  # C. difficile is a slow-growing bacterium. 
  # the generation time of C. difficile can range from 27 minutes to 216 minutes, 
  # depending on the strain and conditions (C. difficile). 
  
  # if we consider rCDdup as the duplication rate per hour, 
  # a reasonable estimate might be between 0.28 (for a 216-minute generation time) 
  # and 2.22 (for a 27-minute generation time) divisions per hour.
  rCDdup = 2.22 # (cell/h)
  
  # normaliation the biomass of the cell between a minimum and maximum value.
  # This could represent the proportion of the cell's current size relative to its maximum size.
  d = (Biom - gDW_CDmin)/(gDW_CDmax - gDW_CDmin)
  # rapresenting the remaining capacity for cell growth or division
  # where 'nBacMax' is the maximum number of cells.
  cap = 1 - (nBac/nBacMax)
  
  # Introduce an exponential factor to the duplication rate calculation
  # The exponential factor increases the rate more rapidly as the biomass and cell number increase
  exp_factor = exp(d * cap)
  dup_rate = d * nBac * rCDdup * cap * exp_factor
  
  return(dup_rate)
}


# Create a sequence of Biom values
Biom = seq(0, gDW_CDmax, length.out = 1000) # (cell)

# Set nBac to 1
nBac <- 1

# Calculate dup_rate for each Biom value
dup_rates <- sapply(Biom, calculate_dup_rate, nBac)

dup_plot = ggplot(data.frame(Biom, dup_rates), aes(x=Biom, y=dup_rates)) +
  geom_line() +
  geom_vline(aes(xintercept=gDW_CDmax), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=gDW_CDmean), color="green", linetype="dashed") +
  geom_vline(aes(xintercept=gDW_CDmin), color="red", linetype="dashed") +
  labs(x = "Biomass (pg)", y = "Duplication Rate (cell/h)", title="Duplication Rate")

####### DeathBac

calculate_rate_death <- function(Biom, nBac) {

  gDW_CDmax = 0.18216 # (pg)
  gDW_CDmin = 0.00777 # (pg)
  gDW_CDmean = 0.05727 # (pg)
  
  # half-life of the bacterial population, which is the time 
  # it takes for half of the bacteria in the population to 
  # undergo the process (e.g., decay).
  half_life = 0.15 # (h)
  
  # This is a constant that determines the steepness of the logistic function.
  k = 100
  
  # representing the difference between the average biomass of a cell 
  # in the bacterial population and the mean biomass for a bacterial cell observed in nature.
  # measure of how much the average biomass of a cell in the population
  # deviates from what is typically observed in nature.
  diff = Biom - gDW_CDmean
  
  # This difference could have various implications depending on the specific biological context. 
  # For example, a positive value could indicate that the bacterial population
  # is well-fed or experiencing favorable conditions, allowing them to grow larger than average. 
  # On the other hand, a negative value could indicate that the bacterial population is under stress 
  # or experiencing unfavorable conditions, causing them to be smaller than average.
  
  # logistic function, such as population decay. 
  # Here, it might represent how the process rate changes with the biomass (Biom), 
  # relative to a mean value (gDW_CDmean). 
  l <- (1/(2 + exp(k * diff)))
  
  rate_death <- half_life * nBac * l
  
  return(rate_death)
}

nBac = nBac = 45700000 # (cell)
# Create a sequence of Biom values
Biom = seq(0, gDW_CDmax, length.out = 1000) # (cell)

rate_death <- calculate_rate_death(Biom, nBac)

death_rate_plot = ggplot(data.frame(Biom, rate_death), aes(x=Biom, y=rate_death)) +
  geom_line() +
  geom_vline(aes(xintercept=gDW_CDmax), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=gDW_CDmean), color="green", linetype="dashed") +
  geom_vline(aes(xintercept=gDW_CDmin), color="red", linetype="dashed") +
  labs(x = "Biomass (pg)", y = "DeathCD Rate (cell/h)", title="DeathCD Rate")

#######

transitions_plot = ((D4T_plot / detox_plot)) | 
  ((dup_plot / death_rate_plot) | (p_starv / (efflux_plot | inflam)))
