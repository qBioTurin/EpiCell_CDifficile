
# UnifiedGreatMod was used to examine the responses of luminal epithelial cells to Clostridium difficile infection

We combined the ODE-based dynamical model with the Genome-Scale Metabolic Model (GSMM) into a unified multiscale hybrid model. _C. difficile_ strain CD196 GSMM was originally coded in the AGORA/MATLAB (1.03 version) format (Magnusdottir et al.,2016) and then elaborated using the R functions collected in [epimod_FBAfunctions](https://github.com/qBioTurin/epimod_FBAfunctions).

## The CDI infection model

The CDI infection Petri Net model is composed by places (graphically represented by circles) corresponding to epithelial cells, bacterial biomass, metabolites and tissue state (i.e. damage to the colonic mucosa) and by transitions (graphically represented by rectangles) corresponding to the interactions among the entities, cellular death, intake or efflux of metabolites, toxin action, intestinal inflammation and drug activity.

## Metabolic model

1. The model was previously published in: Please cite: Magnusdottir, S., Heinken, A., Kutt, L., Ravcheev, D.A., Bauer, E., Noronha, A., Greenhalgh, K., Jager, C., Baginska, J., Wilmes, P., Fleming, R.M.T., Thiele, I., "Generation of genome-scale metabolic reconstructions for 773 members of the human gut microbiota", Nature Biotechnology, 35(1):81-89 (2017). The original model (Clostridium\_difficile\_CD196.mat) was parsed and converted to *.RData.

2. The format conversion is accompanied by the text file (BIGGdata\_CD196\_react.tsv) model's reactions infromation compliant with BiGG database (http://bigg.ucsd.edu/).

## Changes from original model

The original model, for example, include some unbalanced pseudo reactions. These pre-defined boundary reactions fulfill a function for modeling by adding to or removing metabolites from the model system but are not based on real biology. In particular, sink reactions are of great use for those compounds produced or consumed by non-metabolic cellular processes. The following changes add a new reaction. We edit the original model adding a previous lacking sink reaction. After complement the following issue, a new activity of such cellular processes was assessed using FBA.

1. the sink reaction for heme **sink\_pheme(c)** is the boundary reaction for protoheme, an alternative route through non-metabolic processes from $sink\_pheme(c) : [c] : pheme <==> \emptyset$. It is allowed in the model to go in reverse direction. It carries artificially high fluxes (lower bound = -1000; upper bound = 1000), so flux of the final reaction was constrained.

2. **sink\_pheme(c)** do not required new cytosolic  metabolite **pheme\_c** (already originally present)

3. As boundary reaction **sink\_pheme(c)** were associated with a blank gene annotation.

# Model analysis and sensitivity routine

The model's parameterisation is primarily based on values extracted from the literature, except for *Death4Treat*, *Detox*, and *IECsDeath*. These parameters are explored using Partial Rank Correlation Coefficients (PRCC) to observe changes in the overall model dynamics. Then the analysis of the model is performed.

# Github structure

In order to obtain the dynamics showed in Fig.2 (B and C) and Fig.S4 of the [main paper](), the R script *Main.R* must be used. All the files of the model are saved in the *Net* folder, while all the R functions to generate the plot or processing the data and files are saved in the *code* folder.
For the computational cost and timing valuation, the  R script *TimingEvaluation.R* must be used.
