
setwd("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile")
source("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Rfunction/SensitivityAnalysis.R")

seed = 1234
combinations = 15000

#################################################################################################################
#################################################### rnorm ######################################################
#################################################################################################################

load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196/CD196.RData")

SensitivityDataFBAStandard = SensitivityAnalysis(seed = seed, combinations = combinations, algo = "fba",
                                                 model = "CD196", verbose = T,
                                                 sink = F, distr = "rnorm", mean = 100, sd = pi*25,
                                                 minBound = 1e-18, maxBound = 1000)
                                                 
load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196/CD196.RData")

SensitivityDataMTFStandard = SensitivityAnalysis(seed = seed, combinations = combinations, algo = "mtf",
                                                 model = "CD196", verbose = T,
                                                 sink = F, distr = "rnorm", mean = 100, sd = pi*25,
                                                 minBound = 1e-18, maxBound = 1000)

load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")

SensitivityDataFBASink = SensitivityAnalysis(seed = seed, combinations = combinations, algo = "fba",
                                             model = "CD196HemeSink", verbose = T,
                                             sink = T, distr = "rnorm", mean = 100, sd = pi*25,
                                             minBound = 1e-18, maxBound = 1000)
                                             
load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")

SensitivityDataMTFSink = SensitivityAnalysis(seed = seed, combinations = combinations, algo = "mtf",
                                             model = "CD196HemeSink", verbose = T,
                                             sink = T, distr = "rnorm", mean = 100, sd = pi*25,
                                             minBound = 1e-18, maxBound = 1000)                                           

#################################################################################################################
#################################################### runif ######################################################
#################################################################################################################

load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196/CD196.RData")

SensitivityDataFBAStandard = SensitivityAnalysis(seed = seed, combinations = combinations, algo = "fba",
                                                 model = "CD196", verbose = T,
                                                 sink = F, distr = "runif", mean = 100, sd = pi*25,
                                                 minBound = 1e-18, maxBound = 1000)
                                                 
load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196/CD196.RData")

SensitivityDataMTFStandard = SensitivityAnalysis(seed = seed, combinations = combinations, algo = "mtf",
                                                 model = "CD196", verbose = T,
                                                 sink = F, distr = "runif", mean = 100, sd = pi*25,
                                                 minBound = 1e-18, maxBound = 1000)

load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")

SensitivityDataFBASink = SensitivityAnalysis(seed = seed, combinations = combinations, algo = "fba",
                                             model = "CD196HemeSink", verbose = T,
                                             sink = T, distr = "runif", mean = 100, sd = pi*25,
                                             minBound = 1e-18, maxBound = 1000)
                                             
load("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.RData")

SensitivityDataMTFSink = SensitivityAnalysis(seed = seed, combinations = combinations, algo = "mtf",
                                             model = "CD196HemeSink", verbose = T,
                                             sink = T, distr = "runif", mean = 100, sd = pi*25,
                                             minBound = 1e-18, maxBound = 1000)                                           
