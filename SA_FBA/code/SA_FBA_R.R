
# Global Sensitivity Analysis using R
#
# We will investigate the sensitivty of the reaction objective ("biomass205") 
# flux solution by varying the D parameters within the bounds (-10.0, 0.0) 
# defined by the problem shown below.

# libraries required
library(fbar)
library(parallel)
library(dplyr)
library(tidyr)
library(sensobol)
library(tibble)

# setting working directory
setwd("SA_FBA")
wd_sa = "."

# setting results directory
result_dir_sa = paste0(wd_sa, "/res")

# functions and methods
source(paste0(wd_sa, "/code/functions/ParallelFBA_sens.R"))

# FBA model tag
model_file = paste0(wd_sa, "/input/FBAmodel.RData")

# FBA model reaction dataframe file
all_react_file = paste0(wd_sa, "/input/all_react.rds")

# sample matrix file
param_values_file = paste0(result_dir_sa, "/param_values.rds")

# full FBA model parametrization
bounds_set_file = paste0(result_dir_sa, "/bounds_set.rds")

# flux distribution solution template bounds
fbasol_t_file = paste0(result_dir_sa, "/fbasol_template.rds")
# flux distribution solutions samples configuration
fbasol_file = paste0(result_dir_sa, "/fbasol.rds")

# output 1
Y_file = paste0(result_dir_sa, "/Y.rds")
# output 2
Y2_reaction = "EX_h2o_e"
Y2_file = paste0(result_dir_sa, "/Y2.rds")

indices_file = paste0(result_dir_sa, "/df_ind.rds")

# FBA model reactions annotation
all_react = readRDS(all_react_file)
br = all_react[[2]]
all_react = all_react[[1]]

# Define settings
N = 2^13
chunks <<- detectCores()*500
D = length(br$react.id)
params = paste("P", 1:D, sep = ".")
b1 = -10.0
b2 = 0.0

# Create sample matrix using Sobol' Quasi Random Numbers.
# The parameter values are sampled using a uniform distribution (its quantile function)
mat = as.data.frame(apply(sobol_matrices(N = N, params = params, order = "first"), 2, function(x) qunif(x, b1, b2)))
mat_list = split(mat, rep(1:(nrow(mat) %/% chunks + 1), each = chunks, length.out = nrow(mat)))

list_files = list.files(result_dir_sa)

param_values_files = list_files[grepl("param_values", list_files)]
bounds_set_files = list_files[grepl("bounds_set", list_files)]
fbasol_files = list_files[grepl("fbasol", list_files)]

for(i in 1:length(mat_list)) {
  
  print(paste0('chunk: ', i, ' out of ', length(mat_list)))
  saveRDS(mat_list[[i]], param_values_file)
  
  ## FBA runcode
  ParallelFBA_sens(model_file = model_file,
                   param_values_file = param_values_file,
                   model_metadata_file = all_react_file,
                   result_dir = result_dir_sa,
                   reactIndx4secondSol = all_react$ReactionPos[which(
                     all_react$React_ID == Y2_reaction)],
                   cores = detectCores(),
                   wd = wd_sa)
  
  file.rename(bounds_set_file, paste0(result_dir_sa, "/bounds_set", "_", i, ".rds"))
  file.rename(fbasol_file, paste0(result_dir_sa, "/fbasol", "_", i, ".rds"))
  
  Y = readRDS(Y_file)
  Y2 = readRDS(Y2_file)
  
  if(i == 1) {
    Y_t = Y
    Y2_t = Y2
  } else {
    Y_t = c(Y_t, Y)
    Y2_t = c(Y2_t, Y2)
  }
  
  gc()
  
}

# saving data
saveRDS(mat, param_values_file)
saveRDS(Y_t, Y_file)
saveRDS(Y2_t, Y2_file)

param_values = readRDS(param_values_file)
Y = readRDS(Y_file)
Y2 = readRDS(Y2_file)

gc()

# Bootstraping is a statistical technique used to estimate the sampling distribution 
# of an estimator by resampling the original data with replacement. 
# It is used to construct bootstrap confidence intervals on sensitivity indices 
# computed by polynomial chaos expansion

# When parallel = "multicore" is used, each worker process inherits
# the environment of the current session, including the workspace and the loaded namespaces and
# attached packages (but not the random number seed: see below).

# Compute and bootstrap Sobol' indices
sensitivity_Y = sobol_indices(Y = Y, N = N, params = params, type = "norm",
                              boot = TRUE, R = 1000, parallel = "multicore")

sensitivity_Y2 = sobol_indices(Y = Y2, N = N, params = params, 
                               boot = TRUE, R = 1000, parallel = "multicore")

sens.ind = rbind(
  "Y" = rbind(first = data.frame(parameter = params,
                                 value = dplyr::filter(
                                   sensitivity_Y$results, sensitivity == "Si")$original,
                                 low_ci = dplyr::filter(
                                   sensitivity_Y$results, sensitivity == "Si")$low.ci,
                                 high_ci = dplyr::filter(
                                   sensitivity_Y$results, sensitivity == "Si")$high.ci),
              total = data.frame(parameter = params,
                                 value = dplyr::filter(
                                   sensitivity_Y$results, sensitivity == "Ti")$original,
                                 low_ci = dplyr::filter(
                                   sensitivity_Y$results, sensitivity == "Ti")$low.ci,
                                 high_ci = dplyr::filter(
                                   sensitivity_Y$results, sensitivity == "Ti")$high.ci)), 
  "Y2" = rbind(first = data.frame(parameter = params,
                                  value = dplyr::filter(
                                    sensitivity_Y2$results, sensitivity == "Si")$original,
                                  low_ci = dplyr::filter(
                                    sensitivity_Y2$results, sensitivity == "Si")$low.ci,
                                  high_ci = dplyr::filter(
                                    sensitivity_Y2$results, sensitivity == "Si")$high.ci),
               total = data.frame(parameter = params,
                                  value = dplyr::filter(
                                    sensitivity_Y2$results, sensitivity == "Ti")$original,
                                  low_ci = dplyr::filter(
                                    sensitivity_Y2$results, sensitivity == "Ti")$low.ci,
                                  high_ci = dplyr::filter(
                                    sensitivity_Y2$results, sensitivity == "Ti")$high.ci)))

df_ind = sens.ind %>%
  rownames_to_column("sol") %>%
  separate(sol, into = c("sol", "order", "param"), sep="\\.")

saveRDS(df_ind, indices_file)

gc()