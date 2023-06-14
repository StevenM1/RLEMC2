#### RL-EMC2!
rm(list=ls())
library(reshape2)
library(EMC2)
library(emcAdapt)
library(pracma)
source('./emc2_functions_overwrite.R')

## Overwrite to pass along prediction errors
assignInNamespace("log_likelihood_race", log_likelihood_race_new,ns="EMC2")

## Overwrite this function to pass prediction errors from behavioral to neural model
assignInNamespace("log_likelihood_joint", log_likelihood_joint_new, ns="EMC2")

## Overwrite calc_ll_manager to allow for optional arguments
assignInNamespace("calc_ll_manager", calc_ll_manager_new, ns="EMC2")


## load samplers & sample
fn <- paste0('./samples/dataset-trondheim_model-rleam_infinitefactors.RData')
if(file.exists(fn)) print(load(fn)) else stop('Samplers file not found')

run_emc(samplers, fileName=fn, cores_per_chain = 7, cores_for_chains = 3, verbose  = TRUE, verboseProgress = TRUE)
