#### RL-EMC2!
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  decisionModel <- args[1]
  learningModel <- args[2]
  task <- args[3]
} else {
  # manual run
  rm(list=ls())
  decisionModel <- 'ARD'    # ONE OF: RD, ARD
  learningModel <- 'delta'  # vkfbinary'    # ONE OF: delta, vkf, vkfbinary
  task <- 'rlsat'           # rlsat or revl
}

source('./models.R')
library(EMC2)
library(emcAdapt)


## settings
save_fn_samples <- paste0('./samples/dataset-trondheim_task-', task, '_model-', decisionModel, '-', learningModel, '.RData')

## load data
print(load(paste0('./data/dataset-trondheim_task-', task, '.RData')))

## functions
s_left <- function(d) {
  gsub("([[:alpha:]]+)-([0-9]+\\.[0-9]+)_([[:alpha:]]+)-([0-9]+\\.[0-9]+)", '\\1', as.character(d$S))
}
p_left <- function(d) {
  as.numeric(gsub("([[:alpha:]]+)-([0-9]+\\.[0-9]+)_([[:alpha:]]+)-([0-9]+\\.[0-9]+)", '\\2', as.character(d$S)))
}
s_right <- function(d) {
  gsub("([[:alpha:]]+)-([0-9]+\\.[0-9]+)_([[:alpha:]]+)-([0-9]+\\.[0-9]+)", '\\3', as.character(d$S))
}
p_right <- function(d) {
  as.numeric(gsub("([[:alpha:]]+)-([0-9]+\\.[0-9]+)_([[:alpha:]]+)-([0-9]+\\.[0-9]+)", '\\4', as.character(d$S)))
}
p_low <- function(d) {
  pmin(d$p_left, d$p_right)
}
p_high <- function(d) {
  pmax(d$p_left, d$p_right)
}
s_low <- function(d) {
  ifelse(d$p_low==d$p_left, d$s_left, d$s_right)
}
s_high <- function(d) {
  ifelse(d$p_high==d$p_left, d$s_left, d$s_right)
}
lS <- function(d) {
  ## latent STIMULUS that matches an accumulator
  factor(d[cbind(1:nrow(d), match(paste0('s_', d$lR), colnames(d)))])
}
lRS <- function(d) {
  ## latent STIMULUS that was chosen
  ifelse(d$R=='left', d$s_left, d$s_right)
}
## latent DIRECTION that matches an accumulator
correct_direction <- function(d) {
  factor(ifelse(d$p_right>d$p_left, 'right', 'left'), levels=levels(d$R))
}

# data$s_left <- s_left(data)
# data$s_right <- s_right(data)
# data$p_left <- p_left(data)
# data$p_right <- p_right(data)
# data$p_low <- p_low(data)
# data$p_high <- p_high(data)
# data$s_low <- s_low(data)
# data$s_high <- s_high(data)
# data$lRS <- lRS(data)
# data$correct_direction <- correct_direction(data)
head(data)

Ffunctions = list(s_left=s_left, s_right=s_right, p_left=p_left, p_right=p_right,
                  p_low=p_low, p_high=p_high, s_low=s_low, s_high=s_high,
                  lS=lS, lRS=lRS, correct_direction=correct_direction)

matchfun <- function(d){d$lR==d$correct_direction} # "correct" (higher expected reward) response

## get learning model set-up
learningModelParameters <- getLearningModel('delta')


## SAT modelling
if(task == 'rlsat') {
  Emat <- matrix(c(0,-1),nrow=2)
  dimnames(Emat) <- list(NULL,c("a-s"))
  Clist=list(cue=Emat)
} else {
  Emat <- NULL
  Clist <- NULL
}

## get decision model set-up
if(decisionModel == 'ARD') {
  decisionParameters <- c('v0', 'v', 'B', 'A', 't0', 's', 'wd', 'ws')
  decisionConstants <- c(A=log(0), v=log(0), s=log(1))
} else if(decisionModel == 'RD') {
  decisionParameters <- c('v0', 'v', 'B', 'A', 't0', 's', 'w')
  decisionConstants <- c(A=log(0), v=log(0), s=log(1))
}
decisionFlist <- lapply(decisionParameters, function(x) as.formula(paste0(x, '~1')))
names(decisionFlist) <- decisionParameters
if(task == 'rlsat') decisionFlist$B <- B~cue


if(learningModel == 'delta' & decisionModel == 'RD') {
  model = rdmRL
} else if(learningModel == 'delta' & decisionModel == 'ARD') {
  model = rdmRLARD
} else if(learningModel %in% c('vkf', 'vkfbinary') & decisionModel == 'RD') {
  model = rdmRLvkf
} else if(learningModel %in% c('vkf', 'vkfbinary') & decisionModel == 'ARD') {
  model = rdmRLARDvkf
}


# combine -----------------------------------------------------------------
adapt <- list(
  useC=TRUE,
  useSMsApproach=TRUE,
  stimulus=list(
    # 2AFC stimulus components whose value is to be adapted
    targets=as.character(unique(c(data$s_left, data$s_right))),

    # adaptive equation, creating "Q" values for each stimulus component
    init_par=learningModelParameters$init_par,
    adapt_par=learningModelParameters$adapt_par,
    feedback=c("reward"),
    adapt_fun_name=learningModel,
    adapt_fun=learningModelParameters$adapt_fun,

    # Output equation, uses Q values and parameters to update rates for stimulus components
    output_par_names=c("v0", "wd", "ws"),
    output_name = c("v"),
    output_fun=function(output_pars,Q)   # rate based on Q value #NB: Q should be a matrix of trials x accumulator (this/other)
      output_pars[,1] + output_pars[,2]*(Q[,1]-Q[,2]) + output_pars[,3]*(Q[,1]+Q[,2])
  ),
  fixed_pars=c("B","t0","A") # Non-adapted parameters
)

matchfun <- function(d){d$lR=="high"}  # "correct" (higher expected reward) response # do we need this?
design <- make_design(
  Ffactors=list(subjects=levels(data$subjects),S=levels(data$S), cue=levels(data$cue)),
  Ffunctions=Ffunctions,
  Fcovariates = c("reward"),
  Rlevels=levels(data$R),matchfun=matchfun,
  Flist=c(decisionFlist, learningModelParameters$Flist),
  Clist=Clist,
  constants = c(decisionConstants, learningModelParameters$constants),
  adapt=adapt, # Special adapt component
  model=model)

## SAVE DATA & DESIGN FOR JOINT MODELLING
# data_behavior <- data
# design_behavior <- design
# save(data_behavior, design_behavior, file='../data/datadesigntrondheim.RData')
##

#### for debugging
# dadm <- design_model(data,design)
# p_vector <- sampled_p_vector(design,doMap = FALSE)
# # # Average parameter table from Miletic for Exp 1
# p_vector["v0"] <- log(1.92)
# p_vector["B"] <- log(2.16)
# p_vector["t0"] <- log(.1)
# p_vector["wd"] <- log(3.09)
# p_vector["ws"] <- log(1)
# p_vector["alpha"] <- qnorm(0.12)
#
# # microbenchmark(
# # debug(update_pars); dadm2 <- dadm; dadm2$R <- NA
# log_likelihood_race(p_vector, dadm)
# #)


# Fit ----------------------------------------------------------------
if(!file.exists(save_fn_samples)) {
  samplers <- make_samplers(data, design, n_chains = 3, type='standard')
  save(samplers, file=save_fn_samples)
} else {
  load(save_fn_samples)
}


samplers <- run_emc(samplers, fileName=save_fn_samples, iter = 1000, cores_for_chains=3,
                    cores_per_chain=7, verbose = TRUE, verboseProgress = TRUE, max_trys=20) #cores_for_chains=3, cores_per_chain=7, max_trys=10)




#
# plot_chains(samplers, selection='mu')
