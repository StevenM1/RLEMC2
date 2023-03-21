#### RL-EMC2!
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  decisionModel <- args[1]
  learningModel <- args[2]
  experimentN <- args[3]
} else {
  # manual run
  rm(list=ls())
  decisionModel <- 'ARD'     # ONE OF: RD, ARD
  learningModel <- 'delta'   # ONE OF: delta, vkf, vkfbinary
  experimentN <- 'exp3'      # exp1 = 'standard' task, exp2 = reversal learning, exp3 = SAT
}

source('./models.R')         # these are a set of functions that parameterise ARD and learning models
library(EMC2)
library(emcAdapt)            # https://github.com/StevenM1/emcAdapt/


## settings
save_fn_samples <- paste0('./samples/data-', experimentN, '_model-', decisionModel, '-', learningModel, '.RData')

## load data
print(load(paste0('./data/data_', experimentN, '_newformat.RData')))

### Set of functions to expand column 'S' into s_left, p_left, s_right, p_right, p_low, p_high, s_low, s_high, lS, lRS, correct_direction
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


## get model set-up
learningModelParameters <- getLearningModel(learningModel)
decisionModelParameters <- getDecisionModel(decisionModel)


# For experiment 3, add SAT manipulation. For now, just assume B varies with cue, and model the difference ACC-SPD (a-s)
if(experimentN=='exp3') {
  Emat <- matrix(c(0,-1),nrow=2)
  dimnames(Emat) <- list(NULL,c("a-s"))
  decisionModelParameters$Flist$B <- B~cue
  Clist=list(cue=Emat)
  Ffactors = list(subjects=levels(data$subjects),S=levels(data$S), cue=levels(data$cue))
} else {
  Clist=NULL
  Ffactors = list(subjects=levels(data$subjects),S=levels(data$S))
}

# Load correct model
if(learningModel == 'delta' & decisionModel == 'RD') {
  model = rdmRL         # UNTESTED! didn't fit this yet within EMC2
  warning('UNTESTED! didnt fit this yet within EMC2')
} else if(learningModel == 'delta' & decisionModel == 'ARD') {
  model = rdmRLARD
} else if(learningModel %in% c('vkf', 'vkfbinary') & decisionModel == 'RD') {
  model = rdmRLvkf      # UNTESTED! didn't fit this yet within EMC2
  warning('UNTESTED! didnt fit this yet within EMC2')
} else if(learningModel %in% c('vkf', 'vkfbinary') & decisionModel == 'ARD') {
  model = rdmRLARDvkf   # UNTESTED! didn't fit this yet within EMC2
  warning('UNTESTED! didnt fit this yet within EMC2')
}


# combine learning model and decision model parameters -----------------------------------------------------------------
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


matchfun <- function(d){d$lR=="high"}  # "correct" (higher expected reward) response # do we need this? not sure here
design <- make_design(
  Ffactors=Ffactors,
  Ffunctions=Ffunctions,
  Fcovariates = c("reward"),
  Rlevels=levels(data$R),
  matchfun=matchfun,
  Flist=c(decisionModelParameters$Flist, learningModelParameters$Flist),
  constants = c(decisionModelParameters$constants, learningModelParameters$constants),
  adapt=adapt,   # Special adapt component
  model=model)



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
                    cores_per_chain=7, verbose = TRUE, verboseProgress = TRUE, max_trys=20)


