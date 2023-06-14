#### RL-EMC2!
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  ROI <- args[1]
} else {
  # manual run
  rm(list=ls())
  ROI <- 'Str'
}

## First prepare events dataframes for all subjects
library(EMC2)
library(emcAdapt)
source('./utility_funcs_fmri2.R')

## Overwrite to pass along prediction errors
log_likelihood_race_new <- function(p_vector, dadm, min_ll=log(1e-10))
  # Race model summed log likelihood
{

  pars <- EMC2:::get_pars(p_vector,dadm)
  if (is.null(attr(pars,"ok"))){
    ok <- !logical(dim(pars)[1])
  } else {
    ok <- attr(pars,"ok")
  }
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(attr(dadm,"model")()$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <- log(1-attr(dadm,"model")()$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- min_ll
  lds <- lds[attr(dadm,"expand")] # decompress
  if (n_acc>1) {
    winner <- dadm$winner[attr(dadm,"expand")]
    ll <- lds[winner]
    if (n_acc==2) {
      ll <- ll + lds[!winner]
    } else {
      ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
    }
    ll[is.na(ll)] <- min_ll
    out <- sum(pmax(min_ll,ll))
  } else {
    out <- sum(pmax(min_ll,lds))
  }
  if(!is.null(attr(pars, 'predictionErrors'))) {
    attr(out, 'predictionErrors') <- attr(pars, 'predictionErrors')
  }
  return(out)
}
assignInNamespace("log_likelihood_race", log_likelihood_race_new,ns="EMC2")

## Overwrite this function to pass prediction errors from behavioral to neural model
log_likelihood_joint_new <- function(proposals, dadms, component = NULL, useC) #pars, dadms, component = NULL)
{
  parPreFixs <- unique(gsub("[|].*", "", colnames(proposals)))
  i <- 0
  total_ll <- 0
  if(!is.null(component)) dadms <- dadms[component]
  for(dadm in dadms){
    if(is.data.frame(dadm)){
      i <- i + 1
      parPrefix <- parPreFixs[i]
      currentPars <- proposals[,grep(paste0(parPrefix, "|"), colnames(proposals), fixed = T)]
      colnames(currentPars) <- gsub(".*[|]", "", colnames(currentPars))

      if(i == 1) {
        this_ll <- EMC2:::calc_ll_manager(currentPars, dadm, useC, attr(dadm, "model")()$log_likelihood, hasPEs=TRUE)
        predictionErrors <- attr(this_ll, 'predictionErrors')
        row.names(predictionErrors) <- dadm[dadm$winner,'originalTrialNumber']
      } else {   # neural data: bypass calc_ll_manager
        this_ll <- EMC2:::calc_ll_manager(currentPars, dadm, useC, attr(dadm, "model")()$log_likelihood, predictionErrors=predictionErrors)
      }
      total_ll <- total_ll + this_ll
      #  OLD PRE-RCPP
      #   i <- i + 1
      #   parPrefix <- parPreFixs[i]
      #   currentPars <- pars[grep(paste0(parPrefix, "|"), names(pars), fixed = T)]
      #   names(currentPars) <- gsub(".*[|]", "", names(currentPars))
      #
      #   if(i == 1) {  # BEHAVIORAL DATA
      #     this_ll <- attr(dadm, "model")()$log_likelihood(currentPars, dadm)
      #     predictionErrors <- apply(attr(this_ll, 'predictionErrors'), 1, sum, na.rm=TRUE)
      #     predictionErrors <- data.frame(predictionErrors=predictionErrors, trial_nr=dadm[dadm$winner,'originalTrialNumber'])
      #   } else { # NEURAL DATA
      #     this_ll <- attr(dadm, "model")()$log_likelihood(currentPars, dadm, predictionErrors=predictionErrors)
      #   }
      #   total_ll <- total_ll + this_ll
    }
  }
  return(total_ll)
}
assignInNamespace("log_likelihood_joint", log_likelihood_joint_new, ns="EMC2")

## Events reading
read_events <- function(subject, run, stc=-1.38/2, include_events=c('response_left', 'response_right'),
                        path_format='/home/Public/trondheim/derivatives/event_files/sub-%s/ses-rlsat/func/sub-%s_ses-rlsat_task-rlsat_run-%d_events.tsv') {
  if(is.numeric(subject)) subject <- sprintf("%03d", subject)

  events_fn = sprintf(path_format, subject, subject, run)
  events = read.csv(events_fn, sep="\t")

  # apply slice timing correction onset difference
  events$onset <- events$onset+stc

  # overwrite duration, we're not modelling this
  events$duration = 0.001

  # subset
  events <- events[events$trial_type %in% include_events,]

  return(events)
}

read_all_events <- function(timeseries, include_events=c('cue_ACC', 'cue_SPD', 'stimulus', 'feedback', 'response_left', 'response_right')) {
  subject <- timeseries[1,'subjects']
  runs <- unique(timeseries$run)
  all_events <- vector(mode='list', length=length(runs))
  for(run in runs) {
    all_events[[run]] <- read_events(subject=subject, run=run, include_events=include_events)
    all_events[[run]]$run <- run
    all_events[[run]]$subjects <- subject
  }
  return(do.call(rbind, all_events))
}

make_fmri_dm_from_ts <- function(timeseries, run=1, include_events=c('cue_ACC', 'cue_SPD', 'stimulus', 'feedback', 'response_left', 'response_right')) {
  subject <- timeseries[1,'subjects']
  runs <- unique(timeseries$run)
  dms <- vector(mode='list', length=length(runs))
  for(run in runs) {
    events = read_events(subject=subject, run=run, include_events = include_events)
    events$modulation <- 1

    ## generate contrasts of left-right
    these_events <- events[events$trial_type %in% c('response_left', 'response_right'),]
    these_events[these_events$trial_type=='response_left','modulation'] <- 1
    these_events[these_events$trial_type=='response_right','modulation'] <- -1
    events[events$trial_type %in% c('response_left', 'response_right'),'trial_type'] <- 'response'
    these_events$trial_type <- 'response_leftminright'
    events <- rbind(events, these_events)
    events <- events[order(events$onset),]

    ## generate contrasts of acc-spd
    these_events <- events[events$trial_type %in% c('cue_SPD', 'cue_ACC'),]
    these_events[these_events$trial_type=='cue_ACC','modulation'] <- 1
    these_events[these_events$trial_type=='cue_SPD','modulation'] <- -1
    events[events$trial_type %in% c('cue_SPD', 'cue_ACC'),'trial_type'] <- 'cue'
    these_events$trial_type <- 'cue_ACCminSPD'
    events <- rbind(events, these_events)
    events <- events[order(events$onset),]

    dms[[run]] <- make_fmri_design_matrix(timeseries[timeseries$run==run, 'time'],
                                          events=events,
                                          hrf_model='glover + derivative', add_intercept=FALSE)
  }
  dm <- do.call(rbind, dms)
  dm$subjects <- subject
  #  dm$intercept <- 1
  return(dm)
}

# Load all timeseries
all_ROI_timeseries <- read.csv('/home/Public/trondheim/scripts/hierarchical_bayesian_GLM/data/signals.tsv', sep='\t')
colnames(all_ROI_timeseries)[2] <- 'subjects'
include_ROIs <- c(ROI) #c('Amg', 'Cl', 'GPe', 'GPi', 'PAG', 'PPN', 'RN', 'SN', 'STN', 'Str', 'Tha', 'VTA')
lr <- c('.l', '.r')
include_ROIs <- paste(rep(include_ROIs, each=2), lr, sep='')

# exclude subject 9 (fell asleep)
all_ROI_timeseries <- all_ROI_timeseries[all_ROI_timeseries$subjects!=9,]
all_ts_ls <- lapply(unique(all_ROI_timeseries$subject), function(x) all_ROI_timeseries[all_ROI_timeseries$subjects==x, c('subjects', 'run', 'time', include_ROIs)])
all_fmri_ts <- do.call(rbind, all_ts_ls)


####### We split the fMRI design matrix into two parts ########
# 1. Which events are modulated (and depend on the behavioral model)? %FOR MODULATION ONLY%
# We can't create this part of the DM yet, that needs to be done within the likelihood function
# Therefore, we get the event timings
modulated_events <- do.call(rbind, lapply(all_ts_ls, read_all_events, include_events='feedback'))
modulated_events$modulation <- NA  # !!! <- these values depend on the behavioral model, and will be filled within the fMRI-likelihood func !!!
modulated_events$trial_type <- 'feedback_pe'

# and some other info we need to store to speed up convolutions in likelihoods
hkernel = EMC2:::hrf_kernel('glover + derivative', tr=1.38, oversampling=5, fir_delays=NULL)

#regressor = EMC2:::sample_condition_(exp_condition, frame_times, oversampling=oversampling)
regressors = lapply(all_ts_ls, function(x, oversampling=5) {
  subject = x$subjects[1]
  runs = unique(x$run)
  all_regressors <- vector(mode='list', length=length(runs))
  for(run in runs) {
    frame_times = x[x$run==run, 'time']
    these_events = modulated_events[modulated_events$subjects==subject&modulated_events$run==run,]
    exp_condition = cbind(these_events[,c('onset', 'duration')], 1)
    regressor = EMC2:::sample_condition_(exp_condition, frame_times, oversampling=oversampling)
    regressor$frame_times=frame_times
    all_regressors[[run]] <- regressor
  }
  all_regressors
}
)
names(regressors) <- unique(all_fmri_ts$subjects)


# 2. For which events do we only estimate an intercept? Here, we can already create the DM
all_fmri_dms <- lapply(all_ts_ls, make_fmri_dm_from_ts,
                       include_events=c('cue_ACC', 'cue_SPD', #'stimulus',      # let's not estimate stimulus: strong covariance with responses
                                        'feedback', 'response_left', 'response_right'))

# Make fMRI design
design <- make_design_fmri2(events=modulated_events,    # events of modulation
                            data=do.call(rbind, all_ts_ls),
                            design_matrix=do.call(rbind, all_fmri_dms),  # fixed part of the design matrix
                            model=normal2,
                            whiten=FALSE,
                            regressors=regressors, hkernel=hkernel)
data <- all_fmri_ts

# Load behavioral data!
print(load('./datadesigntrondheim_rcpp.RData'))
data_behavior <- droplevels(data_behavior[as.character(data_behavior$subjects)!='9',])   # also exclude subject 9 here
data$subjects <- factor(data$subjects, levels=levels(data_behavior$subjects))
samplers <- make_samplers(list(data_behavior, data), list(design_behavior, design))


# sample
fn <- paste0('./samples/dataset-trondheim_model-rleam_roi-', ROI, '.RData')
if(file.exists(fn)) print(load(fn))
run_emc(samplers, fileName=fn, cores_per_chain = 7, cores_for_chains = 3, verbose  = T, useC = FALSE)
