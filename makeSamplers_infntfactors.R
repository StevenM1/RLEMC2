#### RL-EMC2!
rm(list=ls())

## First prepare events dataframes for all subjects
library(reshape2)
library(EMC2)
library(emcAdapt)
library(pracma)
source('./utility_funcs_fmri.R')


fn <- paste0('./samples/dataset-trondheim_model-rleam_infinitefactors.RData')
# prevent overwrite
if(file.exists(fn)) {
  stop('Samplers file already exist!')
}


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
    if('response_left' %in% include_events) {
      these_events <- events[events$trial_type %in% c('response_left', 'response_right'),]
      these_events[these_events$trial_type=='response_left','modulation'] <- 1
      these_events[these_events$trial_type=='response_right','modulation'] <- -1
      events[events$trial_type %in% c('response_left', 'response_right'),'trial_type'] <- 'response'
      these_events$trial_type <- 'response_leftminright'
      events <- rbind(events, these_events)
      events <- events[order(events$onset),]
    }

    ## generate contrasts of acc-spd
    if('cue_SPD' %in% include_events) {
      these_events <- events[events$trial_type %in% c('cue_SPD', 'cue_ACC'),]
      these_events[these_events$trial_type=='cue_ACC','modulation'] <- 1
      these_events[these_events$trial_type=='cue_SPD','modulation'] <- -1
      events[events$trial_type %in% c('cue_SPD', 'cue_ACC'),'trial_type'] <- 'cue'
      these_events$trial_type <- 'cue_ACCminSPD'
      events <- rbind(events, these_events)
      events <- events[order(events$onset),]
    }

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
include_ROIs <- c('Amg', 'Cl', 'GPe', 'GPi', 'PAG', 'PPN', 'RN', 'SN', 'STN', 'Str', 'Tha', 'VTA')
lr <- c('.l', '.r')
include_ROIs <- paste(rep(include_ROIs, each=2), lr, sep='')

# collapse across left/right
ROI_timeseries_long = melt(all_ROI_timeseries[,c('subjects', 'run', 'time', include_ROIs)], id.vars=c('subjects', 'run', 'time'))
ROI_timeseries_long$ROI <- sapply(as.character(ROI_timeseries_long$variable), function(x) strsplit(x, '.', fixed=TRUE)[[1]][1])
ROI_timeseries_long$hemisphere <- sapply(as.character(ROI_timeseries_long$variable), function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
ROI_timeseries_aggregated <- aggregate(value~time*run*subjects*ROI, ROI_timeseries_long, mean)
ROI_timeseries_aggregated <- ROI_timeseries_aggregated[,c('subjects', 'run', 'time', 'ROI', 'value')]
ROI_timeseries <- reshape(ROI_timeseries_aggregated, direction='wide', idvar=c('subjects', 'run', 'time'), timevar = 'ROI')
colnames(ROI_timeseries) <- gsub('value.', '', colnames(ROI_timeseries), fixed=TRUE)
ROI_timeseries = ROI_timeseries[ROI_timeseries$subjects != 9,]  # exclude subject 9 (fell asleep)

ROI_timeseries_list <- lapply(unique(ROI_timeseries$subjects), function(x) ROI_timeseries[ROI_timeseries$subjects==x,])


####### We split the fMRI design matrix into two parts ########
# 1. Which events are modulated (and depend on the behavioral model)? %FOR MODULATION ONLY%
# We can't create this part of the DM yet, that needs to be done within the likelihood function
# Therefore, we get the event timings
modulated_events <- do.call(rbind, lapply(ROI_timeseries_list, read_all_events, include_events='feedback'))
modulated_events$modulation <- NA  # !!! <- these values depend on the behavioral model, and will be filled within the fMRI-likelihood func !!!
modulated_events$trial_type <- 'feedback_pe'

# and some other info we need to store to speed up convolutions in likelihoods
hkernel = EMC2:::hrf_kernel('glover + derivative', tr=1.38, oversampling=5, fir_delays=NULL)

#regressor = EMC2:::sample_condition_(exp_condition, frame_times, oversampling=oversampling)
regressors = lapply(ROI_timeseries_list, function(x, oversampling=5) {
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
names(regressors) <- unique(ROI_timeseries$subjects)


# 2. For which events do we only estimate an intercept? Here, we can already create the DM
all_fmri_dms <- lapply(ROI_timeseries_list, make_fmri_dm_from_ts,
                       include_events=c('cue_ACC', 'cue_SPD', #'stimulus',      # let's not estimate stimulus: strong covariance with responses
                                        'feedback'))

# Make fMRI design
design <- make_design_fmri2(events=modulated_events,    # events of modulation
                            data=do.call(rbind, ROI_timeseries_list),
                            design_matrix=do.call(rbind, all_fmri_dms),  # fixed part of the design matrix
                            model=normal2,
                            whiten=FALSE,
                            regressors=regressors, hkernel=hkernel)
data <- ROI_timeseries


# Load behavioral data, this was pre-created in fit_RLEAMs_trondheim.R
print(load('./data/datadesigntrondheim.RData'))
data_behavior <- droplevels(data_behavior[as.character(data_behavior$subjects)!='9',])   # also exclude subject 9 here
data$subjects <- factor(data$subjects, levels=levels(data_behavior$subjects))

## bit ugly: make standars samplers first to be able to extract the parameter names, and then re-create samplers with infnt-factor approach
samplers <- make_samplers(list(data_behavior, data), list(design_behavior, design))
pnames <- samplers[[1]]$par_names
## we only model the full var-covar matrix of: t0 wd ws v0 alpha B B_{cuea-s} feedback_pe cue_ACCminSPD
## so we should have 12rois*2 = 24 + 7 = 31 non-nuisance pars
nuisance_pars <- grepl('_derivative', pnames) | grepl('_feedback$', pnames) | grepl('_cue$', pnames) | grepl('_sd$', pnames)
samplers <- make_samplers(list(data_behavior, data), list(design_behavior, design), nuisance=which(nuisance_pars), type = "infnt_factor")

save(samplers, file=fn)
