## RLSAT
rm(list=ls())
library(moments)
# Load data ----------------
dataPath <- '/home/Public/trondheim/derivatives/behavior/sub*/ses-*/func/*task*learning*.tsv'
allFns <- Sys.glob(dataPath)

dat <- NULL

for(fn in allFns) {
  #  path = file.path(dataPath, fn)
  datThisSub <- read.csv(fn, sep='\t')
  #  datThisSub <- datThisSub[datThisSub$event_type == 'response' & !is.na(datThisSub$choice_outcome),]
  datThisSub$pp <- strsplit(strsplit(strsplit(fn, '_')[[1]][1], '-')[[1]][2], '/')[[1]][1]
  datThisSub$task <- strsplit(strsplit(fn, '_')[[1]][2], '-')[[1]][2]
  datThisSub <- datThisSub[,c('pp', 'task', 'trial_nr', 'response', 'rt', 'stimulus_symbol_left', 'stimulus_symbol_right', 'p_win_left', 'p_win_right', 'block_nr', 'correct_response_direction', 'cue', 'choice_direction', 'choice_outcome', 'ease', 'accuracy')]
  dat <- rbind(dat, datThisSub)
}
dat$pp <- factor(as.numeric(as.character(dat$pp)))
dat <- droplevels(dat)
dat <- dat[order(dat$pp),]
dat$choiceIsHighP <- NA
dat$choiceIsHighP[dat$accuracy=='True'] <- 1
dat$choiceIsHighP[dat$accuracy=='False'] <- 0

# some factors
dat$choice_direction <- factor(dat$choice_direction, levels=c('left', 'right'), labels=c(0,1))
dat$correct_response_direction <- factor(dat$correct_response_direction, levels=c('left', 'right'), labels=c(0,1))

# Remove NAs
dat <- dat[!is.na(dat$rt),]

# select trials that should be excluded *AFTER* the updating has occured (ie, do take into account for the RL part, but do not include in likelihood)
# Specifically, we select the extremely fast or slow responses to be removed
dat$removeAfterUpdating <- FALSE
dat$removeAfterUpdating[!is.na(dat$rt) & (dat$rt<.15 | dat$rt>1.75)] <- TRUE

# For NA trials: assume stimulus value is correct, can be inferred from previous trials
# For too low/too late: assume PE is correct, can be calculated, but do not use stimulus value in GLM
# (because not really a normal choice). Hence, calcalate value&PEs for these trials, but do not include in likelihood

# Experiment checks -------------------------------------------------------
# Did the participants get the correct "choice outcomes"?
dat$choiceP <- NA
dat[dat$choice_direction == 0, 'choiceP'] <- dat[dat$choice_direction == 0, 'p_win_left']
dat[dat$choice_direction == 1, 'choiceP'] <- dat[dat$choice_direction == 1, 'p_win_right']
aggregate(choice_outcome ~ choiceP, dat, mean)

# What was the average probability of reward for each stimulus symbol? Should approximate 0.5 after all subs
dat$choiceSymbol <- NA
dat[dat$choice_direction == 0, 'choiceSymbol'] <- dat[dat$choice_direction == 0, 'stimulus_symbol_left']
dat[dat$choice_direction == 1, 'choiceSymbol'] <- dat[dat$choice_direction == 1, 'stimulus_symbol_right']
aggregate(choice_outcome ~ choiceSymbol, dat, mean)


# make stimulus_set column
# dat$stimulus_set = NA
# dat$ease = factor(abs(dat$p_win_left-dat$p_win_right), levels=c(0.6, 0.4, 0.2))
# for(task in unique(dat$task)) {
#   thisTaskIdx <- dat$task==task
#
#   dat$stimulus_set[thisTaskIdx] = as.integer(dat$ease[thisTaskIdx]) + (dat$block_nr[thisTaskIdx]-1) * max(as.integer(dat$ease[thisTaskIdx]))
# }
# these are ordered by block, difficulty; so: block1,diff1 = 1; block1,diff2=2, block1,diff3=3, block2,diff1=4, etc
# dat$reward <- dat$outcome <- dat$choice_outcome

# no full participant exclusions so far
dat$excl <- FALSE

# exclude participant sub-009's second half of first block and entire second block; seems to have fallen asleep. The rest does look OK.
dat <- dat[!((dat$pp==9 & dat$task=='rlsat') & (dat$trial_nr > 57 & dat$trial_nr < 228)),]

# get reversal points (relevant for revl only)
dat$reversal <- 0
######## find reversal point per subject
idx1 = dat$task=='rbrevl'
for(pp in unique(dat[idx1,]$pp)) {
  idx2 <- idx1 & dat$pp == pp

  for(block_nr in unique(dat[idx2,'block_nr'])) {
    idx3 = idx2 & dat$block_nr==block_nr

    orders = aggregate(p_win_left~stimulus_symbol_left, dat[idx3,], unique)
    for(stim_sym_left in orders$stimulus_symbol_left) {
      dat[idx3&dat$stimulus_symbol_left==stim_sym_left & dat$p_win_left==orders[orders$stimulus_symbol_left==stim_sym_left,'p_win_left'][1],'reversal'] = 0
      dat[idx3&dat$stimulus_symbol_left==stim_sym_left & dat$p_win_left==orders[orders$stimulus_symbol_left==stim_sym_left,'p_win_left'][2],'reversal'] = 1
    }
    for(ease in unique(dat[idx3, 'ease'])) {
      idx4 <- idx3 & dat$ease == ease
      dat[idx4, 'switch'] = c(0, diff(dat[idx4,'reversal']))
    }
  }
}

# find switch point, trialN relative to reversal
# What's the switch trial?
#dat$switch <- c(0, diff(dat$reversal))
dat$switch[dat$trial_nr == 0] <- 0
dat$switch[dat$switch == -1] <- 0 ## second reversal, ignore this here
max(dat[dat$switch==1&dat$task=='rbrevl','trial_nr'])

# add bins, trialN relative to (first) reversal
dat$trialBin <- NA
dat$hasReversed <- 0
dat$trialNreversal <- 0
for(pp in unique(dat$pp[dat$task=='rbrevl'])) {
  for(block in unique(dat$block_nr[dat$pp==pp&dat$task=='rbrevl'])) {
    for(ease in unique(dat$ease[dat$block_nr==block&dat$pp==pp&dat$task=='rbrevl'])) {
      idx <- dat$pp == pp & dat$block == block & dat$task=='rbrevl' & dat$ease==ease
      dat[idx, 'trialBin'] <- cut(dat[idx, 'trial_nr'], 10, labels=FALSE)
      dat[idx, 'hasReversed'] <- cumsum(dat[idx, 'switch'])

      idxPreReversal <- idx & dat$hasReversed == 0
      idxPostReversal <- idx & dat$hasReversed > 0 & dat$switch == 0
      dat[idxPreReversal,'trialNreversal'] <- seq(-sum(idxPreReversal), -1, 1)
      dat[idxPostReversal,'trialNreversal'] <- seq(1, sum(idxPostReversal), 1)
    }
  }
}
dat[dat$task=='rlsat','reversal'] <- 0
dat[dat$task=='rlsat','switch'] <- 0


##
data <- dat[!dat$excl,]
data$pp <- as.factor(dat$pp)
data <- data[!is.na(data$rt),]

data$s_left <- data$stimulus_symbol_left
data$s_right <- data$stimulus_symbol_right

data$p_left <- data$p_win_left
data$p_right <- data$p_win_right

data$choice_symbol <- NA
data$choice_symbol[data$choice_direction == 0] <- as.character(data$stimulus_symbol_left[data$choice_direction == 0])
data$choice_symbol[data$choice_direction == 1] <- as.character(data$stimulus_symbol_right[data$choice_direction == 1])

## ugly hack to get all relevant information
data$S <- factor(paste(paste(data$s_left, data$p_left, sep='-'), paste(data$s_right, data$p_right, sep='-'), sep='_'))

##
data <- data[,c("pp","trial_nr","S", 'cue', "choice_outcome", "choiceIsHighP", "rt", "choice_direction", "choice_symbol",
                's_left', 's_right', 'p_left', 'p_right', 'trialNreversal')]
if(max(data$reward)>1) { data$reward <- data$reward/100 }
names(data) <- c("subjects","originalTrialNumber","S", 'cue', "reward","Racc","rt", "R", "Rs", 's_left', 's_right', 'p_left', 'p_right', 'trialNreversal')
data$R <- factor(ifelse(data$R==0, 'left', 'right'), levels=c('left', 'right'))
data$cue <- factor(data$cue, levels=c('ACC', 'SPD'))
allData <- data

## Subset on experiment, save
for(task in c('rlsat', 'revl')) {
  if(task == 'rlsat') {
    data <- allData[!is.na(allData$cue),]
  } else if(task == 'revl') {
    data <- allData[is.na(allData$cue),]
  }
  data <- EMC2:::add_trials(data)
  data <- droplevels(data)
  levels(data$cue) <- c('ACC', 'SPD')
  save(data, file=paste0('../data/dataset-trondheim_task-', task, '.RData'))
}



# Experiment 3 (SAT) ----------------------------------------
# experimentN <- 'exp3'
# print(load(paste0('~/Projects/EMC_latest/RLEMC/data/data_', experimentN, '.RData')))
#
# data <- dat[!dat$excl,] # Excluded participants
# data$pp <- factor(as.character(data$pp))
# data <- data[!is.na(data$rt),]
#
# data$s_left <- data$stimulus_symbol_left
# data$s_right <- data$stimulus_symbol_right
#
# data$p_left <- data$p_win_left
# data$p_right <- data$p_win_right

#data$choiceIsHighP <- data$choiceIsHighPpreRev
# data$choice_symbol <- NA
# data$choice_symbol[data$choice_direction == 0] <- as.character(data$stimulus_symbol_left[data$choice_direction == 0])
# data$choice_symbol[data$choice_direction == 1] <- as.character(data$stimulus_symbol_right[data$choice_direction == 1])


## ugly hack to get all relevant information
# data$S <- factor(paste(paste(data$s_left, data$p_left, sep='-'), paste(data$s_right, data$p_right, sep='-'), sep='_'))
#
# data <- data[,c("pp","trial_nr","S", 'cue', "outcome","choiceIsHighP","rt", "choice_direction", "choice_symbol",
#                 's_left', 's_right', 'p_left', 'p_right')]
# if(max(data$reward)>1) { data$reward <- data$reward/100 }
# names(data) <- c("subjects","trials","S", 'cue', "reward","Racc","rt", "R", "Rs", 's_left', 's_right', 'p_left', 'p_right')
# data$R <- factor(ifelse(data$R==0, 'left', 'right'), levels=c('left', 'right'))
#
# # overwrite trial numbers
# data <- EMC2:::add_trials(data)
# head(data)
# save(data, file='~/Projects/EMC_latest/RLEMC/data/data_exp3_newformat.RData')
