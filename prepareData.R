## prepare data of experiments in Miletic et al (2021) Elife.
# here, we reformat the data to create the S column that contains information about
# which symbols and corresponding probability is presented on the left & right side on each trial

# experiment 1 (classical) ------------------------------------------------
experimentN <- 'exp1'
print(load(paste0('./data/data_', experimentN, '.RData')))

data <- dat[!dat$excl,] # Excluded participants
data$pp <- factor(as.character(data$pp))
data <- data[!is.na(data$rt),]

data$s_left <- data$appear_left
data$s_right <- data$appear_right
data$p_left <- data$p_win_left
data$p_right <- data$p_win_right
data$choice_symbol <- data[cbind(1:nrow(data), match(paste0('s_',data$choiceDirection), colnames(data)))]

data$S <- factor(paste(paste(data$s_left, data$p_left, sep='-'), paste(data$s_right, data$p_right, sep='-'), sep='_'))
data <- data[,c("pp","TrialNumber","S","reward","choiceIsHighP","rt", "choiceDirection", "choice_symbol",
                's_left', 's_right', 'p_left', 'p_right')]

if(max(data$reward)>1) { data$reward <- data$reward/100 }
names(data) <- c("subjects","trials","S","reward","Racc","rt", "R", "Rs", 's_left', 's_right', 'p_left', 'p_right')
data$R <- factor(data$R)

# overwrite trial numbers
data <- add_trials(data)
head(data)
save(data, file='./data/data_exp1_newformat.RData')



# Experiment 3 (SAT) ----------------------------------------
experimentN <- 'exp3'
print(load(paste0('./data/data_', experimentN, '.RData')))

data <- dat[!dat$excl,] # Excluded participants
data$pp <- factor(as.character(data$pp))
data <- data[!is.na(data$rt),]

data$s_left <- data$stimulus_symbol_left
data$s_right <- data$stimulus_symbol_right

data$p_left <- data$p_win_left
data$p_right <- data$p_win_right

#data$choiceIsHighP <- data$choiceIsHighPpreRev
data$choice_symbol <- NA
data$choice_symbol[data$choice_direction == 0] <- as.character(data$stimulus_symbol_left[data$choice_direction == 0])
data$choice_symbol[data$choice_direction == 1] <- as.character(data$stimulus_symbol_right[data$choice_direction == 1])


## ugly hack to get all relevant information
data$S <- factor(paste(paste(data$s_left, data$p_left, sep='-'), paste(data$s_right, data$p_right, sep='-'), sep='_'))

data <- data[,c("pp","trial_nr","S", 'cue', "outcome","choiceIsHighP","rt", "choice_direction", "choice_symbol",
                's_left', 's_right', 'p_left', 'p_right')]
if(max(data$reward)>1) { data$reward <- data$reward/100 }
names(data) <- c("subjects","trials","S", 'cue', "reward","Racc","rt", "R", "Rs", 's_left', 's_right', 'p_left', 'p_right')
data$R <- factor(ifelse(data$R==0, 'left', 'right'), levels=c('left', 'right'))

# overwrite trial numbers
data <- EMC2:::add_trials(data)
head(data)
save(data, file='./data/data_exp3_newformat.RData')



