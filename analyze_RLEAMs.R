#### RL-EMC!
rm(list=ls())
library(EMC2)
library(emcAdapt)

## settings
experimentN <- 'exp3'
decisionModel <- 'ARD'
learningModel <- 'delta'
save_fn_samples <- paste0('/home/stevenm/Projects/RLEMC2/samples/data-', experimentN, '_model-', decisionModel, '-', learningModel, '.RData')

## load data
print(load(paste0('~/Projects/EMC_latest/RLEMC/data/data_', experimentN, '_newformat.RData')))

## load samples
print(load(save_fn_samples)); chain_n(samplers); filter <- ifelse(chain_n(samplers)[1,3]>0, 'sample', 'burn')
plot_chains(samplers, selection='mu', filter=filter)

# posteriors? look good to me
samples_combined <- do.call(cbind, lapply(samplers, function(x) x$samples$theta_mu[,x$samples$stage=='sample']))
samples_combined <- rdmRLARD()$Ntransform(t(samples_combined))
round(apply(samples_combined, 2, mean), 3)

## other checks
##
EMC2:::iat_pmwg(samplers)
EMC2:::gd_pmwg(samplers, return_summary = TRUE)


## posteriors
# Plotting ----------------------------------------------------------------
calculateByBin <- function(pp, data, n_cores=15, byColumn=NULL) {
  ## calculate bin numbers by S
  if(!is.null(byColumn))  {
    pp$bin <- pp[,byColumn]
    data$bin <- data[,byColumn]
  } else {
    for(S in unique(pp$S)) {
      idx <- pp$S==S
      pp[idx,'bin'] <- as.numeric(cut(pp[idx,'trials'], breaks=10, ordered_result=TRUE))

      idx <- data$S==S
      data[idx,'bin'] <- as.numeric(cut(data[idx,'trials'], breaks=10, ordered_result=TRUE))
    }
  }

  if(n_cores > 1) {
    library(snowfall)
    cl = sfInit(parallel=TRUE, cpus=n_cores)
    #  ppDf <- data.frame(pp)
    #  sfExport('ppDf')
    ppBySub <- lapply(unique(pp$subjects), function(x) pp[pp$subjects==x,])
    #  out = sfLapply(ppBySub, function(x) aggregate(rt~bin*S*postn*subjects, x, quantile, seq(.1, .9, .4)))
    rtByBinBySByPostnBySubjects <- do.call(rbind, sfLapply(ppBySub, function(x) aggregate(rt~bin*S*postn*subjects, x, quantile, seq(.1, .9, .4))))
    accByBinBySByPostnBySubjects <- do.call(rbind, sfLapply(ppBySub, function(x) aggregate(acc~bin*S*postn*subjects, x, mean)))

    sfStop()
  } else {
    rtByBinBySByPostnBySubjects <- aggregate(rt~bin*S*postn*subjects, data.frame(pp), quantile, probs=seq(.1,.9,.4))
    accByBinBySByPostnBySubjects <- aggregate(acc~S*postn*bin*subjects, pp, mean)
  }

  # Response times per trial bin in PP
  rtByBinPP <- aggregate(cbind(`10%`, `50%`, `90%`)~bin*S,                                                 # THIRD STEP: aggregate across posterior predictives (quantile .025, 0.975)
                         aggregate(rt~bin*S*postn,                                                         # SECOND STEP: aggregate across subjects (mean)
                                   rtByBinBySByPostnBySubjects, #aggregate(rt~bin*S*postn*subjects, pp, quantile, probs=seq(.1,.9,.4)),  # FIRST STEP: calculate by bin (aggregate over trials)
                                   mean),
                         quantile, c(.025, 0.975))

  accByBinPP <- aggregate(acc~bin*S,                                                          # THIRD STEP: aggregate across predictives (get .025 and .0975 quantils)
                          aggregate(acc~S*postn*bin,                                          # SECOND STEP: aggregate across subjects (mean)
                                    accByBinBySByPostnBySubjects, #aggregate(acc~S*postn*bin*subjects, pp, mean),       # FIRST STEP: calculate by bin (aggregate over trials)
                                    mean),
                          quantile, c(.025, .975))

  # Response times per trial bin in data
  # data$acc <- data$R=='high'
  rtByBinData <- aggregate(rt~bin*S,
                           aggregate(rt~bin*S*subjects, data, quantile, c(.1, .5, .9)),
                           mean)
  accByBinData <- aggregate(acc~bin*S,
                            aggregate(acc~bin*S*subjects, data, mean),
                            mean)

  return(list('rtByBinPP'=rtByBinPP, 'accByBinPP'=accByBinPP, 'rtByBinData'=rtByBinData, 'accByBinData'=accByBinData))
}

plotPosteriors <- function(rtByBinPP, accByBinPP, rtByBinData, accByBinData) {
  xminmax <- range(rtByBinPP$bin)
  # par(mfrow=c(2,4))
  # accuracy
  yminmax <- range(accByBinPP$acc)
  for(S in unique(rtByBinData$S)) {
    plot(accByBinData$bin[accByBinData$S==S], accByBinData$acc[accByBinData$S==S], type='l', xlab='Trial bin', ylab='Accuracy', ylim=yminmax, lwd=2, main=S)
    polygon(c(accByBinPP$bin[accByBinPP$S==S], rev(accByBinPP$bin[accByBinPP$S==S])),
            c(accByBinPP$acc[,2][accByBinPP$S==S], rev(accByBinPP$acc[,1][accByBinPP$S==S])), col=rgb(t(col2rgb( "red" )), alpha=50, maxColorValue=255), border = FALSE)
  }

  #RT
  for(S in unique(rtByBinData$S)) {
    plot(0,0,type='n',xlab='Trial bin', ylab='RT', ylim=c(0.4, 1.1), xlim=xminmax, main=S)
    for(q in 1:3) {
      lines(rtByBinData$bin[rtByBinData$S==S], rtByBinData[,q+2][rtByBinData$S==S], lwd=2)
      polygon(c(rtByBinPP$bin[rtByBinPP$S==S], rev(rtByBinPP$bin[rtByBinPP$S==S])),
              c(rtByBinPP[rtByBinPP$S==S,q+2][,2], rev(rtByBinPP[rtByBinPP$S==S,q+2][,1])), col=rgb(t(col2rgb( "red" )), alpha=50, maxColorValue=255), border = FALSE)
    }
  }
}


pp <- EMC2:::post_predict(samplers, n_cores = 20)

if(experimentN == 'exp1') {
  # make S column with difficulty
  if(!'Sorig' %in% colnames(data)) data$Sorig <- data$S
  data$S <- as.factor(abs(data$p_left-data$p_right))
  if(!'Sorig' %in% colnames(pp)) pp$Sorig <- pp$S
  pp$S <- as.factor(abs(pp$p_left-pp$p_right))
  pp$acc <- as.character(pp$R)==as.character(pp$correct_direction)
  data$acc <- data$Racc

  # get descriptives per bin
  descriptives <- calculateByBin(pp, data, n_cores=20)
  for(varname in names(descriptives)) {
    assign(varname, descriptives[[varname]])
  }

  # plot
  par(mfrow=c(2,4))
  plotPosteriors(rtByBinPP=rtByBinPP, accByBinPP=accByBinPP, rtByBinData=rtByBinData, accByBinData=accByBinData)
}

if(experimentN == 'exp2') {
  ## reversal learning
  # or collapse across all stimuli
  accuracy <- function(x) {
    x$acc <- 0
    preReversal <- x$trialNreversal < 0
    x[x$R=='right' & (x$p_right > x$p_left) & preReversal, 'acc'] <- 1
    x[x$R=='left' & (x$p_right < x$p_left) & preReversal, 'acc'] <- 1

    # FLIP!
    x[x$R=='right' & (x$p_right < x$p_left) & !preReversal, 'acc'] <- 1
    x[x$R=='left' & (x$p_right > x$p_left) & !preReversal, 'acc'] <- 1
    return(x$acc)
  }
  pp2 <- pp[!is.na(pp$R),]
  data$acc <- accuracy(data)
  pp2$acc <- accuracy(pp2)

  if(!'Sorig' %in% colnames(data)) data$Sorig <- data$S
  data$S <- 'A'
  if(!'Sorig' %in% colnames(pp2)) pp2$Sorig <- pp$S
  pp2$S <- 'A'

  descriptives <- calculateByBin(pp2, data, n_cores=20, byColumn='trialNreversal')
  for(varname in names(descriptives)) {
    assign(varname, descriptives[[varname]])
  }

  par(mfrow=c(2,1))
  plotPosteriors(rtByBinPP=rtByBinPP, accByBinPP=accByBinPP, rtByBinData=rtByBinData, accByBinData=accByBinData)
}


if(experimentN == 'exp3') {
  # make S column with difficulty
  if(!'Sorig' %in% colnames(data)) data$Sorig <- data$S
  data$S <- as.factor(data$cue)
  if(!'Sorig' %in% colnames(pp)) pp$Sorig <- pp$S
  pp$S <- as.factor(pp$cue)
  pp$acc <- as.character(pp$R)==as.character(pp$correct_direction)
  data$acc <- data$Racc

  ## tmp
  data$trialsOrig <- data$trials
  pp$trialsOrig <- pp$trials
  data$trials <- data$trials %% (max(data$trials)/3)
  pp$trials <- pp$trials %% (max(pp$trials)/3)

  # get descriptives per bin
  descriptives <- calculateByBin(pp, data, n_cores=20)
  for(varname in names(descriptives)) {
    assign(varname, descriptives[[varname]])
  }

  # plot
  par(mfrow=c(2,2))
  plotPosteriors(rtByBinPP=rtByBinPP, accByBinPP=accByBinPP, rtByBinData=rtByBinData, accByBinData=accByBinData)
}


# tmp <- pp
# tmp <- tmp[tmp$postn==1,]
# tmp <- tmp[order(tmp$subjects, tmp$trials),]
# cols <- c('subjects', 'trials', 'rt', 'R', 'Rs', 's_left', 's_right', 'p_left', 'p_right', 'acc')
# tmp[,cols]

# tmp$difficulty <- abs(tmp$p_left-tmp$p_right)
# aggregate(acc~difficulty, tmp, mean)
#
#
# debugonce(update_pars); post_predict(samplers, n_cores=1)


tmp <- aggregate(acc~trials, data, mean)
plot(tmp$trials, tmp$acc)

tmp <- aggregate(acc~trials, pp2, mean)
plot(tmp$trials, tmp$acc)

tmp2 <- pp[order(pp$subjects, pp$trials, pp$postn),]
