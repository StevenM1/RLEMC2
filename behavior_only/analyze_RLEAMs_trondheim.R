rm(list=ls())
library(EMC2)
library(emcAdapt)
library(parallel)
source('./behavior_only/utility_funcs_RLEAMs.R')

## settings
task <- 'rlsat'
decisionModel <- 'ARD'
learningModel <- 'delta'
save_fn_samples <- paste0('./samples/dataset-trondheim_task-', task, '_model-', decisionModel, '-', learningModel, '.RData')

## load data
print(load(paste0('./data/dataset-trondheim_task-', task, '.RData')))

## load samples
print(load(save_fn_samples)); chain_n(samplers); filter <- ifelse(chain_n(samplers)[1,4]>0, 'sample', 'burn')
plot_chains(samplers, selection='mu', filter=filter)

# posteriors? look good to me
samples_combined <- do.call(cbind, lapply(samplers, function(x) x$samples$theta_mu[,x$samples$stage=='sample']))
samples_combined <- rdmRLARD()$Ntransform(t(samples_combined))
round(apply(samples_combined, 2, mean), 3)

## other checks
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
    ppBySub <- lapply(unique(pp$subjects), function(x) pp[pp$subjects==x,])
    rtByBinBySByPostnBySubjects <- do.call(rbind, mclapply(ppBySub, function(x) aggregate(rt~bin*S*postn*subjects, x, quantile, seq(.1, .9, .4)), mc.cores = 15))
    accByBinBySByPostnBySubjects <- do.call(rbind, mclapply(ppBySub, function(x) aggregate(acc~bin*S*postn*subjects, x, mean), mc.cores = 15))

  } else {
    rtByBinBySByPostnBySubjects <- aggregate(rt~bin*S*postn*subjects, data.frame(pp), quantile, probs=seq(.1,.9,.4))
    accByBinBySByPostnBySubjects <- aggregate(acc~bin*S*postn*subjects, pp, mean)
  }

  # Response times per trial bin in PP
  rtByBinPP <- aggregate(cbind(`10%`, `50%`, `90%`)~bin*S,                                                 # THIRD STEP: aggregate across posterior predictives (quantile .025, 0.975)
                         aggregate(rt~bin*S*postn,                                                         # SECOND STEP: aggregate across subjects (mean)
                                   rtByBinBySByPostnBySubjects, #aggregate(rt~bin*S*postn*subjects, pp, quantile, probs=seq(.1,.9,.4)),  # FIRST STEP: calculate by bin (aggregate over trials)
                                   mean),
                         quantile, c(.025, 0.975))

  accByBinPP <- aggregate(acc~bin*S,                                                          # THIRD STEP: aggregate across predictives (get .025 and .0975 quantils)
                          aggregate(acc~bin*S*postn,                                          # SECOND STEP: aggregate across subjects (mean)
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

calculatePEByBin <- function(pp, n_cores=15, byColumn=NULL) {
  if(!is.null(byColumn))  {
    pp$bin <- pp[,byColumn]
  } else {
    for(S in unique(pp$S)) {
      idx <- pp$S==S
      pp[idx,'bin'] <- as.numeric(cut(pp[idx,'trials'], breaks=10, ordered_result=TRUE))
    }
  }

  if(n_cores > 1) {
    # library(snowfall)
    # cl = sfInit(parallel=TRUE, cpus=n_cores)
    #  ppDf <- data.frame(pp)
    #  sfExport('ppDf')
    ppBySub <- lapply(unique(pp$subjects), function(x) pp[pp$subjects==x,])
    #  out = sfLapply(ppBySub, function(x) aggregate(rt~bin*S*postn*subjects, x, quantile, seq(.1, .9, .4)))
    PEByBinBySByPostnBySubjects <- do.call(rbind, mclapply(ppBySub, function(x) aggregate(PE~bin*S*postn*subjects, x, mean), mc.cores = 15))

    sfStop()
  } else {
    PEByBinBySByPostnBySubjects <- aggregate(PE~bin*S*postn*subjects, pp, mean)
  }

  PEByBinPP <- aggregate(PE~bin*S,                                                          # THIRD STEP: aggregate across predictives (get .025 and .0975 quantils)
                          aggregate(PE~bin*S*postn,                                          # SECOND STEP: aggregate across subjects (mean)
                                    PEByBinBySByPostnBySubjects, #aggregate(acc~S*postn*bin*subjects, pp, mean),       # FIRST STEP: calculate by bin (aggregate over trials)
                                    mean),
                          quantile, c(.025, .975))

  # Response times per trial bin in data

  return(list('PEByBinPP'=PEByBinPP))
}

pp <- EMC2:::post_predict(samplers, n_cores = 20)
pp$s_left <- s_left(pp)
pp$s_right <- s_right(pp)
pp$p_left <- p_left(pp)
pp$p_right <- p_right(pp)
pp$p_low <- p_low(pp)
pp$p_high <- p_high(pp)
pp$s_low <- s_low(pp)
pp$s_high <- s_high(pp)
pp$lRS <- lRS(pp)
pp$correct_direction <- correct_direction(pp)

# get_learn <- function(samplers, pars, n_cores=1) {
#   data <- attr(samplers,"data_list")
#   design <- attr(samplers,"design_list")
#   model <- attr(samplers,"model_list")
#   n_post <- length(pars)
#
#   for(j in 1:length(data)){
#     subjects <- levels(data[[j]]$subjects)
#     ## TMP OVERWRITE
# #    model[[j]] <- rdmRLARD
#     ##
#
#     data2 <- design_model(
#       EMC2:::add_accumulators(data[[1]],design[[j]]$matchfun,simulate=TRUE,
#                               type=model[[j]]()$type,Fcovariates=design[[j]]$Fcovariates),
#       design[[j]],model[[j]],add_acc=FALSE,compress=FALSE,verbose=FALSE,
#       rt_check=FALSE)
#     if(n_cores == 1) {
#       PEs <- Qvalues <- vector(mode="list",length=n_post)
#       for(i in 1:n_post) {
#         cat('.', sep='')
#         pars2 <- model[[j]]()$Ttransform(model[[j]]()$Ntransform(EMC2:::map_p(
#           model[[j]]()$transform(add_constants(pars[[i]],design[[j]]$constants)),data2
#         )),data2)
#
#         PEs[[i]] <- attr(pars2, 'predictionErrors')
#         Qvalues[[i]] <- attr(pars2, 'Qvalues')
#       }
#     }
#   }
#   return(list(predictionErrors=PEs, Qvalues=Qvalues))
# }

# out <- get_learn(samplers, attr(pp, 'pars'))
# PEs <- lapply(out$predictionErrors, function(x) apply(x, 1, sum, na.rm=TRUE))
# PEs2 <- unlist(PEs) # do.call(rbind, PEs)
# # plot(apply(PEs2, 1, mean))
# pp <- pp[order(pp$postn, pp$subjects, pp$trials),]
# pp$PE <- PEs2  # apply(PEs2, 1, mean)
#
# tmp3 <- aggregate(PE~trialNreversal, pp, mean)
# plot(tmp3$trialNreversal, tmp3$PE, type='l')

# plot(0,0, type='n', xlab='trialN', xlim=c(-40, 40), ylim=c(-.5, .6), ylab='PE')
# doPlot1 <- function(PE, trialN) {
#   tmp <- data.frame(PE=PEs, trialN=trialN)
#   tmp2 <- aggregate(PE~trialN, tmp, mean)
#   lines(tmp2$trialN, tmp2$PE)
# }
# for(i in 1:length(PEs)) {
#   doPlot1(PEs[[i]], data$trialNreversal)
# }



# plot(PE~trials, data)
# data$absPE <- data$PE
# aggregated <- aggregate(absPE~trials, data, mean)
# plot(aggregated$trials, aggregated$absPE)

if(task == 'revl') {
  ## reversal learning
  # or collapse across all stimuli
  accuracy <- function(x) {
    x$acc <- 0
    preReversal <- x$trialNreversal < 0
    x[x$R=='right' & (x$p_right > x$p_left) & preReversal, 'acc'] <- 1
    x[x$R=='left' & (x$p_right < x$p_left) & preReversal, 'acc'] <- 1

    # FLIP!
#    x[x$R=='right' & (x$p_right < x$p_left) & !preReversal, 'acc'] <- 1
#    x[x$R=='left' & (x$p_right > x$p_left) & !preReversal, 'acc'] <- 1
    return(x$acc)
  }
  pp2 <- pp[!is.na(pp$R),]
  data$acc <- accuracy(data)
  pp2$acc <- accuracy(pp2)

  if(!'Sorig' %in% colnames(data)) data$Sorig <- data$S
  data$S <- 'A'
  if(!'Sorig' %in% colnames(pp2)) pp2$Sorig <- pp$S
  pp2$S <- 'A'

  descriptivesPE <- calculatePEByBin(pp2, n_cores=1, byColumn='trialNreversal')
  plot(descriptivesPE$PEByBinPP$bin, descriptivesPE$PEByBinPP$PE[,1], type='l')
  lines(descriptivesPE$PEByBinPP$bin, descriptivesPE$PEByBinPP$PE[,2])

  descriptives <- calculateByBin(pp2, data, n_cores=1, byColumn='trialNreversal')
  for(varname in names(descriptives)) {
    assign(varname, descriptives[[varname]])
  }

  par(mfrow=c(2,1))
  plotPosteriors(rtByBinPP=rtByBinPP, accByBinPP=accByBinPP, rtByBinData=rtByBinData, accByBinData=accByBinData)
}


if(task == 'rlsat') {
  # make S column with difficulty
  if(!'Sorig' %in% colnames(data)) data$Sorig <- data$S
  data$S <- as.factor(data$cue)
  if(!'Sorig' %in% colnames(pp)) pp$Sorig <- pp$S
  pp$S <- as.factor(pp$cue)
  pp$acc <- as.character(pp$R)==as.character(pp$correct_direction)
  data$acc <- data$Racc

  ## tmp
  if(!'trialsOrig' %in% colnames(data)) {
    data$trialsOrig <- data$trials
    data$trials <- data$trials %% (max(data$trials)/3)
  }
  if(!'trialsOrig' %in% colnames(pp)) {
    pp$trialsOrig <- pp$trials
    pp$trials <- pp$trials %% (max(pp$trials)/3)
  }

  # get descriptives per bin
  descriptives <- calculateByBin(pp, data, n_cores=1)
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


### For anne: save prediction errors
meanPP <- aggregate(PE~trials*subjects, pp, mean)
meanPP$S <- pp[pp$postn==1,'S']
meanPP$trialNreversal <- pp[pp$postn==1,'trialNreversal']
meanPP$originalTrialNumber <- pp[pp$postn==1,'originalTrialNumber']


write.csv(meanPP, file='../RLEMC2/prediction_errors_for_anne.csv')
#tmp5 <- aggregate(PE~trialNreversal, meanPP, mean)
#plot(tmp5$trialNreversal, tmp5$PE, type='l')
