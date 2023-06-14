rm(list=ls())
library(EMC2)
library(emcAdapt)
library(ggplot2)
library(reshape2)
library(corrplot)

source('./utility_funcs_fmri2.R')

makeBarPlot <- function(x, use_sd=FALSE) {
  means <- apply(x, 2, mean)
  if(use_sd){
    sds <- apply(x, 2, sd)
  } else {
    ## use upper/lower
    quants <- apply(x, 2, quantile, c(.025, .975))
  }
  mu_summary <- data.frame(ROI=names(means), mean=as.numeric(means), lower=quants[1,], upper=quants[2,])

  ggplot(mu_summary) +
    geom_bar( aes(x=ROI, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=ROI, ymin=lower, ymax=upper), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

makeBarPlot2 <- function(x) {
  means <- aggregate(value~parameter*ROI, x, mean)
  quants <- aggregate(value~parameter*ROI, x, quantile, c(.025, .975))

  mu_summary <- data.frame(parameter=means$parameter, ROI=means$ROI, mean=means$value, q025=quants[,3][,1], q975=quants[,3][,2])
  mu_summary$parameter <- factor(mu_summary$parameter)
  mu_summary$ROI <- factor(mu_summary$ROI)

  ggplot(mu_summary) +
    geom_bar( aes(x=parameter, y=mean, fill=ROI), stat="identity",  alpha=0.7, position='dodge') +
    geom_errorbar( aes(x=parameter, group=ROI, ymin=q025, ymax=q975), position=position_dodge(width=0.90), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


savePlots <- TRUE
ROI <- 'RN'

# now what ----------------------------------------------------------------
fn <- paste0('./samples/dataset-trondheim_model-rleam_roi-', ROI, '.RData')
# fn <- paste0('./samples/first_joint_model8.RData')
print(load(fn)); chain_n(samplers)

if(savePlots) pdf(file=paste0('./figures/dataset-trondheim_model-rleam_roi-', ROI, '.pdf'), width=10, height=10)
plot_chains(samplers, filter='sample', selection = 'mu')

# ## Extract main 'interest' here
stage <- 'sample'
allsamples = t(do.call(cbind, lapply(samplers, function(x) x$samples$theta_mu[,x$samples$stage==stage])))
allsamples[,grepl('_sd', colnames(allsamples))] <- exp(allsamples[,grepl('_sd', colnames(allsamples))])


makeBarPlot(allsamples[,1:7])
# makeBarPlot(allsamples[,8:ncol(allsamples)], use_sd = FALSE)



# Calculate peak of BOLD response -----------------------------------------
## OLD stuff - we now have contrasts
# getPeak <- function(beta_dg, beta_deriv, frame_times=seq(0,20,.1)) {
#   events_ <- data.frame(onset=c(0), duration=c(0.001), trial_type='event')
#   X <- as.matrix(make_fmri_design_matrix(frame_times=frame_times, events=events_, hrf_model='glover + derivative', add_intercept=FALSE))
#   y_hat <- X %*% as.numeric(c(beta_dg, beta_deriv))
#   return(max(y_hat))
# }



# left/right hemisphere plotting --------------------------------------------------------------
samples_long <- melt(allsamples[,grepl('2\\|', colnames(allsamples))], varnames=c('X', 'parameter'))[,2:3]
regex <- "([0-9])\\|([[:alpha:]]+)\\.([[:alpha:]])_(.*)"
#gsub("([0-9])\\|([[:alpha:]]+)\\.([[:alpha:]])_(.*)", '\\3', "2|Str.l_response_left")
samples_long$component <- sapply(as.character(samples_long$parameter), function(x) gsub(regex, '\\1', x))
samples_long$ROI <- sapply(as.character(samples_long$parameter), function(x) gsub(regex, '\\2', x))
samples_long$ROI <- sapply(as.character(samples_long$parameter), function(x) gsub(regex, '\\3', x))
samples_long$parameter <- sapply(as.character(samples_long$parameter), function(x) gsub(regex, '\\4', x))



makeBarPlot2(samples_long)


# Correlations -----------------------------------------------------------
par(mfrow=c(1,1))
sgma <- apply(samplers[[1]]$samples$theta_var, 1:2, mean)
rownames(sgma) <- colnames(sgma) <- sub('derivative', 'd', colnames(sgma))
corrplot::corrplot(cov2cor(sgma))

if(savePlots) dev.off()


# pdf(file='./figures/')
# plot_chains(samplers, selection='correlation')
