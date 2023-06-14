## joint barplots
rm(list=ls())
library(EMC2)
library(emcAdapt)
library(ggplot2)
library(reshape2)
library(corrplot)
library(abind)

source('./utility_funcs_fmri2.R')

to_long <- function(samples) {
  samples_long <- melt(samples[,grepl('2\\|', colnames(samples))], varnames=c('X', 'parameter'))[,2:3]
  regex <- "([0-9])\\|([[:alpha:]]+)\\.([[:alpha:]])_(.*)"
  #gsub("([0-9])\\|([[:alpha:]]+)\\.([[:alpha:]])_(.*)", '\\3', "2|Str.l_response_left")
  samples_long$component <- gsub(regex, '\\1', samples_long$parameter)
  samples_long$ROI <- gsub(regex, '\\2', samples_long$parameter)
  samples_long$hemisphere <- gsub(regex, '\\3', samples_long$parameter)
  samples_long$parameter <- gsub(regex, '\\4', samples_long$parameter)
  samples_long
}

savePlots <- TRUE
ROIs <- c('Amg', 'Cl', 'GPe', 'GPi', 'PAG', 'PPN', 'RN', 'SN', 'STN', 'Str', 'Tha', 'VTA')
allSamples <- list()
for(ROI in ROIs) {
  print(ROI)
  # if(ROI == 'Str')  {
  #   fn <- paste0('./samples/first_joint_model8.RData')
  # } else {
    fn <- paste0('./samples/dataset-trondheim_model-rleam_roi-', ROI, '.RData')
  # }
  load(fn)
  stage <- 'sample'
  tmp = t(do.call(cbind, lapply(samplers, function(x) x$samples$theta_mu[,x$samples$stage==stage])))
  tmp[,grepl('_sd', colnames(tmp))] <- exp(tmp[,grepl('_sd', colnames(tmp))])

  idx <- which(samplers[[1]]$samples$stage==stage)
  correlations <- lapply(samplers, function(x) abind(lapply(idx, function(y) cov2cor(x$samples$theta_var[,,y])), along=3))

  cor1 <- t(do.call(cbind, lapply(correlations, function(x) x['1|B_cuea-s',
                                                                c(paste0('2|', ROI, '.l_cue_ACCminSPD'), paste0('2|', ROI, '.r_cue_ACCminSPD')),])))
  cor2 <- t(do.call(cbind, lapply(correlations, function(x) x['1|t0',
                                                                c(paste0('2|', ROI, '.l_cue_ACCminSPD'), paste0('2|', ROI, '.r_cue_ACCminSPD')),])))
  cor3 <- t(do.call(cbind, lapply(correlations, function(x) x['1|B',
                                                              c(paste0('2|', ROI, '.l_cue_ACCminSPD'), paste0('2|', ROI, '.r_cue_ACCminSPD')),])))

  # cor4 <- t(do.call(cbind, lapply(correlations, function(x) x[paste0('2|', ROI, '.l_cue_ACCminSPD'),
  #                                                             paste0('2|', ROI, '.r_cue_ACCminSPD'),])))
  # covar2 <- t(do.call(cbind, lapply(samplers, function(x) x$samples$theta_var['1|t0',
  #                                                                             c(paste0('2|', ROI, '.l_cue_ACCminSPD'), paste0('2|', ROI, '.r_cue_ACCminSPD')),
  #                                                                             x$samples$stage==stage])))
  dimnames(cor1)[[2]] <- paste0(dimnames(cor1)[[2]], '_covar_with_B_cuea-s')
  dimnames(cor2)[[2]] <- paste0(dimnames(cor2)[[2]], '_covar_with_t0')
  dimnames(cor3)[[2]] <- paste0(dimnames(cor3)[[2]], '_covar_with_B')
  # dimnames(cor4)[[2]] <- 'cor_left_right'

  tmp <- cbind(tmp, cor1, cor2, cor3) #, cor4)

  allSamples[[ROI]] <- to_long(tmp)
}


# barplot ----------------------------------------------------------------
allSamplesLong <- do.call(rbind, allSamples)
makeBarPlot <- function(x, plot_order=NULL) {
  means <- aggregate(value~parameter*ROI*hemisphere, x, mean)
  quants <- aggregate(value~parameter*ROI*hemisphere, x, quantile, c(.025, .975))

  mu_summary <- data.frame(parameter=means$parameter, ROI=means$ROI, hemisphere=means$hemisphere, mean=means$value, q025=quants[,4][,1], q975=quants[,4][,2])
  if(is.null(plot_order)) {
    mu_summary$parameter <- factor(mu_summary$parameter)
  } else {
    mu_summary$parameter <- factor(mu_summary$parameter, levels=plot_order)
  }
  mu_summary$ROI <- factor(mu_summary$ROI)
  mu_summary$hemisphere <- factor(mu_summary$hemisphere)

  # labeller <- sapply(unique(x$parameter), function(y) return(x[x$parameter==y,'unit'][1]))
  # names(labeller) <- unique(x$parameter)
  #


  #gridExtra::grid.arrange(p1 , p2, p3)

  params <- plot_order #levels(mu_summary$parameter)

  for(i in 1:4) assign(paste0('p', i),
                       ggplot(mu_summary[as.numeric(mu_summary$parameter)==i,]) +
                         geom_bar( aes(x=ROI, y=mean, fill=hemisphere), stat="identity",  alpha=0.7, position='dodge') +
                         geom_errorbar( aes(x=ROI, group=hemisphere, ymin=q025, ymax=q975),
                                        position=position_dodge(width=0.90), width=0.4, colour="orange", alpha=0.9, size=1.3) +
                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                         ylab(ifelse(i==4, 'Correlation (r)', '% Signal change')) +
                         xlab(NULL) +
                         ggtitle(levels(mu_summary$parameter)[i])
                       )
  # p1 <- ggplot(mu_summary[as.numeric(mu_summary$parameter)==1,]) +
  #     geom_bar( aes(x=ROI, y=mean, fill=hemisphere), stat="identity",  alpha=0.7, position='dodge') +
  #     geom_errorbar( aes(x=ROI, group=hemisphere, ymin=q025, ymax=q975),
  #                    position=position_dodge(width=0.90), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #   ylab('% Signal change') +
  #   xlab(NULL) +
  #   ggtitle(levels(mu_summary$parameter)[1])
  #
  # p2 <- ggplot(mu_summary[as.numeric(mu_summary$parameter)==2,]) +
  #   geom_bar( aes(x=ROI, y=mean, fill=hemisphere), stat="identity",  alpha=0.7, position='dodge') +
  #   geom_errorbar( aes(x=ROI, group=hemisphere, ymin=q025, ymax=q975),
  #                  position=position_dodge(width=0.90), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('% Signal change') + xlab(NULL) +
  #   ggtitle(levels(levels(mu_summary$parameter))[2])
  #
  # p3 <- ggplot(mu_summary[as.numeric(mu_summary$parameter)==3,]) +
  #   geom_bar( aes(x=ROI, y=mean, fill=hemisphere), stat="identity",  alpha=0.7, position='dodge') +
  #   geom_errorbar( aes(x=ROI, group=hemisphere, ymin=q025, ymax=q975),
  #                  position=position_dodge(width=0.90), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('% Signal change') + xlab(NULL) +
  #   ggtitle(levels(levels(mu_summary$parameter))[3])
  #
  # p4 <- ggplot(mu_summary[as.numeric(mu_summary$parameter)==4,]) +
  #   geom_bar( aes(x=ROI, y=mean, fill=hemisphere), stat="identity",  alpha=0.7, position='dodge') +
  #   geom_errorbar( aes(x=ROI, group=hemisphere, ymin=q025, ymax=q975),
  #                  position=position_dodge(width=0.90), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('Correlation (r)') + xlab(NULL)  +
  #   ggtitle(levels(levels(mu_summary$parameter))[4])

  combined <- p1 + p2 + p3 + p4 & theme(legend.position = "right")
  combined + plot_layout(guides = "collect")


  # ggplot(mu_summary) +
  #   geom_bar( aes(x=ROI, y=mean, fill=hemisphere), stat="identity",  alpha=0.7, position='dodge') +
  #   geom_errorbar( aes(x=ROI, group=hemisphere, ymin=q025, ymax=q975),
  #                  position=position_dodge(width=0.90), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #   facet_wrap(. ~ parameter, ncol=2, scales = 'free', labeller = as_labeller(labeller)) #, strip.position = "left") + ylab(NULL) +
    # theme(strip.background = element_blank(),
    #       strip.placement = "outside")
  # ggplot(mu_summary) +
  #   geom_bar( aes(x=ROI, y=mean, fill=hemisphere), stat="identity",  alpha=0.7, position='dodge') +
  #   geom_errorbar( aes(x=ROI, group=hemisphere, ymin=q025, ymax=q975),
  #                  position=position_dodge(width=0.90), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #   facet_wrap(. ~ parameter, ncol=2, scales = 'free', labeller = as_labeller(labeller))
}

allSamplesLong$unit <- '% Signal change'
allSamplesLong$unit[grepl('_covar_', allSamplesLong$parameter)] <- 'Correlation (r)'
allSamplesLong[allSamplesLong$hemisphere=='l','hemisphere'] <- 'Left'
allSamplesLong[allSamplesLong$hemisphere=='r','hemisphere'] <- 'Right'


allSamplesLong2 <- allSamplesLong
## Swap left&right to contra/ipsi
allSamplesLong2[allSamplesLong2$parameter=='response_leftminright' & allSamplesLong2$hemisphere == 'l', 'value'] <- -allSamplesLong2[allSamplesLong2$parameter=='response_leftminright' & allSamplesLong2$hemisphere == 'l', 'value']
allSamplesLong2[allSamplesLong2$parameter=='response_leftminright', 'parameter'] <- 'response_contra-ipsi'

# better names
new_names <- c('Response contra-ipsi', 'ACC - SPD cue (intercept)', 'Reward prediction errors', 'ACC - SPD cue (correlation B~BOLD)')
allSamplesLong2[allSamplesLong2$parameter=='response_contra-ipsi', 'parameter'] <- new_names[1]
allSamplesLong2[allSamplesLong2$parameter=='cue_ACCminSPD', 'parameter'] <- new_names[2]
allSamplesLong2[allSamplesLong2$parameter=='cue_ACCminSPD_covar_with_B_cuea-s', 'parameter'] <- new_names[4]
allSamplesLong2[allSamplesLong2$parameter=='feedback_pe', 'parameter'] <- new_names[3]

pdf('./figures/barplot_joint_models-threshold.pdf', width=12, height=8)
makeBarPlot(allSamplesLong2[allSamplesLong2$parameter %in% new_names,], plot_order=new_names)
dev.off()


# Barplots for t0 ---------------------------------------------------------
new_names <- c('Response contra-ipsi', 'ACC - SPD cue (intercept)', 'Reward prediction errors', bquote('ACC - SPD cue (correlation t0[intercept] ~ BOLD[acc-spd])'))
allSamplesLong3 <- allSamplesLong2
allSamplesLong3[allSamplesLong3$parameter=='cue_ACCminSPD_covar_with_t0', 'parameter'] <- new_names[4]

pdf('./figures/barplot_joint_models-t0.pdf', width=12, height=8)
makeBarPlot(allSamplesLong3[allSamplesLong3$parameter %in% new_names,], plot_order=new_names)
dev.off()


# new_names <- c('Response contra-ipsi', 'ACC - SPD cue', 'Reward prediction errors', 'cor(B, cue_acc-spd)')
# allSamplesLong4 <- allSamplesLong2
# allSamplesLong4[allSamplesLong4$parameter=='cue_ACCminSPD_covar_with_B', 'parameter'] <- new_names[4]
# makeBarPlot(allSamplesLong4[allSamplesLong4$parameter %in% new_names,], plot_order=new_names)

# ##
# #dev.off()
# pdf('./figures/wtf.pdf', width=12, height=12)
# par(mfrow=c(1,1))
# sgma <- apply(samplers[[1]]$samples$theta_var, 1:2, mean)
# rownames(sgma) <- colnames(sgma) <- sub('derivative', 'd', colnames(sgma))
# corrplot::corrplot(cov2cor(sgma))
# dev.off()
