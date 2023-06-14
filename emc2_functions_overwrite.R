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


log_likelihood_joint_new <- function(proposals, dadms, component = NULL) #pars, dadms, component = NULL)
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
        this_ll <- EMC2:::calc_ll_manager(currentPars, dadm, attr(dadm, "model")()$log_likelihood, hasPEs=TRUE)
        predictionErrors <- attr(this_ll, 'predictionErrors')
        row.names(predictionErrors) <- dadm[dadm$winner,'originalTrialNumber']
      } else {   # neural data: bypass calc_ll_manager
        this_ll <- EMC2:::calc_ll_manager(currentPars, dadm, attr(dadm, "model")()$log_likelihood, predictionErrors=predictionErrors)
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
calc_ll_manager_new <- function(proposals, dadm, ll_func, component = NULL, ...){
  if(!is.data.frame(dadm)){
    lls <- EMC2:::log_likelihood_joint(proposals, dadm, component)
  } else{
    c_name <- attr(dadm,"model")()$c_name
    if(is.null(c_name)){ # use the R implementation
      #### SM ADDITIONS
      dots <- list(...)
      if('predictionErrors' %in% names(dots)) {
        lls <- vector(mode='numeric', length=nrow(proposals))
        for(p in 1:nrow(proposals)) {
          lls[p] <- ll_func(proposals[p,], dadm, predictionErrors=dots$predictionErrors[,p])
        }
        # lls <- apply(proposals,1, ll_func, dadm = dadm, predictionErrors = dots$predictionErrors)
      } else if('hasPEs' %in% names(dots)) {
        PEs <- matrix(NA, ncol=nrow(proposals), nrow=length(unique(dadm$trials)))
        lls <- vector(mode='numeric', length=nrow(proposals))
        for(p in 1:nrow(proposals)) {
          this_ll <- ll_func(proposals[p,], dadm)
          lls[p] <- this_ll
          PEs[,p] <- attr(this_ll, 'predictionErrors')
        }
        attr(lls, 'predictionErrors') <- PEs
      } else {
        #### END SM ADDITIONS
        lls <- apply(proposals,1, ll_func, dadm = dadm)
      }

    } else{
      p_types <- attr(dadm,"model")()$p_types
      designs <- list()
      for(p in p_types){
        designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
      }
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA
      n_trials = nrow(dadm)
      if(c_name == "DDM"){
        levels(dadm$R) <- c(0,1)
        pars <- EMC2:::get_pars(proposals[1,],dadm)
        pars <- cbind(pars, dadm$R)
        parameter_char <- apply(pars, 1, paste0, collapse = "\t")
        parameter_factor <- factor(parameter_char, levels = unique(parameter_char))
        parameter_indices <- split(seq_len(nrow(pars)), f = parameter_factor)
        names(parameter_indices) <- 1:length(parameter_indices)
      } else{
        parameter_indices <- list()
      }
      lls <- EMC2:::calc_llcalc_ll(proposals, dadm, constants = constants, n_trials = n_trials, designs = designs, type = c_name, p_types = p_types,
                                   min_ll = log(1e-10), winner = dadm$winner, expand = attr(dadm, "expand"), group_idx = parameter_indices)
    }
  }
  return(lls)
}
