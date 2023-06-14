
make_design_fmri2 <- function(data,
                              model,
                              events=NULL,
                              design_matrix=NULL,
                              hrf_model='glover + derivative',
                              whiten=FALSE,
                              add_intercept=FALSE,
                              ...) {
  data_names <- colnames(data)[!colnames(data) %in% c('subjects', 'run', 'time')]

  # generate parameter names
  if(!is.null(design_matrix)) {
    par_names <- colnames(design_matrix)[colnames(design_matrix) != "subjects"]
  } else {
    par_names <- c()
  }
  if(!is.null(events)) {
    if(hrf_model == 'glover + derivative') {
      par_names_from_events <- unique(events$trial_type)  # events contains all events
      par_names <- c(par_names, paste(rep(par_names_from_events, each=2), c('', '_derivative'), sep=''))
    }
  }
  if(!('intercept' %in% par_names) & add_intercept) {
    par_names <- c(par_names, 'intercept')
  }

  if(!is.null(events)) events_list <- split(events, f=events$subjects)

  dms_list <- split(design_matrix, f=design_matrix$subjects)
  model <- model()
  if(!is.null(events)) {
    model$events <- lapply(events_list, FUN = function(x){
      y <- x[,colnames(x) != "subjects"]
      y
    })
  }
  model$design_matrix <- lapply(dms_list, FUN=function(x) {
    y <- x[,colnames(x) != 'subjects']
    data.matrix(y)
  })

  dots <- list(...)
  if('hkernel' %in% names(dots)) {
    model$hkernel <- hkernel
  }
  if('regressors' %in% names(dots)) {
    model$regressors <- regressors
    model$fast_convolve <- TRUE
  }


  df_par_names <- expand.grid(c(par_names, "sd"), data_names)
  par_names <- paste0(df_par_names[,2], "_", df_par_names[,1])

  n_pars <- length(par_names)
  Flist <- vector("list", n_pars)
  for(i in 1:n_pars){
    Flist[[i]] <- as.formula(paste0(par_names[i], "~1"))
  }
  model_function <- function() {return(model)}
  design <- list(model = model_function, Flist = Flist)
  attr(design, "p_vector") <- par_names
  return(design)
}



normal2 <- function(){
  return(
    list(
      type="MRI",
      p_types=c("sd"), #This is a bit hacky for now
      # Transform to natural scale
      Ntransform=function(x) {

        if(is.null(dim(x))) {
          is_sd <- grepl("_sd", names(x))
          x[is_sd] <- exp(x[is_sd])+.001
        } else {
          is_sd <- grepl("_sd", dimnames(x)[[2]])
          x[,is_sd] <- exp(x[,is_sd])+.001
        }
        return(x)
      },
      # Trial dependent parameter transform
      Ttransform = function(pars,dadm=NULL)
      {
        pars
      },
      # p_vector transform
      transform = function(x) x,
      # Random function for racing accumulators
      rfun=function(lR,pars) rNORMAL(lR,pars),
      # Density function (PDF) for single accumulator
      dfun=function(rt,pars) dNORMAL(rt,pars),
      # Probability function (CDF) for single accumulator
      pfun=function(rt,pars) pNORMAL(rt,pars),
      # Race likelihood combining pfun and dfun
      log_likelihood=function(p_vector, dadm, predictionErrors=NULL, min_ll=log(1e-10)){
        # data
        y <- as.matrix(dadm[,!colnames(dadm) %in% c("subjects", 'run', 'time', "trials")])

        # first part of the design matrix is already generated, and fixed across parameters
        X <- attr(dadm, 'model')()$design_matrix[[as.character(dadm$subjects[1])]]

        # transform parameters
        p_vector <- normal2()$Ntransform(p_vector)

        # get events once
        if(!is.null(predictionErrors)) {
          if(!is.null(attr(dadm, 'model')()$fast_convolve)) {
            events <- attr(dadm, "model")()$events[[as.character(dadm$subjects[1])]]
            newPEvector <- rep(0, nrow(events))
            newPEvector[as.numeric(names(predictionErrors))+1] <- predictionErrors  # shift by 1: trial_nr = 0-indexed!!
            regressor <- attr(dadm,'model')()$regressors[[as.character(dadm$subjects[1])]]
            hkernel <- attr(dadm,'model')()$hkernel
            modulator <- (newPEvector-mean(newPEvector))/sd(newPEvector)
            X2 <- do.call(rbind, lapply(1:length(regressor), function(run) {
              EMC2:::quick_convolve(regressor[[run]], modulator[events$run==run], hkernel, frame_times=regressor[[run]]$frame_times)
            }))
          } else {
            eventsOrig <- attr(dadm, "model")()$events[[as.character(dadm$subjects[1])]]
            events <- merge(eventsOrig, data.frame(trial_nr=names(predictionErrors), predictionErrors=predictionErrors), on='trial_nr')
            events$modulation <- (events$predictionErrors - mean(events$predictionErrors)) / sd(events$predictionErrors)
            runs <- unique(events$run)
            X2 <- as.matrix(do.call(rbind, lapply(runs, function(x) {
              make_fmri_design_matrix(frame_times=dadm[dadm$run==x,'time'],
                                      events=events[events$run==x,],
                                      hrf_model = 'glover + derivative',
                                      add_intercept=FALSE, oversampling = 5)
            })))
          }


          # combine 'fixed' and modulation parts of the DM
          X <- cbind(X, X2)
        }

        # grab the right parameters
        is_sd <- grepl("sd", names(p_vector))
        sigma <- p_vector[is_sd]
        betas <- p_vector[!is_sd]
        betas <- matrix(betas, ncol = length(sigma))

        # get rid of intercept
        y_hat <- X %*% betas
        y_hat <- y_hat - mean(y_hat)

        ## SM style
        if(ncol(betas > 1)) {
          total_sum <- sum(pmax(dnorm(as.matrix(y), mean = y_hat, sd = rep(sigma, each=nrow(X)), log = T)))
        } else {
          total_sum <- sum(pmax(dnorm(y, mean = y_hat, sd = sigma, log = T), min_ll))
        }

        return(total_sum)
      }
    )
  )
}

# make_design_fmri <- function(design_matrix, data, model, whiten = FALSE){
#   if(is.null(design_matrix$subjects)) stop("Design matrix must have a subjects column")
#   if(is.null(data$subjects)) stop("data must have a subjects column")
#   par_names <- colnames(design_matrix)[colnames(design_matrix) != "subjects"]
#   data_names <- colnames(data)[colnames(data) != "subjects"]
#
#   dm_list <- split(design_matrix, f = design_matrix$subjects)
#   model <- model()
#   model$design_matrix <- lapply(dm_list, FUN = function(x){
#     y <- x[,colnames(x) != "subjects"]
#     data.matrix(y)
#   })
#
#   if(!whiten){
#     df_par_names <- expand.grid(c(par_names, "sd"), data_names)
#     par_names <- paste0(df_par_names[,2], "_", df_par_names[,1])
#   } else{
#     df_par_names <- expand.grid(c(par_names, "sd", "rho"), data_names)
#     par_names <- paste0(df_par_names[,2], "_", df_par_names[,1])
#     model$toeplitz <- lapply(model$design_matrix, FUN = function(x) return(toeplitz(0:(nrow(x) - 1))))
#   }
#   n_pars <- length(par_names)
#   Flist <- vector("list", n_pars)
#   for(i in 1:n_pars){
#     Flist[[i]] <- as.formula(paste0(par_names[i], "~1"))
#   }
#   model_function <- function() {return(model)}
#   design <- list(model = model_function, Flist = Flist)
#   attr(design, "p_vector") <- par_names
#   return(design)
# }
#
# make_beta_matrix <- function(betas){
#   return(matrix(betas))
# }
#
#
# normal <- function(){
#   return(
#     list(
#       type="MRI",
#       p_types=c("sd"), #This is a bit hacky for now
#       # Transform to natural scale
#       Ntransform=function(x) {
#
#         x[,dimnames(x)[[2]] == "sd"] <- exp(x[,dimnames(x)[[2]] == "sd"])
#         x
#       },
#       # Trial dependent parameter transform
#       Ttransform = function(pars,dadm=NULL)
#       {
#         pars
#       },
#       # p_vector transform
#       transform = function(x) x,
#       # Random function for racing accumulators
#       rfun=function(lR,pars) rNORMAL(lR,pars),
#       # Density function (PDF) for single accumulator
#       dfun=function(rt,pars) dNORMAL(rt,pars),
#       # Probability function (CDF) for single accumulator
#       pfun=function(rt,pars) pNORMAL(rt,pars),
#       # Race likelihood combining pfun and dfun
#       log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
#         y <- dadm[,!colnames(dadm) %in% c("subjects", "trials")]
#         X <- attr(dadm, "model")()$design_matrix[[as.character(dadm$subjects[1])]]  ## tada
#         is_sd <- grepl("sd", names(p_vector))
#
#         #        betas[betas < -10] <- -10
#         #        betas[betas > 10] <- 10
#         sigma <- exp(p_vector[is_sd])+.001
#         betas <- p_vector[!is_sd]
#         betas <- matrix(betas, ncol = length(sigma))
#
#         #        writeLines(names(p_vector), fileConn1, sep='\t')
#         #        writeLines(as.character(p_vector), fileConn2, sep='\t')
#         #        writeLines(as.character(X), fileConn3, sep='\t')
#         #print(p_vector)
#         # print('')
#         # print(p_vector)
#         # print(dadm[1,])
#
#         ## get rid of intercept
#         y_hat <- X%*%betas
#         y_hat <- y_hat-mean(y_hat)
#
#         ## SM style
#         if(ncol(betas > 1)) {
#           total_sum <- sum(pmax(dnorm(as.matrix(y), mean = y_hat, sd = rep(sigma, each=nrow(X)), log = T)))
#         } else {
#           total_sum <- sum(pmax(dnorm(y, mean = y_hat, sd = sigma, log = T), min_ll))
#         }
#
#         # if(ncol(betas) == 1){
#         #   total_sum <- sum(pmax(dnorm(y, mean = X %*% betas, sd = sigma, log = T), min_ll))
#         # } else{
#         #   total_sum <- 0
#         #   for(i in 1:ncol(betas)){
#         #     total_sum <- total_sum + sum(pmax(dnorm(y[,i], mean = X %*% betas[,i], sd = sigma[i], log = T), min_ll))
#         #   }
#         # }
#
#         # print(lls)
#         #save(y,X,p_vector,lls,dadm,file='./samples-last_iteration_values.RData')
#         return(total_sum)
#         #        return(max(sum(dnorm(y, mean=X %*% betas, sd = sigma, log = T)), min_ll*length(y)))
#       }
#     )
#   )
# }
#
#
# normal_whitened <- function(){
#   return(
#
#     list(
#       type="MRI",
#       p_types=c("sd", "rho"), #This is a bit hacky for now
#       # Transform to natural scale
#       Ntransform=function(x) {
#
#         get_p_types <- function(nams)
#           unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
#
#         nams <- get_p_types(dimnames(x)[[2]])
#         x[,nams == "sd"] <- exp(x[,nams == "sd"])
#         x
#       },
#       # Trial dependent parameter transform
#       Ttransform = function(pars,dadm=NULL)
#       {
#         pars
#       },
#       # p_vector transform
#       transform = function(x) x,
#       # Random function for racing accumulators
#       rfun=function(lR,pars) rNORMAL(lR,pars),
#       # Density function (PDF) for single accumulator
#       dfun=function(rt,pars) dNORMAL(rt,pars),
#       # Probability function (CDF) for single accumulator
#       pfun=function(rt,pars) pNORMAL(rt,pars),
#       # Race likelihood combining pfun and dfun
#       log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
#         y <- dadm[,!colnames(dadm) %in% c("subjects", "trials")]
#         X <- attr(dadm, "model")()$design_matrix[[dadm$subjects[1]]]
#         order <- attr(dadm, "model")()$toeplitz[[dadm$subjects[1]]]
#         is_sd <- grepl("sd", names(p_vector))
#         is_rho <- grepl("rho", names(p_vector))
#
#         sigma <- exp(p_vector[is_sd])
#         betas <- p_vector[!is_sd & !is_rho]
#         betas <- matrix(betas, ncol = length(sigma))
#         rho <- p_vector[is_rho]
#         V <- rho ** order * (sigma**2)
#         return(max(mvtnorm::dmvnorm(as.matrix(y), mean=X %*% betas, sigma = V, log = T, checkSymmetry = F), min_ll*length(y)))
#
#       }
#     )
#   )
# }
#
#
#
#
#
#
#
#
# make_design_fmri_HIT <- function(events=NULL, data, model,
#                                  design_matrix=NULL, HIT_parametrize=TRUE, add_intercept=FALSE,
#                                  hrf_model=NULL, whiten=FALSE) {
#   data_names <- colnames(data)[!colnames(data) %in% c('subjects', 'run', 'time')]
#
#   # generate parameter names
#   if(!is.null(design_matrix)) {
#     par_names <- colnames(design_matrix)[colnames(design_matrix) != "subjects"]
#   } else {
#     par_names <- c()
#   }
#   if(!is.null(hrf_model)) {
#     if(hrf_model == 'glover + derivative') {
#       par_names_from_events <- unique(events$trial_type)  # events contains all events
#       par_names <- c(par_names, paste(rep(par_names_from_events, each=2), c('', '_derivative'), sep=''))
#     }
#   }
#   if(!'intercept' %in% par_names & add_intercept) {
#     par_names <- c(par_names, 'intercept')
#   }
#
#   if(!is.null(events)) {
#     events_list <- split(events, f=events$subjects)
#   }
#   dms_list <- split(design_matrix, f=design_matrix$subjects)
#
#   ## normalize as required for HIT
#   if(HIT_parametrize) {
#     dms_list = lapply(dms_list, function(x) apply(x, 2, function(y) y/sqrt(sum(y^2))))
#
#     is_derivative <- grepl('_derivative', par_names)
#     alphas <- !is.na(match(paste0(par_names, '_derivative'), par_names))
#     par_names[alphas] <- paste0(par_names[alphas], '_alpha')
#   }
#
#   model <- model()
#   if(!is.null(events)) {
#     model$events <- lapply(events_list, FUN = function(x){
#       y <- x[,colnames(x) != "subjects"]
#       y
#     })
#   }
#   model$design_matrix <- lapply(dms_list, FUN=function(x) {
#     y <- x[,colnames(x) != 'subjects']
#     data.matrix(y)
#   })
#
#   df_par_names <- expand.grid(c(par_names, "sd"), data_names)
#   par_names <- paste0(df_par_names[,2], "_", df_par_names[,1])
#
#   n_pars <- length(par_names)
#   Flist <- vector("list", n_pars)
#   for(i in 1:n_pars){
#     Flist[[i]] <- as.formula(paste0(par_names[i], "~1"))
#   }
#   model_function <- function() {return(model)}
#   design <- list(model = model_function, Flist = Flist)
#   attr(design, "p_vector") <- par_names
#   return(design)
# }
#
# normal_HIT <- function(){
#   return(
#     list(
#       type="MRI",
#       p_types=c("sd"), #This is a bit hacky for now
#       # Transform to natural scale
#       Ntransform=function(x) {
#
#         if(is.matrix(x)) {
#           pnames <- dimnames(x)[[2]]
#           # sigma is estimated on log-scale, add a minor shift to prevent SD=0
#           is_sd <- grepl('_sd', pnames)
#           x[,is_sd] <- exp(x[,is_sd])+.001
#
#           # HIT parametrisation: we estimate two parameters:
#           # 1. alpha: the signed total effect size
#           # 2. beta2: the parameter of the *derivative*, estimated as a signed proportion of alpha
#           is_alpha <- grepl('_alpha', pnames)
#           alpha <- x[,is_alpha]
#           #          alpha <- exp(alpha) ## alpha is on log scale
#
#           ## beta2 from real to signed proportion of alpha
#           is_beta2 <- grepl('_derivative', pnames)
#           beta2 <- x[,is_beta2]
#
#           # beta2 on probit scale?
#           #          beta2 <- 2*(pnorm(beta2)-.5)*alpha
#           beta2 <- pnorm(beta2)*alpha
#           beta1 <- sign(alpha)*(sqrt(alpha^2-beta2^2))
#
#           x[,is_alpha] <- beta1
#           x[,is_beta2] <- beta2
#           dimnames(x)[[2]] <- sub('_alpha', '', dimnames(x)[[2]])  # back to 'normal' beta1/beta2 parametrisation
#         } else {
#           # sigma is estimated on log-scale, add a minor shift to prevent SD=0
#           is_sd <- grepl('_sd', names(x))
#           x[is_sd] <- exp(x[is_sd])+.001
#
#           # HIT parametrisation: we estimate two parameters:
#           # 1. alpha: the signed total effect size
#           # 2. beta2: the parameter of the *derivative*, estimated as a signed proportion of alpha
#           is_alpha <- grepl('_alpha', names(x))
#           alpha <- x[is_alpha]
#           #          alpha <- exp(alpha) ## alpha is on log scale
#           is_beta2 <- grepl('_derivative', names(x))
#           beta2 <- x[is_beta2]
#
#           ## beta2 from real to signed proportion of alpha
#           #          beta2 <- 2*(pnorm(beta2)-.5)*alpha
#           beta2 <- pnorm(beta2)*alpha
#           beta1 <- sign(alpha)*(sqrt(alpha^2-beta2^2))
#           #        beta1 <- sign(alpha)*sqrt(alpha^2-beta2^2)
#
#           x[is_alpha] <- beta1
#           x[is_beta2] <- beta2
#         }
#         return(x)
#       },
#       # Trial dependent parameter transform
#       Ttransform = function(pars,dadm=NULL)
#       {
#         pars
#       },
#       # p_vector transform
#       transform = function(x) x,
#       # Random function for racing accumulators
#       rfun=function(lR,pars) rNORMAL(lR,pars),
#       # Density function (PDF) for single accumulator
#       dfun=function(rt,pars) dNORMAL(rt,pars),
#       # Probability function (CDF) for single accumulator
#       pfun=function(rt,pars) pNORMAL(rt,pars),
#       # Race likelihood combining pfun and dfun
#       log_likelihood=function(p_vector, dadm, predictionErrors=NULL, min_ll=log(1e-10)){
#         y <- as.matrix(dadm[,!colnames(dadm) %in% c("subjects", 'run', 'time', "trials")])
#
#         # first part of the design matrix is already generated
#         X <- attr(dadm, 'model')()$design_matrix[[as.character(dadm$subjects[1])]]
#
#         # now it gets tricky: add prediction errors first
#         if(!is.null(predictionErrors)) {
#           events <- attr(dadm, "model")()$events[[as.character(dadm$subjects[1])]]
#           events <- merge(events, predictionErrors, on='trial_nr')
#           events$modulation <- (events$predictionErrors - mean(events$predictionErrors)) / sd(events$predictionErrors)
#           runs <- unique(events$run)
#           X2 <- as.matrix(do.call(rbind, lapply(runs, function(x) {
#             make_fmri_design_matrix(frame_times=dadm[dadm$run==x,'time'],
#                                     events=events[events$run==x,],
#                                     hrf_model = 'glover + derivative',
#                                     add_intercept=FALSE, oversampling = 5)
#           })))
#
#           # combine 'fixed' and modulation parts of the DM
#           X <- cbind(X, X2)
#         }
#         normfunc <- normal_HIT()$Ntransform
#
#         p_vector <- normfunc(p_vector)
#         is_sd <- grepl("sd", names(p_vector))
#         sigma <- p_vector[is_sd]
#         betas <- p_vector[!is_sd]
#         betas <- matrix(betas, ncol = length(sigma))
#
#         ## Re-calculate betas to peak X=1
#         #        betas <- betas*apply(X, 2, max)
#
#         ## demean to get rid of intercept?
#         y_hat <- X %*% betas
#         #        y_hat <- y_hat-mean(y_hat)
#
#         ## SM style
#         if(ncol(betas > 1)) {
#           total_sum <- sum(pmax(dnorm(as.matrix(y), mean = y_hat, sd = rep(sigma, each=nrow(X)), log = T)))
#         } else {
#           total_sum <- sum(pmax(dnorm(y, mean = y_hat, sd = sigma, log = T), min_ll))
#         }
#
#         return(total_sum)
#       }
#     )
#   )
# }
#
#
#
