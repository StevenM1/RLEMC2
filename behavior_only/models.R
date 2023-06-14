## Model setup functions

# Learning models ---------------------------------------------------------
getLearningModel <- function(learningModel='delta') {
  if(learningModel == 'delta') {
    init_par = c('q0')
    adapt_par <- c('alpha')
    learningParameters <- c(init_par, adapt_par)
    learningFlist <- lapply(learningParameters, function(x) as.formula(paste0(x, '~1')))
    learningConstants <- c(q0=log(0.0001))
    adapt_fun <- function(lastValues, parameters, reward) {
      Q <- lastValues[1]
      alpha <- parameters
      Qnew <- Q+alpha*(reward-Q)
      return(Qnew)
    }
  } else if(learningModel == 'vkfbinary') {
    init_par = c('q0', 'volatility0', 'w0')     # predictions0, volatility0, uncertainty0
    learningParameters <- c(init_par, 'alpha')  # alpha here = volatilitylearningrate!
    learningFlist <- lapply(learningParameters, function(x) as.formula(paste0(x, '~1')))
    learningConstants <- c(q0=log(0.0001), alpha=log(.5), volatility0=log(1))  # basically fix everything BUT w0
    adapt_par <- c('alpha', 'w0')
    adapt_fun <- function(lastValues, parameters, reward) {
      ## lastvalues should be a VECTOR of (prediction, volatility, uncertainty)
      prediction <- lastValues[1]  # m
      volatility <- lastValues[2]  # v
      uncertainty <- lastValues[3] # w
      # learning rate = alpha, kalman gain = k

      volatilityLearningRate <- parameters[1]  # lambda
      omega <- parameters[2]  # w0

      # updating
      predictionError <- reward - (1/(1+exp(-prediction)))
      kalmanGain <- (uncertainty+volatility)/(uncertainty+volatility+omega)
      learningRate <- sqrt(uncertainty + volatility)
      newPrediction <- prediction + learningRate*predictionError
      newUncertainty <- (1-kalmanGain) * (uncertainty+volatility)

      wcov <- (1-kalmanGain)*uncertainty
      volatilityError <- (newPrediction-prediction)^2 + newUncertainty + uncertainty - 2*wcov - volatility
      newVolatility <- volatility + volatilityLearningRate * volatilityError

      return(c(newPrediction, newVolatility, newUncertainty))
    }
  } else if(learningModel == 'vkf') {
    init_par = c('q0', 'volatility0', 'w0') # predictions0, volatility0, uncertainty0
    learningParameters <- c(init_par, 'alpha')  # alpha = volatilitylearningrate!
    learningFlist <- lapply(learningParameters, function(x) as.formula(paste0(x, '~1')))
    learningConstants <- c(q0=log(0.0001), alpha=log(.5), volatility0=log(1))
    adapt_par <- c('alpha', 'w0')
    adapt_fun <- function(lastValues, parameters, reward) {
      ## lastvalues should be a VECTOR of (prediction, volatility, uncertainty)
      prediction <- lastValues[1]  # m
      volatility <- lastValues[2]  # v
      uncertainty <- lastValues[3] # w
      # learning rate = k

      volatilityLearningRate <- parameters[1]  # lambda
      sigma2 <- parameters[2]  # w0

      # updating
      predictionError <- reward - prediction
      learningRate <- (uncertainty+volatility)/(uncertainty+volatility+sigma2)
      newPrediction <- prediction + learningRate*predictionError
      newUncertainty <- (1-learningRate) * (uncertainty+volatility)

      wcov <- (1-learningRate)*uncertainty
      volatilityError <- (newPrediction-prediction)^2 + newUncertainty + uncertainty - 2*wcov - volatility
      newVolatility <- volatility + volatilityLearningRate * volatilityError

      return(c(newPrediction, newVolatility, newUncertainty))
    }
  }

  return(list(init_par=init_par,
              parameters=learningParameters,
              Flist=learningFlist,
              constants=learningConstants,
              adapt_par=adapt_par,
              adapt_fun=adapt_fun))
}


# Decision models ---------------------------------------------------------
getDecisionModel <- function(decisionModel='ARD') {
  if(decisionModel == 'ARD') {
    decisionParameters <- c('v0', 'v', 'B', 'A', 't0', 's', 'wd', 'ws')
    decisionFlist <- lapply(decisionParameters, function(x) as.formula(paste0(x, '~1')))
    names(decisionFlist) <- decisionParameters
    decisionConstants <- c(A=log(0), v=log(0), s=log(1))
  } else if(decisionModel == 'RD') {
    decisionParameters <- c('v0', 'v', 'B', 'A', 't0', 's', 'w')
    decisionFlist <- lapply(decisionParameters, function(x) as.formula(paste0(x, '~1')))
    names(decisionFlist) <- decisionParameters
    decisionConstants <- c(A=log(0), v=log(0), s=log(1))
  }

  return(list(parameters=decisionParameters,
              Flist=decisionFlist,
              constants=decisionConstants))
}
