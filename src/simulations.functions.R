source("./functions.R")
library(doParallel)
library(caret)
library(reshape2)


# Optimal policy functions ------------------------------------------------

GetOptimalDecisionsForFirstSimulation <- function(covariates) {
  stopifnot(is.data.frame(covariates))
  stopifnot("V1" %in% names(covariates))
  return(as.matrix(1 + 0.5 * covariates$V1 + 0.5 * covariates$V2))
}


# Q functions  ------------------------------------------------------------

GetQfunctionValuesForFirstSimulation <- function(covariates, given.treatement, optimal.treatment) {
  linear.part <-  8 + 4 * covariates$V1 - 2 * covariates$V2 - 2 * covariates$V3
  non.linear.part <- - 25 * (optimal.treatment - given.treatement) ** 2
  return (linear.part + non.linear.part)
}


# RewardFunctions ---------------------------------------------------------

GetRewardGivenQfunctionValuesAsMeanVec <- function(q.function.values, variance=1,
                                                   only.positive=TRUE, eps=0.01) {
  raw.reward <- rnorm(NROW(q.function.values), q.function.values, variance)
  reward <-  raw.reward - min(raw.reward) + eps
  return (list(raw.reward=as.matrix(raw.reward), 
               reward=as.matrix(reward)))
}


# Get Data List -----------------------------------------------------------

GetSimulationData <- function(sample.size, number.of.covariates, add.intercept=T) {
  covariates <- matrix(runif(sample.size * number.of.covariates, min=-1, max=1),
                       ncol=number.of.covariates)
  if (isTRUE(add.intercept)) {
    covariates <- model.matrix( ~ ., as.data.frame(covariates))
  }
  treatment <- matrix(runif(sample.size, min=0, max=2), ncol=1)
  optimal.treatment <- GetOptimalDecisionsForFirstSimulation(as.data.frame(covariates))
  q.function.values <- GetQfunctionValuesForFirstSimulation(
      covariates = as.data.frame(covariates),
      given.treatement =  treatment,
      optimal.treatment = optimal.treatment)
  reward.list <- GetRewardGivenQfunctionValuesAsMeanVec(q.function.values)
  stopifnot(is.list(reward.list))
  # FIXME: Find out how to propely estimate propensity scores
  # prop.scores <- GetPropensityScores(
  #     data.frame("treatment"=treatment, "covars"=covariates),
  #     balance.formula = formula (treatment ~ . - covars..Intercept. - 1),
  #     two.step = T)
  prop.scores <- rep(1, nrow(covariates))
   return (
     c(list(covariates=as.matrix(covariates), 
            treatment=as.matrix(treatment),
            prop.scores=as.matrix(prop.scores)), 
       reward.list))
}


# Optimization routines ---------------------------------------------------

# covariates should already contain or not intercept
OptimizeParamsOfPolicyFunction <- function(obs.data, offset, policy.function, lambda,
                                           hyperparams=list(), opt.hyperparams=list()) {
  init.pars <- opt.hyperparams$init.pars
  if (isTRUE("obj.func"  %in% names(opt.hyperparams))) {
    number.of.params <- ncol(obs.data$covariates)
    lower.params.threshold <- rep(-100, number.of.params)
    upper.params.threshold <- rep(100, number.of.params)
    optimized <- GenSA(par = initial.params, fn = opt.hyperparams$obj.func,
        lower.params.threshold, upper.params.threshold,
        control=list(smooth=T, verbose=TRUE, nb.stop.improvement=1000,
                     maxit=6000, temperature = 6000),
        # additional arguments which goes directly to ObjectiveFunction
        obs.data, offset, policy.function, lambda)
    params  <- optimized$par
  } else if (isTRUE("opt.func"  %in% names(opt.hyperparams))) {
    params <- opt.hyperparams$opt.func(params=init.pars, obs.data,  
                                       offset, policy.function, lambda, 
                                       hyperparams=hyperparams,
                                       opt.hyperparams=opt.hyperparams)
  } else {
    stop("Uknown optimization technique! Please provide 'opt.func' key in hyperparams!")
  }
  return (params)
}


# Get DTR value  ----------------------------------------------------------

GetDtrValuesOnFolds <- function(folds, obs.data, offset, control.offset,  
                                policy.function, lambda, opt.hyperparams=list()) {
  dtrs <- numeric(length = 0)
  for (control.fold in folds) {
    tune.obs.data <- list(
      covariates  = as.matrix(obs.data$covariates [-control.fold, ]), 
      treatment   = as.matrix(obs.data$treatment  [-control.fold, ]),
      reward      = as.matrix(obs.data$reward     [-control.fold, ]),
      prop.scores = as.matrix(obs.data$prop.scores[-control.fold, ])
    )
    params <- OptimizeParamsOfPolicyFunction(tune.obs.data, offset,
        policy.function, lambda, opt.hyperparams)
    control.obs.data <- list(
      covariates  = as.matrix(obs.data$covariates [control.fold, ]), 
      treatment   = as.matrix(obs.data$treatment  [control.fold, ]),
      raw.reward      = as.matrix(obs.data$raw.reward [control.fold, ]),
      prop.scores = as.matrix(obs.data$prop.scores[control.fold, ])
    )
    dtr.value <- ValueFunction(params = params,  obs.data = control.obs.data,  
                               offset = control.offset, 
                               policy.function = policy.function)
    dtrs <- c(dtrs, dtr.value)
  }
  return(dtrs)
}



#  Plot decisions versus observed treatment  ------------------------------

PlotDecsionsVersusObserved <- function(obs.data, policy.function, params, offset,
                                       title="DCA decisions versus observed"){
  decision.values <- policy.function(params, obs.data$covariates)
  rewards.scaled.0.1 <- (obs.data$reward - min(obs.data$reward) ) / 
                        (max(obs.data$reward) - min(obs.data$reward))
  plot(decision.values, obs.data$treatment,
       col=rgb(1 - rewards.scaled.0.1, rewards.scaled.0.1, 0),
       pch=20, main=title)
  abline(0, 1)  
  abline(offset, 1, col = "blue", lty = 2)  
  abline(-offset, 1, col = "blue", lty = 2)  
}



# GetSimulationInfo -------------------------------------------------------

GetMetricsForParams  <- function(params, datasets, offset, policy.function, lambda) {
  stat.list  <- list()
  for(data.name in names(datasets)) {
    d = datasets[[data.name]]
    for (name  in names(params)) {
      p = params[[name]] 
      vf = ValueFunction(params = p, obs.data = d, offset = offset,  policy.function)
      objf = ObjectiveFunction(params = p, obs.data = d, offset, policy.function, lambda)

      pred = pmin(pmax(policy.function(p, d$covariates), 0),2)
      optimal = GetOptimalDecisionsForFirstSimulation(as.data.frame(d$covariates))
      value = mean(GetQfunctionValuesForFirstSimulation(
                   covariates = as.data.frame(d$covariates),
                   given.treatement = pred,
                   optimal.treatment = optimal))

      stat.list[[paste("VF", name, toupper(data.name), sep=".")]] = vf
      stat.list[[paste("OBJF", name, toupper(data.name), sep=".")]] = objf
      stat.list[[paste("Qfun", name, toupper(data.name), sep=".")]] = value
    }
  }
return (t(as.matrix(stat.list)))
}



