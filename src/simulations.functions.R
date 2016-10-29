setwd("~/yandexDisk/DIPLOMA/CODE/src")
source("./functions.R")
library(doParallel)
library(caret)
library(reshape2)
library(data.table)
library(tgp)
library(gridExtra)
library(ggplot2)
# install.packages("truncnorm")
source("../../OWL/O_learning_functions.r")
# hack to overcome   "do_rtruncnorm" not resolved from current namespace (truncnorm)  bug 
rtruncnorm <- function (n, a = -Inf, b = Inf, mean = 0, sd = 1)  {
  if (length(n) > 1) 
    n <- length(n)
  if (length(n) > 1) 
    n <- length(n)
  else if (!is.numeric(n)) 
    stop("non-numeric argument n.")
  .Call("do_rtruncnorm", as.integer(n), a, b, mean, sd, PACKAGE="truncnorm")
}


# Optimal policy functions ------------------------------------------------

GetOptimalDecisionsForFirstSimulation <- function(covariates) {
  stopifnot(is.data.frame(covariates))
  stopifnot("V1" %in% names(covariates))
  return(as.matrix(1 + 0.5 * covariates$V1 + 0.5 * covariates$V2))
}

GetOptimalDecisionsForSecondSimulation <- function(covariates) {
  stopifnot(is.data.frame(covariates))
  stopifnot("V1" %in% names(covariates))
  return (with(data = covariates, expr = {  
      Iv1 <-  (-0.5 < V1)  & (V1 < 0.5) 
      as.matrix(0.6 * Iv1 + 1.2 * (!Iv1)   + V4 ** 2 + 0.5 * log(abs(V7) + 1) - 0.6)
    })
  )
}




# Q functions  ------------------------------------------------------------

GetQfunctionValuesForFirstSimulation <- function(covariates, given.treatement, optimal.treatment) {
  with(as.data.frame(covariates), {
    linear.part <-  8 + 4 * V1 - 2 * V2 - 2 * V3
    non.linear.part <- - 25 * (optimal.treatment - given.treatement) ** 2
    return (linear.part + non.linear.part)
  })
}

GetQfunctionValuesForSecondSimulation <- function(covariates, given.treatement, optimal.treatment) {
  return(with(as.data.frame(covariates), {
    8 + 4 * cos(2 * pi * V2) - 2 * V4 - 8 * V5 ^ 3  - 15 * abs(optimal.treatment - given.treatement)
  }))
}

GetQValue <- function(params, data, policy.function) {
  pred = pmin(pmax(policy.function(params, data$covariates), 0),2)
  return (with(data, { mean(GetQFunctionValues(covariates, pred, optimal.treatment)) }))
}



# RewardFunctions ---------------------------------------------------------

GetRewardGivenQfunctionValuesAsMeanVec <- function(q.function.values, sd=1,
                                                   only.positive=TRUE, eps=0.01) {
  raw.reward <- rnorm(NROW(q.function.values), q.function.values, sd)
  reward <-  raw.reward - min(raw.reward) + eps
  return (list(raw.reward=as.matrix(raw.reward), 
               reward=as.matrix(reward)))
}


#  GetSimulationData ----------------------------------------------------------


# A <- rbinom(sample.size, 1, 0.5)
# 
# sample.size <- 20 
# n.of.covars <- 5
# covariates <- matrix(runif(sample.size * n.of.covars, min=-1, max=1), ncol=n.of.covars)
# treatment <- runif(sample.size, -1, 1)
# mu <- with(as.data.frame(covariates), 1 + V1  + V2 + 2 * V3 +  0.5 * V4)
# delta <- function(X) { with(as.data.frame(X), 0.4 *  (V2 - 0.25 * V1 ^ 2  - 1)) } 
# Q.values <- rnorm(sample.size, mu)
# 
# ggg <- list(x = 1:5, y = 1:5, z = delta(grid))
# 
# image(ggg)
# 
# grid <- expand.grid(V1 = seq(-1, 1, length.out = 100), 
#                     V2 = seq(-1, 1, length.out = 100))



ReplaceExtremeWithUnif <- function(values, min.val, max.val){
  violatros.index <- values < min.val | values > max.val
  values[violatros.index] <- runif(sum(violatros.index), min = min.val, max=max.val)
  return(values)
}

Shvechikov.2.fopt <- function(x) {
  return (((x - .5) / .5) ** 2)
} 


GetDataForSchvechikov.2 <- function(sample.size,  sd, min.A=0, max.A=1, min.X=0, max.X=1){
  GetQFunctionValues <- function(covars, given.A, optimal.A=NULL) {
    if (is.null(optimal.A)) {
      optimal.A <- Shvechikov.2.fopt(covars)
    }
    return (-(given.A - optimal.A) ** 2)
  }
  X <- ReplaceExtremeWithUnif(rgamma(n=sample.size,  shape=3, rate=7.3), min.X, max.X)
  A <- ReplaceExtremeWithUnif(rgamma(n=sample.size,  shape=3, rate=7.3), min.A, max.A)
  A.opt <-  Shvechikov.2.fopt(X)
  Q.vals <-GetQFunctionValues(X, A, A.opt)
  R.list <- GetRewardGivenQfunctionValuesAsMeanVec(Q.vals, sd = sd)
  data <- list(covariates = X, treatment=A, optimal.treatment=A.opt, 
               prop.scores = rep(1, sample.size))
  data$GetQFunctionValues = GetQFunctionValues
  data <- c(data, R.list)
  data <- lapply(data, function(x)  {if (is.numeric(x)) as.matrix(x) else x} )
  data$GetOptimalTreatment <- Shvechikov.2.fopt
  return (data)
}





Shvechikov.1.fopt <- function(x) {
  x = x * 100  # initially function was written to support x from 0 to 100, not from 0 to 1
  step.start = 15
  breakage.start = 40
  first.part <- x**2  *  exp(- x**(.8)) +  (x - step.start) / 50 + .5
  first.and.changed.part <- ( - x**1.9 / 900 + 6 *  sin(x/6) + 1.5 * (sin(x**1.01 ) + 1) * x**(1/7) )* 1/20 
  second.part <-  (- first.part + first.and.changed.part) * (x > breakage.start)
  whole <- first.part + second.part
  f <- function(x) { # scaled Runge function
    x <- x / 60   - 0.3
    cos(3*pi*x)/(1+25*(x-0.25)^2) 
  }
  return ((whole + f(x)) / 2 + 0.35)
} 


GetDataForSchvechikov.1 <- function(sample.size, sd, min.A=0, max.A=1, min.X=0, max.X=100) {
  GetQFunctionValues <- function(covars, given.A, optimal.A=NULL) {
    if (is.null(optimal.A)) {
      optimal.A <- Shvechikov.1.fopt(covars)
    }
    return (-(given.A - optimal.A) ** 2)
  }
  X <- ReplaceExtremeWithUnif(rgamma(n=sample.size,  shape=13, rate=0.3), min.X, max.X) 
  X <- X / 100 # hack to bring X in [0,1] range
  A <- ReplaceExtremeWithUnif(rgamma(n=sample.size,  shape=2, rate=4.21), min.A, max.A)
  A.opt <-  Shvechikov.1.fopt(X)
  Q.vals <-GetQFunctionValues(X, A, A.opt)
  R.list <- GetRewardGivenQfunctionValuesAsMeanVec(Q.vals, sd = sd)
  data <- list(covariates = X, treatment=A, optimal.treatment=A.opt, 
               prop.scores = rep(1, sample.size))
  data$GetQFunctionValues = GetQFunctionValues
  data <- c(data, R.list)
  data <- lapply(data, function(x)  {if (is.numeric(x)) as.matrix(x) else x} )
  data$GetOptimalTreatment <- Shvechikov.1.fopt
  return (data)
}






## COMPATABLIITY WITH OLD CODE IS BROKEN BY INTRODUCING NEW GetQFunctionValue
GetSimulationData <- function(sample.size,  scenario="chen.1", add.intercept=T, ...) {
  
  if (grepl("shvechikov", scenario, ignore.case = T)) { 
    return (switch (scenario,
      "shvechikov.1" = GetDataForSchvechikov.1(sample.size = sample.size, ...), 
      "shvechikov.2" = GetDataForSchvechikov.2(sample.size = sample.size, ...)
    ))
  }
  n.of.covars <- switch(scenario, "chen.1" = 30, "chen.2" = 10)
  GetOptimalDecisions <- switch (scenario,
    "chen.1" = GetOptimalDecisionsForFirstSimulation,
    "chen.2" = GetOptimalDecisionsForSecondSimulation
  )
  GetQFunctionValues <- switch (scenario,
    "chen.1" = GetQfunctionValuesForFirstSimulation,
    "chen.2" = GetQfunctionValuesForSecondSimulation
  )
  
  covariates <- matrix(runif(sample.size * n.of.covars, min=-1, max=1), 
                       ncol=n.of.covars)
  if (isTRUE(add.intercept)) {
    covariates <- model.matrix( ~ ., as.data.frame(covariates))
  }
  treatment <- matrix(runif(sample.size, min=0, max=2), ncol=1)
  
  
  optimal.treatment <- GetOptimalDecisions(as.data.frame(covariates))
  q.function.values <- GetQFunctionValues(
      covariates = as.data.frame(covariates),
      given.treatement =  treatment,
      optimal.treatment = optimal.treatment
  )
  
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
           prop.scores=as.matrix(prop.scores),
           optimal.treatment=optimal.treatment, 
           GetQFunctionValues=GetQFunctionValues), 
     reward.list)
  )
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
      raw.reward  = as.matrix(obs.data$raw.reward [control.fold, ]),
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

ScaleToZeroOne <- function(data) {
  return ((data - min(data))  / (max(data) - min(data)))
}



PlotDecsionsVersusObserved <- function(obs.data, policy.function, params, offset,
                                       title="DCA decisions versus observed", highlihght=NULL){
  decision.values <- policy.function(params, obs.data$covariates)
  rewards.scaled.0.1 <- (obs.data$reward - min(obs.data$reward) ) / 
                        (max(obs.data$reward) - min(obs.data$reward))
  plot(decision.values, obs.data$treatment,
       col=rgb(1 - rewards.scaled.0.1, rewards.scaled.0.1, 0),
       pch=19, main=title,  xlim = c(0,2), xlab="Predicted treatment", 
       ylab="Observed treatment")
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
      qf = GetQValue(p, d, policy.function)

      stat.list[[paste("VF", name, toupper(data.name), sep=".")]] = vf
      stat.list[[paste("OBJF", name, toupper(data.name), sep=".")]] = objf
      stat.list[[paste("Qfun", name, toupper(data.name), sep=".")]] = qf
    }
  }
return  (unlist(stat.list)) 
}



GetScoreForEachParamSet <- function(hyp.grid, opt.fname, opt.fparams, 
  score.func=function(p) QFunctionFirstScenario(p, data=test, policy.function=PolicyFunLinearKernel)) {
  # Replaces necessary params in opt.fparams with row of hyp.grid; 
  # Computes metrics using the score function with two pars: 
  #   opt.result  and  replaced opt.fparams
  # 
  # Args:
  #  score.func - Function to call with result of opt.fname( ) invocation with updated params from hyp.grid.
  #               Should have signature  score.func(opt.result, opt.modified.params)
  opt.results = list()
  stat.list = list()
  for (i in seq_len(NROW(hyp.grid))) {
    opt.modified.params = as.list(modifyList(opt.fparams, hyp.grid[i, , drop=F]))
    param.name = paste(paste(names(hyperparams.grid),  hyp.grid[i, ], sep="="), collapse = ".")
    opt.results[[i]] <- do.call(opt.fname, opt.modified.params)
    stat.list[[i]] <-  score.func(opt.results[[i]])
  } 
  return (unlist(stat.list))
}





# Quantile regression from quantreg hacked to support weights -------------


rq.with.weights <- function (formula, tau = 0.5, data, subset, weights, na.action, 
                             method = "br", model = TRUE, contrasts = NULL, ...) 
{
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), 
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  if (method == "model.frame") 
    return(mf)
  mt <- attr(mf, "terms")
  weights <- as.vector(model.weights(mf))
  Y <- model.response(mf)
  X <- model.matrix(mt, mf, contrasts)
  eps <- .Machine$double.eps^(2/3)
  Rho <- function(u, tau) u * (tau - (u < 0))
  if (length(tau) > 1) {
    print("assuming many taus!")
    if (any(tau < 0) || any(tau > 1)) 
      stop("invalid tau:  taus should be >= 0 and <= 1")
    if (any(tau == 0)) 
      tau[tau == 0] <- eps
    if (any(tau == 1)) 
      tau[tau == 1] <- 1 - eps
    coef <- matrix(0, ncol(X), length(tau))
    rho <- rep(0, length(tau))
    fitted <- resid <- matrix(0, nrow(X), length(tau))
    for (i in 1:length(tau)) {
      z <- {
        if (length(weights)) 
          rq.wfit.with.weights(X, Y, tau = tau[i], weights, method, 
                               ...)
        else rq.fit(X, Y, tau = tau[i], method, ...)
      }
      coef[, i] <- z$coefficients
      resid[, i] <- z$residuals
      rho[i] <- sum(Rho(z$residuals, tau[i]))
      fitted[, i] <- Y - z$residuals
    }
    taulabs <- paste("tau=", format(round(tau, 3)))
    dimnames(coef) <- list(dimnames(X)[[2]], taulabs)
    dimnames(resid) <- list(dimnames(X)[[1]], taulabs)
    fit <- z
    fit$coefficients <- coef
    fit$residuals <- resid
    fit$fitted.values <- fitted
    if (method == "lasso") 
      class(fit) <- c("lassorqs", "rqs")
    else if (method == "scad") 
      class(fit) <- c("scadrqs", "rqs")
    else class(fit) <- "rqs"
  }
  else {
    process <- (tau < 0 || tau > 1)
    if (tau == 0) 
      tau <- eps
    if (tau == 1) 
      tau <- 1 - eps
    fit <- {
      if (length(weights)) 
        rq.wfit.with.weights(X, Y, tau = tau, weights, method, ...)
      else rq.fit(X, Y, tau = tau, method, ...)
    }
    if (process) 
      rho <- list(x = fit$sol[1, ], y = fit$sol[3, ])
    else {
      dimnames(fit$residuals) <- list(dimnames(X)[[1]], 
                                      NULL)
      rho <- sum(Rho(fit$residuals, tau))
    }
    if (method == "lasso") 
      class(fit) <- c("lassorq", "rq")
    else if (method == "scad") 
      class(fit) <- c("scadrq", "rq")
    else class(fit) <- ifelse(process, "rq.process", "rq")
  }
  fit$na.action <- attr(mf, "na.action")
  fit$formula <- formula
  fit$terms <- mt
  fit$xlevels <- .getXlevels(mt, mf)
  fit$call <- call
  fit$tau <- tau
  fit$weights <- weights
  fit$residuals <- drop(fit$residuals)
  fit$rho <- rho
  fit$method <- method
  fit$fitted.values <- drop(fit$fitted.values)
  attr(fit, "na.message") <- attr(m, "na.message")
  if (model) 
    fit$model <- mf
  fit
}



rq.wfit.with.weights <- function (x, y, tau = 0.5, weights, method = "br", ...)  {
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  if (length(tau) > 1) {
    tau <- tau[1]
    warning("Multiple taus not allowed in rq.wfit: solution restricted to first element")
  }
  contr <- attr(x, "contrasts")
  wx <- x * weights
  wy <- y * weights
  fit <- switch(method, fn = rq.fit.fnb(wx, wy, tau = tau,  ...), 
                fnb = rq.fit.fnb(wx, wy, tau = tau, ...), 
                br = rq.fit.br(wx,  wy, tau = tau, ...), 
                fnc = rq.fit.fnc(wx, wy, tau = tau, ...), 
                pfn = rq.fit.pfn(wx, wy, tau = tau, ...), { 
                  what <- paste("rq.fit.", method, sep = "") 
                    if (exists(what, mode = "function")) 
                      (get(what, mode = "function"))(wx,  wy, tau, ...) 
                    else stop(paste("unimplemented method:",  method)) 
                  })
  if (length(fit$sol)) 
    fit$fitted.values <- x %*% fit$sol[-(1:3), ]
  else fit$fitted.values <- x %*% fit$coef
  fit$residuals <- y - fit$fitted.values
  fit$contrasts <- attr(x, "contrasts")
  fit$weights <- weights
  fit
}



# GP ERA  ------------------------------------------------------------



GetBestPredictions <- function(gp_model, s = 2) {
  means <- gp_model$ZZ.mean
  stds <- sqrt(gp_model$ZZ.s2)
  dt <- data.table(gp_model$XX)
  dt[, LB:= means - s * stds]
  col_names <- c(grep('^C', names(dt), value = T))
  return(dt[, .(A_pred=A[which.max(LB)], LB_max=max(LB)), by=col_names])
}

GetValueOfPredictedA <- function(best.estim.dt, data.obj)  {
  best.estim.dt[, mean(data.obj$GetQFunctionValues(C, A_pred))]
}

PlotDecisionSurface <- function(models, s=2) {
  for(m_name in names(models)) {
    m <-  models[[m_name]]
    surf <- matrix(m$ZZ.mean - s * sqrt(m$ZZ.s2), nrow=length(unique(m$XX$A)))
    plt1 <- levelplot(surf, col.regions = gray(0:100/100),  xlab="C",  
                      ylab="A", main=paste("Decision Surface, s = ", s))
    plt2 <- wireframe(surf, xlab="C", ylab="A", zlab="decision surf",  main=m_name, 
                      par.settings = list(axis.line = list(col = "transparent")))
    grid.arrange(plt1, plt2, ncol=2)
  }
}

FitAndPlotAllModels <- function(noise_sd=NULL, n_samples=100, scenario=NULL, s=2, 
                                model_names=c("bgp", "bgpllm", "btgp", "btgpllm"), train=NULL, test=NULL) {
  # model_names=c("blm", "btlm", "bcart", "bgp", "bgpllm", "btgp", "btgpllm")) {
  stopifnot(!is.null(scenario))
  # stop if we do not know noise.sd and there are no explicit data for either train or test
  stopifnot(xor(is.null(noise.sd), is.null(train) || is.null(test)))
  if (is.null(train)) train <- GetSimulationData(n_samples, scenario = scenario, sd=noise_sd)
  if (is.null(test))  test <- GetSimulationData(n_samples, scenario = scenario, sd=noise_sd)
  
  A_grid <- seq(0, 1, length.out = min(n_samples, 80)) 
  X <- with(train, data.frame(C=covariates, A=treatment))
  Y <- train$reward
  ZZ <- expand.grid(seq(0,1,length.out = n_samples), A_grid)  
  
  # fit selected models - may take time
  models <- lapply(model_names, function(f_name) do.call(f_name, list(X, Y, ZZ)))
  names(models) <- model_names
  
  predictions <- list()
  best_Q <- list()
  for (m in models) {
    res_dt <- GetBestPredictions(m, s=s)
    predictions[[1]]  <- res_dt$C
    predictions[[length(predictions) + 1]]  <- res_dt$A_pred
    best_Q[[length(best_Q) + 1]] <- GetValueOfPredictedA(res_dt, train)
  }
  dt <- as.data.table(predictions)
  formatted_names <- paste(model_names, paste(", Q =",  round(unlist(best_Q), 5)), sep="")
  names(dt) <- c("C", formatted_names)
  dt_melted <- melt(dt, id.vars = "C", variable.name = "Algo",  value.name = "est_A")
  gg <- ggplot(dt_melted, aes(C, est_A)) + geom_point() + geom_smooth() + facet_wrap(~ Algo, nrow = 2)
  print(gg)
  gg <- ggplot(dt_melted, aes(C, est_A, col=Algo)) + geom_smooth()
  print(gg)
  dt_melted[, A_opt := test$GetOptimalTreatment(C)]
  gg <- ggplot(dt_melted, aes(x=C, y=est_A)) + geom_point() + geom_smooth() + 
    geom_line(aes(C, A_opt, col="red")) + facet_wrap(~ Algo, nrow = 2)
  print(gg)
  
  for (m_name in names(models)){
    plot(models[[m_name]], main=paste(m_name, ", s = ", s))
  }
  PlotDecisionSurface(models, s = s)
}


ChangeFormatFromOurToChenEnriched <- function(our_data) {
  with(our_data, list(X=covariates, A=treatment, R=raw.reward, D_opt=optimal.treatment, 
                      mu=GetQFunctionValues(covariates, treatment, optimal.treatment),
                      GetQFunctionValues=GetQFunctionValues, 
                      GetOptimalTreatment=GetOptimalTreatment))
}



# EnrhichChenDataWithQfunction <- function(chen_data) {
#   stopifnot(!is.null(chen_data$GetQFunctionValues))
#   stopifnot(!is.null(chen_data$GetOptimalTreatment))
#   
#   if (!is.null(chen_data$GetQFunctionValues)) {
#     chen_data$GetQFunctionValues <- chen_data$GetQFunctionValues
#   } else {
#     warning("This function works only for scenario.4! Don't use it for anything else!!!")
#     chen_data$GetQFunctionValues <- function(X, A, A_opt=chen_data$D_opt){
#       stopifnot(is.matrix(X))
#       8 + 4*cos(2*pi*X[,2]) - 2*X[,4] - 8*X[,5]^3 - 15*abs(A_opt-A)
#     }
#   }
#   if (!is.null(chen_data$GetOptimalTreatment)) {
#     chen_data$GetOptimalTreatment <- chen_data$GetOptimalTreatment
#   } else {
#     warning("This function works only for scenario.4! Don't use it for anything else!!!")
#     chen_data$GetOptimalTreatment <- function(X) {
#       stopifnot(is.matrix(X))
#       I(X[,1] > -0.5)*I(X[,1] < 0.5)*0.6 + 1.2*I(X[,1] > 0.5) + 1.2*I(X[,1] < -0.5) + 
#         X[,4]^2 + 0.5*log(abs(X[,7])+1) - 0.6
#     }
#   }
#   return (chen_data)
# }


ChangeFormatFromChenEnrichedToOur <- function(chen_data, eps=0.01) {
  with(chen_data, list(covariates=X, treatment=A, raw.reward=R, 
                       optimal.treatment=D_opt, reward=R - min(R) + eps, 
                       GetQFunctionValues=GetQFunctionValues,
                       GetOptimalTreatment=GetOptimalTreatment))
}


# slightly rewritten function from KO-learning pred_s2 to achieve flexibility
pred_s4 <- function(model,test) {
  A_pred <- pmin(pmax(predict(model,test$X), 0), 2)
  pred_value <- with(test, mean(GetQFunctionValues(X, A_pred, D_opt)))
  return(list(A_pred=A_pred, Q=pred_value))
}


GetKOLearningValueAndPredictedDose <- function(train, test, q = 0.6) {
  train$weight = with(train, R - min(quantile(R, q),0))
  index = with(train, which(R > quantile(R,q)))
  model = with(train, svm(x = X[index,], y = A[index], w=weight[index],   
                          type="eps-regression", epsilon = 0.15, scale=FALSE))
  return(pred_s4(model,test))
}


GetGPValueAndPredictedDose <- function(train, test, model_name=NULL, s=2, eps=0.1) {
  stopifnot(is.character(model_name))
  n_samples <- length(train$reward)
  X <- with(train, data.frame(C=covariates, A=treatment))
  Y <- train$reward
  granularity <- min(n_samples, 80)
  A_grid <- with(train, seq(min(treatment)-eps, max(treatment)+eps, length.out = granularity))
  ZZ <-  data.table(C = test$covariates)
  A_grid  <- data.table(rep(A_grid, nrow(ZZ)))
  ZZ <- ZZ[rep(seq.int(1, nrow(ZZ)), each=granularity), ]
  ZZ[, A:=A_grid]
  model <- do.call(model_name, list(X,Y,ZZ)) 
  res_dt <- GetBestPredictions(model, s=s)
  col_names <- c(grep('^C', names(res_dt), value = T))
  stopifnot(length(col_names) > 0)
  C_matrix <-  as.matrix(res_dt[, col_names, with=F]) 
  pred_value <- mean(test$GetQFunctionValues(C_matrix, res_dt$A_pred))
  # why doesn't this work?
  # pred_value <- with(test, res_dt[, mean(GetQFunctionValues(as.matrix(.SD), A_pred)), .SDcols=col_names])
  return(list(A_pred=res_dt$A_pred, Q=pred_value))
}













# Why this is not working ? 
# res_dt[, .(.SD, A_pred),  .SDcols=col_names ]
# res_dt[, f(.SD, A_pred),  .SDcols=col_names ]
# Why this doesn't return a matrix ? 
# res_dt[, as.matrix(col_names), with=F]











