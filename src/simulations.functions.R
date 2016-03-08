source("./functions.R")
library(doParallel)
library(caret)


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
  prop.scores <- GetPropensityScores(
      data.frame("treatment"=treatment, "covars"=covariates),
      balance.formula = formula (treatment ~ . - covars..Intercept. - 1),
      two.step = T)
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
                                           opt.hyperparams=list()) {
  regression.model <- lm(obs.data$treatment ~ obs.data$covariates - 1, model=T, x=T)
  # TODO: Find out what to do with NA after regression
  initial.params <- regression.model$coefficients
  initial.params <- replace(initial.params, is.na(initial.params), 0)
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
    params <- opt.hyperparams$opt.func(params=initial.params, obs.data,  
                                       offset, policy.function, lambda)
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

PlotDecsionsVersusObserved <- function(obs.data, policy.function, params, 
                                       title="DCA decisions versus observed"){
  decision.values <- policy.function(params, obs.data$covariates)
  rewards.scaled.0.1 <- (obs.data$reward - min(obs.data$reward) ) / 
                        (max(obs.data$reward) - min(obs.data$reward))
  plot(decision.values, obs.data$treatment,
       col=rgb(1 - rewards.scaled.0.1, rewards.scaled.0.1, 0),
       pch=20, main=title)
  abline(0, 1)  
}














# # Do grid search ----------------------------------------------------------
#
# dtrs.mean.on.test <-  list()
# for (i in seq(1, number.of.simulations)) {
#   global.result <- list()
#   for (offset in  offsets) {
#     print(paste("Optimization offset: ", offset))
#     result <- foreach (lambda = lambdas) %dopar% {
#       folds <- createFolds(train.treatment)
#       mean.dtr.value.gen.sa <- GetDtrValuesOnFolds(folds, train.treatment,
#           train.covariates, train.prop.scores, train.reward, offset,
#           control.offset, PolicyFunLinearKernel, lambda,
#           list("obj.func"=ObjectiveFunction))
#       mean.dtr.value.dca <- GetDtrValuesOnFolds(folds, train.treatment,
#           train.covariates, train.prop.scores, train.reward, offset,
#           control.offset, PolicyFunLinearKernel, lambda)
#       return (list("offset"=offset, "lambda"=lambda,
#                    "value.function.gen.sa"=mean.dtr.value.gen.sa,
#                    "value.function.dca"=mean.dtr.value.dca, ))
#     }
#     global.result <- c(global.result, result)
#   }
#   global.result.df <- as.data.frame(t(sapply(global.result, cbind)))
#   colnames(global.result.df)  <- names(global.result[[1]])
#   pars <-  global.result.df[ which.max(global.result.df$value.function),  ]
#   best.dtr.value.on.train <- max(unlist(global.result.df$value.function))
#   best.offset  <- unlist(pars$offset)
#   best.lambda  <- unlist(pars$lambda)
#
#   opt.result <- OptimizeParamsOfPolicyFunction(train.treatment,
#       train.covariates, train.prop.scores, train.reward, best.offset,
#       PolicyFunLinearKernel, best.lambda, ObjectiveFunction)
#   dtr.value.on.test <- ValueFunction(opt.result$par, test.treatment,
#       test.covariates,  test.prop.scores,  test.reward,
#       best.offset, PolicyFunLinearKernel)
#
#   dtrs.mean.on.test <- c(dtrs.mean.on.test,
#      list(list("dtr.value.on.test"=dtr.value.on.test,
#                "best.offset"=best.offset, "best.lambda"=best.lambda,
#                "best.dtr.value.on.train"=best.dtr.value.on.train,
#                "opt.params"=opt.result$par)))
# }
#
#
#
#
#
#
#
#
# # Single shot comparison of different optimization techniques -------------
#
# global.result <- list()
# for (offset in  offsets) {
#   print(paste("Optimization offset: ", offset))
#   result <- foreach (lambda = lambdas) %dopar% {
#     folds <- createFolds(train.treatment)
#     mean.dtr.value.gen.sa <- GetDtrValuesOnFolds(folds, train.treatment,
#         train.covariates, train.prop.scores, train.reward, offset,
#         control.offset, PolicyFunLinearKernel, lambda,
#         list("obj.func"=ObjectiveFunction))
#     mean.dtr.value.dca <- GetDtrValuesOnFolds(folds, train.treatment,
#         train.covariates, train.prop.scores, train.reward, offset,
#         control.offset, PolicyFunLinearKernel, lambda)
#     mean.dtr.value.dca.l2 <- GetDtrValuesOnFolds(folds, train.treatment,
#         train.covariates, train.prop.scores, train.reward, offset,
#         control.offset, PolicyFunLinearKernel, lambda,
#         list(""))
#     return (list("offset"=offset, "lambda"=lambda,
#                  "value.function.gen.sa"=mean.dtr.value.gen.sa,
#                  "value.function.dca"=mean.dtr.value.dca,
#                  "value.function.dca.l2"=mean.dtr.value.dca.l2))
#   }
#   global.result <- c(global.result, result)
# }
#
# global.result.df <- as.data.frame(t(sapply(global.result, cbind)))
# colnames(global.result.df)  <- names(global.result[[1]])
# result.to.plot <- global.result.df[, names(global.result.df) %in%
#                                      c("offset", "lambda",
#                                        "value.function.dca",
#                                        "value.function.gen.sa")]
# result.to.plot <- as.data.frame(lapply(result.to.plot, as.numeric))
# result.to.plot$lambda  <- as.factor(result.to.plot$lambda)
#
# par(mfrow=c(2,2))
# par(mfrow=c(1,2))
#
# # DCA plotting ------------------------------------------------------------
#
# best.dtr.value.on.train.dca <- max(unlist(global.result.df$value.function.dca))
# pars.dca <-  global.result.df[ which.max(global.result.df$value.function.dca),  ]
# best.offset.dca <- unlist(pars.dca$offset)
# best.lambda.dca <- unlist(pars.dca$lambda)
# opt.params.dca <- OptimizeParamsOfPolicyFunction(train.treatment,
#     train.covariates, train.prop.scores, train.reward, best.offset.dca,
#     PolicyFunLinearKernel, best.lambda.dca)
# dtr.value.on.test.dca <- ValueFunction(opt.params.dca, test.treatment,
#     test.covariates,  test.prop.scores,  test.reward,
#     best.offset.dca, PolicyFunLinearKernel)
#
#
# # Plot DCA tuning results  ------------------------------------------------
#
# ggplot(result.to.plot,
#        aes(x = offset, y = value.function.dca, colour = lambda)) +
#   geom_line() +
#   ggtitle("Cross validation results for Difference of Convex functions optimization")
#
# # Plot DCA treatment assignment comparison --------------------------------
#
# decision.values <- PolicyFunLinearKernel(opt.params.dca, test.covariates)
# rewards.scaled.0.1 <- (test.reward - min(test.reward) ) / (max(test.reward) - min(test.reward))
# plot(decision.values, test.treatment,
#      col=rgb(1 - rewards.scaled.0.1, rewards.scaled.0.1, 0),
#      pch=20, main="DCA decisions versus observed")
# abline(0, 1)
#
#
# # Plot DCA treatment assignment density -----------------------------------
#
# plot(density(decision.values), col="green", lwd=4,  ylim=c(0, 0.8),
#      main = "Treatment density")
# lines(density(test.treatment), col="blue", lwd=4)
# legend("topleft", c("Recommended by DTR", "Observed"),
#        lwd=c(2.5,2.5),col=c("green","blue"), bty = "n")
#
#
#
# # Simulated Annealing plotting --------------------------------------------
#
# best.dtr.value.on.train.gen.sa <- max(unlist(global.result.df$value.function.gen.sa))
# pars.gen.sa <-  global.result.df[ which.max(global.result.df$value.function.gen.sa),  ]
# best.offset.gen.sa  <- unlist(pars.gen.sa$offset)
# best.lambda.gen.sa  <- unlist(pars.gen.sa$lambda)
# opt.params.gen.sa <- OptimizeParamsOfPolicyFunction(train.treatment,
#     train.covariates, train.prop.scores, train.reward, best.offset.gen.sa,
#     PolicyFunLinearKernel, best.lambda.gen.sa, list("obj.func"=ObjectiveFunction))
# dtr.value.on.test.gen.sa <- ValueFunction(opt.params.gen.sa, test.treatment,
#     test.covariates,  test.prop.scores,  test.reward,
#     best.offset.gen.sa, PolicyFunLinearKernel)
#
#
#
# # Plot tuning results  ----------------------------------------------------
#
# ggplot(result.to.plot,
#        aes(x = offset, y = value.function.gen.sa, colour = lambda)) +
#   geom_line() +
#   ggtitle("Cross validation results for Simulated Annealing optimization")
#
#
# # Plot treatment assignment comparison ------------------------------------
#
# decision.values <- PolicyFunLinearKernel(opt.params.gen.sa, test.covariates)
# rewards.scaled.0.1 <- (test.reward - min(test.reward) ) /
#   (max(test.reward) - min(test.reward))
# plot(decision.values, test.treatment,
#      col=rgb(1 - rewards.scaled.0.1, rewards.scaled.0.1, 0),
#      pch=20,  main="Sim. Anneal decisions versus observed")
# abline(0, 1)
#
#
# # Plot treatment assignment density ---------------------------------------
#
# plot(density(decision.values), col="green", lwd=4,  ylim=c(0, 0.9),
#      main = "Treatment density")
# lines(density(test.treatment), col="blue", lwd=4)
# legend("topleft", c("Recommended by DTR", "Observed"),
#        lwd=c(2.5,2.5),col=c("green","blue"), bty = "n")
#
#
# hist(decision.values)
# hist(test.reward)
#
#
# Sys.sleep(100)
# save.image("./src/Workspace.Rdata")
