library(doParallel)
library(caret)
source("./src/functions.R")
# for accidental sourcisource("./src/functions.R")
Sys.sleep(10)


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

GetRewardGivenQfunctionValuesAsMeanVec <- function(q.function.values,
                                                   variance=1, eps=0.01) {
  norm.sample <- rnorm(NROW(q.function.values), q.function.values, variance)
  return(as.matrix(norm.sample - min(norm.sample) + eps))
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
  reward <- GetRewardGivenQfunctionValuesAsMeanVec(q.function.values)
  prop.scores <- as.matrix(rep(1, nrow(covariates)))
  # train.prop.scores <- GetPropensityScores(
  #     as.data.frame(cbind(train.treatment, train.covariates)),
  #     balance.formula = formula (V1 ~ .), two.step = T)
  return(list(covariates=covariates, treatment=treatment, reward=reward,
              prop.scores=prop.scores))
}

GetMeanDtrValueOnFolds <- function(folds, treatment, covariates, prop.scores,
     reward, offset, control.offset, policy.function, lambda, obj.func) {
  dtrs <- numeric(length = 0)
  for (control.fold in folds) {
    tune.covariates  <-  covariates [-control.fold, ]
    tune.treatment   <-  treatment  [-control.fold, ]
    tune.reward      <-  reward     [-control.fold, ]
    tune.prop.scores <-  prop.scores[-control.fold, ]
    opt.result <- OptimizeParamsOfPolicyFunction(treatment = tune.treatment,
        covariates = tune.covariates, prop.scores = tune.prop.scores,
        reward = tune.reward, offset = offset,
        policy.function = PolicyFunLinearKernel, lambda = lambda,
        obj.func = ObjectiveFunction, regress.init.params=T)
    control.covariates  <-  covariates [control.fold, ]
    control.treatment   <-  treatment  [control.fold, ]
    control.reward      <-  reward     [control.fold, ]
    control.prop.scores <-  prop.scores[control.fold, ]
    dtr.value <- ValueFunction(opt.result$par, control.treatment, control.covariates,
        control.prop.scores,  control.reward, control.offset, PolicyFunLinearKernel)
    dtrs <- c(dtrs, dtr.value)
    cat("DTRs: ", dtrs)
  }
  return(mean(dtrs))
}




# Optimization routines ---------------------------------------------------

# covariates should already contain or not intercept
OptimizeParamsOfPolicyFunction <- function(treatment, covariates,
    prop.scores, reward, offset, policy.function, lambda, obj.func,
    regress.init.params=T) {
  number.of.params <- ncol(covariates)
  lower.params.threshold <- rep(-1000, number.of.params)
  upper.params.threshold <- rep(1000, number.of.params)
  if (isTRUE(regress.init.params)) {
    regression.model <- lm(treatment ~ covariates - 1, model=T, x=T)
    initial.params <-  regression.model$coefficients # TODO: Find out what to do with NA after regression
    initial.params <- replace(initial.params, is.na(initial.params), 0)
  } else {
    initial.params <- NULL
  }
  optimized <- GenSA(par = initial.params, fn = obj.func,
      lower.params.threshold, upper.params.threshold,
      control=list(smooth=T, verbose=TRUE, maxit=6000, temperature = 6200),
      # additional arguments which goes directly to ObjectiveFunction
      treatment, covariates, prop.scores, reward, offset, policy.function, lambda)
  return(optimized)
}




# Defining control execution constants ------------------------------------

registerDoParallel(cores = 4)
train.data.sample.sizes <- c(50, 100, 200, 400, 800)
test.data.sample.size <- 10000
number.of.covariates = 30
sample.size <-  100
offsets = seq(0.01, 1, length = 10)
control.offset = min(offsets) / 2 
lambdas = seq(1, 10,  length=5) * 0.01


# Prepare data  -----------------------------------------------------------

train.data.list <- GetSimulationData(sample.size, number.of.covariates)
train.covariates  = as.matrix(train.data.list$covariates)
train.treatment   = as.matrix(train.data.list$treatment)
train.reward      = as.matrix(train.data.list$reward)
train.prop.scores = as.matrix(train.data.list$prop.scores)

test.data.list <- GetSimulationData(test.data.sample.size, number.of.covariates)
test.covariates  = test.data.list$covariates
test.treatment   = test.data.list$treatment
test.reward      = test.data.list$reward
test.prop.scores = test.data.list$prop.scores


# Do grid search ----------------------------------------------------------

global.result <- list()

for (offset in  offsets) {
  print(paste("Optimization offset: ", offset))
  result <- foreach (lambda = lambdas) %dopar% {
    folds <- createFolds(train.treatment)
    mean.dtr.value <- GetMeanDtrValueOnFolds(folds, train.treatment, 
        train.covariates, train.prop.scores, train.reward, offset, 
        control.offset, PolicyFunLinearKernel, lambda, obj.func)
    return (list("offset"=offset, "lambda"=lambda,
                 "value.function"=mean.dtr.value))
  }
  global.result <- c(global.result, result)
}


# dtr.value.on.test <- ValueFunction(opt.result$par, test.treatment, 
#     test.covariates, test.prop.scores,  test.reward, 
#     control.offset, PolicyFunLinearKernel)


global.result.df <- as.data.frame(t(sapply(global.result, cbind)))
colnames(global.result.df)  <- names(global.result[[1]])
result.to.plot <- global.result.df[, names(global.result.df) %in%
                                     c("offset", "lambda", "value.function")]
result.to.plot <- as.data.frame(lapply(result.to.plot, as.numeric))
result.to.plot$lambda  <- as.factor(result.to.plot$lambda)
ggplot(result.to.plot, aes(x = offset, y = value.function, colour = lambda)) + geom_line()

pars = global.result.df[ which.max(global.result.df$value.function),  ]
best.offset = unlist(pars$offset)
best.lambda = unlist(pars$lambda)
opt.result <- OptimizeParamsOfPolicyFunction(train.treatment,
    train.covariates, train.prop.scores, train.reward, best.offset,
    PolicyFunLinearKernel, best.lambda, ObjectiveFunction)
dtr.value.on.test <- ValueFunction(opt.result$par, test.treatment, 
    test.covariates,  test.prop.scores,  test.reward, 
    best.offset, PolicyFunLinearKernel)
decision.values <- PolicyFunLinearKernel(opt.result$par, test.covariates)
rewards.scaled.0.1 <- (test.reward - min(test.reward) ) / (max(test.reward) - min(test.reward))
plot(decision.values, test.treatment,
     col=rgb(1 - rewards.scaled.0.1, rewards.scaled.0.1, 0),
     pch=20)
abline(0, 1)

plot(density(decision.values), col="green", lwd=4,  ylim=c(0, 0.5))
lines(density(test.treatment), col="blue", lwd=4)
