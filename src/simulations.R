source("./src/functions.R")
# for accidental sourcisource("./src/functions.R")
Sys.sleep(10)


# Optimal policy functions ------------------------------------------------

GetOptimalDecisionsForFirstSimulation <- function(covariates) {
  return(as.matrix(1 + 0.5 * covariates[,1] + 0.5 * covariates[,2]))
}


# Q functions  ------------------------------------------------------------

GetQfunctionValueForFirstSimulation <- function(covariates, given.treatement, optimal.treatment) {
  linear.part <-  8 + 4 * covariates[, 1] - 2 * covariates[, 2] - 2 * covariates[, 3]
  non.linear.part <- - 25 * (optimal.treatment - given.treatement) ** 2
  return (linear.part + non.linear.part)
}


# RewardFunctions ---------------------------------------------------------

GetRewardGivenQfunctionValuesAsMeanVec <- function(q.function.values,
                                                   variance=1, eps=0.01) {
  norm.sample <- rnorm(NROW(q.function.values), q.function.values, variance)
  return(norm.sample - min(norm.sample) + eps)
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
      control=list(smooth=T, verbose=TRUE, maxit=6000, temperature = 6230),
      # additional arguments which goes directly to ObjectiveFunction
      treatment, covariates, prop.scores, reward, offset, policy.function, lambda)
}



# Defining control execution constants ------------------------------------

train.data.sample.sizes <- c(50, 100, 200, 400, 800)
test.data.sample.size <- 10000
number.of.covariates = 30
sample.size <-  50
offsets = seq(0.0001, 1, length = 2)
lambdas = seq(1, 10,  length=2) * 0.001



# Prepare data  -----------------------------------------------------------

train.treatment <- matrix(runif(sample.size, min=0, max=2), ncol=1)
train.covariates <- matrix(runif(sample.size * number.of.covariates,
                                 min=-1, max=1),
                           ncol=number.of.covariates)
train.optimal.treatment <- GetOptimalDecisionsForFirstSimulation(train.covariates)
q.function.values <- GetQfunctionValueForFirstSimulation(
    covariates = train.covariates,
    given.treatement =  train.treatment,
    optimal.treatment = train.optimal.treatment)
train.reward <- GetRewardGivenQfunctionValuesAsMeanVec(q.function.values)

train.prop.scores <- GetPropensityScores(
    as.data.frame(cbind(train.treatment, train.covariates)),
    balance.formula = formula (V1 ~ .), two.step = T)
train.prop.scores <- rep(1, nrow(train.covariates))




# Do grid search ----------------------------------------------------------

library(doParallel)
registerDoParallel(cores = 4)

global.result <- list()
                 
for (offset in  offsets) {
  result <- foreach (lambda = lambdas, .combine = rbind) %dopar% {
    opt.result  <- OptimizeParamsOfPolicyFunction(treatment = train.treatment, 
        covariates = train.covariates, prop.scores = train.prop.scores, 
        reward = train.reward, offset = offset, 
        policy.function = PolicyFunLinearKernel, lambda = lambda, 
        obj.func = ObjectiveFunction, regress.init.params=T)
    dtr.value <- ValueFunction(opt.result$par, train.treatment, train.covariates,
        train.prop.scores,  train.reward, offset, PolicyFunLinearKernel)
    return (list("offset"=offset, "lambda"=lambda, "params"=opt.result$par,
            "value.function"=dtr.value))
  }
  global.result <- c(global.result, result)
}


res <- foreach(i = 1:10, .combine = rbind)  %dopar% {
  PolicyFunGaussKernel(opt.result$par, covariates)
 return (list("i" = i, "somestuff"=list("value"=i+10)))
}


decision.values <- PolicyFunLinearKernel(opt.result$par, train.covariates)
rewards.scaled.0.1 <- (train.reward - min(train.reward) ) / (max(train.reward) - min(train.reward))
plot(decision.values, train.treatment,
     col=rgb(1 - rewards.scaled.0.1, rewards.scaled.0.1, 0),
     pch=19)
abline(0, 1)





library(ggplot2)
library(reshape)
data <- data.frame(time = seq(0, 23), noob = rnorm(24), plus = runif(24),
                   extra = rpois(24, lambda = 1))
Molten <- melt(data, id.vars = "time")
ggplot(Molten, aes(x = time, y = value, colour = variable)) + geom_line()
