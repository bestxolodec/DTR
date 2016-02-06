source("./code/functions.R")

# for accidental sourcing the whole R file
Sys.sleep(10)




# Prepare Data ------------------------------------------------------------

data <- GetDataToWorkWith()
# propensity.scores <- GetPropensityScores(observations.with.metadata=data, two.step = T)
propensity.scores <- rep(1, nrow(data))
reward.values  <- GetReward(observations.with.metadata=data)
treatment.values <- data$Therapeutic.Dose.of.Warfarin
exclude.from.model.colnames  <- c("PharmGKB.Subject.ID",
                                  "INR.on.Reported.Therapeutic.Dose.of.Warfarin",
                                  "Therapeutic.Dose.of.Warfarin")
# remove colums that are not necessary in learning process rescale and add intercept
patient.covariates.with.factors <- data[, -which(names(data) %in% exclude.from.model.colnames)]
patient.covariates <- model.matrix( ~ . -1 , data = patient.covariates.with.factors)
patient.covariates <- reshape::rescaler(patient.covariates, type = "range")
patient.covariates <- model.matrix( ~ ., data = as.data.frame(patient.covariates))

number.of.params <- ncol(patient.covariates)




# Optimization of linear model --------------------------------------------
lower.params.threshold <- rep(-1000, number.of.params)
upper.params.threshold <- rep(1000, number.of.params)
offset <- 2
lambda <- 0.000001
lambda <- 0
regression.model <- lm(treatment.values ~ patient.covariates - 1, model=T, x=T)
initial.params <-  regression.model$coefficients # TODO: Find out what to do with NA after regression
initial.params <- replace(initial.params, is.na(initial.params), 0)
# initial.params <-  runif(length(colnames(patient.covariates)))

optimized <- GenSA(par = initial.params, fn = ObjectiveFunction, lower.params.threshold,
      upper.params.threshold,
      control=list(smooth=T, verbose=TRUE, maxit=6000, temperature = 6230),
      # additional arguments which goes directly to ObjectiveFunction
      patient.covariates, propensity.scores, offset, PolicyFunLinearKernel,
      lambda, reward.values, treatment.values)
optimized[c("value",  "par", "counts")]
plot(optimized$trace.mat[, "function.value"])
plot(optimized$trace.mat[, "current.minimum"])
PolicyFunLinearKernel(optimized$par, patient.covariates =  patient.covariates)
ObjectiveFunction(optimized$par, patient.covariates, propensity.scores, offset,
                  PolicyFunLinearKernel, lambda, reward.values, treatment.values)
PolicyFunLinearKernel(initial.params, patient.covariates =  patient.covariates)
ObjectiveFunction(initial.params, patient.covariates, propensity.scores, offset,
                  PolicyFunLinearKernel, lambda, reward.values, treatment.values)

PolicyFunLinearKernel(params, observations.with.metadata = patient.covariates)


decision.values <- PolicyFunLinearKernel(optimized$par, patient.covariates =  patient.covariates)
rewards.scaled.0.1 <- (reward.values - min(reward.values) ) / (max(reward.values) - min(reward.values)) 
plot(decision.values, treatment.values, 
     col=rgb(1 - rewards.scaled.0.1, rewards.scaled.0.1, 0),
     pch=19)
abline(0, 1)

plot(density(decision.values), col="green", lwd=4)
lines(density(treatment.values), col="blue", lwd=4)





# Optimization for gauss kernel -------------------------------------------

data <- GetDataToWorkWith()
propensity.scores <- GetPropensityScores(observations.with.metadata=data, two.step = T)
# propensity.scores <- rep(1, nrow(data))
reward.values  <- GetReward(observations.with.metadata=data)
treatment.values <- data$Therapeutic.Dose.of.Warfarin
exclude.from.model.colnames  <- c("PharmGKB.Subject.ID",
                                  "INR.on.Reported.Therapeutic.Dose.of.Warfarin",
                                  "Therapeutic.Dose.of.Warfarin")
patient.covariates.with.factors <- data[, -which(names(data) %in% exclude.from.model.colnames)]
patient.covariates <- model.matrix( ~ . , data = patient.covariates.with.factors)
patient.covariates <- reshape::rescaler(patient.covariates, type = "range")

number.of.params <- nrow(patient.covariates)
lower.params.threshold <- rep(-100000, number.of.params)
upper.params.threshold <- rep(100000, number.of.params)
offset <- 1
lambda <- 0.1
gamma <- 1
optimized <- GenSA(par = NULL, fn = ObjectiveFunction, lower.params.threshold,
      upper.params.threshold,
      control=list(smooth=TRUE, verbose=TRUE),
      # additional arguments which goes directly to ObjectiveFunction
      patient.covariates, propensity.scores, offset, PolicyFunGaussKernel,
      lambda, reward.values, treatment.values, hyperparams=list(gamma=gamma))
optimized[c("value",  "par", "counts")]
plot(optimized$trace.mat[, "function.value"])
plot(optimized$trace.mat[, "current.minimum"])
PolicyFunLinearKernel(optimized$par, patient.covariates =  patient.covariates)
ObjectiveFunction(optimized$par,
      patient.covariates, propensity.scores, offset, PolicyFunGaussKernel,
      lambda, reward.values, treatment.values, hyperparams=list(gamma=gamma))

PolicyFunGaussKernel(optimized$par, patient.covariates =  patient.covariates,
                     hyperparams = list(gamma=gamma))
ObjectiveFunction(initial.params,
      patient.covariates, propensity.scores, offset, PolicyFunGaussKernel,
      lambda, reward.values, treatment.values, hyperparams=list(gamma=gamma))

# Empirical Value Function  -----------------------------------------------

# NOT YET READY
ValueFunction <- function(data, offset, propensity.scores, policy.function,
                           observed.treatment.colname = "Therapeutic.dose.of.warfarin") {
  loss.value <- 1 - (abs(data[, observed.treatment.colname]) - policy.function(data)) / offset
  multiplier <- max(loss.value, 0)
  return(mean(get.reward(data) / propensity.scores / offset  * multiplier))
}

# Split data on test and train  -------------------------------------------

NUMBER.OF.SPLITS = 10
partition <- createDataPartition(data$Therapeutic.Dose.of.Warfarin,
                                 times = NUMBER.OF.SPLITS, p=0.46)

all.indexies = 1:NROW(data$Therapeutic.Dose.of.Warfarin)
for (train.indexies in partition) {
  test.indexies = setdiff(all.indexies,   train.indexies)
  test.data  <- sapply(data[test.indexies, ], as.numeric)
  train.data <- sapply(data[train.indexies, ], as.numeric)
}


test.indexies = setdiff(all.indexies,   train.indexies)
test.data  <- sapply(data[test.indexies, ], as.numeric)
test.prop.scores <- GetPropensityScores(train.data, formula = balance.formula )
train.data <- sapply(data[train.indexies, ], as.numeric)
train.prop.scores <- GetPropensityScores(train.data, formula = balance.formula )

value.function(data=train.data,
               offset = 0.001,
               propensity.scores = train.prop.scores,
               policy.function =  policy.function,
               observed.outcome.colname = "Therapeutic.Dose.of.Warfarin")









cbind(reward.values, treatment.values, PolicyFunLinearKernel(initial.params, patient.covariates = patient.covariates))













# Test with toy example  --------------------------------------------------

patient.covariates.test <- data.frame(c(1,1), c(2,10), c(3,6)) 
treatment.values.test <- c(20,35)
reward.values.test <- c(-0.5, -5)
propensity.scores.test <- rep(1, nrow(patient.covariates.test)) 

number.of.params <- ncol(patient.covariates.test)
lower.params.threshold <- rep(-100, number.of.params)
upper.params.threshold <- rep(100, number.of.params)
offset <- 5
lambda <- 0
initial.params.test <- c(1, 3, 1)
optimized <- GenSA(par = initial.params.test, fn = ObjectiveFunction, lower.params.threshold,
      upper.params.threshold,
      control=list(smooth=TRUE, verbose=TRUE),
      # additional arguments which goes directly to ObjectiveFunction
      patient.covariates.test, propensity.scores.test, offset, PolicyFunLinearKernel,
      lambda, reward.values.test, treatment.values.test, hyperparams=list())
optimized[c("value",  "par", "counts")]
plot(optimized$trace.mat[, "function.value"])
plot(optimized$trace.mat[, "current.minimum"])
PolicyFunLinearKernel(optimized$par, patient.covariates =  patient.covariates.test)
ObjectiveFunction(optimized$par,
      patient.covariates.test, propensity.scores.test, offset, PolicyFunLinearKernel,
      lambda, reward.values.test, treatment.values.test, hyperparams=list())


ObjectiveFunction(initial.params.test,
      patient.covariates.test, propensity.scores.test, offset, PolicyFunLinearKernel,
      lambda, reward.values.test, treatment.values.test, hyperparams=list())
ObjectiveFunction(c(0,3,1),
      patient.covariates.test, propensity.scores.test, offset, PolicyFunLinearKernel,
      lambda, reward.values.test, treatment.values.test, hyperparams=list())

PolicyFunLinearKernel(c(0,10,0), patient.covariates =  patient.covariates.test)
ObjectiveFunction(c(0,10,0),
      patient.covariates.test, propensity.scores.test, offset, PolicyFunLinearKernel,
      lambda, reward.values.test, treatment.values.test, hyperparams=list())





plot(PolicyFunLinearKernel(optimized$par, patient.covariates =  patient.covariates), 
     data$Therapeutic.Dose.of.Warfarin, 
     col= rgb(1 -  (reward.values - min(reward.values) ) / (max(reward.values) - min(reward.values)), 
              (reward.values - min(reward.values) ) / (max(reward.values) - min(reward.values)), 
              0), 
     pch=19)

