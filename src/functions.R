library(caret)
library(GenSA)
library(CBPS)
library(reshape)
library(quantreg)


# Functions to work with data ---------------------------------------------

GetDataToWorkWith <- function(file.path, choosen.columns.names, columns.to.float) {
  # Read, clean, and prepare data to work with
  #
  # Args:
  #   file.path: Path to csv file with data
  #   choosen.columns.names: Which columns whould be used as a subset of the data
  #   columns.to.float: List of colnames to convert from float representation with
  #                     ',' as a separator
  #
  # Returns:
  #   The data frame with fully available observations
  if (missing(file.path)) {
    file.path <- "./data/ipwc_data.csv"
  }
  if (missing(choosen.columns.names)) {
    choosen.columns.names <- c("PharmGKB.Subject.ID", "Age", "Height..cm.",
                               "Weight..kg.", "Race..Reported.", "VKORC1.497.consensus",
                               "Cyp2C9.genotypes", "Amiodarone..Cordarone.",
                               "INR.on.Reported.Therapeutic.Dose.of.Warfarin",
                               "Therapeutic.Dose.of.Warfarin",
                               "VKORC1.genotype...1639.G.A..3673...chr16.31015190..rs9923231..C.T",
                               "Macrolide.Antibiotics")
  }
  if (missing(columns.to.float)) {
    columns.to.float  <- c("Height..cm.", "Therapeutic.Dose.of.Warfarin",
                           "INR.on.Reported.Therapeutic.Dose.of.Warfarin", "Weight..kg.")
  }
  iwpc <- read.csv(file.path, na.strings = c("NA", ""))
  complete.cases.indexies <- complete.cases(iwpc[, choosen.columns.names])
  data <- iwpc[complete.cases.indexies, choosen.columns.names]
  # column names with "," as a separator
  for (column in columns.to.float) {
    data[, column] <-  as.numeric(gsub(",",".", as.character(data[, column])))
  }
  return (as.data.frame(droplevels(data)))
}


# Get propensity scores  --------------------------------------------------

GetPropensityScores <- function(observations.with.metadata, balance.formula,
                                two.step=F, return.cbps.object=F) {
  # Get propensity scores of row observations in observed.data with
  # supplied balance.formula.
  #
  # Args:
  #   observations.with.metadata: data for balancing
  #   balance.formula:  of the from: treatment ~ cov1 + cov2  + ...
  #
  # Returns:
  #   data.frame with columns being identifiers of rows of the observed.data and
  #   values (one row) being propensity scores of that rows
  if(missing(observations.with.metadata)) {
    stop("No observations.with.metadata data supplied to GetPropensityScores function")
  }
  if (missing(balance.formula)) {
    balance.formula <- formula(Therapeutic.Dose.of.Warfarin ~ Height..cm.
                               + Age + Weight..kg. + Race..Reported. +  VKORC1.497.consensus
                               + Cyp2C9.genotypes  + Amiodarone..Cordarone.
                               + VKORC1.genotype...1639.G.A..3673...chr16.31015190..rs9923231..C.T
                               + Macrolide.Antibiotics)
  }
  x.balanced <- CBPS(balance.formula, data=observations.with.metadata, twostep=two.step)
  if (isTRUE(return.cbps.object)) {
    return (x.balanced)
  } else {
    return (x.balanced$fitted.values)
  }
}


# Reward function ---------------------------------------------------------

GetReward <- function(observations.with.metadata, ideal.outcome.value = 2.5,
                      observed.outcome.colname = "INR.on.Reported.Therapeutic.Dose.of.Warfarin",
                      eps=0.01) {
  # Get reward value for each object (row) in observations.with.metadata as minus
  # the absolute deviation of observed.outcome.colname from ideal.value
  #
  # Args:
  #   observations.with.metadata: Data frame with colname observed.outcome.colname.
  #   ideal.outcome.value: Ideal value of an outcome.
  #   observed.outcome.colname: Name of the colname which contains an outcome for
  #                             paricular objects, being rows of observations.with.metadata.
  #
  # Returns:
  #   List of Reward values, with indexies being the row
  #   indexies of observations.with.metadata.
  observed.outcome <- observations.with.metadata[, observed.outcome.colname]
  minus.abs.outcome.deviance <-  - abs(observed.outcome - ideal.outcome.value)
  return (minus.abs.outcome.deviance - min(minus.abs.outcome.deviance) + eps)
  # return  (- abs(value - ideal.outcome.value) + eps)
}


# Policy function ---------------------------------------------------------

PolicyFunLinearKernel <- function(params, patient.covariates,
                                  hyperparams=list(), ...) {
  # Make linear decision for patient treatment.
  #
  # Args:
  #   params: Parameters of the linear model.
  #   patient.covariates:  Data.frame with rows being equal to patients covariates.
  #   hyperparams: list of parameters for policy function, empty for linearkernel
  #
  # Returns:
  #   List of treatment's values, with indexies being the row
  #   indexies of patient.covariates.
  return (as.matrix(patient.covariates) %*%  as.matrix(params))
}


GaussKernel <- function(X, gamma) {
  return(exp(-1 * as.matrix(dist(X)^2) / gamma))
}


PolicyFunGaussKernel <- function(params, patient.covariates,
                                 hyperparams=list(gamma=0.1), ...) {
  # Make decision for patient treatment based on Gauss kernels
  #
  # Args:
  #   params: Parameters of the linear model.
  #   patient.covariates:  Data.frame with rows being equal to patients covariates.
  #   hyperparams: list of parameters for policy function, empty for linearkernel
  #
  # Returns:
  #   List of treatment's values, with indexies being the row
  #   indexies of patient.covariates.
  return (GaussKernel(patient.covariates, hyperparams$gamma) %*% params)
}



# Objective function ------------------------------------------------------

ObjectiveFunction <- function(params, treatment, covariates, prop.scores,
                              reward, offset, policy.function, lambda,
                              hyperparams=list()) {
  # Objective (risk) function for minimization problem. Includes regularization.
  #
  # Args:
  #   params: Model parameters for policy.function -- are being changed
  #           during minimization procedure.
  #   treatment: Observed treatment values for patients.
  #   covariates: Obseravation data of patients.
  #   prop.scores: Correspond to observations in data.
  #   reward: Precomputed reward values of treatment treatment.
  #           The greater reward, the better.
  #   offset: Smallest value of difference between
  #           observational treatment value and the predicted one,
  #           which is interpeted as full "1" loss
  #   policy.function: Decision function of the form function(params, data, ...)
  #                    which returns treatment given covariates.
  #   lambda: Regularization parameter.
  #   hyperparams:  list of hyperparameters for policy function
  #
  # Returns:
  #   Regularized risk function of DTR specified with policy.function.
  prediction <- policy.function(params, covariates, hyperparams)
  multiplier <- pmin(abs(treatment - prediction) / offset, 1)
  risk.function.value <- mean(reward / prop.scores / (2 * offset) * multiplier)
  # tail is because we dont penalize intercept parameter
  regularization.value <- lambda * sum(tail(params, -1) ** 2)
  return(risk.function.value + regularization.value)
}




# Optimization with DC functions ------------------------------------------


DifferenceConvexOptimize <- function(params=NULL, treatment, covariates,
    propensity.scores, reward, offset, policy.function, lambda, 
     hyperparams=list(), tolerance = 0.00001) {
  # browser()
  stopifnot(is.matrix(covariates))
  if (is.null(params)) {
    params <- runif(ncol(covariates), min=-1, max=1)
  }
  data.with.target <- c(data.frame(treatment), as.data.frame(covariates))
  t.params <- params
  # simply for first loop iteration
  t.next.params <- t.params
  iteration <- 0
  repeat{
    cat("Iteration", iteration, "\n")
    t.params <- t.next.params
    t.prediction <- policy.function(t.params, covariates, hyperparams)
    t.abs.deviance.from.treatment  <- abs(treatment - t.prediction)
    t.Q.values  <-  ifelse(t.abs.deviance.from.treatment <= offset, 0, 1)
    t.weights  <- reward * t.Q.values / offset ** 2
    # TODO: think of how to use rq.fit.lasso() here
    # -1 as an intercept term is already present in data matirx
    t.next.model <- rq(treatment ~ . - 1, tau=.5,
                       data = as.data.frame(data.with.target),
                       weights = t.weights,  method="lasso")
    t.next.params <-  t.next.model$coefficients
    if(sum((t.next.params - t.params) ** 2) < tolerance){
      break
    }
    iteration = iteration + 1
  }
  return(t.next.params)
}


s <-   DifferenceConvexOptimize(params=NULL, train.treatment, train.covariates, 
    train.prop.scores,  train.reward, offset, PolicyFunLinearKernel, lambda, 
     hyperparams=list())

opt.with.sa <-  OptimizeParamsOfPolicyFunction(train.treatment, train.covariates, 
    train.prop.scores,  train.reward, offset, PolicyFunLinearKernel, lambda, ObjectiveFunction)
opt.with.sa$par  
s
ValueFunction(params=opt.with.sa$par, train.treatment, train.covariates, 
    train.prop.scores,  train.reward, offset, PolicyFunLinearKernel)
ValueFunction(params=s, train.treatment, train.covariates, 
    train.prop.scores,  train.reward, offset, PolicyFunLinearKernel)
ValueFunction(params=opt.with.sa$par, test.treatment, test.covariates, 
    test.prop.scores,  test.reward, control.offset, PolicyFunLinearKernel)
ValueFunction(params=s, test.treatment, test.covariates, 
    test.prop.scores,  test.reward, control.offset, PolicyFunLinearKernel)


# Empirical Value Function  -----------------------------------------------

ValueFunction <- function(params, treatment, covariates, prop.scores, reward, offset,
                          policy.function) {
  loss <- pmax(1 - abs(treatment - policy.function(params, covariates)) / offset, 0)
  return(mean(reward * loss / prop.scores / offset))
}