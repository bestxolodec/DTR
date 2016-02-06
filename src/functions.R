library(caret)
library(GenSA)
library(CBPS)
library(reshape)


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
    file.path <- "~/yandexDisk/DIPLOMA/data/ipwc_data.csv"
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
  value <- observations.with.metadata[, observed.outcome.colname]
  tmp <-   - abs(value - ideal.outcome.value) 
  return (tmp - min(tmp) + eps)
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

ObjectiveFunction <- function(params, patient.covariates, propensity.scores,
                              offset, policy.function, lambda, reward.values,
                              treatment.values, hyperparams=list()) {
  # Objective (risk) function for minimization problem. Includes regularization.
  #
  # Args:
  #   params: Model parameters for policy.function -- are being changed
  #           during minimization procedure.
  #   patient.covariates: Obseravation data of patients.
  #   propensity.scores: Correspond to observations in data.
  #   offset: Smallest value of difference between
  #           observational treatment value and the predicted one,
  #           which is interpeted as full "1" loss
  #   policy.function: Decision function of the form function(params, data, ...)
  #                    which returns treatment given covariates.
  #   lambda: Regularization parameter.
  #   reward.values: Precomputed reward values of treatment.values treatment.
  #                  The greater reward, the better.
  #   treatment.values: Observed treatment values of patients.
  #   hyperparams:  list of hyperparameters for policy function
  #
  # Returns:
  #   Regularized risk function of DTR wich is specified with policy.function.
  # browser()
  prediction <- policy.function(params, patient.covariates, hyperparams)
  multiplier <- pmin(abs(treatment.values - prediction) / offset, 1)
  risk.function.value <- mean(reward.values / propensity.scores / (2 * offset) * multiplier)
  # tail is because we dont normilize intercept parameter
  regularization.value <- lambda * norm(as.matrix(tail(params, -1)), type="F") ** 2
  return(risk.function.value + regularization.value)
}




# Optimization with DC functions ------------------------------------------


DifferenceConvexOptimize <- function(initial.params, patient.covariates,
                                     propensity.scores, offset,
                                     policy.function, lambda, reward.values,
                                     treatment.values, hyperparams=list(),
                                     next.gen.obj.fun, tolerance = 0.01) {
  data.with.target <- cbind(treatment.values, patient.covariates)
  t.params <- initial.params
  # simply for first loop iteration
  t.next.params <- t.params
  while(sum((t.next.params-t.params)^2) > tolerance) {
    t.params <- t.next.params
    t.prediction <- policy.function(t.params, patient.covariates, hyperparams)
    t.abs.deviance.from.treatment  <- abs(treatment.values - t.prediction)
    t.Q.values  <-  ifelse(t.abs.deviance.from.treatment <= offset, 0, 1)
    t.weights  <- reward.values * t.Q.values / offset**2
    # TODO: think of how to use rq.fit.lasso() here
    # -1 as an intercept term is already present ind data matirx
    t.next.model <-  rq(treatment.values ~ . , tau=.5, 
                        data = as.data.frame(data.with.target), 
                        weights = t.weights,  method="lasso")
  }
}




# Experiment linear model  ------------------------------------------------

OLSCostFunction <- function(params, patient.covariates) {
  
}

