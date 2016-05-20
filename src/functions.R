library(caret)
library(GenSA)
library(CBPS)
library(reshape)
library(quantreg)
library(maxLik)


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
  x.balanced <- CBPS(balance.formula, data=observations.with.metadata,
                     twostep=two.step, ATT=0)
  if (isTRUE(return.cbps.object)) {
    return (x.balanced)
  } else {
    return (x.balanced$fitted.values)
  }
}

# x.balanced <- CBPS(treatment ~ . - covars..Intercept.,
#    data=data.frame("treatment"=test.treatment,  "covars"=test.covariates),
#    twostep=T, ATT=0)
#
# prop.scores <- GetPropensityScores(
#   data.frame("treatment"=train.treatment, "covars"=train.covariates),
#   balance.formula = formula (treatment ~ . - covars..Intercept.),
#   two.step = T)




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







# Optimization with DC functions ------------------------------------------

GetNextStepParams  <- function(covars, treatment, relative.weights, lambda) {
  # remove intercept
  covars  <- covars[, ! grepl("intercept", tolower(colnames(covars)))]
  # check there is no more intercept term
  stopifnot(apply(covars, 2, sd) != 0)

  sum.rel.weights  <- sum(relative.weights)
  adj.vector <- colSums(as.vector(relative.weights) * covars / sum.rel.weights)
  weighted.covars <- t(as.vector(relative.weights) * covars)
  value.matrix  <- weighted.covars %*% sweep(covars, 2, adj.vector)
  reg.matrix <- lambda * diag(x = ncol(value.matrix))
  inv.matrix  <- solve(reg.matrix  + value.matrix)

  subst.coef <-  sum(relative.weights * treatment) / sum.rel.weights
  right.vector <- as.matrix(colSums(as.vector(relative.weights) * covars *
                                    as.vector(treatment - subst.coef)))
  w.next = inv.matrix %*% right.vector
  b.next = sum(as.vector(relative.weights) * 
               (treatment - covars %*% w.next)) /  sum.rel.weights
  return (as.numeric(c(b.next, w.next)))
}



## TODO: Possible improvements in first coefs gessing with an lm model
DCOptimizeWithMML2Penalized <- function(params=NULL, obs.data, offset,
    policy.function, lambda, hyperparams=list(), 
    opt.hyperparams=list()) {
  default.list <- list(debug.file=NULL, approximation.eps=NULL,  tolerance=1e-6)
  opt.hyperparams <- modifyList(default.list, opt.hyperparams)
  stopifnot(is.matrix(obs.data$covariates))
  data.with.target <- c(data.frame(obs.data$treatment),
                        as.data.frame(obs.data$covariates))
  ## TODO: remove this !
  if (is.null(params)) {
    t.params <- runif(ncol(obs.data$covariates), min=-1, max=1)
  } else {
    t.params <- params
  }
  subsequent.converged.iters <- 5
  converged.iters <- 0
  iteration <- 0
  iter.info = list()
  repeat {
    if (iteration > 1000) {
      save(iter.info, file = opt.hyperparams$debug.file)
      stop("Infinite iterations of MM algorithm! Saved iter.info to ", 
           opt.hyperparams$debug.file)
    }
    iteration = iteration + 1
    t.prediction <- policy.function(t.params, obs.data$covariates, hyperparams)
    t.abs.deviance.from.treatment <- abs(obs.data$treatment - t.prediction)
    t.Q.values <- ifelse(t.abs.deviance.from.treatment <= offset, 1, 0)
    t.weights <- obs.data$reward * t.Q.values / (offset ** 2) / obs.data$prop.scores
    if (is.null(opt.hyperparams$approximation.eps)) {
      # cat("OLD METHOD!\n")
      t.relative.weights <- t.weights / t.abs.deviance.from.treatment
    } else {
      t.relative.weights <- 2 * t.weights ** 2 / 
        (opt.hyperparams$approximation.eps + 2 * t.weights * t.abs.deviance.from.treatment)
    }
    t.next.params <- GetNextStepParams(obs.data$covariates, obs.data$treatment, 
                                       t.relative.weights, lambda)
    discrepancy <- sum((t.next.params - t.params) ** 2)
    
    
    if (! is.null(opt.hyperparams$debug.file)) {
      vf.test = ValueFunction(t.next.params, test, offset, policy.function)
      vf.train = ValueFunction(t.next.params, train, offset, policy.function)
      objf.test = ObjectiveFunction(t.next.params, obs.data = test, 
                               policy.function = policy.function, 
                               offset=offset, lambda=lambda)
      objf.train = ObjectiveFunction(t.next.params, obs.data = train, 
                               policy.function = policy.function, 
                               offset=offset, lambda=lambda)
      next.dctfun.val = DCTargetFunction(t.next.params, t.params, train, offset, 
                                    policy.function, lambda)  
      prev.dctfun.val = DCTargetFunction(t.params, t.params, train, offset, 
                                    policy.function, lambda)  
      info <- list("discrepancy" =  discrepancy, "valuefunc.train" = vf.train, 
                   "valuefunc.test" = vf.test, "obj.func.train" = objf.train,  
                   "obj.func.test" = objf.test, "dct" = prev.dctfun.val, 
                   "deltadct" =   prev.dctfun.val - next.dctfun.val, 
                   "params" = t.next.params, "iteration"=iteration)
      iter.info[[length(iter.info) + 1]] <- info
    }
    if(discrepancy < opt.hyperparams$tolerance) {
      if (converged.iters == subsequent.converged.iters) {
        if (! is.null(opt.hyperparams$debug.file)) {
          save(iter.info, file = opt.hyperparams$debug.file)
          cat("Converged after ", iteration, " iterations\n")
        }
        break 
      } else {
        converged.iters <-  converged.iters + 1
      }
    } else {
      converged.iters <- 0
    }
    t.params  <- t.next.params
  }
  return(t.next.params)
}



GetFuncMinWithMM <- function(init.pars, covars, treatment, weights, lambda, eps, tol) {
  params <- init.pars
  sequential.converged.iters <- 3
  converged.iters <- 0
  iteration <- 0
  repeat {
    iteration <- iteration + 1
    prediction <- PolicyFunLinearKernel(params, covars)
    abs.deviance <- abs(treatment - prediction)
    stopifnot(sum(abs.deviance < 1e-17) == 0)
    if (is.null(eps)) {
      relative.weights <- weights / abs.deviance
    } else {
      relative.weights <- 2 * weights ** 2 / (eps + 2 * weights * abs.deviance)
    }
    next.params <- GetNextStepParams(covars, treatment, relative.weights, lambda)
    discrepancy <- sum((next.params - params) ** 2)
    if(discrepancy < tol) {
      if (converged.iters == sequential.converged.iters) {
        break 
      } else {
        converged.iters <-  converged.iters + 1
      }
    } else {
      converged.iters <- 0
    }
    params <- next.params
  }
  # cat("MM converged after ", iteration - sequential.converged.iters, " iterations\n")
  return (next.params)
}

## TODO: Possible improvements in first coefs gessing with an lm model
DCOptimizeWithMML2PenalizedProperIters <- function(params=NULL, obs.data, offset,
    policy.function, lambda, hyperparams=list(), 
    opt.hyperparams=list()) {
  default.list <- list(debug.file=NULL, approximation.eps=NULL, tolerance=1e-6, q=0.65)
  opt.hyperparams <- modifyList(default.list, opt.hyperparams)
  stopifnot(is.matrix(obs.data$covariates))
  stopifnot(! is.null(params))
  t.params <- params
  # if (is.null(params)) {
  #   
  #   index = obs.data$reward > quantile(obs.data$reward, opt.hyperparams$q) 
  #   rqmodel = rq(obs.data$treatment[index] ~ obs.data$covariates[index,] - 1, tau=.5, 
  #                method="lasso", weights=obs.data$reward[index], lambda = 1)
  #   t.params = coef(rqmodel) 
  #   cat("GENERATE initial approximation of decision!\n")
  # } else {
  #   cat("Setting supplied params!\n")
  #   t.params <- params
  # }
  iteration <- 0
  t.prev.Q.values <- rep(0, length(obs.data$treatment))
  repeat {
    iteration <- iteration + 1
    if (iteration > 1000) {
      stop("Infinite iterations (more than 1000) of DC algorithm!")
    }
    t.prediction <- policy.function(t.params, obs.data$covariates, hyperparams)
    t.abs.deviance.from.treatment <- abs(obs.data$treatment - t.prediction)
    t.Q.values <- ifelse(t.abs.deviance.from.treatment <= offset, 1, 0)
    if (sum(t.prev.Q.values !=  t.Q.values) == 0)  {
     break 
    }
    t.weights <- obs.data$reward * t.Q.values / (offset ** 2) / obs.data$prop.scores
    t.params <- GetFuncMinWithMM(t.params, obs.data$covariates, 
                                 obs.data$treatment, t.weights,  
                                 lambda, opt.hyperparams$approximation.eps, 
                                 opt.hyperparams$tolerance)
    t.prev.Q.values <- t.Q.values
  }
  cat("MML2PenalizedProperIters Converged after ", iteration, " iterations\n")
  return(t.params)
}





DCOptimizeL1Penalized <- function(params=NULL, obs.data, offset,
    policy.function, lambda, hyperparams=list(),
    opt.hyperparams=list()) {
  opt.hyperparams <- modifyList(opt.hyperparams, list(tolerance=1e-8))
  stopifnot(is.matrix(obs.data$covariates))
  data.with.target <- c(data.frame(obs.data$treatment),
                        as.data.frame(obs.data$covariates))
  if (is.null(params)) {
    t.params <- runif(ncol(obs.data$covariates), min=-1, max=1)  # default params
  } else {
    t.params <- params
  }
  t.next.params <- t.params  # simply for first loop iteration
  iteration <- 0
  subsequent.converged.iters <- 3
  converged.iters <- 0
  repeat {
    t.params <- t.next.params
    t.prediction <- policy.function(t.params, obs.data$covariates, hyperparams)
    t.abs.deviance.from.treatment  <- abs(obs.data$treatment - t.prediction)
    t.Q.values  <-  ifelse(t.abs.deviance.from.treatment <= offset, 1, 0)
    t.weights  <- obs.data$reward * t.Q.values / (offset ** 2)  / obs.data$prop.scores
    t.next.model <- rq(obs.data.treatment ~ . - 1, tau=.5,
                       data = as.data.frame(data.with.target),
                       weights = t.weights,  method="lasso", lambda=lambda)
    t.next.params <-  coef(t.next.model)
    discrepancy <- sum((t.next.params - t.params) ** 2)
    if(discrepancy < opt.hyperparams$tolerance) {
      if (converged.iters == subsequent.converged.iters) {
        break 
      } else {
        converged.iters <-  converged.iters + 1
      }
    } else {
      converged.iters <- 0
    }
    iteration = iteration + 1
  }
  cat("DCOptimizeL1Penalized converged after ", iteration, " iterations\n")
  return(t.next.params)
}

  


GetOwlParams <- function (data, lambda, weights=F, q = 0.6, t = 0.5) {
  constant = min(quantile(data$raw.reward, q), 0)
  data$weight = data$raw.reward - constant
  index = which(data$raw.reward > quantile(data$raw.reward, q))
  if (weights) {
    rqmodel = rq.with.weights(data$treatment[index] ~ data$covariates[index, ] - 1, tau = t, 
                              method = "lasso", weights = data$weight[index], lambda = lambda)
  } else {
    rqmodel = rq(data$treatment[index] ~ data$covariates[index, ] - 1, tau = t, method = "lasso", 
                 weights = data$weight[index], lambda = lambda)
  }
  coefs = coef(rqmodel)
  return(matrix(coefs))
}




# s <- DifferenceConvexOptimize(params=NULL, train.treatment, train.covariates,
#      train.prop.scores,  train.reward, offset, PolicyFunLinearKernel, lambda)
# opt.with.sa.par <-  OptimizeParamsOfPolicyFunction(train.treatment, train.covariates,
#     train.prop.scores,  train.reward, offset, PolicyFunLinearKernel, lambda,
#     ObjectiveFunction)
# ValueFunction(params=opt.with.sa.par, test.treatment, test.covariates,
#     test.prop.scores,  test.reward, control.offset, PolicyFunLinearKernel)
# ValueFunction(params=s, test.treatment, test.covariates,
#     test.prop.scores,  test.reward, control.offset, PolicyFunLinearKernel)



# Empirical Value Function  -----------------------------------------------

ValueFunction <- function(params, obs.data,  offset, policy.function, 
                          hyperparams=list()) {
  prediction <- policy.function(params, obs.data$covariates, hyperparams)
  abs.deviance <- abs(obs.data$treatment - prediction)
  gain <- pmax(1 - abs.deviance / offset, 0)
  return(mean(obs.data$raw.reward * gain / obs.data$prop.scores / offset))
}

# Objective function ------------------------------------------------------

ObjectiveFunction <- function(params, obs.data, offset, policy.function, lambda,
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
  prediction <- policy.function(params, obs.data$covariates, hyperparams)
  loss <- pmin(abs(obs.data$treatment - prediction) / offset, 1)
  denominator <- obs.data$prop.scores * (2 * offset)
  risk.function.value <- mean(obs.data$raw.reward * loss / denominator)
  # tail is because we dont penalize intercept parameter
  regularization.value <- lambda * sum(tail(params, -1) ** 2)
  return(risk.function.value + regularization.value)
}



# DC Target function ------------------------------------------------------
DCTargetFunction <- function(next.params, prev.params, obs.data, offset, 
                             policy.function, lambda, hyperparams=list()) {
  t.prediction <- policy.function(prev.params, obs.data$covariates, hyperparams)
  t.abs.deviance.from.treatment = abs(obs.data$treatment  -  t.prediction)
  t.Q.R.values <- ifelse(t.abs.deviance.from.treatment <= offset, obs.data$reward, 0)
  
  t.next.prediction <- policy.function(next.params, obs.data$covariates, hyperparams)
  t.next.abs.deviance.from.treatment <- abs(obs.data$treatment - t.next.prediction)
  risk.func.value <- sum(t.Q.R.values * t.next.abs.deviance.from.treatment) / (offset ** 2)
  # tail is because we dont penalize intercept parameter
  regularization.value <- 0.5 * lambda * sum(tail(next.params, -1) ** 2)
  return(risk.func.value + regularization.value) 
}



# Newton Raphson optimization of Meshalkin function -----------------------

TriangularLoss <- function (u, offset) {
  return ((abs(u) <= offset) * abs(u) / offset  + (abs(u) > offset) * 1)
}

MeshalkinLoss <- function(u, offset=1) {
  offset = offset / 5
  return  (1 - exp(-u**2 / offset ))
}
MeshalkinLoss.grad  <- function(u, offset=1) {
  offset = offset / 5
  return (2 * u / offset * exp(-u**2 / offset))
}
MeshalkinLoss.hess <- function(u, offset=1) {
  offset = offset / 5
  return (2 / offset * exp(-u**2 / offset) * (1 - 2 * u ** 2 / offset))
}

# WARN: Using only 
# TODO: Decide what to use here: reward or raw reward ? 
MeshalkinObjFunc <- function(params, obs.data,  offset, policy.function, lambda_reg, hyperparams=list()) { 
  prediction <- policy.function(params, obs.data$covariates, hyperparams)
  deviance <- obs.data$treatment - prediction
  loss <- MeshalkinLoss(deviance, offset = offset)
  return(mean(obs.data$reward * loss / obs.data$prop.scores / offset)  + lambda_reg * params %*% params )
}

MeshalkinObjFunc.grad <- function(params, obs.data,  offset, policy.function, lambda_reg, hyperparams=list()) { 
  prediction <- policy.function(params, obs.data$covariates, hyperparams)
  deviance <- obs.data$treatment - prediction
  loss <- MeshalkinLoss.grad(deviance, offset = offset)
  obj.func.coefs <- obs.data$reward * loss / obs.data$prop.scores / offset 
  obj.func.value <-  colMeans(Diagonal(n=NROW(obj.func.coefs), x=obj.func.coefs) %*% (-obs.data$covariates))
  return(obj.func.value  + 2 * lambda_reg * params)
}


MeshalkinObjFunc.hess <- function(params, obs.data,  offset, policy.function, lambda_reg, hyperparams=list()) { 
  prediction <- policy.function(params, obs.data$covariates, hyperparams)
  deviance <- obs.data$treatment - prediction
  loss <- MeshalkinLoss.hess(deviance, offset = offset)
  obj.func.matrix <- obs.data$reward * loss / obs.data$prop.scores / offset
  n <- length(obs.data$reward)
  acc <- 0 
  for (i in seq(1, n)) {
    acc <- acc + obj.func.matrix[i] * obs.data$covariates[i, ] %o% obs.data$covariates[i, ] 
  }
  return(acc / n  + 2 * lambda_reg * diag(NCOL(obs.data$covariates)))
}

MeshalkinObjFunc.2max <- function(...) {- MeshalkinObjFunc(...)}
MeshalkinObjFunc.grad.2max <- function(...) {- MeshalkinObjFunc.grad(...)}
MeshalkinObjFunc.hess.2max <- function(...) {- MeshalkinObjFunc.hess(...)}
















# Cauchy loss function ----------------------------------------------------

CauchyLoss <- function(u, c=1) {
  return (c / 2 * log(1  + (u / c) ^ 2)) 
}    

CauchyLoss.grad <- function(u, c=1) {
  return (u  /  (1 + (u / c) ^ 2))
}    

CauchyLoss.hess <- function(u, c=1) {
  return  ((1 - (u / c) ^ 2)  /  (1 + (u / c)^2)^2)
}  

# WARN: Using only 
# TODO: Decide what to use here: reward or raw reward ? 
CauchyObjFunc <- function(params, obs.data,  offset, policy.function, lambda_reg, hyperparams=list()) { 
  prediction <- policy.function(params, obs.data$covariates, hyperparams)
  deviance <- obs.data$treatment - prediction
  loss <- CauchyLoss(deviance, offset)
  return(mean(obs.data$reward * loss / obs.data$prop.scores / offset)  + lambda_reg * params %*% params )
}

CauchyObjFunc.grad <- function(params, obs.data,  offset, policy.function, lambda_reg, hyperparams=list()) { 
  prediction <- policy.function(params, obs.data$covariates, hyperparams)
  deviance <- obs.data$treatment - prediction
  loss <- CauchyLoss.grad(deviance, offset)
  obj.func.coefs <- obs.data$reward * loss / obs.data$prop.scores / offset 
  obj.func.value <-  colMeans(
    Diagonal(n=NROW(obj.func.coefs), x=obj.func.coefs) %*% (-obs.data$covariates)
  )
  return(obj.func.value  + 2 * lambda_reg * params)
}


CauchyObjFunc.hess <- function(params, obs.data,  offset, policy.function, lambda_reg, hyperparams=list()) { 
  prediction <- policy.function(params, obs.data$covariates, hyperparams)
  deviance <- obs.data$treatment - prediction
  loss <- CauchyLoss.hess(deviance, offset)
  obj.func.matrix <- obs.data$reward * loss / obs.data$prop.scores / offset
  n <- length(obs.data$reward)
  acc <- 0 
  for (i in seq(1, n)) {
    acc <- acc + obj.func.matrix[i] * obs.data$covariates[i, ] %o% obs.data$covariates[i, ] 
  }
  return(acc / n  + 2 * lambda_reg * diag(NCOL(obs.data$covariates)))
}


CauchyObjFunc.2max <- function(...) {- CauchyObjFunc(...)}
CauchyObjFunc.grad.2max <- function(...) {- CauchyObjFunc.grad(...)}
CauchyObjFunc.hess.2max <- function(...) {- CauchyObjFunc.hess(...)}




# Newton Raphson optimization ---------------------------------------------



NROptimize <- function (params, obs.data,  offset, policy.function, lambda,  
                        hyperparams=list(), opt.hyperparams=list()) {
  res <- maxNR(opt.hyperparams$func, grad=opt.hyperparams$func.grad,  
               hess=opt.hyperparams$hess, start=opt.hyperparams$init.pars, 
               obs.data=train,offset=offset, policy.function=policy.function, lambda_reg=lambda) 
  cat("Code: ", res$code, " Message: ", res$message, "\n")
  return (res$estimate)
} 




