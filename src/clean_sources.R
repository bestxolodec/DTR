
#source("~/yandexDisk/DIPLOMA//OWL/O_learning_functions.r")

library("truncnorm")
library(devtools)

dirname <- dirname(sys.frame(1)$ofile)
load_all(file.path(dirname, "SVMW"))
library(SVMW)  

rtruncnorm <- function (n, a = -Inf, b = Inf, mean = 0, sd = 1)  {
  if (length(n) > 1)
    n <- length(n)
  if (length(n) > 1)
    n <- length(n)
  else if (!is.numeric(n))
    stop("non-numeric argument n.")
  .Call("do_rtruncnorm", as.integer(n), a, b, mean, sd, PACKAGE="truncnorm")
}



GetRewardGivenQfunctionValuesAsMeanVec <- function(q.function.values, sd=1,
                                                   only.positive=TRUE, eps=0.01) {
  raw.reward <- rnorm(NROW(q.function.values), q.function.values, sd)
  reward <-  raw.reward - min(raw.reward) + eps
  return (list(raw.reward=as.matrix(raw.reward),
               reward=as.matrix(reward)))
}

ReplaceExtremeWithUnif <- function(values, min.val, max.val){
  violatros.index <- values < min.val | values > max.val
  values[violatros.index] <- runif(sum(violatros.index), min = min.val, max=max.val)
  return(values)
}




# Our scenarios -----------------------------------------------------------



Shvechikov.1.fopt <- function(x) {
  return (((x - .5) / .5) ** 2)
}



GetDataForShvechikov.1 <- function(sample.size, seed, sd=0.05, min.A=0, max.A=1, min.X=0, max.X=1){
  GetQFunctionValues <- function(covars, given.A, optimal.A=NULL) {
    if (is.null(optimal.A)) {
      optimal.A <- Shvechikov.1.fopt(covars)
    }
    return (-(given.A - optimal.A) ** 2)
  }
  set.seed(seed)
  X <- ReplaceExtremeWithUnif(rgamma(n=sample.size,  shape=3, rate=7.3), min.X, max.X)
  A <- ReplaceExtremeWithUnif(rgamma(n=sample.size,  shape=3, rate=7.3), min.A, max.A)
  A.opt <- Shvechikov.1.fopt(X)
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


Shvechikov.2.fopt <- function(x) {
  x = x * 100  # initially function was written to support x from 0 to 100, not from 0 to 1
  step.start = 15
  breakage.start = 40
  first.part <- x**2 * exp(-x**(.8)) +  (x - step.start) / 50 + .5
  first.and.changed.part <- ( - x**1.9 / 900 + 6 *  sin(x/6) + 1.5 * (sin(x**1.01 ) + 1) * x**(1/7) )* 1/20
  second.part <-  (- first.part + first.and.changed.part) * (x > breakage.start)
  whole <- first.part + second.part
  f <- function(x) { # scaled Runge function
    x <- x / 60   - 0.3
    cos(3*pi*x)/(1+25*(x-0.25)^2)
  }
  return ((whole + f(x)) / 2 + 0.35)
}


GetDataForShvechikov.2 <- function(sample.size, seed, sd=0.05, min.A=0, max.A=1, min.X=0, max.X=100) {
  GetQFunctionValues <- function(covars, given.A, optimal.A=NULL) {
    if (is.null(optimal.A)) {
      optimal.A <- Shvechikov.2.fopt(covars)
    }
    return (-(given.A - optimal.A) ** 2)
  }
  set.seed(seed)
  X <- ReplaceExtremeWithUnif(rgamma(n=sample.size,  shape=13, rate=0.3), min.X, max.X)
  X <- X / 100 # hack to bring X in [0,1] range
  A <- ReplaceExtremeWithUnif(rgamma(n=sample.size,  shape=2, rate=4.21), min.A, max.A)
  A.opt <- Shvechikov.2.fopt(X)
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





# Wrappers to allow using Chen code ---------------------------------------

ChangeFormatFromOurToChenEnriched <- function(our_data) {
  with(our_data, list(X=covariates, A=treatment, R=raw.reward, D_opt=optimal.treatment,
                      mu=GetQFunctionValues(covariates, treatment, optimal.treatment),
                      GetQFunctionValues=GetQFunctionValues,
                      GetOptimalTreatment=GetOptimalTreatment))
}


ChangeFormatFromChenEnrichedToOur <- function(chen_data, eps=0.01) {
  with(chen_data, list(covariates=X, treatment=A, raw.reward=R,
                       optimal.treatment=D_opt, reward=R - min(R) + eps,
                       GetQFunctionValues=GetQFunctionValues,
                       GetOptimalTreatment=GetOptimalTreatment))
}



# Predictions with Chen method --------------------------------------------

PredValueGeneral <- function(A_pred, test) {
  if ("covariates" %in% names(test)) {
    with(test, mean(GetQFunctionValues(covariates, A_pred, optimal.treatment)))
  } else { # Chen notation
    with(test, mean(GetQFunctionValues(X, A_pred, D_opt)))
  }
}

# slightly rewritten function from KO-learning to allow code reuse
pred_ko <- function(model, test) {
  # if (any((o > 2) | (o < 0)))
  warning("Clipping predicted treatment to [0,2]!")
  A_pred <- pmin(pmax(predict(model,test$X), 0), 2)
  return(list(A_pred=A_pred, Q=PredValueGeneral(A_pred, test)))
}

GetKOLearningValueAndPredictedDose <- function(train, test, q = 0.6) {
  train$weight = with(train, R - min(quantile(R, q),0))
  index = with(train, which(R > quantile(R,q)))
  model = with(train, svm(x = X[index,], y = A[index], w=weight[index],
                          type="eps-regression", epsilon = 0.15, scale=FALSE))
  return(model)
  return(pred_ko(model,test))
}



# Scenarios from Chen, 2016 -----------------------------------------------

Scenario1Enriched <- function(size, seed, ncov=30){
  GetOptimalTreatment <- function(X) {
    1 + 0.5*X[,2] + 0.5*X[,1]
  }
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  A = runif(size,0,2)
  D_opt =  GetOptimalTreatment(X)
  GetQFunctionValues <- function(X, A, A_opt=D_opt) {
    8 + 4*X[,1] - 2*X[,2] - 2*X[,3] - 25*((A_opt-A)^2)
  }
  mu <- GetQFunctionValues(X, A)
  R = rnorm(length(mu),mu,1)
  datainfo = list(X=X, A=A, R=R, D_opt=D_opt, mu=mu,
                  GetQFunctionValues=GetQFunctionValues,
                  GetOptimalTreatment=GetOptimalTreatment)
  return(datainfo)
}



Scenario2Enriched <- function(size, seed, ncov=10){
  GetOptimalTreatment <- function(X) {
    I(X[,1] > -0.5)*I(X[,1] < 0.5)*0.6 + 1.2*I(X[,1] > 0.5) + 1.2*I(X[,1] < -0.5) +
      X[,4]^2 + 0.5*log(abs(X[,7])+1) - 0.6
  }
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  A = runif(size,0,2)
  D_opt = GetOptimalTreatment(X)
  GetQFunctionValues <- function(X, A, A_opt=D_opt) {
    8 + 4*cos(2*pi*X[,2]) - 2*X[,4] - 8*X[,5]^3 - 15*abs(D_opt-A)
  }
  mu = GetQFunctionValues(X, A)
  R = rnorm(length(mu),mu,1)
  datainfo = list(X=X, A=A, R=R, D_opt=D_opt, mu=mu,
                  GetQFunctionValues=GetQFunctionValues,
                  GetOptimalTreatment=GetOptimalTreatment)
  return(datainfo)
}


Scenario4Enriched <- function(size, seed, ncov=10){
  GetOptimalTreatment <- function(X) {
    I(X[,1] > -0.5)*I(X[,1] < 0.5)*0.6 + 1.2*I(X[,1] > 0.5) + 1.2*I(X[,1] < -0.5) +
      X[,4]^2 + 0.5*log(abs(X[,7])+1) - 0.6
  }
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  D_opt = GetOptimalTreatment(X)
  A = rtruncnorm(size,a=0,b=2,mean=D_opt,sd=0.5)
  GetQFunctionValues <- function(X, A, A_opt=D_opt) {
    8 + 4*cos(2*pi*X[,2]) - 2*X[,4] - 8*X[,5]^3 - 15*abs(A_opt-A)
  }
  mu = GetQFunctionValues(X, A)
  R = rnorm(length(mu),mu,1)
  datainfo = list(X=X, A=A, R=R, D_opt=D_opt, mu=mu,
                  GetQFunctionValues=GetQFunctionValues,
                  GetOptimalTreatment=GetOptimalTreatment)
  return(datainfo)

}

## Sanity checks
# check_function <- function(ko, o) {
#   print(all.equal(o$X, ko$X))
#   print(all.equal(o$A, ko$A))
#   print(all.equal(o$R, ko$R))
#   print(all.equal(o$D_opt, ko$D_opt))
#   print(all.equal(o$mu, ko$mu))
# }
# check_function(Scenario1(100, 30, 0), Scenario1Enriched(100,30,0))
# check_function(Scenario2(100, 30, 0), Scenario2Enriched(100,30,0))
# check_function(Scenario4(100, 30, 0), Scenario4Enriched(100,30,0))



# HOW TO USE :
#
# n_samples <- 100
# n_test_samples <- 100
# noise_sd <- 0
#
# GetData  <- GetDataForShvechikov.2
# GetData  <- GetDataForShvechikov.1
# train <- GetData(n_samples, noise_sd)
# test <- GetData(n_test_samples, noise_sd)
# ko_train <- ChangeFormatFromOurToChenEnriched(train)
# ko_test <- ChangeFormatFromOurToChenEnriched(test)
#
#
#
#
#
# n_samples <- 100
# n_test_samples <- 200
# n_covariates <- 10
#
# GetData <- Scenario1Enriched
# GetData <- Scenario2Enriched
# GetData <- Scenario4Enriched
# ko_train <- GetData(n_samples, ncov = n_covariates, seed = 0)
# ko_test <- GetData(n_test_samples, ncov = n_covariates, seed = 1)
# train <- ChangeFormatFromChenEnrichedToOur(ko_train)
# test <- ChangeFormatFromChenEnrichedToOur(ko_test)







