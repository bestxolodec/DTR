library(devtools)
load_all("../MIDAs/")
source("./src/simulations.functions.R")

train <- GetSimulationData(500, 30)
m <- MIDAs(train$covariates, train$reward, train$treatment, discrete = F, 
           depth = 4, minBucket = 10, bandwidth = 0, lambda = 0.1, kernel = 0, 
           nBoot = 3, eps = 1e-6, randomSeed = -1, 
           nSeed = 0, smallData = T, look = FALSE, nquantile = NULL) 
m
plot(m)

