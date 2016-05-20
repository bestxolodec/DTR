source('../../OWL/O_learning_functions.r')

##### simulation Examples ##########
##### L-O-L for Scenario 1.  #######

ind=5
allloop =0
datagen = Scenario1
covsize=30
epsilon = 0.1
pred = pred_s1

results_lol_s1 <- list()
## It is fine to run all the simulation in one loop
# for(samplesize in c(50,100,200,400,800)){
for(samplesize in c(800)){
  allloop = allloop + 1
  ### simsize is number of simulation runs, in the paper we use 200
  simsize = 200
  lol_value = rep(0,simsize)
  pred_dose <- list()
  test = datagen(5000,covsize,seed=ind+30000)
  for (loop in 1:simsize){
    print(loop)
    train = datagen(samplesize,covsize,seed=loop+20002+samplesize)
    constant = min(quantile(train$R,0.6),0)
    train$weight = train$R - constant
    index = which(train$R > quantile(train$R,0.6))
    ### We use the penalized quantile regression to enhance the model fitting.
    rqmodel = rq(train$A[index] ~ train$X[index,], .5, method="lasso", weights=train$weight[index], lambda = 1)
    tmpresults = pred_s1(rqmodel,test)
    lol_value[loop] = tmpresults$pred_value
    pred_dose[[loop]] = tmpresults$pred_dose
  }
  results=list(lol=lol_value,dose=pred_dose,size=samplesize)
  results_lol_s1[[allloop]] <- results
}

sapply(results_lol_s1, function(x) mean(x$lol))
sapply(results_lol_s1, function(x) sd(x$lol))



covsize=30
epsilon = 0.1
lambda = c(0.1, 0.2, 0.5, 1, 3) 
pred = pred_s1
samplesize = 200
datagen = Scenario1


train = datagen(samplesize,covsize,seed=loop+20002+samplesize)
train$R = train$R - min(train$R)
dose.cv(train, epsilon, lambda, cvfold=5, tunefunc=loss_O_learning, optimfunc=dc_solution,seed=0)






















# Defining control execution constants 
train.data.sample.sizes <- c(50, 100, 200, 400, 800)
test.data.sample.size <- 10000
number.of.covariates = 30
sample.size <- 800
offsets = seq(0.2, 1, length = 10)
control.offset = min(offsets) / 2
lambdas = seq(0,  10,  length=12)

registerDoParallel(cores = 4)
number.of.iters <- 10

offset = 0.1
lambda = 10


GetInitPars <- function(train, q=0.6) {
  index = train$reward > quantile(train$reward, q)
  rqmodel = rq(train$treatment[index] ~ train$covariates[index,] - 1, tau=.5, 
               method="lasso", weights=train$reward[index], lambda = lambda)
  return(coef(rqmodel))
}

train <- GetSimulationData(sample.size, number.of.covariates)
opt.params = list("opt.func"=DCOptimizeWithMML2PenalizedProperIters, "tolerance"=1e-7, "init.pars"=GetInitPars(train, q=q))
cs[[length(cs) + 1]] <- OptimizeParamsOfPolicyFunction(train, offset, PolicyFunLinearKernel, lambda=lambda, opt.hyperparams=opt.params)











data(stackloss)
rq(stack.loss ~ stack.x,.5)  #median (l1) regression  fit for the stackloss data. 
rq(stack.loss ~ stack.x,.25)  #the 1st quartile, 
#note that 8 of the 21 points lie exactly on this plane in 4-space! 
rq(stack.loss ~ stack.x, tau=-1)   #this returns the full rq process
rq(rnorm(50) ~ 1, ci=FALSE)    #ordinary sample median --no rank inversion ci
rq(rnorm(50) ~ 1, weights=runif(50),ci=FALSE)  #weighted sample median 
#plot of engel data and some rq lines see KB(1982) for references to data
data(engel)
attach(engel)
plot(income,foodexp,xlab="Household Income",ylab="Food Expenditure",type = "n", cex=.5)
points(income,foodexp,cex=.5,col="blue")
taus <- c(.05,.1,.25,.75,.9,.95)
xx <- seq(min(income),max(income),100)
f <- coef(rq((foodexp)~(income),tau=taus))
yy <- cbind(1,xx)%*%f
for(i in 1:length(taus)){
  lines(xx,yy[,i],col = "gray")
}
abline(lm(foodexp ~ income),col="red",lty = 2)
abline(rq(foodexp ~ income), col="blue")
legend(3000,500,c("mean (LSE) fit", "median (LAE) fit"),
       col = c("red","blue"),lty = c(2,1))
#Example of plotting of coefficients and their confidence bands
plot(summary(rq(foodexp~income,tau = 1:49/50,data=engel)))
#Example to illustrate inequality constrained fitting
n <- 100
p <- 5
X <- matrix(rnorm(n*p),n,p)
y <- .95*apply(X,1,sum)+rnorm(n)
#constrain slope coefficients to lie between zero and one
R <- cbind(0,rbind(diag(p),-diag(p)))
r <- c(rep(0,p),-rep(1,p))
rq(y~X,R=R,r=r,method="fnc")






data <- list()
object.index = 1:4
data$X <- as.matrix(data.frame(intercept=rep(1,4), x1=c(1,1, 5, 10), x2=c(1, 6, 4, 1)))
data$R <- as.matrix(data.frame(a1=c(0, 49, 101, 0), a0=c(101, 100, 0, 50)))
R.index <- c(1, 2 , 1, 2)
R.observed <- data$R[cbind(object.index, R.index)]
R.observed.01.scaled <- ScaleToZeroOne(R.observed)
plot(data$X[, 2:3], col=rgb(1 - R.observed.01.scaled, R.observed.01.scaled ,  0), pch=20, cex=3)
grid()
residuals <- ScaleToZeroOne(resid(lm(R.observed ~ data$X -1)))
plot(data$X[, 2:3], col=rgb(1 - residuals, residuals,  0), pch=20, cex=3)
plot(object.index, residuals, col=rgb(1 - residuals, residuals,  0),  pch=20, cex=2)
text(abs_losses, percent_losses, labels=namebank, cex= 0.7, pos=3)





wireframe(Height ~ x*y, data = elevation.fit,
          xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
          main = "Surface elevation data",
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -60, x = -60)
)






xspace <- 
expand.grid()

persp(seq(10, 300, 5), seq(10, 300, 5), z, phi = 45, theta = 45,
      xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
      main = "Surface elevation data"
)







wireframe(volcano, shade = TRUE,
          aspect = c(61/87, 0.4),
          light.source = c(10,0,10))

g <- expand.grid(x = 1:10, y = 5:15, gr = 1:2)
g$z <- log((g$x^g$gr + g$y^2) * g$gr)
wireframe(z ~ x * y, data = g, 
          scales = list(arrows = FALSE),
          drape = TRUE, colorkey = TRUE,
          screen = list(z = 30, x = -60))













library(maxLik)

xspace <-  seq(-10, 10, 0.01)


lossf <- function(u, offset=1) {
  return (log(1 + exp(-u)))
}

lossf <- function(u, offset=1) {
  return (1 - exp(-u**2 / offset))
}

lossf_diff <- function(u, offset=1) {
  return (2 * u * exp(-u**2 / offset))
}

lossf_hess <- function(u, offset=1) {
  return (2 / offset * exp(-u**2 / offset) * (1 - 2 * u ** 2 / offset))
}



plot(xspace, lossf(xspace, 10), type="l")
plot(xspace, lossf_diff(xspace), type="l", ylim=c(-2.5,2.5))
plot(xspace, lossf_hess(xspace), type="l")


neg_lossf <- function(p) {- lossf(p)} 
neg_lossf_diff <- function(p) {- lossf_diff(p)} 
neg_lossf_hess <- function(p) {- lossf_hess(p)} 

par(mfrow=c(3,1))
init.pars <- seq(-5, 5, 0.1)
max_vals <- sapply(init.pars, function(p) (maxNR(neg_lossf, grad=neg_lossf_diff, hess=neg_lossf_hess,  start=p)$maximum))
plot(init.pars, max_vals)
lines(init.pars, neg_lossf(init.pars), type="l")
plot(init.pars, neg_lossf_diff(init.pars))
plot(init.pars, neg_lossf_hess(init.pars))




for (i in seq(-10, 10, 0.1)) {
  res <- maxNR(neg_lossf, grad=neg_lossf_diff, hess=neg_lossf_hess,  start=i, control=list(steptol=1e-20))
  cat("For initial ", i , " Maximum = ", res$maximum, "\n")
}

maxNR(neg_lossf, grad=neg_lossf_diff, hess=neg_lossf_hess,  start=0.5, control=list(steptol=1e-20))






plot(xspace, lossf(xspace,  offset=1), type="l", ylim=c(0,1.1))
lines(xspace, lossf(xspace, offset=3), type="l")
lines(xspace, lossf(xspace, offset=1), type="l")
lines(xspace, lossf(xspace, offset=0.1), type="l")
lines(xspace, lossf(xspace, offset=0.001), type="l")
grid()
# 
plot(xspace,  lossf_diff(xspace, offset=5), type="l", ylim=c(-2.5,2.5))
lines(xspace, lossf_diff(xspace, offset=1), type="l")
lines(xspace, lossf_diff(xspace, offset=0.1), type="l")
grid()

init.par = -100
optim(init.par, lossf, lossf_diff, method = "L-BFGS-B")
optim(init.par, lossf,  method = "SANN")

optim(init.par, lossf, lossf_diff, method = "BFGS")
optim(init.par, lossf, lossf_diff, method = "CG", control = list(type = 1))
optim(init.par, lossf, lossf_diff, method = "CG", control = list(type = 2))
optim(init.par, lossf, lossf_diff, method = "CG", control = list(type = 3))
optim(init.par, lossf, lossf_diff, method = "Brent")




fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
  x1 <- x[1]
  x2 <- x[2]
  c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 *      (x2 - x1 * x1))
}
optim(c(-1.2,1), fr)
(res <- optim(c(-1.2,1), fr, grr, method = "BFGS"))
optimHess(res$par, fr, grr)
optim(c(-1.2,1), fr, NULL, method = "BFGS", hessian = TRUE)
## These do not converge in the default number of steps

xspace [lossf(xspace, b=0.5) >= 1]

maxNR(function(p) {- lossf(p)}, grad=lossf_diff , start=1)











## estimate the exponential distribution parameter by ML
t <- rexp(100, 2)
loglik <- function(theta) sum(log(theta) - theta*t)
## Note the log-likelihood and gradient are summed over observations
gradlik <- function(theta) sum(1/theta - t)
hesslik <- function(theta) -100/theta^2
## Estimate with finite-difference gradient and Hessian
a <- maxNR(loglik, start=1, control=list(printLevel=2))
summary(a)
## You would probably prefer 1/mean(t) instead ;-)
## Estimate with analytic gradient and Hessian
a <- maxNR(loglik, gradlik, hesslik, start=1)
summary(a)

## BFGS estimation with finite-difference gradient
a <- maxBFGSR( loglik, start=1 )
summary(a)

## For the BHHH method we need likelihood values and gradients
## of individual observations
loglikInd <- function(theta) log(theta) - theta*t
gradlikInd <- function(theta) 1/theta - t
## Estimate with analytic gradient
a <- maxBHHH(loglikInd, gradlikInd, start=1)
summary(a)

##
## Example with a vector argument:  Estimate the mean and
## variance of a random normal sample by maximum likelihood
## Note: you might want to use maxLik instead
##
loglik <- function(param) {
  mu <- param[1]
  sigma <- param[2]
  ll <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
  ll
}
x <- rnorm(100, 1, 2) # use mean=1, stdd=2
N <- length(x)
res <- maxNR(loglik, start=c(0,1)) # use 'wrong' start values
summary(res)
##
## The previous example with named parameters and fixed values
##
resFix <- maxNR(loglik, start=c(mu=0, sigma=1), fixed="sigma")
summary(resFix)  # 'sigma' is exactly 1.000 now.
###
### Constrained optimization
###
## We maximize exp(-x^2 - y^2) where x+y = 1
hatf <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## Note: you may prefer exp(- theta %*% theta) instead
}
## use constraints: x + y = 1
A <- matrix(c(1, 1), 1, 2)
B <- -1
res <- maxNR(hatf, start=c(0,0), constraints=list(eqA=A, eqB=B),
             control=list(printLevel=1))
print(summary(res))
      
      
      



xspace <- seq(-10, 10, length.out = 1000)
c = seq(0.001, 20, length.out = 20)      
CauchyGain <- function(x, c) {
  return (c / 2 * log(1  + (x / c) ^ 2 )) 
}      
plot(xspace, CauchyLoss(xspace, tail(c, 1)), type="l")
for (ci  in c) {
  lines(xspace, CauchyLoss(xspace, ci))
}


plot(xspace, CauchyGain(xspace, head(tail(c, 3), 1)), type="l")
plot(xspace, CauchyGain(xspace, .1), type="l")
plot(xspace, CauchyGain(xspace, .03), type="l")
plot(xspace, CauchyGain(xspace, .004), type="l")
plot(xspace, CauchyGain(xspace, .100), type="l")










      
      