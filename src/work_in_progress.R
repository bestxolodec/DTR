

xspace > 0 
xpos = xspace[xspace > 0]
plot(xpos, log(1 + xpos^2) + exp(xpos))



diff(log(1 + xpos^2) + exp(xpos)) 



X <- read.csv("/tmp/bla.csv")
condition1 <- c(30, 20, 50) 
condition2 <- c(35, 30, 35)
X <- cbind( condition1, condition2 )
rownames(X) <- c( 'choice1', 'choice2', 'choice3' )
print(X)


N = 5000
n=12
s <- lapply(seq(12, 1000, 2), function(n) 
            sapply(seq_len(N),  function(...) { (sum(runif(n)) -  n / 2) / sqrt(n / 12) }))
p.vals <- sapply(s, function(x) shapiro.test(x)$p.value)
qplot(seq_along(s), p.vals)
hist(p.vals, breaks=100)
image(kde2d(seq_along(s), p.vals))
hist(s, breaks=100)
shapiro.test(s)




library(ggplot2)
qplot(displ, hwy, data = mpg, facets = drv ~ fl) 
qplot(hwy, data=mpg, fill=drv)
qplot(drv, hwy, data=mpg, geom="boxplot", color=manufacturer)
qplot(displ, hwy, data = mpg, color = drv, geom=c("point", "smooth"))
qplot(displ, hwy, data=mpg, geom=c("point", "smooth"), 
      facets = . ~ drv)



Obviously, theres a DATA FRAME which contains the data youre trying to plot. Then the AESTHETIC MAPPINGS determine how
| data are mapped to color, size, etc. The GEOMS (geometric objects) are what you see in the plot (points, lines, shapes)
| and FACETS are the panels used in conditional plots. Youve used these or seen them used in the first ggplot2 (qplot)
| lesson.
| There are 3 more. STATS are statistical transformations such as binning, quantiles, and smoothing which ggplot2 applies
| to the data. SCALES show what coding an aesthetic map uses (for example, male = red, female = blue). Finally, the plots
| are depicted on a COORDINATE SYSTEM. When you use qplot these were taken care of for you.




gg  <- ggplot(mpg, aes(displ,hwy)) 
gg + geom_point(size=4,alpha=1/2, aes(color="pink"))   
gg + geom_point(color="pink",size=4,alpha=1/2)


g + geom_point(aes(color = drv)) + 
  labs(title="Swirl Rules!", x="Displacement",  y="Hwy Mileage") 



ggplot(mpg, aes(x=displ, y=hwy, color=factor(year))) + 
  geom_point() + 
  facet_grid(drv~cyl, margins = F)
  

qplot(price,data=diamonds,geom="density", color=cut)  
qplot(price, data=diamonds, binwidth=10, fill=cut)


# Как посмотреть на kde в каждом из facet ? 
g + geom_point(alpha=1/3) + facet_grid(cut~car2) + geom_smooth(method="lm",size=1,color="pink")
g + geom_density(alpha=1/7) + facet_grid(cut~car2)  + geom_smooth(method="density")
g + geom_density()

ggplot(diamonds,aes(depth,price)) + geom_point(alpha=1/3) + 
  facet_grid(cut~car2, margins = F) + 
  geom_smooth(method="lm",size=1,color="pink")


ggplot(diamonds,aes(depth,price, fill=..density..)) + geom_point(alpha=1/3) + 
  facet_grid(cut~car2, margins = F) + 
  geom_smooth(method="lm",size=1,color="pink")


d <- ggplot(diamonds, aes(carat, price, fill = ..density..)) +
  xlim(0, 2) + stat_binhex(na.rm = TRUE) + theme(aspect.ratio = 1)
d + facet_wrap(~ color)
d


ggplot(diamonds,aes(depth,price)) + geom_density2d() + 
  facet_grid(cut~car2, margins = T)

ggplot(diamonds, aes(carat, count, fill = cut)) +
  geom_density(position = "fill")


qplot(displ, hwy, data=mpg,color=factor(year)) + facet_grid(cyl~class)



p <- ggplot(diamonds, aes(x = price)) + 
  geom_density(aes(fill = "epanechnikov"), kernel = "rectangular") + 
  facet_grid(~cut) + 
  ggtitle("Kernel density estimate with Facets")
print(p)



m <- ggplot(faithful, aes(x = eruptions, y = waiting)) +
  geom_point() +
  xlim(0.5, 6) +
  ylim(40, 110)
m + geom_density2d()




get.book.prob <- function(N, n, k) {
  sum = 0 
  for (l in seq(k, N)) {
    sum = sum + (-1)^(l-k) * nCm(l, k) *  nCm(N, l) * (1 - l / N)^n
    print(sum)
    }
  return (sum)
}


get.book.prob(10, 7, 5)
get.my.sol.prob(10, 7, 5)
get.my.sol.prob(3, 3, 1)
get.book.prob(3, 3, 1)



N = 3
n = 3


# стр 16 в Notes Plus, записи про Лагутина
get.my.sol.prob <- function(N, n, k) {
  return (nCm(N, k) * nCm(n, N-k) * (N-k)  * factorial(n - 1) /  N^n)
}


qplot(x, y)
qplot(x, y, geom = c("point", "smooth"))
qplot(x, y) + geom_smooth(method = "lm")

qplot(carat, price, data = diamonds, geom = c("point", "smooth"), method = "lm")

cut(y, breaks=5)
table(y, cut(y, breaks=5))

c(by(y, cut(y, breaks=5), mean))



qplot(newx, newy) + geom_smooth(method = "lm") 
qplot(x, y) + geom_smooth(method = "lm")

qplot( geom = c("point", "smooth"), method="lm")


lm(y ~ x)
lm(newy ~ newx)




size = 100
breaks = 10
x = 10 * runif(size)
y = x + rnorm(size)
plot(isoreg(x, y))
newy = tapply(y, cut(x, breaks=breaks), mean)
newx = tapply(x, cut(x, breaks=breaks), mean)
lines(isoreg(newx, newy), col="green")
plot(isoreg(newx, newy), col="green")
s = rbinom(size, 1, prob = seq(0,1, length.out = size) )
plot(s)
plot(isoreg(s ~ x))
news =  tapply(s, cut(x, breaks=breaks), mean)
lines(isoreg(news ~ newx), col="green")





rqmodel = rq.lasso.fit(train$X[index,],train$A[index], .5,method="lasso",weights=train$weight[index],lambda = 1)


cor(c(1,2), c(2,3), method = "spearman")








# Warfarin  ---------------------------------------------------------------





wf <- read.csv("./data/ipwc_data.csv", dec = ",")
wf$Age <- as.numeric(wf$Age)
wf$Penalty <- abs(wf$INR.on.Reported.Therapeutic.Dose.of.Warfarin - 2.5)

# breaks = grep("Project", wf$Comments.regarding.Project.Site.Dataset)
# wf$Project <- as.numeric(cut(1:nrow(wf), breaks))

f <- function(x, dose) {
  if (is.factor(x)){
    NA
  } else {
    spear <- cor(dose, x, method="spearman", use = "pairwise.complete.obs" )
    pears <- cor(dose, x, use = "pairwise.complete.obs" )
    # list(diff=spear-pears, s=spear, p=pears)
    spear - pears
  }
}
cors <- sapply(warfarin, function(x) f(x, warfarin$Therapeutic.Dose.of.Warfarin))
sort(sapply(cors[!is.na(cors)], function(x) x$diff))

by(wf, wf$Project.Site, dim)
sort(table(wf$Project.Site))




g <- function(dfs) {
  sapply(dfs, function(col) f(col, dfs$Therapeutic.Dose.of.Warfarin))
} 
res <- as.data.frame(do.call(rbind, by(wf, wf$Project.Site, g,  simplify = F)))

res.filtered <- res[, sapply(res, function(x) any(x, na.rm = T))]
# age has the biggest difference in spearmen - pearson 
which.max(sapply(res.filtered, function(x) max(x, na.rm = T)))

which.max(res.filtered$Age)


library(ggplot2)

table(cut(wf$INR.on.Reported.Therapeutic.Dose.of.Warfarin[wf$Project.Site == 3], 
          breaks = c(2, 3)), exclude = NULL)




ggplot(data=wf, aes(Age, Therapeutic.Dose.of.Warfarin)) + 
  facet_wrap(~Project.Site,  ncol = 5) + 
  geom_point() 

ggplot(data=wf, aes(Weight..kg., log(Therapeutic.Dose.of.Warfarin))) + 
  facet_wrap(~Project.Site,  ncol = 4) + 
  geom_point()


ggplot(data=wf, aes(Weight..kg., log(Therapeutic.Dose.of.Warfarin), col=as.factor(Age))) + 
  facet_wrap(~Project.Site,  ncol = 4) + 
  geom_point()

ggplot(data=wf, aes(Weight..kg., INR.on.Reported.Therapeutic.Dose.of.Warfarin, 
                    col=as.factor(Age))) + 
  facet_wrap(~Project.Site,  ncol = 4) + 
  geom_point()

ggplot(data=wf, aes(INR.on.Reported.Therapeutic.Dose.of.Warfarin, 
                    Therapeutic.Dose.of.Warfarin, col=as.factor(Age))) + 
  facet_wrap(~Project.Site,  ncol = 4) +
  coord_cartesian(ylim=c(0, 75)) +
  geom_point()

ggplot(data=wf, aes(Penalty, 
                    Therapeutic.Dose.of.Warfarin, col=Weight..kg., alpha=0.8)) + 
  facet_wrap(~Project.Site,  ncol = 7) +
  coord_cartesian(ylim=c(0, 75)) +
  geom_point()

ggplot(data=wf, aes(Penalty,  
                    Therapeutic.Dose.of.Warfarin, col=Height..cm.)) + 
  facet_wrap(~Project.Site,  ncol = 7) +
  coord_cartesian(ylim=c(0, 75)) +
geom_point(aes(alpha = 0.7))



ggplot(data=wf, aes(INR.on.Reported.Therapeutic.Dose.of.Warfarin, log(Therapeutic.Dose.of.Warfarin))) + 
  facet_wrap(~Project.Site,  ncol = 4) +
  geom_point()


ggplot(data=wf, aes(Weight..kg., Therapeutic.Dose.of.Warfarin)) + 
  facet_wrap(~Project.Site,  ncol = 7) +
  coord_cartesian(ylim=c(0, 100)) +
  geom_point(aes(col=factor(Age)))


plot(density(wf$Therapeutic.Dose.of.Warfarin, na.rm = T))
lines(density(wf$Weight..kg., na.rm = T))
lines(density(wf$Height..cm., na.rm = T))










# Importance info for dose prediction -------------------------------------



library(caret)
library(gbm)
gbm.fit <- train(Therapeutic.Dose.of.Warfarin ~ . 
                 - INR.on.Reported.Therapeutic.Dose.of.Warfarin 
                 - PharmGKB.Subject.ID -PharmGKB.Sample.ID 
                 - Project.Site, 
                 subset=!is.na(wf$Therapeutic.Dose.of.Warfarin), data=wf, method="gbm", na.action=na.pass)
gbm.imp <- varImp(gbm.fit)
gbm.imp

# only 20 most important variables shown (out of 3843)
# Overall
# Weight..kg.                                                          100.000
# Age                                                                   89.281
# VKORC1.1173.consensusT/T                                              69.425
# VKORC1..1639.consensusG/G                                             66.480
# Height..cm.                                                           36.177
# VKORC1.genotype...1639.G.A..3673...chr16.31015190..rs9923231..C.TG/G  31.355
# VKORC1.genotype..1173.C.T.6484...chr16.31012379..rs9934438..A.GC/C    31.276
# VKORC1.1542.consensusG/G                                              23.306
# Indication.for.Warfarin.Treatment1; 2                                 20.388
# Amiodarone..Cordarone.                                                17.212
# Ethnicity..Reported.North Africa                                      16.797
# Current.Smoker                                                        15.388
# CYP2C9.consensus*1/*3                                                 14.720
# VKORC1.genotype..1542G.C..6853...chr16.31012010..rs8050894..C.GG/G    14.297
# Carbamazepine..Tegretol.                                              12.965
# Cyp2C9.genotypes*1/*3                                                 12.181
# Cyp2C9.genotypes*1/*2                                                  9.845
# VKORC1.genotype..2255C.T..7566...chr16.31011297..rs2359612..A.GT/T     9.106
# VKORC1.2255.consensusT/T                                               8.701
# VKORC1.genotype..1173.C.T.6484...chr16.31012379..rs9934438..A.GT/T     8.495



wf3 <- subset(wf, subset = Project.Site==3)
wf3$Project.Site <- NULL
wf3$Comments.regarding.Project.Site.Dataset <- NULL


wf1 <- subset(wf, subset = Project.Site==1)
wf1$Project.Site <- NULL
wf1$Comments.regarding.Project.Site.Dataset <- NULL


ggplot(data=wf3, aes(Weight..kg., Therapeutic.Dose.of.Warfarin)) + 
  facet_wrap(~Age,  ncol = 7) +
  coord_cartesian(ylim=c(0, 100)) +
  geom_point(aes(col=factor(VKORC1.1173.consensus)))

ggplot(data=wf, aes(Weight..kg., Therapeutic.Dose.of.Warfarin)) + 
  facet_wrap(~VKORC1.1173.consensus,  ncol = 7) +
  coord_cartesian(ylim=c(0, 150)) +
  geom_point(aes(col=Age))

ggplot(data=wf, aes(Weight..kg., Therapeutic.Dose.of.Warfarin)) + 
  facet_wrap(~Project.Site,  ncol = 7) +
  coord_cartesian(ylim=c(0, 150)) +
  geom_point(aes(col=VKORC1.1173.consensus)) + 
  geom_smooth()


ggplot(data=wf, aes(Weight..kg., Therapeutic.Dose.of.Warfarin)) + 
  facet_wrap(~Project.Site,  ncol = 7) +
  coord_cartesian(ylim=c(0, 150)) +
  geom_point(aes(col=VKORC1..1639.consensus)) + 
  geom_smooth()



# VKORC1.1173.consensus 
tapply(wf$Therapeutic.Dose.of.Warfarin, wf$VKORC1.1173.consensus, 
       function(x) mean(x, na.rm = T))
# C/C      C/T      T/T 
# 42.93582 31.47737 20.56822 

lm.model <- lm(Therapeutic.Dose.of.Warfarin ~ Weight..kg. + VKORC1.1173.consensus + Age, data=wf)
lm.model


library(MASS)

weight <- na.omit(wf$Weight..kg.)
plot(density(weight))
x <- seq(min(weight)-10, max(weight)+10)
# est.m <- mean(weight, na.rm = T)
# est.sd <- sd(weight., na.rm = T)
# lines(x, dnorm(x, est.m, est.sd), col=3)
fit <- fitdistr(weight, "gamma")
# shape           rate    
# 13.879328750    0.178275988 
# ( 0.263566085) ( 0.003447252)
lines(x, dgamma(x, shape = fit$estimate[1], rate = fit$estimate[2]), col=2)

reward <- na.omit(wf$INR.on.Reported.Therapeutic.Dose.of.Warfarin)
plot(density(reward))

plot(wf$Therapeutic.Dose.of.Warfarin, wf$INR.on.Reported.Therapeutic.Dose.of.Warfarin)



dose <- na.omit(wf$Therapeutic.Dose.of.Warfarin)
plot(density(dose))
x <- seq(min(dose), max(dose))
fit <- fitdistr(dose, "gamma")
# shape         rate    
# 3.818609297   0.123268957 
# (0.069682506) (0.002403988)
lines(x, dgamma(x, shape = fit$estimate[1], rate = fit$estimate[2]), col=2)



gbm.reward.model <- train(INR.on.Reported.Therapeutic.Dose.of.Warfarin  
                            ~ Therapeutic.Dose.of.Warfarin + Weight..kg., 
                 subset=!is.na(wf$INR.on.Reported.Therapeutic.Dose.of.Warfarin), 
                 data=wf, method="gbm", na.action=na.pass)
plot(gbm.reward.model)







# Determining fopt(X) -----------------------------------------------------

which.optimal <- (wf$Penalty < 0.5)   & (!is.na(wf$Penalty))
wf.subset  <- wf[which.optimal, ]
wf.subset <- wf

wf.subset$Penalty.bins <- cut(wf.subset$Penalty, 
                              c(0, 0.1, 0.2, 0.3, 0.5, 0.9,  5), exclude = NULL)
ggplot(data=wf.subset, aes((Weight..kg.), (Therapeutic.Dose.of.Warfarin))) +
  geom_point(aes(col=(Penalty))) + 
  facet_wrap( ~ Penalty.bins)

ggplot(data=wf.subset, aes((Age), (Therapeutic.Dose.of.Warfarin))) +
  geom_point(aes(col=(Penalty))) + 
  facet_wrap( ~ Penalty.bins)


plot(wf$Weight..kg., wf$Therapeutic.Dose.of.Warfarin,






# Predicting reward  ------------------------------------------------------


wf.without.na.in.reward <- wf[!is.na(wf$Penalty), ]
  
gbm.fit.reward <- train(Penalty ~ . 
                           - INR.on.Reported.Therapeutic.Dose.of.Warfarin 
                           - PharmGKB.Subject.ID - PharmGKB.Sample.ID 
                           - Project.Site
                           - Comments.regarding.Project.Site.Dataset ,  
                        data=wf.without.na.in.reward, method="gbm", 
                        na.action=na.pass, tuneLength = 2)

(gbm.imp.reward <- varImp(gbm.fit.reward))

# Overall
# Race..Reported.Japanese                                                                     100.000
# Subject.Reached.Stable.Dose.of.Warfarin                                                      24.202
# Estimated.Target.INR.Range.Based.on.Indication1.7-3.3                                        19.802
# Valve.Replacement                                                                            15.848
# Indication.for.Warfarin.Treatment4                                                           14.713
# Target.INR                                                                                   10.815
# Therapeutic.Dose.of.Warfarin                                                                  8.397
# Race..Reported.Caucasian                                                                      7.173
# ComorbiditiesValve                                                                            5.863
# ComorbiditiesAtrial Fibrillation/Flutter                                                      4.961
# Ethnicity..Reported.White                                                                     4.242
# Gendermale                                                                                    4.137
# Medicationsnot cordarone; not pacerone; not amiodarone; not celebrex; not vioxx; not bextra   3.799
# Age                                                                                           3.736
# Indication.for.Warfarin.Treatment8                                                            3.202
# Race..Reported.Han Chinese                                                                    2.655
# Weight..kg.                                                                                   2.550
# Estimated.Target.INR.Range.Based.on.Indication2-3.5                                           1.913
# Race..OMB.White                                                                               1.309
# Race..Reported.Chinese                                                                        1.110










# Deisgn X  ---------------------------------------------------------------


train <- GetSimulationData(10000, scenario = "shvechikov.2", sd=0.1, noise=F)

curve(Shvechikov.1.fopt, from=0, to=100)
curve(Shvechikov.2.fopt, from=0, to=100)
with(train, plot(covariates, optimal.treatment, 
                 pch=19, cex=1.2, col=rgb(0,0.5,0.5, 0.03)))
with(train, plot(covariates, treatment - optimal.treatment, 
                 pch=19, cex=1.2, col=rgb(0,0.5,0.5, .03)))
with(train, plot(covariates, abs(treatment - optimal.treatment)**2, 
                 pch=19, cex=1.2, col=rgb(0,0.5,0.5, 0.03)))
with(train, plot(covariates, GetQFunctionValues(covariates, treatment, optimal.treatment), 
                 pch=19, cex=1.2, col=rgb(0,0.5,0.5, 0.03)))
# then for each Q.value there applied gauss noise with mean = Q.value and sd=var


library(ggplot2)

train <- GetSimulationData(1000, scenario = "shvechikov.2", sd=0.0)
granularity <- 4
levels <- 1:granularity / granularity
d <- with(train, data.frame(reward=reward, treatment=treatment, covariates=covariates, 
                            treat.bins=cut(treatment, c(0, quantile(treatment, levels)), include.lowest = T), 
                            rew.bins=cut(reward, c(0, quantile(reward, levels)), include.lowest = T), 
                            cov.bins=cut(covariates, c(0, quantile(covariates, levels)), include.lowest = T)))

ggplot(d, aes(covariates, treatment, col=reward)) +  geom_point() + geom_smooth() + facet_wrap(~rew.bins) +
   scale_color_gradient(low="red",  high="green")


ggplot(d, aes(treatment, reward, col=covariates)) + geom_point() + geom_smooth() + facet_wrap(~cov.bins)
ggplot(d, aes(covariates, reward, col=reward)) + geom_point() + geom_smooth() + facet_wrap(~treat.bins)
ggplot(d, aes(covariates, treatment, col=reward)) +  geom_point() + geom_smooth() + facet_wrap(~rew.bins) + geom_density2d(col="black")




train <- GetSimulationData(200, scenario = "shvechikov.1")
d <- with(train, data.frame(reward=reward, treatment=treatment, covariates=covariates, 
                            treat.bins=cut(treatment, c(0, quantile(treatment, levels)), include.lowest = T), 
                            rew.bins=cut(reward, c(0, quantile(reward, levels)), include.lowest = T), 
                            cov.bins=cut(covariates, c(0, quantile(covariates, levels)), include.lowest = T)))
library(lattice)
for (i in seq(-180, 180, 3)) {
  print(i)
  # x - treatment;  y - covariates; z - reward
cl <- cloud(reward ~  treatment * covariates, data=d, screen = list(z = i, x = -60, y = 0))
  print(cl)
  Sys.sleep(0.2)
}












# GP regression  ----------------------------------------------------------

library(tgp)
library(ggplot2)
library(data.table)
library(gridExtra)
library(lattice)



test <- GetSimulationData(n_samples, scenario = "shvechikov.2", sd=0.02)
A.grid <- seq(0, 1, length.out = 100) 
ZZ <- expand.grid(seq(0,1,length.out = 100), A.grid)

gp_bgp <- bgp(X, Y)
par(mfrow=c(1,1))


m <- bgp(X, Y, XX=ZZ)
m <- models[["bgp"]]
plot(m)
levelplot(matrix(m$ZZ.q, nrow=100),  xlab="C",   ylab="A") # q.95 - q.05
levelplot(matrix(m$ZZ.s2, nrow=100), col.regions = gray(0:100/100), xlab="C",   ylab="A") # predictive variance

levelplot(matrix)


GetBestPredictions <- function(gp_model, s = 1) {
  means <- gp_model$ZZ.mean
  stds <- sqrt(gp_model$ZZ.s2)
  dt <- data.table(gp_model$XX) 
  dt[, LB:= means - s * stds]
  return(dt[, .(est_A=A[which.max(LB)], est_LB=max(LB)), keyby=C])
}
list_of_figures <- list()
for (s in 0:8) {
  preds <- GetBestPredictions(m, s)
  preds[, optimal_A := test$GetOptimalTreatment(C*100)]
  gg <- ggplot(preds, aes(x=C, y=est_A)) + geom_point() + geom_smooth() +
          geom_line(aes(C, optimal_A, col="red")) + ggtitle(paste("Substract ", s, "deviations"))
  print(gg)
  Sys.sleep(1)
  list_of_figures[[length(list_of_figures) + 1]] <- gg 
}
s <- arrangeGrob(grobs=list_of_figures)
grid.arrange(s)




# TGP approach  -----------------------------------------------------------


# No noise, 100 observations
n_samples = 10
train <- GetSimulationData(n_samples, scenario = "shvechikov.2", sd=0)
X <- with(train, data.frame(C=covariates/100, A=treatment))
Y <- train$reward

test <- GetSimulationData(n_samples, scenario = "shvechikov.2", sd=0)
A.grid <- seq(0, 1, length.out = 100) 
ZZ <- with(train, expand.grid(covariates / 100, A.grid))
# ZZ <- expand.grid(seq(0,1,length.out = 100), A.grid)



####  MEAN WITH S*\SQRT(VAR)  APPROACH  #######

func_names <- c("blm", "btlm", "bcart", "bgp", "bgpllm", "btgp", "btgpllm")
models <- lapply(func_names, function(f_name) do.call(f_name, list(X, Y, ZZ)))
names(models) <- func_names

GetBestPredictions <- function(gp_model, s = 1) {
  means <- gp_model$ZZ.mean
  stds <- sqrt(gp_model$ZZ.s2)
  dt <- data.table(gp_model$XX) 
  dt[, LB:= means - s * stds]
  return(dt[, .(est_A=A[which.max(LB)], est_LB=max(LB)), keyby=C])
}

GetValueOfPredictedA <- function(best.estim.dt, test.obj)  
  best.estim.dt[, mean(test.obj$GetQFunctionValues(C, est_A))]

GetTunedBestPredictions <- function(model, test.obj=test) {
  seq.of.s <- seq(0, 1000, length.out = 100)
  seq.of.Q.vals <- sapply(seq.of.s, function(s) {
    result.dt <- GetBestPredictions(model, s)
    GetValueOfPredictedA(result.dt, test.obj)
  }) 
  plot(seq.of.Q.vals)
  best.s <- seq.of.s[which.max(seq.of.Q.vals)]
  best.result.dt <-  GetBestPredictions(model, s=best.s)
  return (list(m_pred=best.result.dt, best.s=best.s))
}

m <- bgp(X, Y, XX=X)
m$ZZ.s2
m$ZZ.ks2

p <- GetTunedBestPredictions(m, test)
p$best.s
ggplot(p$m_pred, aes(C, est_A)) + geom_smooth() + geom_point() + geom_line(aes(C, test$GetOptimalTreatment(C*100)), col="red")


m_pred <- GetBestPredictions(m, s=10)
m_pred <- GetBestPredictions(m, s=3)
ggplot(m_pred, aes(C, est_A)) + geom_smooth() + geom_point() + geom_line(aes(C, test$GetOptimalTreatment(C*100)))

plot(m)


predictions <- list()
best.s <- list()
for (m_name in names(models)) {
  tuned.res.dt <- GetTunedBestPredictions(models[[m_name]])
  predictions[[1]]  <- tuned.res.dt$m_pred$C
  predictions[[length(predictions) + 1]]  <- m_pred$est_A
  best.s[[length(best.s) + 1]]  <- tuned.res.dt$best.s
}
dt <- as.data.table(predictions)
names(dt) <- c("C", func_names)
dt.melted <- melt(dt, id.vars = "C", variable.name = "Algo",  value.name = "est_A")
ggplot(dt.melted, aes(C, est_A)) + geom_point() + geom_smooth() + facet_wrap(~ Algo, nrow = 2)
ggplot(dt.melted, aes(C, est_A, col=Algo)) + geom_smooth()
dt.melted[, optimal_A := test$GetOptimalTreatment(C*100)]
ggplot(dt.melted, aes(x=C, y=est_A)) + geom_point() + geom_smooth() + 
  geom_line(aes(C, optimal_A, col="red")) + facet_wrap(~ Algo, nrow = 2)
dt[, c("C", "est_A") := .(C, bgp)]

Q.vals  <- sapply(dt[, !"C", with=F],  f)
barplot(Q.vals)


plot(models$blm)
plot(models$btlm)
plot(models$bcart)
plot(models$bgp)
plot(models$bgpllm)
plot(models$btgp)
plot(models$btgpllm)





####  QUANTILE   APPROACH #######


par(mfrow=c(1,1))
head(gp_bgp_with_test$ZZ.mean)
head(gp_bgp_with_test$ZZ.km)
hist(gp_bgp_with_test$ZZ.mean  - gp_bgp_with_test$ZZ.km, breaks = 100)

# ZZ.mean	 Vector of mean predictive estimates at XX locations


gp_r_xa <- GP_fit(X, Y, corr=list(type="exponential",power=2))



































































 





library(kernlab) # gausspr
gp_gausspr <- gausspr(X, Y)
str(gp_gausspr, max.level = 2)
gp_gausspr@alpha






library(gptk)





# GPfit -------------------------------------------------------------------



# 1. Провести быстрое моделирование разнеымми методами (всеми)
# 2. Численно посмотреть, что работает лучше. 




library(GPfit)


library(GPfit)
# No noise, 100 observations
train <- GetSimulationData(100, scenario = "shvechikov.1", sd=0)
test <- GetSimulationData(100, scenario = "shvechikov.1", sd=0)
X <- with(train, data.frame(covariates/100, treatment))
Y <- train$reward
gp_r_xa <- GP_fit(X, Y, corr=list(type="exponential",power=2))
plot(gp_r_xa)
plot(gp_r_xa,  surf_check = TRUE)
plot(gp_r_xa,  surf_check = TRUE, screen=c(x=-50))



heatmap(GPfit::corr_matrix(X, gp_r_xa$beta, corr=list(type="exponential",power=2)))


K <- corr_matrix(gp_obj$X, gp_obj$beta, corr=gp_obj$correlation_param)
  
GetParamsFromGPobject <- function(gp_obj)  {
  K <- corr_matrix(gp_obj$X, gp_obj$beta, corr=gp_obj$correlation_param)
  A <- solve(K + gp_obj$sig2 * diag(NROW(K)))
  A <- solve(K)
  y <- gp_obj$Y
  stopifnot(gp_obj$correlation_param$type == "exponential")
  GetCovCorr <- function(c_test, ci) {
    beta_index = pmatch("covariat",  colnames(gp_obj$X))
    stopifnot(length(beta_index) == 1)
    stopifnot(!is.na(beta_index))
    beta <- gp_obj$beta[beta_index]
    return(exp(-10**beta * abs(c_test - ci) ** gp_obj$correlation_param$power))
  }
  GetTreatCorr <- function(a_test, a) {
    beta_index = pmatch("treat",  colnames(gp_obj$X))
    stopifnot(length(beta_index) == 1)
    stopifnot(!is.na(beta_index))
    beta <- gp_obj$beta[beta_index]
    return (exp(-10**beta * abs(a_test - a) ** gp_obj$correlation_param$power))
  }
  return (list(K=K, A=A, y=y, GetCovCorr=GetCovCorr, GetTreatCorr=GetTreatCorr))
}

gp_params <- GetParamsFromGPobject(gp_r_xa)
list2env(gp_params, globalenv())

a_test <- 0.9
c_test <- test$covariates[1, 1] / 100
gp_obj <- gp_r_xa
x_test = cbind(c_test, a_test)

with(gp_params, {
  KC <-  diag(sapply(gp_obj$X[, "covariates.100"] , function(c) GetCovCorr(c_test, c)))  
  KT <-  matrix(sapply(gp_obj$X[, "treatment"] , function(a) GetTreatCorr(a_test, a)), nrow=1)
  as.numeric(KT %*% KC %*% A %*% y)
})

predict(gp_r_xa, x_test)


# test if correlation was correctly found 
corr_matrix(rbind(x_test, gp_r_xa$X[1, ]), gp_r_xa$beta, corr = gp_r_xa$correlation_param)
GetCovCorr(c_test, gp_r_xa$X[1, 1]) *  GetTreatCorr(a_test, gp_r_xa$X[1, 2])

x_test











heatmap(A %*% (K + sigma2 * diag(NROW(K))), Rowv = NA, Colv = NA, revC = T)













# Noise with 100 observations
train <- GetSimulationData(100, scenario = "shvechikov.1", sd=0.05)
X <- with(train, data.frame(covariates/100, treatment))
names(X)
Y <- train$reward
gp_r_xa <- GP_fit(X, Y)
plot(gp_r_xa)
plot(gp_r_xa,  surf_check = TRUE)
plot(gp_r_xa,  surf_check = TRUE, screen=list(z=90, x=-60, y=-70))


# More extreme noise with 100 observations
train <- GetSimulationData(100, scenario = "shvechikov.1", sd=0.1)
X <- with(train, data.frame(covariates/100, treatment))
names(X)
Y <- train$reward
gp_r_xa <- GP_fit(X, Y)
plot(gp_r_xa)
plot(gp_r_xa,  surf_check = TRUE)
plot(gp_r_xa,  surf_check = TRUE, screen=list(z=90, x=-60, y=-70))


x_test <- as.matrix(seq(0,1, length.out = 100))
a_test <-

predict(gp_r_xa, x_test)






curve(dnorm(x, mean=1) + dnorm(x, mean=0), from=0, to=1) 






centers <- train$treatment * 100
# centers <- c(20, 40)
f <- function(lambdas, a) {
  sum(lambdas * exp(-(centers - a)**2 / 2))
}
lambdas <- runif(n=100, max=20) 
seq.treatments <- seq(0, 100, length.out = 1000)
values <- sapply(seq.treatments,  function(x) f(lambdas, x))
# values <- sapply(seq.treatments,  function(x) f(c(1, 5), x))
par(mfrow=c(1,1))
plot(seq.treatments, values, type="l")
for (t in train$treatment) {
  abline(v=t * 100, col="red")
}
guess <- centers[which.max(lambdas)]
# abline(v=guess, col="green", lwd=4) 




  

order(lambdas * centers, decreasing = T)[1]

sort(lambdas * centers)
(lambdas * centers)[98]


sort(lambdas)

alphas
treatments

plot(treatments, values / 8000, type="l")

lines(density(train$treatment * 100)) 


fit <- fitdistr(, "gamma")



a <- c(1, 1, 2, 4)
outer(a,a)
sum(outer(a,a))








install.packages("plotly")
library(plotly)

plot_ly(df, x = ~x, y = ~y, z = ~z) %>% add_markers()



























#################        NEW ERA          ###################

#  KO-Learning on OUR data ------------------------------------------------ 

n_samples <- 100
n_test_samples <- 100
noise.sd <- 0
scenario <- "shvechikov.1"


train <- GetSimulationData(sd=noise.sd, sample.size = n_samples, scenario = scenario )
test <- GetSimulationData(sd=noise.sd, sample.size = n_test_samples, scenario = scenario)
ko_train <- ChangeFormatFromOurToChenEnriched(train)
ko_test <- ChangeFormatFromOurToChenEnriched(test)


res <- list()
# possible problem here - we are making max(0, min(pred, 2))
res$ko <- GetKOLearningValueAndPredictedDose(ko_train, ko_test)

# res$bgp <- GetGPValueAndPredictedDose(train, test, model_name = "bgp")

# data_list <- GetListOfTrainTestData(train, test, eps)
# str(data_list)



model_name <- "bgp"
model_joint <- do.call(model_name, c(data_list, pred.n=FALSE))
plot(model_joint)
plot(model_joint, center="km", as="ks2")

model_map <- do.call(model_name,  with(data_list, list(X, Z, pred.n=FALSE, verb=4)))
out.kp <- predict(model, data_list$XX, krige=use_krige)
str(model_map$X)
str(data_list$XX)

plot(out.kp)
plot(out.kp, center="km", as="ks2")



with(test, plot(test$covariates, test$optimal.treatment))










if (use_MAP) {
  model <- do.call(model_name,  with(data_list, list(X, Z, pred.n=FALSE, verb=4)))
} else {
  model <- do.call(model_name, c(data_list, pred.n=FALSE))
}
if (use_krige) {
  return (list(means=model$ZZ.km, vars=model$ZZ.ks2, model=model))
} else {
  return (list(means=model$ZZ.mean, vars=model$ZZ.s2, model=model))
}
}

preds_on_test <- FitGPAndPredictOnTest(model_name, data_list, use_MAP=F, use_krige=F)
mesh_grid_with_pred <- GetArgmaxTreatmentsAndMaxLowerBounds(preds_on_test, s=s)
pred_value <- GetPredValue(mesh_grid_with_pred, test)
return(list(A_pred=mesh_grid_with_pred$A_pred, Q=pred_value, model=preds_on_test$model))





data.list <- GetCovarsTreatsMeshGrid(train, test, 0.1)
model_name <- "bgp"
system.time ({
  model_mcmcm <-  do.call(model_name, list(X, Y, ZZ))
})
# user      system  elapsed 
# 159.670   5.616   171.408 

system.time ({
  model_mcmcm <-  do.call(model_name, list(X, Y, ZZ, pred.n=F, verb=4))
})

system.time ({
  model_map <-  do.call(model_name, data.list[-3])
  model_map_prediction <-  predict(model_map, data.list[[3]])
})
# user     system elapsed 
# 18.978   1.155  20.899 

model_map_prediction <-  predict(model_map, data.list$XX, R = 0)
model_map_prediction$ZZ.km
model_map_prediction$ZZ.mean

r <-  GetBestPredictions(model_map_prediction)
gp_model <- model_map_prediction
dim(ZZ)

plot((gp_model$ZZ.km - gp_model$ZZ.mean))

plot(gp_model$ZZ.km)
plot(gp_model$ZZ.mean)


means <- gp_model$ZZ.mean
length(means)
stds <- sqrt(gp_model$ZZ.ks2)
length(stds)
dt <- data.table(gp_model$XX)
dt[, LB:= means - s * stds]
col_names <- c(grep('^C', names(dt), value = T))
dt[, .(A_pred=A[which.max(LB)], LB_max=max(LB)), by=col_names]



  
with(test, {
 plot(covariates, res$bgp$A_pred, pch=19)
 plot(covariates, , pch=19, col=2)
})

par(mfrow=c(1,2))
plot(test$covariates,  res$ko$A_pred, main=res$ko$Q)
plot(test$covariates,  res$bgp$A_pred, main=res$bgp$Q)
par(mfrow=c(1,1))


with(test, plot(res$ko$A_pred, covariates))
with(test, plot(res$bgp$A_pred, covariates))




# OUR model on Kosorok data -----------------------------------------------
  


Scenario1Enriched <- function(size,ncov,seed){
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


Scenario2Enriched <- function(size,ncov,seed){
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
  mu =   GetQFunctionValues(X, A)
  R = rnorm(length(mu),mu,1)
  datainfo = list(X=X, A=A, R=R, D_opt=D_opt, mu=mu, 
                  GetQFunctionValues=GetQFunctionValues, 
                  GetOptimalTreatment=GetOptimalTreatment)
  return(datainfo)
}


Scenario4Enriched <- function(size,ncov,seed){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  GetOptimalTreatment <-function(X) {
    I(X[,1] > -0.5)*I(X[,1] < 0.5)*0.6 + 1.2*I(X[,1] > 0.5) + 1.2*I(X[,1] < -0.5) +
    X[,4]^2 + 0.5*log(abs(X[,7])+1) - 0.6
  }
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

n_samples <- 100
n_test_samples <- 200
n_covariates <- 10


GenData <- Scenario4Enriched
ko_train <- GenData(size = n_samples, ncov = n_covariates, seed = 0)
ko_test <- GenData(size = n_test_samples, ncov = n_covariates, seed = 1)
train <- ChangeFormatFromChenEnrichedToOur(ko_train)
test <- ChangeFormatFromChenEnrichedToOur(ko_test)

res_ko_data <- list()
res_ko_data$ko <- GetKOLearningValueAndPredictedDose(ko_train, ko_test)

use_MAP = T; use_krige = T
res_ko_data$bgp_map_not_krige <- GetGPValueAndPredictedDose(train, test, model_name = "bgp", use_MAP = T, use_krige = F)
res_ko_data$bgp_map_krige <- GetGPValueAndPredictedDose(train, test, model_name = "bgp", use_MAP = use_MAP, use_krige = use_krige)
res_ko_data$bgp <- GetGPValueAndPredictedDose(train, test, model_name = "bgp", use_MAP = use_MAP, use_krige = use_krige)

res_ko_data$bgpllm <- GetGPValueAndPredictedDose(train, test,  model_name = "bgpllm", use_MAP = use_MAP, use_krige = use_krige)
res_ko_data$btgp <- GetGPValueAndPredictedDose(train, test, model_name = "btgp", use_MAP = use_MAP, use_krige = use_krige)
res_ko_data$btgpllm <- GetGPValueAndPredictedDose(train, test, model_name = "btgpllm", use_MAP = use_MAP, use_krige = use_krige)

for (n in names(res_ko_data)){
  Q = res_ko_data[[n]]$Q
  print(paste(n, ":  Q = ", Q))
}


# saveRDS(res_ko_data, file = "../data/scenario_1_train_50_teset_200_covars_30")
# saveRDS(res_ko_data, file = "../data/scenario_1_train_50_teset_200_covars_10")
# saveRDS(res_ko_data, file = "../data/scenario_1_train_100_teset_200_covars_10")
# saveRDS(res_ko_data, file = "../data/scenario_1_train_100_teset_200_covars_30")
# saveRDS(res_ko_data, file = "../data/scenario_1_train_200_teset_200_covars_30")
# saveRDS(res_ko_data, file = "../data/scenario_1_train_300_teset_200_covars_30")
# saveRDS(res_ko_data, file = "../data/scenario_1_train_800_teset_200_covars_30")

# saveRDS(res_ko_data, file = "../data/scenario_2_train_50_teset_200_covars_10") 
# saveRDS(res_ko_data, file = "../data/scenario_2_train_100_teset_200_covars_10")  # we are better !
# saveRDS(res_ko_data, file = "../data/scenario_2_train_200_teset_200_covars_10")
# saveRDS(res_ko_data, file = "../data/scenario_2_train_300_teset_200_covars_10")

# saveRDS(res_ko_data, file = "../data/scenario_4_train_100_teset_200_covars_10")  # we loose
# saveRDS(res_ko_data, file = "../data/scenario_4_train_200_teset_200_covars_10")  # we win
# saveRDS(res_ko_data, file = "../data/scenario_4_train_300_teset_200_covars_10")
# saveRDS(res_ko_data, file = "../data/scenario_4_train_800_teset_200_covars_10")



library(gtools) # for mixedsort

res <- list()
for (sc in c("scenario_1, scenario_2, scenario_4")) {
  files <- grep(sc, list.files("../data/"), value = T)
  counter <- 1
  for (f in mixedsort(files)) {
    s <-  readRDS(paste("../data/", f, sep = ""))
    for (n in names(s)){
      res[[sc]][[]]s[[n]]$Q
      print(paste(n, ":  Q = ", s[[n]]$Q))
    }
    for (n in names(s)){
      print(paste(n, ":  Q = ", s[[n]]$Q))
    }
    cat("\n\n")
    counter <- counter
  }
}


files <- grep("scenario", list.files("../data/"), value = T)
for (f in mixedsort(files)) {
  cat(f, "\n")
  s <-  readRDS(paste("../data/", f, sep = ""))
  for (n in names(s)){
    print(paste(n, ":  Q = ", s[[n]]$Q))
  }
  cat("\n\n")
}



res_ko_data$bgp_s0 <- GetGPValueAndPredictedDose(train, test, s=0, model_name = "bgp")
res_ko_data$bgpllm_s0 <- GetGPValueAndPredictedDose(train, test, s=0, model_name = "bgpllm")
res_ko_data$btgp_s0 <- GetGPValueAndPredictedDose(train, test, s=0, model_name = "btgp")
res_ko_data$btgpllm_s0 <- GetGPValueAndPredictedDose(train, test, s=0, model_name = "btgpllm")


plot(ko_test$D_opt - res_ko_data$ko$A_pred)
plot(ko_test$D_opt - res_ko_data$bgp$A_pred)


sum((test$optimal.treatment - res_ko_data$ko$A_pred)**2)
sum((test$optimal.treatment - res_ko_data$bgp$A_pred)**2)

with(res_ko_data, plot(bgp$A_pred, ko$A_pred))


res_ko_data$bgp_map_not_krige <- GetGPValueAndPredictedDose(train, test, model_name = "bgp", use_MAP = T, use_krige = F)

model <- do.call(model_name, with(data_list, list(X, Z, trace=T)))




xspace = seq(from=0,to=1,length.out = 50)
plot(xspace, dbeta(xspace, 1, 10), type="l")
plot(xspace, dbeta(xspace, 2, 1000), type="l")

curve(dbeta(x, 1, .1))








