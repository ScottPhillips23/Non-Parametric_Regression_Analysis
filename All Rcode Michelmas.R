#ECDF calculations

set.seed(1)
n = 10
x <- runif(n)
xseq <- seq(0, 1, length = 100)

Fhat <- ecdf(x)
plot(Fhat, verticals = TRUE, main = 'Uniform n = 10')
lines(xseq, punif(xseq), col = 'blue', lwd = 2)

set.seed(1)
n <- 10
x <- rnorm(n)
xseq <- seq(-3, 3, length = 100)
Fhat = ecdf(x)
plot(Fhat, verticals  = TRUE, main = 'Normal n = 10')
lines(xseq, pnorm(xseq), col = 'red', lwd = 2)

#Bootstrap

bootstrap_variance_median <- function(X, B){
  n <- length(X)
  mhat <- numeric(B)
  for(j in 1:B){
    X.boot <- sample(X, n, replace = TRUE)
    mhat[j] <- median(X.boot)
  }
  var = var(mhat)
  return(list(var = var, mhat = mhat))
}

#Now we implement bootstrap to find estimators

set.seed(1)
n = 500
X <- rgamma(n, shape = 5, rate = 2)
mhat = median(X)
B = 10000
results = bootstrap_variance_median(X, B)
var.boot = results$var
var.boot
mhat.boot = results$mhat
hist(mhat.boot, xlab = 'mhat', main = '', col = 'blue', breaks = 20)
points(mhat, 0, pch = 22, col = 'red', bg = 'red')


#Implemetation of bootstrap to confidence intervals before corrections(
#non normality and bias), try implementing these yourself

#Normal Confidence intervals

alpha = 0.05
se = sqrt(var(mhat.boot))
se
mhat_normalCI = c(mhat - qnorm(1 - alpha/2)*se, mhat + qnorm(1 - alpha/2)*se)
mhat_normalCI

#Percentile Confidence intervals

mhat_percentileCI = quantile(mhat.boot, probs = c(0.025, 0.975))
mhat_percentileCI


#Bias estimation via bootstrap

require(boot)
head(bigcity)
barplot(t(bigcity), beside = TRUE,
        ylab = 'Population (thousands)', 
        xlab = 'cities',
        legend.text = c('1920', '1930'), 
        axisnames = FALSE, 
        col = c('lightblue2', 'dodgerblue3'),
        border = NA,
        args.legend = list(bty = n, border = NA)
        )

pop1920 <- bigcity$u
pop1930 <- bigcity$x
ratiohat <- mean(pop1930)/mean(pop1920)
ratiohat

#Bootstrap for target being ratio of 2 empirical means
set.seed(1)
n = nrow(bigcity)
B = 10000
ratiohat.boot = rep(0, B)
for(j in 1:B){
  ind = sample(n, replace = TRUE)
  pop1920mean.boot <- mean(pop1920[ind])
  pop1930mean.boot <- mean(pop1930[ind])
  ratiohat.boot[j] <- pop1930mean.boot / pop1920mean.boot
}
bias.boot <- mean(ratiohat.boot) - ratiohat
bias.boot
debiased.ratiohat <- ratiohat - bias.boot
debiased.ratiohat




##################
#Linear Regression with simulated data

n = 50
b0 = 1.5
b1 = 0.9
set.seed(1)
epsilon <- rt(n, df = 4)
x <- rnorm(n, 0, 1)
X <- cbind(1,x)
y = b0 + b1*x + epsilon
y


NPBbootstrap_betas <- function(y, X, B){
  n.boot <- length(y)
  beta.hat <- matrix(NA, nrow = B, ncol= 2)
  for(j in 1:B){
    indices = sample(1:n.boot, replace = TRUE)
    X.boot <- X[indices,]
    y.boot = y[indices]
    beta.hat[j,] <- solve(t(X.boot) %*% X.boot) %*% t(X.boot) %*% y.boot
  }
  return(beta.hat)
}

set.seed(1)
beta.boot = NPBbootstrap_betas(y, X, 1000)
head(beta.boot, n = 5)

var.boot <- apply(beta.boot, 2, var)
var.boot

#sample ols estimates
beta.ols <- c(solve(t(X) %*% X) %*% t(X) %*% y)
beta.ols

alpha = 0.05

#normal confidence intervals
se.boot = sqrt(var.boot)
se.boot
normalCI.beta0 = c(beta.ols[1] - qnorm(1-alpha/2) * se.boot[1],
                   beta.ols[1] + qnorm(1-alpha/2) * se.boot[1])
normalCI.beta1 = c(beta.ols[2] - qnorm(1-alpha/2) * se.boot[2],
                   beta.ols[2] + qnorm(1-alpha/2) * se.boot[2])

normalCI <- rbind(normalCI.beta0, normalCI.beta1)
normalCI

#percentile confidence intervals

percentileCI.beta0 <- quantile(beta.boot[,1], probs = c(0.025, 0.975))
percentileCI.beta1 <- quantile(beta.boot[,2], probs = c(0.025, 0.975))

percentileCI <- rbind(percentileCI.beta0, percentileCI.beta1)
percentileCI

#Wald confidence intervals

waldCI = confint(lm(y ~ x))
waldCI

#Plots

par(mfrow = c(1,2))
hist(beta.boot[,1], main = expression(beta[0]), xlab = '')
abline(v = normalCI[1,], lwd = 2.5)
abline(v = percentileCI[1,], col = 2, lwd = 2.5)
abline(v = waldCI[1,], col = 3, lwd = 2.5)
legend('topright', bty = 'n', col = c(1,2,3),
       legend = c('Normal CI', 'Percentile CI', 'Wald CI'), fill = c(1,2,3))
hist(beta.boot[,2], main = expression(beta[0]), xlab = '')
abline(v = normalCI[2,], lwd = 2.5)
abline(v = percentileCI[2,], col = 2, lwd = 2.5)
abline(v = waldCI[2,], col = 3, lwd = 2.5)
legend('topright', bty = 'n', col = c(1,2,3),
       legend = c('Normal CI', 'Percentile CI', 'Wald CI'), fill = c(1,2,3))
par(mfrow = c(1,1))




#Semiparamteric residual bootstrap

n = 20
b0= 1.5
b1 = 0.4
set.seed(1)
epsilon = rt(n, df = 4)
hist(epsilon)
x = rnorm(n, 0, 1)
X = cbind(1, x)
y = b0 + b1*x^3 + epsilon

plot(x, y, pch = 16)
abline(lm(y~x), col = 'blue')

NPBbootstrap_betas <- function(y, X, B){
  n.boot <- length(y) 
  beta.hat <- matrix(NA, nrow = B, ncol= 2)
  for(j in 1:B){
    indices = sample(1:n.boot + 1, replace = TRUE)
    X.boot <- X[indices,]
    y.boot = y[indices]
    beta.hat[j,] <- solve(t(X.boot) %*% X.boot) %*% t(X.boot) %*% y.boot
  }
  return(beta.hat)
}

SRBbootstrap_betas <- function(y, X, B){
  beta.OLS = solve(t(X) %*% X) %*% t(X) %*% y
  e.hat = y - (X %*% beta.OLS)
  n.boot = length(y)
  beta.hat = matrix(NA, nrow = B, ncol = 2)
  for(j in 1:B){
    indices = sample(1:n.boot, replace = TRUE)
    e.boot = e.hat[indices]
    y.boot = X %*% beta.OLS + e.boot
    beta.hat[j,] = solve(t(X) %*% X) %*% t(X) %*% y.boot
  }
  return(beta.hat)
}

set.seed(1)
NPBbeta.boot = NPBbootstrap_betas(y, X, B = 1000)
SRBbeta.boot = SRBbootstrap_betas(y, X, B = 1000)

NPBvar.boot = apply(NPBbeta.boot, 2, var)
SRBvar.boot = apply(SRBbeta.boot, 2, var)

NPBvar.boot
SRBvar.boot

#Sample OLS estimates
beta.ols = c(solve(t(X) %*% X) %*% t(X) %*% y)
beta.ols

#Normal CI's

alpha = 0.05
NPB.se.boot = sqrt(NPBvar.boot)
NPB.se.boot
NPB.normalCI.beta0 = c(beta.ols[1] - qnorm(1-alpha/2)*NPB.se.boot[1],
                       beta.ols[1] + qnorm(1-alpha/2)*NPB.se.boot[1])
NPB.normalCI.beta1 = c(beta.ols[2] - qnorm(1-alpha/2)*NPB.se.boot[2],
                       beta.ols[2] + qnorm(1-alpha/2)*NPB.se.boot[2])

NPB.normalCI = rbind(NPB.normalCI.beta0, NPB.normalCI.beta1)

SRB.se.boot = sqrt(SRBvar.boot)
SRB.se.boot
SRB.normalCI.beta0 = c(beta.ols[1] - qnorm(1-alpha/2)*SRB.se.boot[1],
                       beta.ols[1] + qnorm(1-alpha/2)*SRB.se.boot[1])
SRB.normalCI.beta1 = c(beta.ols[2] - qnorm(1-alpha/2)*SRB.se.boot[2],
                       beta.ols[2] + qnorm(1-alpha/2)*SRB.se.boot[2])

SRB.normalCI = rbind(SRB.normalCI.beta0, SRB.normalCI.beta1)


#Percentile confidence intervals
NPB.percentileCI.beta0 = quantile(NPBbeta.boot[,1], probs = c(0.025, 0.975))
NPB.percentileCI.beta1 = quantile(NPBbeta.boot[,2], probs = c(0.025, 0.975))

NPB.percentileCI = rbind(NPB.percentileCI.beta0, NPB.percentileCI.beta1)

SRB.percentileCI.beta0 = quantile(SRBbeta.boot[,1], probs = c(0.025, 0.975))
SRB.percentileCI.beta1 = quantile(SRBbeta.boot[,2], probs = c(0.025, 0.975))

SRB.percentileCI = rbind(SRB.percentileCI.beta0, SRB.percentileCI.beta1)

waldCI <- confint(lm(y~x))

round(rbind(waldCI, NPB.normalCI, SRB.normalCI, 
            NPB.percentileCI, SRB.percentileCI), 3)




#Bayesian Bootstrap

install.packages('MCMCpack')
library(MCMCpack)

sample(x, length(x), replace = TRUE)

#Conventional Bootstrap
mean.cb <- function(data, B){
  boot.mean = numeric()
  for(j in 1:B){
    boot.mean[j] = mean(sample(data, length(data), replace = TRUE))
  }
  return(boot.mean)
}

#Bayesian Bootstrap

mean.bb = function(data, B, alpha){
  pi_mat = rdirichlet(B, rep(1, length(data)))
  bayes.boot.mean = apply(pi_mat, 1, weighted.mean, x = data)
  return(bayes.boot.mean)
}

set.seed(1)
n1 = 10
x1 = rnorm(n1, 0, 1)
n2 = 100
x2 = rnorm(n2, 0, 1)

cb1 <- mean.cb(data = x1, B = 1000)
bb1_flat <- mean.bb(data = x1, B = 1000, alpha = 1)
bb1_zero <- mean.bb(data = x1, B = 1000, alpha = 0.001)

cb2 <- mean.cb(data = x2, B = 1000)
bb2_flat <- mean.bb(data = x2, B = 1000, alpha = 1)
bb2_zero <- mean.bb(data = x2, B = 1000, alpha = 0.001)


par(mfrow = c(2,3))



#Dirichlet process

set.seed(1)
alpha = 1
G_0 <- function(n){
  return(rnorm(n, 0, 1))
}
n = 10
b = rbeta(n, 1, alpha)
p = numeric(n)
p[1] = b[1]
p[2:n] = sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
m = G_0(n)
theta = sample(m, prob = p, replace = TRUE)
round(theta, 3)

#Dirichlet process's are so fucking cool!!!!! They make so much sense now!!!

#Application to synthetic data now

set.seed(1)
n = 50
y = rpois(n, 1.5)
observed.dist = table(y)
observed.dist

#Now using observed.dis as data, pretend we dont know its poisson

alpha = 1
G_0 <- function(n) rpois(n, 1)
p <- numeric(n)
b <- rbeta(n, 1, alpha)

p[1] <- b[1]
p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
m = G_0(n)
theta <- sample(m, prob = p, replace = TRUE)

table(theta)

row.names(observed.dist)

indicator <- as.numeric(row.names(observed.dist))
indicator

prior.dist <- table(row.names(observed.dist)) - 1
prior.dist

for(i in indicator + 1){
  prior.dist[[i]] = sum(theta==i-1)
}

prior.dist

plot(observed.dist, col = 4, lwd = 4, bty = 'n',
     ylab = 'Counts', ylim = c(0, max(observed.dist, prior.dist)))
lines(prior.dist, type = 'l', lty = 3, col = 1, lwd = 3)


avg.theta <- function(repetitions, n, alpha){
  for(j in 1:repetitions){
    p = numeric(n)
    beta <- rbeta(n, 1, alpha)
    p[1] <- beta[1]
    p[2:n] <- sapply(n-1, function(i) b[i] * prod(1 - b[1:(i-1)]))
    m = G_0(n)
    theta <- rbind(theta, sample(m, prob = p, replace = TRUE))
  }
  return(atheta)
}
set.seed(1)
avg.theta(10, 50, 1)




#Sampling from a DP posterior


#First part is as before

set.seed(1)
n = 50
y = rpois(n, 1.5)
observed.dist = table(y)
observed.dist

alpha.prior <- 1

G_0.prior <- function(n) rpois(n, 1)
p.prior <- numeric(n)
b.prior <- rbeta(n, 1, alpha.prior)
p.prior[1] <- b.prior[1]
p.prior[2:n] <- sapply(2:n, function(i) b.prior[i] * prod( 1 - b.prior[1:(i-1)]))
m = G_0.prior(n)
theta.prior <- sample(m, prob = p.prior, replace = TRUE)
table(theta.prior)

#Sampling from posterior
alpha.post = alpha.prior + n

p.post <- numeric(n)
b.post <- rbeta(n, 1, alpha.post)
p.post[1] <- b.post[1]
p.post[2:n] <- sapply(2:n, function(i) b.post[i] * prod(1 - b.post[1:(i-1)]))

y.star <- sample(y, prob = rep(1/n, n), replace = TRUE)

w1 <- alpha.prior/alpha.post
w2 <- n/alpha.post

M <- matrix(c(theta.prior, y.star), byrow = TRUE, nrow = 2)
m.star <- numeric(n)
for(i in 1:length(m.star)){
  m.star[i] <- sample(M[,i], size = 1, prob = c(w1,w2))
}
theta.post <- sample(m.star, prob = p.post, replace = TRUE)


indicator <- as.numeric(row.names(observed.dist))
post.dist <- table(row.names(observed.dist)) - 1
for(i in indicator+1){
  post.dist[[i]] = sum(theta.post == i-1)
}

plot(observed.dist, col = 4, lwd = 4, bty = 'n', ylab = 'counts', 
     ylim = c(0, max(observed.dist, post.dist)))
lines(post.dist, type = 'l', lty = 3, col = 1, lwd = 3)


#The Dirichlet process mixture model

set.seed(1)
n1 = 40
n2 = 20
y1 = rnorm(n1, 0, 1)
y2 = rnorm(n2, 4, 1)
y = c(y1, y2)


####################################################################################################

#Practical 1

#Dataset

data()
data(package = 'MASS')

data('iris')
class(iris)

dim(iris)

head(iris, n = 5)

y = iris[, 1:4]
n = nrow(y)

par(mfrow = c(2,2))
for(i in 1:4){
  hist(y[,i], col = i+1, main = '', xlab = names(y)[i])
}


#Bootstrap implementation

sepal.length = y[,1]
n

median.hat = median(sepal.length)
median.hat

bootstrap_variance_median <- function(X, B){
  n = length(X)
  mhat = numeric(B)
  for(i in 1:B){
    X.boot <- sample(X, n, replace = TRUE)
    mhat[i] <- median(X.boot)
  }
  var <- var(mhat)
  return(list(var = var, mhat = mhat))
}

B = 10000
set.seed(1)
results = bootstrap_variance_median(sepal.length, B)
var.boot <- results$var
var.boot

mhat.boot = results$mhat

hist(mhat.boot, xlab = 'Bootstrap samples of the median.', 
     main = '', col ='lightblue2', breaks = 20)

points(median.hat, 0, pch = 22, col = 'red', bg = 'red')


#Computing manual bootstrap confidence intervals

#3.1 Constructing Normal confidence intervals

normal.CI <- c(5.8 - qnorm(1-0.01/2) * sqrt(var.boot), 5.8 + qnorm(1-0.01/2) * sqrt(var.boot))

round(normal.CI, 1)

#3.2 Percentile confidence intervals

percentile.CI <- quantile(results$mhat, probs = c(0.005, 0.995))
percentile.CI

#3.3

set.seed(1)
B1 = 1000
mhat.se.boot <- numeric(B1)
for(j in 1:B1){
  mhat.double.boot <- sample(mhat.boot, B, replace = TRUE)
  mhat.se.boot[j] <- sd(mhat.double.boot)
}
mhat.se <- mean(mhat.se.boot)
mhat.se

#3.4

xi = (mhat.boot - median.hat)/mhat.se

#3.5

zeta.percentiles <- quantile(xi, probs = c(0.005, 0.995))
zeta.percentiles

#3.6
bootstrap_t.CI <- c(median.hat - zeta.percentiles[2]*mhat.se, 
                    median.hat - zeta.percentiles[1]*mhat.se)

bootstrap_t.CI

#3.7
alpha = 0.01
b0 <- qnorm(1/B*sum(mhat.boot <= median.hat))
p1 <- pnorm(qnorm(alpha/2) + 2*b0)
p2 <- pnorm(qnorm(1-alpha/2) + 2*b0)

bc.BC <- quantile(mhat.boot, prob = c(p1,p2))
bc.BC



#Part4

library(boot)

foo <- function(data, indices, cor.type){
  dt <- data[indices,]
  c(cor(dt[,1], dt[,3], method = cor.type),
    median(dt[,1]),
    median(dt[,3]))
}

set.seed(1)
myBootstrap <- boot(iris, foo, R = 10000, cor.type = 's')

head(myBootstrap$t, n = 10)

myBootstrap$t0

myBootstrap

colMeans(myBootstrap$t) - myBootstrap$to

apply(myBootstrap$t, 2, sd)

plot(myBootstrap, index = 1)

#Task 9

plot(myBootstrap, index = 2)

boot.ci(myBootstrap, index = 1, conf = 0.99)

#Task 10
boot.ci(myBootstrap, index = 2)

boot.ci(myBootstrap, index = 3)

###############################################################################################

#Applied Lecture 1


#Analysis of dataset 


library('MASS')
data(Cars93)
class(Cars93)

dim(Cars93)
head(Cars93)


sum(is.na(Cars93))

Cars93 = na.omit(Cars93)
n = nrow(Cars93)
n

attach(Cars93)

plot(Horsepower, Price, col = 2, bty = 'l', pch = 16)
linear.model <- lm(Price ~ Horsepower, data = Cars93)
abline(linear.model, col = 'blue')


#Non-regression analysis

library(boot)

boot.cor.f <- function(y, x, data, indices){
  dt <- data[indices,]
  r <- cor(dt[,y], dt[,x])
  return(r)
}

boot.cor.f <- function(y, x, data, indices){
  dt <- data[indices,]
  r <- cor(dt[y], dt[x])
  return(r)
}

set.seed(1)
bootstrap.cor <- boot(y = 'Price', x = 'Horsepower', data = Cars93, 
                      statistic = boot.cor.f, R = 5000)
head(bootstrap.cor$t, n = 10)

bootstrap.cor

plot(bootstrap.cor)

boot.ci(bootstrap.cor, conf = 0.95)

boot.cor.nested.f <- function(y, x, data, indices, iter){
  dt <- data[indices,]
  r <- cor(dt[y], dt[x])
  nested = boot(
    data = data,
    y = y, 
    x = x, 
    R = iter,
    statistic = boot.cor.f)
  v = var(nested$t, na.rm = TRUE)
  return(c(r, v))
}

bootstrap.cor <- boot(y = 'Price', x = 'Horsepower', data = Cars93,
                      statistic = boot.cor.nested.f, R = 1000, iter = 100)

boot.ci(bootstrap.cor, conf = 0.95, type = 'all')

#############################################################################################

#Semi manual bootstrap linear regression with one predictor

summary(linear.model)

coef(linear.model)

confint(linear.model)

sigma(linear.model)

summary(linear.model)$r.squared

cor(Horsepower, Price)^2

#Nonparametric paired bootstrap

NPB.reg <- function(y, x, data, B){
  n = nrow(data)
  result = matrix(NA, B, 4)
  colnames(result) = c('b0', 'b1', 'sigma', 'Rsquare')
  for(j in 1:B){
    indices = sample(1:n, n, replace = TRUE)
    reg = lm(y[indices] ~ x[indices])
    betas = coef(reg)
    sigma = sigma(reg)
    R2 = summary(reg)$r.squared
    result[j,] = c(betas, sigma, R2)
  }
  result
}

set.seed(1)
NPB.reg.results <- NPB.reg(y = Price, x = Horsepower, data = Cars93, B = 1000)
head(NPB.reg.results)

par(mfrow = c(2,2))
hist(NPB.reg.results[,1], main = expression(hat(beta)[0]), xlab = '', col = 2)
hist(NPB.reg.results[,2], main = expression(hat(beta)[1]), xlab = '', col = 3)
hist(NPB.reg.results[,3], main = expression(hat(sigma)), xlab = '', col = 4)
hist(NPB.reg.results[,4], main = expression(Rˆ2), xlab = '', col = 5)

dev.off()

apply(NPB.reg.results, 2, sd)

apply(NPB.reg.results, 2, quantile, probs = c(0.025, 0.975))

plot(1, bty = 'l', xlab = 'Horsepower', ylab = 'Price', xlim = range(Horsepower),
     ylim = range(Price))
for(j in 1:1000){
  abline(a = NPB.reg.results[j,1], b = NPB.reg.results[j,2], col = 'pink', lwd = 2)
}
abline(linear.model, col = 1, lwd = 2)
points(Horsepower, Price, bty = 'l', pch = 16)

#Semiparametric residual bootstrap

SRB.reg <- function(y, x, data, B){
  n = nrow(data)
  result = matrix(NA, B, 4)
  colnames(result) = c('b0', 'b1', 'sigma', 'Rsquare')
  res = lm(y~x)$residuals
  beta.ols <- coef(lm(y~x))
  X = cbind(1,x)
  for(j in 1:B){
    indices = sample(1:n, n, replace= TRUE)
    res.boot = res[indices]
    y.boot = X %*% beta.ols + res.boot
    reg = lm(y.boot ~ x)
    betas = coef(reg)
    sigma = sigma(reg)
    R2 = summary(reg)$r.squared
    result[j, ] = c(betas, sigma, R2)
  }
  return(result)
}

set.seed(1)
SRB.reg.results = SRB.reg(y = Price, x = Horsepower, data = Cars93, B = 1000)

apply(SRB.reg.results, 2, sd)

apply(SRB.reg.results, 2, quantile, probs = c(0.025, 0.975))


#Plots
par(mfrow = c(1,2))
# NPB
plot(1, bty="l", xlab="Horsepower", ylab="Price", main = 'Nonparametric paired',
     xlim = range(Horsepower), ylim = range(Price))
for (j in 1:1000) {
  abline(a = NPB.reg.results[j,1], b = NPB.reg.results[j,2], col='pink', lwd = 2)
}
abline(linear.model, col=1, lwd = 2)
points(Horsepower, Price, bty = 'l', pch = 16)
#SRB
plot(1, bty="l", xlab="Horsepower", ylab="Price", main = 'Semiparametric residual',
     xlim = range(Horsepower), ylim = range(Price))
for (j in 1:1000) {
  abline(a = SRB.reg.results[j,1], b = SRB.reg.results[j,2], col='lightgreen', lwd = 2)
}
abline(linear.model, col=1, lwd = 2)
points(Horsepower, Price, bty = 'l', pch = 16)





#Multiple Linear Regression
install.packages('car')
library(car)

linear.model1 <- lm(Price~Horsepower + Fuel.tank.capacity + Passengers, data = Cars93)
summary(linear.model)

set.seed(1)
NPB.results1 <- Boot(linear.model1, R = 1000)
summary(NPB.results1, high.moments = TRUE)


set.seed(1)
NPB.results2 <- Boot(linear.model1, method = 'case',
                     f = function(model){c(coef(model), sigma(model), summary(model)$r.squared)},
                     labels = c(names(coef(linear.model1)), 'sigmaHat', 'Rsquared'), R = 1000)

summary(NPB.results2, high.moments = TRUE)

confint(NPB.results2, level = 0.95, type = 'perc')

hist(NPB.results2, legend = 'top')

head(NPB.results2$t)



set.seed(1)
SRB.results2 = Boot(linear.model1, method = 'residual',
                    f = function(model){c(coef(model), sigma(model), summary(model)$r.squared)},
                    labels = c(names(coef(linear.model1)), 'sigmaHat', 'Rsquared'),
                    R = 1000)



###################################################################################################

#Practical 2

library(MASS)
nrow(Cars93)

install.packages('bayesboot')
library(bayesboot)


samples = c(1.865, 3.053, 1.401, 0.569, 4.132)

library(boot)
foo = function(data, indices){
  dt <- data[indices]
  return(c(mean(dt)))
}

cb <- boot(data = samples, statistic = foo, R = 10000)

bb <- bayesboot(data = samples, statistic = mean, R = 10000)

head(cb$t)
head(bb$V1)

#Fancy graphs

par(mfrow=c(1,2))
hist(cb$t, main = 'Classical bootstrap', col = 'lightgreen', xlab = '')
hist(bb$V1, main = 'Bayesian bootstrap', col = 'lightblue', xlab = '')

par(mfrow = c(1,1))
hist(bb$V1, col = 'lightblue', xlab = '', main = 'Bootstrap estimates of the mean')
hist(cb$t, col = 'lightgreen', xlab = '', add = TRUE, main = '')
legend('topright', legend = c('Classical bootstrap','Bayesian bootstrap'),
       col = c('lightgreen','lightblue'), pch = 15, bty = 'n')

hist(bb$V1, col = rgb(0,0,1,1/4), xlab = '', main = 'Bootstrap estimates of the mean')
hist(cb$t, col = rgb(1,0,0,1/4), xlab = '', add = TRUE, main = '')
legend('topright', legend = c('Classical bootstrap','Bayesian bootstrap'),
       col = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), pch = 15, bty = 'n')


######USING BAYESBOOT

#means and variances
set.seed(1)
bb.mean1 <- bayesboot(data = Price, statistic = mean, R = 5000, use.weights = FALSE)

set.seed(1)
bb.mean2 <- bayesboot(data = Price, statistic = weighted.mean, R = 5000, use.weights = TRUE)

summary(bb.mean1)

hist(bb.mean1$V1, main= 'Bootstrap means of Price', xlab = '', col = 'pink')

plot(bb.mean1)

plot(bb.mean1, main = 'Bootstrap means of Price', cex.main = 1, col = 'pink')

#Task 1

hist(bb.mean2$V1, main = 'Boostrap means of Price', xlab = '', col = 'yellow')

plot(bb.mean2, main = 'Bootstrap means of Price', cex.main = 1, col = 'yellow')

#Task 2

set.seed(1)
bb.var <- bayesboot(data = Horsepower, statistic = var, R = 5000, use.weights = FALSE)

hist(bb.var$V1, main = 'Boostrap variances of Horsepower', xlab = '', col = 'green')

plot(bb.var, main = 'Bootstrap variances of Horsepower', cex.main = 1, col = 'green')


#correlations

PrHo <- data.frame(Price, Horsepower)
cor(PrHo)

cor.fun <- function(data){
  cor(data)[1,2]
}

set.seed(1)
bb.cor1 <- bayesboot(data = PrHo, statistic = cor.fun, R = 5000, use.weights = FALSE)
plot(bb.cor1, main = 'Bootstrap correlations of Price and Horsepower',
     cex.main = 1, xlab = '', col = 'magenta')

corr(PrHo)

#Task 3

set.seed(1)
bb.cor2 <- bayesboot(data = PrHo, statistic = corr, R = 5000, use.weights = TRUE)
plot(bb.cor2, main = 'Bootstrap correlations of Price and Horsepower', 
     cex.main = 1, xlab = '', col = 'magenta', showMode = TRUE, showCurve = TRUE, cred.mass = 0.99)

#4.3 Multiple statistics

foo1 <- function(data){
  return(c(mean(data), sd(data)))
}

set.seed(1)
bb.mult1 <- bayesboot(data = Price, statistic = foo1, R = 5000, use.weights = FALSE)

names(bb.mult1) <- c('Mean', 'Standard deviation')

plot(bb.mult1, cex = 1, cex.lab = 1)
plot(bb.mult1, cex = 1, cex.lab = 1, showCurve = TRUE)

summary(bb.mult1)

head(bb.mult1)

bb.mean <- bb.mean1$mean
bb.sd <- bb.mult1$`Standard deviation`


#Task 4
install.packages('moments')
library(moments)

foo2 <- function(data){
  return(c(skewness(data), kurtosis(data)))
}

set.seed(1)
bb.mult2 <- bayesboot(data = Horsepower, statistic = foo2, R = 5000, use.weights = FALSE)

names(bb.mult2) <- c('Skew', 'Kurtosis')
plot(bb.mult2, col = 'magenta', cex.lab = 1, cex = 1, showCurve = TRUE, cred.mass = 0.95)



#############Using Bayesboot for regression

#One predictor

lm.coefs1 <- function(data){
  coef(lm(Price ~ Horsepower, data = data))
}
lm.coefs1(Cars93)

set.seed(1)
bb.reg1 <- bayesboot(data = Cars93, statistic = lm.coefs1, R = 1000, use.weights = FALSE)
plot(bb.reg1, cex = 1, cex.lab = 1)
summary(bb.reg1)

plot(1, bty="l", xlab="Horsepower", ylab="Price", xlim = range(Horsepower), ylim = range(Price))
for (j in 1:1000) {
  abline(coef = bb.reg1[j,], col='pink', lwd = 2)
}
abline(coef = coef(lm(Price ~ Horsepower, data = Cars93)), col=1, lwd = 2)
points(Horsepower, Price, bty = 'l', pch = 16)

b0 <- bb.reg1$`(Intercept)`
b1 <- bb.reg1$Horsepower



#Task 5


lm.coefs2 <- function(data, w){
  xlm <- coef(lm(Price ~ Horsepower, , data = data , weights = w))
  return(xlm)
}

set.seed(1)
bb.reg2 <- bayesboot(Cars93, statistic = lm.coefs2, R = 1000, use.weights = TRUE)

plot(bb.reg2, cex.lab = 1, cex = 1, col = 'magenta', cred.mass = 0.95, showCurve = TRUE)


#With multiple predictors

lm.coefs3 <- function(data){
  coef(lm(Price ~ Horsepower + Fuel.tank.capacity + Passengers, data = data))
}

lm.coefs3(Cars93)

bb.reg3 <- bayesboot(data = Cars93, statistic = lm.coefs3, R = 1000, use.weights = FALSE)

plot(bb.reg3, cex.lab = 1, cex = 1, col = 'lightgreen')



###################################################################################################

#Practical 3

data(iris)

dim(iris)

y = iris[,1:4]
n = nrow(y)

par(mfrow = c(2,2))
for(i in 1:4){
  hist(y[,i], col = i+1, main = '', xlab = names(y)[i])
}
pairs(y)

############BNPmix

install.packages('BNPmix')
library(BNPmix)

#Univariate application

#Task 1

PYcalibrate(Ek = 3, n = n)
mcmc <- list(niter = 11000, nburn = 1000)
prior <- list(strength = 1, discount = 0, hyper = 'FALSE')
output <- list(grid = seq(min(y[,3]), max(y[,3]), length.out = 100), out_param = TRUE)
set.seed(1)
fit1 <- PYdensity(y = y[,3], mcmc = mcmc, prior = prior, output = output)
summary(fit1)

plot(fit1, band = TRUE, show_clust = TRUE, xlab = "Petal length", ylab = 'Density')

#Finding MAP number of clusters

names(fit1)

cluster.labels <- fit1$clust
class(cluster.labels)
dim(cluster.labels)


n.cluster = apply(cluster.labels, 1, unique)
n.cluster = lapply(n.cluster, length)
n.cluster = unlist(n.cluster)
counts = table(n.cluster)
counts

plot(counts/sum(counts), col = 1:length(counts), bty = 'l', xlab = 'Number of clusters',
     ylab = 'Posterior probability', main = 'alpha = 1')


#Task 2

DPrior <- PYcalibrate(Ek = 3, n = n, discount = 0)
DPrior$strength

PYcalibrate(Ek = 3, n = n)
mcmc <- list(niter = 11000, nburn = 1000)
prior <- list(strength = DPrior$strength, discount = 0, hyper = 'FALSE')
output <- list(grid = seq(min(y[,3]), max(y[,3]), length.out = 100), out_param = TRUE)
set.seed(1)
fit2 <- PYdensity(y = y[,3], mcmc = mcmc, prior = prior, output = output)
summary(fit1)
plot(fit2, band = TRUE, show_clust = TRUE, xlab = "Petal length", ylab = 'Density')

cluster.labels <- fit2$clust

n.cluster = apply(cluster.labels, 1, unique)
n.cluster = lapply(n.cluster, length)
n.cluster = unlist(n.cluster)
counts = table(n.cluster)
counts

plot(counts/sum(counts), col = 1:length(counts), bty = 'l', xlab = 'Number of clusters',
     ylab = 'Posterior probability', main = 'alpha = 1')



#Posterior inference for parameters

means <- fit1$mean


post.means <- matrix(NA, nrow = 10000, ncol = 2)
for(i in 1:10000){
  post.means[i,1] <- means[[i]][1]
  post.means[i,2] <- means[[i]][2]
}

summary(post.means)


#Task 7

sigmas <- fit1$sigma2

post.sigmas <- matrix(NA, nrow = 10000, ncol = 2)
for(i in 1:10000){
  post.sigmas[i,1] <- sigmas[[i]][1]
  post.sigmas[i,2] <- sigmas[[i]][2]
}

plot(density(post.sigmas[,1]), main = '', lwd = 2, ylab = 'posterior density',
     xlab = expression(sigma[1]), bty  = 'l', col = 1, cex.lab = 1)

plot(density(post.sigmas[,2]), main = '', lwd = 2, ylab = 'posterior density',
     xlab = expression(sigma[2]), bty  = 'l', col = 1, cex.lab = 1, xlim = c(0,0.2))


names(fit1)

probs <- fit1$probs

post.probs <- matrix(NA, nrow = 10000, ncol = 2)
for(i in 1:10000){
  post.probs[i,1] <- probs[[i]][1]
  post.probs[i,2] <- probs[[i]][2]
}

plot(density(post.probs[,1]), main = '', lwd = 2, ylab = 'posterior density',
     xlab = expression(probs[1]), bty  = 'l', col = 1, cex.lab = 1)

plot(density(post.probs[,2]), main = '', lwd = 2, ylab = 'posterior density',
     xlab = expression(probs[2]), bty  = 'l', col = 1, cex.lab = 1)




#Controlling the specification of prior hyper parameters

set.seed(1)
mean(1/rgamma(1000, shape = 1, rate = 1000))
var(1/rgamma(1000, shape = 1, rate = 1000))


#Task 3


PYcalibrate(Ek = 3, n = n)
mcmc <- list(niter = 11000, nburn = 1000)
prior1 <- list(strength = DPrior$strength, discount = 0, hpyer = 'FALSE', m0 = 0, k0 = 1, a0 = 1, b0 = 100)
output <- list(grid = seq(min(y[,3]), max(y[,3]), length.out = 100), out_param = TRUE)

set.seed(1)
fit3 <- PYdensity(y = y[,3], mcmc = mcmc, prior = prior1, output = output)
summary(fit3)

plot(fit3, band = TRUE, show_hist = TRUE)
plot(fit3, band = TRUE, show_clust = TRUE)

cluster.labels3 <- fit3$clust

n.cluster3 = apply(cluster.labels3, 1, unique)
n.cluster3 = lapply(n.cluster3, length)
n.cluster3 = unlist(n.cluster3)
counts3 = table(n.cluster3)
counts3

plot(counts3)

means <- fit3$mean

post.means <- matrix(NA, nrow = 10000, ncol = 2)
for(i in 1:10000){
  post.means[i,1] <- means[[i]][1]
  post.means[i,2] <- means[[i]][2]
}

plot(density(post.means[,1]))
plot(density(post.means[,2]))


#Empirical observations

DPrior <- PYcalibrate(Ek = 3, n = n, discount = 0)
mcmc <- list(niter = 11000, nburn = 1000)
prior <- list(strength = DPrior$strength, discount = 0, hyper = 'FALSE')
output = list(grid = seq(min(y[,4]), max(y[,4]), length.out = 100),
              out_param = TRUE)

set.seed(1)
fit3 = PYdensity(y = y[,4], mcmc = mcmc, prior = prior, output = output)

y.scaled <- scale(y)
apply(y.scaled, 2, mean); apply(y.scaled, 2, var)

output = list(grid = seq(min(y.scaled[,4]), max(y.scaled[,4]), length.out = 100),
              out_param = TRUE)
set.seed(1)
fit3 = PYdensity(y = y.scaled[,4], mcmc = mcmc, prior = prior, output = output)
summary(fit3)

plot(fit3, band = TRUE, show_hist = TRUE, xlab = "Petal width",
     ylab = 'Density', col = 3)




mcmc = list(niter = 11000, nburn = 1000, method = 'MAR')
output = list(grid = seq(min(y[,4]), max(y[,4]), length.out = 100),
              out_param = TRUE)
set.seed(1)
fit3 = PYdensity(y = y[,4], mcmc = mcmc, prior = prior, output = output)

summary(fit3)

plot(fit3, band = TRUE, show_hist = TRUE, xlab = "Petal width",
     ylab = 'Density', col = 3)

#################################################################################################

#Applied Lecture 2

data('iris')
class(iris)

dim(iris)

head(iris, n=5)

y = iris[, 1:4]
n = nrow(y)
library(BNPmix)

#HPDM model for variable petal length in iris dataset

DPprior <- PYcalibrate(Ek = 2, n = n, discount = 0)
DPprior

mcmc <- list(niter = 11000, nburn = 1000)
prior <- list(strength = DPprior$strength, discount = 0, hyper = TRUE)
output <- list(grid = seq(min(y[,3]), max(y[,3]), length.out = 100),
               out_param = TRUE)

set.seed(1)
fit4 = PYdensity(y = y[,3], mcmc = mcmc, prior = prior, output = output)

summary(fit4)

plot(fit4, band = TRUE, show_hist = TRUE, show_clust = TRUE,
     xlab = "Petal length", ylab = 'Density')



prior <- list(strength = DPprior$strength, discount = 0, hyper = 'TRUE', 
              m1 = 0, s21 = 1000, k0 = 1, a0 = 3, a1 = 0.001, b1 = 0.001)
mcmc <- list(niter = 11000, nburn = 1000)

set.seed(1)
fit5 <- PYdensity(y = y[,3], mcmc = mcmc, prior = prior, output = output)

summary(fit5)

plot(fit5, band = TRUE, show_hist = TRUE, show_clust = TRUE,
     xlab = "Petal length", ylab = 'Density', col = 7)

cluster.labels = fit5$clust
n.cluster = apply(cluster.labels, 1, unique)
n.cluster = lapply(n.cluster, length)
n.cluster = unlist(n.cluster)
counts = table(n.cluster)
counts

plot(counts/sum(counts), col = 1:length(counts), bty = 'l',
     xlab = 'Number of clusters', ylab = 'Distribution of clusters')


#Multivariate DPM 

grid.length = 20
grid1 = seq(min(y[,1]), max(y[,1]), length.out = grid.length)
grid2 = seq(min(y[,2]), max(y[,2]), length.out = grid.length)
grid3 = seq(min(y[,3]), max(y[,3]), length.out = grid.length)
grid4 = seq(min(y[,4]), max(y[,4]), length.out = grid.length)

grid <- expand.grid(grid1, grid2, grid3, grid4)
dim(grid)

mcmc.total.n <- 2000
mcmc.burnin <- 1000


#To increase numbers of groups, play around with Ek and strength parameters
DPprior = PYcalibrate(Ek = 4, n = n, discount = 0)
mcmc = list(niter = mcmc.total.n, nburn = mcmc.burnin, method = 'ICS')
prior = list(strength = DPprior$strength, discount = 0, hyper = FALSE)
output = list(grid = grid, out_param = TRUE)
set.seed(1)
mult.fit1 = PYdensity(y = y, mcmc = mcmc, prior = prior, output = output)
summary(mult.fit1)

cluster.labels = mult.fit1$clust
n.cluster = apply(cluster.labels, 1, unique)
n.cluster = lapply(n.cluster, length)
n.cluster = unlist(n.cluster)
counts = table(n.cluster)
counts

prior = list(strength = DPprior$strength, discount = 0, hyper = FALSE,
             m0 = rep(0,4), Sigma0 = diag(1000,4))
set.seed(1)
mult.fit2 = PYdensity(y = y, mcmc = mcmc, prior = prior, output = output)
summary(mult.fit2)



#Multivariate realisations

p1 = plot(mult.fit2, dim = c(1, 1), xlab = "Sepal length", ylab = "Density")
p12 = plot(mult.fit2, dim = c(1, 2), show_clust = TRUE,
           xlab = "Sepal length", ylab = "Sepal width")
p1
p12

p1 = plot(mult.fit2, dim = c(1, 1), xlab = "Sepal length", ylab = "Density",
          col = 2)
p2 = plot(mult.fit2, dim = c(2, 2), xlab = "Sepal width", ylab = "Density",
          col = 3)
p3 = plot(mult.fit2, dim = c(3, 3), xlab = "Petal length", ylab = "Density",
          col = 4)
p4 = plot(mult.fit2, dim = c(4, 4), xlab = "Petal width", ylab = "Density",
          col = 5)
p12 = plot(mult.fit2, dim = c(1, 2), show_clust = TRUE, xlab = "Sepal length",
           ylab = "Sepal width", col = 2)
p13 = plot(mult.fit2, dim = c(1, 3), show_clust = TRUE, xlab = "Sepal length",
           ylab = "Petal length", col = 2)
p14 = plot(mult.fit2, dim = c(1, 4), show_clust = TRUE, xlab = "Sepal length",
           ylab = "Petal width", col = 2)
p21 = plot(mult.fit2, dim = c(2, 1), show_clust = TRUE, xlab = "Sepal width",
           ylab = "Sepal length", col = 3)
p23 = plot(mult.fit2, dim = c(2, 3), show_clust = TRUE, xlab = "Sepal width",
           ylab = "Petal length", col = 3)
p24 = plot(mult.fit2, dim = c(2, 4), show_clust = TRUE, xlab = "Sepal width",
           ylab = "Petal width", col = 3)
p31 = plot(mult.fit2, dim = c(3, 1), show_clust = TRUE, xlab = "Petal length",
           ylab = "Sepal length", col = 4)
p32 = plot(mult.fit2, dim = c(3, 2), show_clust = TRUE, xlab = "Petal length",
           ylab = "Sepal width", col = 4)
p34 = plot(mult.fit2, dim = c(3, 4), show_clust = TRUE, xlab = "Petal length",
           ylab = "Petal width", col = 4)
p41 = plot(mult.fit2, dim = c(4, 1), show_clust = TRUE, xlab = "Petal width",
           ylab = "Sepal length", col = 5)
p42 = plot(mult.fit2, dim = c(4, 2), show_clust = TRUE, xlab = "Petal width",
           ylab = "Sepal width", col = 5)
p43 = plot(mult.fit2, dim = c(4, 3), show_clust = TRUE, xlab = "Petal width",
           ylab = "Petal length", col = 5)


library(gridExtra)
grid.arrange(p1, p12, p13, p14, p21, p2, p23, p24, p31, p32, p3, p34, p41, p42,
             p43, p4, layout_matrix = matrix(1:16, 4, 4))

grid.arrange(p1, p2, p3, p4, layout_matrix = matrix(1:4, 1, 4))
grid.arrange(p12, p13, p14, p21, p23, p24, p31, p32, p34, p41, p42, p43,
             layout_matrix = matrix(1:12, 3, 4))

cluster.labels = mult.fit2$clust
n.cluster = apply(cluster.labels, 1, unique)
n.cluster = lapply(n.cluster, length)
n.cluster = unlist(n.cluster)
counts = table(n.cluster)
counts/sum(counts)

names(mult.fit2)

means <- mult.fit2$mean
length(means)

post.means1 = matrix(NA, nrow = mcmc.total.n-mcmc.burnin, ncol = 2)
post.means2 = matrix(NA, nrow = mcmc.total.n-mcmc.burnin, ncol = 2)
post.means3 = matrix(NA, nrow = mcmc.total.n-mcmc.burnin, ncol = 2)
post.means4 = matrix(NA, nrow = mcmc.total.n-mcmc.burnin, ncol = 2)
for(i in 1:(mcmc.total.n-mcmc.burnin)){
  post.means1[i,] = means[[i]][,1][1:2]
  post.means2[i,] = means[[i]][,2][1:2]
  post.means3[i,] = means[[i]][,2][1:2]
  post.means4[i,] = means[[i]][,4][1:2]
}
summary(post.means1)


apply(post.means1,2,var)

par(mfrow = c(1,2))
plot(density(post.means1[,1]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(mu[1]), bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.means1[,2]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(mu[2]), bty = 'l', col = 2, cex.lab = 1.3)
mtext('Group means for Sepal length', side = 3, line = -2, outer = TRUE,
      cex = 1.5)

library(ggplot2)
all.means1 = data.frame(cluster = factor(rep(c('Cluster1','Cluster2'), 
                                             each = mcmc.total.n-mcmc.burnin)), 
                        samples = c(post.means1[,1],post.means1[,2]))
head(all.means1)

plot.m1 = ggplot(all.means1, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), 
        axis.text=element_text(size=10)) +
  xlab(expression(mu)) + ylab('Posterior density') + 
  ggtitle('Sepal length') 
plot.m1

all.means2 = data.frame(cluster = factor(rep(c('Cluster1','Cluster2'), each = mcmc.total.n-mcmc.burnin)), samples = c(post.means2[,1],post.means2[,2]))
all.means3 = data.frame(cluster = factor(rep(c('Cluster1','Cluster2'), each = mcmc.total.n-mcmc.burnin)), samples = c(post.means3[,1],post.means3[,2]))
all.means4 = data.frame(cluster = factor(rep(c('Cluster1','Cluster2'), each = mcmc.total.n-mcmc.burnin)), samples = c(post.means4[,1],post.means4[,2]))

plot.m2 = ggplot(all.means2, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab(expression(mu)) + ylab('Posterior density') + 
  ggtitle('Sepal width') 

plot.m3 = ggplot(all.means3, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab(expression(mu)) + ylab('Posterior density') + 
  ggtitle('Petal length') 

plot.m4 = ggplot(all.means4, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab(expression(mu)) + ylab('Posterior density') + 
  ggtitle('Petal width')

grid.arrange(plot.m1, plot.m2, plot.m3, plot.m4, 
             layout_matrix = matrix(1:4, 2, 2))


covariances = mult.fit2$sigma2
covariances[[1]] 

Covs1 = list()
Covs2 = list()
for(i in 1 : c(mcmc.total.n - mcmc.burnin)){
  Covs1[[i]] = covariances[[i]][,,1]
  Covs2[[i]] = covariances[[i]][,,2]
}

post.mean.cov1 = Reduce('+', Covs1)/(mcmc.total.n - mcmc.burnin)
post.mean.cov2 = Reduce('+', Covs2)/(mcmc.total.n - mcmc.burnin)
colnames(post.mean.cov1) = c('Sepal length', 'Sepal width', 'Petal length',
                             'Petal width')
colnames(post.mean.cov2) = c('Sepal length', 'Sepal width', 'Petal length',
                             'Petal width')
rownames(post.mean.cov1) = colnames(post.mean.cov1)
rownames(post.mean.cov2) = colnames(post.mean.cov2)
post.mean.cov1
post.mean.cov2

post.mean.cor1 = cov2cor(post.mean.cov1)
post.mean.cor2 = cov2cor(post.mean.cov2)
post.mean.cor1
post.mean.cor2

install.packages('pheatmap') 
library(pheatmap)
pheatmap(post.mean.cor1, display_numbers = TRUE, cluster_rows = FALSE,
         cluster_cols = FALSE, fontsize_number = 15,
         main = 'Average posterior correlations of means in Cluster 1')

pheatmap(post.mean.cor2, display_numbers = TRUE, cluster_rows = FALSE,
         cluster_cols = FALSE, fontsize_number = 15,
         main = 'Average posterior correlations of means in Cluster 2')

#install.packages('corrplot')
library(corrplot)
corrplot.mixed(post.mean.cor1, lower = "number", upper = "ellipse")
mtext('Average posterior correlations of means in Cluster 1', side = 2,
      line = -4, outer = TRUE)

corrplot.mixed(post.mean.cor2, lower = "number", upper = "ellipse")
mtext('Average posterior correlations of means in Cluster 2', side = 2,
      line = -4, outer = TRUE)


pi = mult.fit2$probs
post.pi = matrix(NA, nrow = mcmc.total.n - mcmc.burnin, ncol = 2)
for(i in 1:(mcmc.total.n - mcmc.burnin)){
  post.pi[i,1] = pi[[i]][1]
  post.pi[i,2] = pi[[i]][2]
}
summary(post.pi)

par(mfrow = c(1,2))
plot(density(post.pi[,1]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(pi[1]), bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.pi[,2]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(pi[2]), bty = 'l', col = 2, cex.lab = 1.3)
mtext('2-component mixture weights', side = 3, line = -2, outer = TRUE,
      cex = 1.5)



#Last considerations and conclusions

y.setosa = y[iris[,5] == 'setosa', ]
n.setosa = nrow(y.setosa)
n.setosa

grid1.setosa = seq(min(y.setosa[,1]), max(y.setosa[,1]),
                   length.out = grid.length)
grid2.setosa = seq(min(y.setosa[,2]), max(y.setosa[,2]),
                   length.out = grid.length)
grid3.setosa = seq(min(y.setosa[,3]), max(y.setosa[,3]),
                   length.out = grid.length)
grid4.setosa = seq(min(y.setosa[,4]), max(y.setosa[,4]),
                   length.out = grid.length)
grid.setosa = expand.grid(grid1.setosa, grid2.setosa, grid3.setosa,
                          grid4.setosa)

DPprior = PYcalibrate(Ek = 4, n = n.setosa, discount = 0)
DPprior

mcmc = list(niter = mcmc.total.n, nburn = mcmc.burnin, method = 'MAR')
prior = list(strength =  DPprior$strength, discount = 0, hyper = FALSE,
             m0 = rep(0,4), Sigma0 = diag(1000,4))
output = list(grid = grid.setosa, out_param = TRUE)
set.seed(1)
mult.fit1.setosa = PYdensity(y = y.setosa, mcmc = mcmc, prior = prior,
                             output = output)
summary(mult.fit1.setosa)

cluster.labels = mult.fit1.setosa$clust
n.cluster = apply(cluster.labels, 1, unique)
n.cluster = lapply(n.cluster, length)
n.cluster = unlist(n.cluster)
counts = table(n.cluster)
counts

p1 = plot(mult.fit1.setosa, dim = c(1, 1), xlab = "Sepal length",
          ylab = "Density", col = 2)
p2 = plot(mult.fit1.setosa, dim = c(2, 2), xlab = "Sepal width",
          ylab = "Density", col = 3)
p3 = plot(mult.fit1.setosa, dim = c(3, 3), xlab = "Petal length",
          ylab = "Density", col = 4)
p4 = plot(mult.fit1.setosa, dim = c(4, 4), xlab = "Petal width",
          ylab = "Density", col = 5)

p12 = plot(mult.fit1.setosa, dim = c(1, 2), show_clust = TRUE,
           xlab = "Sepal length", ylab = "Sepal width", col = 2)
p13 = plot(mult.fit1.setosa, dim = c(1, 3), show_clust = TRUE,
           xlab = "Sepal length", ylab = "Petal length", col = 2)
p14 = plot(mult.fit1.setosa, dim = c(1, 4), show_clust = TRUE,
           xlab = "Sepal length", ylab = "Petal width", col = 2)

p21 = plot(mult.fit1.setosa, dim = c(2, 1), show_clust = TRUE,
           xlab = "Sepal width", ylab = "Sepal length", col = 3)
p23 = plot(mult.fit1.setosa, dim = c(2, 3), show_clust = TRUE,
           xlab = "Sepal width", ylab = "Petal length", col = 3)
p24 = plot(mult.fit1.setosa, dim = c(2, 4), show_clust = TRUE,
           xlab = "Sepal width", ylab = "Petal width", col = 3)

p31 = plot(mult.fit1.setosa, dim = c(3, 1), show_clust = TRUE,
           xlab = "Petal length", ylab = "Sepal length", col = 4)
p32 = plot(mult.fit1.setosa, dim = c(3, 2), show_clust = TRUE,
           xlab = "Petal length", ylab = "Sepal width", col = 4)
p34 = plot(mult.fit1.setosa, dim = c(3, 4), show_clust = TRUE,
           xlab = "Petal length", ylab = "Petal width", col = 4)

p41 = plot(mult.fit1.setosa, dim = c(4, 1), show_clust = TRUE,
           xlab = "Petal width", ylab = "Sepal length", col = 5)
p42 = plot(mult.fit1.setosa, dim = c(4, 2), show_clust = TRUE,
           xlab = "Petal width", ylab = "Sepal width", col = 5)
p43 = plot(mult.fit1.setosa, dim = c(4, 3), show_clust = TRUE,
           xlab = "Petal width", ylab = "Petal length", col = 5)

grid.arrange(p1, p12, p13, p14, p21, p2, p23, p24, p31, p32, p3, p34, p41, p42,
             p43, p4, layout_matrix = matrix(1:16, 4, 4))


y.versicolor = y[iris[,5] == 'versicolor', ]
y.virginica = y[iris[,5] == 'virginica', ]

par(mfrow = c(1,3))
plot(as.table(apply(y.setosa,2, mean)), ylim = c(0,7), col = c(2:5),
     lwd = 5, ylab = 'Centimeters', main = 'Setosa')
plot(as.table(apply(y.versicolor,2, mean)), ylim = c(0,7), col = c(2:5),
     lwd = 5, ylab = 'Centimeters', main = 'Versicolor')
plot(as.table(apply(y.virginica,2, mean)), ylim = c(0,7), col = c(2:5),
     lwd = 5, ylab = 'Centimeters', main = 'Virginica')


###################################################################################
#Practical 4

#install.packages('AppliedPredictiveModeling') 
library(AppliedPredictiveModeling)
library(BNPmix)
library(ggplot2)
library(gridExtra)

data(package = 'AppliedPredictiveModeling')

data(abalone)
class(abalone)
head(abalone)
dim(abalone)
sum(is.na(abalone))

attach(abalone)
x1 = LongestShell
x2 = Diameter
x3 = Height
x4 = WholeWeight
x5 = ShuckedWeight
x6 = VisceraWeight
x7 = ShellWeight

Age = Rings + 1.5
y = Age
hist(y, xlab = 'Abalone age', col = 51, main = '')


par(mfrow = c(3,3))
plot(LongestShell, Age, col = Type, bty = 'l')
plot(Diameter,Age, col = Type, bty = 'l')
plot(Height,Age, col = Type, bty = 'l')
plot(WholeWeight,Age, col = Type, bty = 'l')
plot(ShuckedWeight,Age, col = Type, bty = 'l')
plot(VisceraWeight,Age, col = Type, bty = 'l')
plot(ShellWeight,Age, col = Type, bty = 'l')

levels(Type)

pairs(abalone[ , -c(1, 9)], col = Type)


#Applications of DPM regression

DPprior <- PYcalibrate(Ek = 3, n = length(y), discount = 0)

prior = list(strength = DPprior$strength, discount = 0, hyper = FALSE)
mcmc <- list(niter = 2000, nburn = 1000)
grid_y <- seq(min(y), max(y), length.out = 100)
grid_x3 <- quantile(x3, probs = c(0.01, 0.25, 0.50, 0.75, 0.99))
gird_x3 <- round(grid_x3, 3)

output <- list(grid_y = grid_y, grid_x = grid_x3, out_type = 'full', out_param = TRUE)
set.seed(1)
fit.reg <- PYregression(y = y, x = x3, prior = prior, mcmc = mcmc, output = output)

summary(fit.reg)

#Task 1

y <- log(Age)

par(mfrow = c(3,3))
plot(LongestShell, y, col = Type, bty = 'l')
plot(Diameter,y, col = Type, bty = 'l')
plot(Height,y, col = Type, bty = 'l')
plot(WholeWeight,y, col = Type, bty = 'l')
plot(ShuckedWeight,y, col = Type, bty = 'l')
plot(VisceraWeight,y, col = Type, bty = 'l')
plot(ShellWeight,y, col = Type, bty = 'l')

#Task 2
#Same as before lol

grid_y <- seq(min(y), max(y), length.out = 100)
output <- list(grid_y = grid_y, grid_x = grid_x3, out_type = 'full', out_param = TRUE)
fit.reg <- PYregression(y = y, x = x3, prior = prior, mcmc = mcmc, output = output)

#Task 3

DPprior <- PYcalibrate(Ek = 2, n = length(y) ,discount = 0)
DPprior
prior = list(strength = DPprior$strength, discount = 0, hyper = FALSE)
fit.reg <- PYregression(y = y, x = x3, prior = prior, mcmc = mcmc, output = output)
summary(fit.reg)

#Task 4

prior = list(strength = DPprior$strength, discount = 0, hyper = 'FALSE', 
             m0 = c(0,0), a0 = 3, b0 = 100)
set.seed(1)
fit.reg <- PYregression(y = y, x = x3, prior = prior, mcmc = mcmc, output = output)
summary(fit.reg)

cluster.labels = fit.reg$clust
n.cluster = apply(cluster.labels, 1, unique)
n.cluster = lapply(n.cluster, length)
n.cluster = unlist(n.cluster)
counts = table(n.cluster)
counts

##

regplot = data.frame(
  dens = as.vector(apply(fit.reg$density, c(1, 2), mean)),
  qlow = as.vector(apply(fit.reg$density, c(1, 2),
                         quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.reg$density, c(1, 2),
                         quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x3)),
  label = factor(rep(paste('Height = ', grid_x3), each = length(grid_y)),
                 level = rep(paste('Height = ', grid_x3)))
)

ggplot(regplot) + theme_bw() +
  geom_line(data = regplot, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot, map = aes(x = grid, ymin = qlow,
                                        ymax = qupp), fill = 'blue', alpha = 0.3) +
  facet_wrap(~label, ncol = 5, nrow = 1) +
  labs(x = 'log(Age)', y = 'Density')


#Posterior summaries and plots

names(fit.reg)

betas = fit.reg$beta
betas[[1]]

post.betas1 = matrix(NA, nrow = 1000, ncol = 2)
post.betas2 = matrix(NA, nrow = 1000, ncol = 2)
post.betas3 = matrix(NA, nrow = 1000, ncol = 2)
for(i in 1:1000){
  post.betas1[i,] = betas[[i]][1,]
  post.betas2[i,] = betas[[i]][2,]
  post.betas3[i,] = betas[[i]][3,]
}
summary(post.betas1)

par(mfrow = c(3,2))
plot(density(post.betas1[,1]), main = '', lwd = 2, ylab = '',
     xlab = expression(beta[0]), bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.betas1[,2]), main = '', lwd = 2, ylab = '',
     xlab = expression(beta[1]), bty = 'l', col = 2, cex.lab = 1.3)
mtext('Cluster 1', side = 3, line = -2, outer = TRUE, cex = 1.2)
plot(density(post.betas2[,1]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(beta[0]), bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.betas2[,2]), main = '', lwd = 2, ylab = '',
     xlab = expression(beta[1]), bty = 'l', col = 2, cex.lab = 1.3)
mtext('Cluster 2', side = 3, line = -14, outer = TRUE, cex = 1.2)
plot(density(post.betas3[,1]), main = '', lwd = 2, ylab = '',
     xlab = expression(beta[0]), bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.betas3[,2]), main = '', lwd = 2, ylab = '',
     xlab = expression(beta[1]), bty = 'l', col = 2, cex.lab = 1.3)
mtext('Cluster 3', side = 3, line = -26, outer = TRUE, cex = 1.2)


#Task 5

sigmas <- fit.reg$sigma2

post.sigmas <- matrix(NA, nrow = 1000, ncol = 3)
for(i in 1:1000){
  post.sigmas[i,] <- sigmas[[i]][1:3]
}

summary(post.sigmas)
par(mfrow = c(1,3))
plot(density(post.sigmas[,1]), main = '', cex = 1, cex.lab = 1, col = 1,
     xlab = expression(sigma[1]^2))
plot(density(post.sigmas[,2]), main = '', cex = 1, cex.lab = 1, col = 1,
     xlab = expression(sigma[2]^2))
plot(density(post.sigmas[,3]), main = '', cex = 1, cex.lab = 1, col = 1,
     xlab = expression(sigma[3]^2))

#Task 6

pi <- fit.reg$probs

post.probs <- matrix(NA, nrow = 1000, ncol = 3)
for(i in 1:1000){
  post.probs[i,] <- pi[[i]][1:3]
}
summary(post.probs)

plot(density(post.probs[,1]), main = '', cex = 1, cex.lab = 1, col = 1,
     xlab = expression(pi[1]))
plot(density(post.probs[,2]), main = '', cex = 1, cex.lab = 1, col = 1,
     xlab = expression(pi[2]))
plot(density(post.probs[,3]), main = '', cex = 1, cex.lab = 1, col = 1,
     xlab = expression(pi[3]))

table(Type)/n


#####Regression with two predictors

X = cbind(x3, x2)

grid_x2 <- quantile(x2, probs = c(0.5))
grid_x2

grid_X <- expand.grid(grid_x3, grid_x2)
grid_X

DPprior

prior = list(strength = DPprior$strength, discount = 0, hyper = FALSE,
             m0 = c(0,0,0), a0 = 3, b0 = 100)
output <- list(grid_x = grid_X, grid_y = grid_y, out_type = 'FULL', out_param = TRUE)
set.seed(1)
fit.reg1 = PYregression(y = y, x = X, prior = prior, mcmc = mcmc, output = output)
summary(fit.reg1)


regplot <- data.frame(
  dens = as.vector(apply(fit.reg1$density, c(1, 2), mean)),
  qlow = as.vector(apply(fit.reg1$density, c(1, 2),
                         quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.reg1$density, c(1, 2),
                         quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x3)),
  label = factor(rep(paste("Height = ", grid_x3), each = length(grid_y)),
                 level = rep(paste("Height = ", grid_x3)))
)

ggplot(regplot) + theme_bw() +
  geom_line(data = regplot, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot, map = aes(x = grid, ymin = qlow,
                                        ymax = qupp), fill = "red", alpha = 0.3) +
  facet_wrap(~label, ncol = 5, nrow = 1) +
  labs(x = "log(Age)", y = "Density given Diameter = 0.425")

grid_x2 = quantile(x2, probs = c(0, 0.50, 1))
grid_x2

grid_X = expand.grid(grid_x3, grid_x2)
output <- list(grid_x = grid_X, grid_y = grid_y, out_type = 'FULL', out_param = TRUE)
set.seed(1)
fit.reg1 = PYregression(y = y, x = X, prior = prior, mcmc = mcmc, output = output)
dim(fit.reg1$density)

regplot1 <- data.frame(
  dens = as.vector(apply(fit.reg1$density[, 1:5,], c(1, 2), mean)),
  qlow = as.vector(apply(fit.reg1$density[, 1:5,], c(1, 2),
                         quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.reg1$density[, 1:5,], c(1, 2),
                         quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x3)),
  label = factor(rep(paste("Height = ", grid_x3), each = length(grid_y)),
                 level = rep(paste("Height = ", grid_x3)))
)
regplot2 <- data.frame(
  dens = as.vector(apply(fit.reg1$density[, 6:10,], c(1, 2), mean)),
  qlow = as.vector(apply(fit.reg1$density[, 6:10,], c(1, 2),
                         quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.reg1$density[, 6:10,], c(1, 2),
                         quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x3)),
  label = factor(rep(paste("Height = ", grid_x3), each = length(grid_y)),
                 level = rep(paste("Height = ", grid_x3)))
)
regplot3 <- data.frame(
  dens = as.vector(apply(fit.reg1$density[, 11:15,], c(1, 2), mean)),
  qlow = as.vector(apply(fit.reg1$density[, 11:15,], c(1, 2),
                         quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.reg1$density[, 11:15,], c(1, 2),
                         quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x3)),
  label = factor(rep(paste("Height = ", grid_x3), each = length(grid_y)),
                 level = rep(paste("Height = ", grid_x3)))
)
plot1 = ggplot(regplot1) + theme_bw() +
  geom_line(data = regplot1, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot1, map = aes(x = grid, ymin = qlow,
                                         ymax = qupp), fill = "green", alpha = 0.3) +
  facet_wrap(~label, ncol = 5, nrow = 1) +
  labs(x = "log(Age)", y = "Diameter = 0.055")
plot2 = ggplot(regplot2) + theme_bw() +
  geom_line(data = regplot2, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot2, map = aes(x = grid, ymin = qlow,
                                         ymax = qupp), fill = "red", alpha = 0.3) +
  facet_wrap(~label, ncol = 5, nrow = 1) +
  labs(x = "log(Age)", y = "Diameter = 0.425")
plot3 = ggplot(regplot3) + theme_bw() +
  geom_line(data = regplot3, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot3, map = aes(x = grid, ymin = qlow,
                                         ymax = qupp), fill = "blue", alpha = 0.3) +
  facet_wrap(~label, ncol = 5, nrow = 1) +
  labs(x = "log(Age)", y = "Diameter = 0.650")

grid.arrange(plot1, plot2, plot3, layout_matrix = matrix(1:3, 3, 1))


#Task 7

betas <- fit.reg1$beta


post.betas1 <- matrix(NA, nrow = 1000, ncol = 3)
post.betas2 <- matrix(NA, nrow = 1000, ncol = 3)
post.betas3 <- matrix(NA, nrow = 1000, ncol = 3)
for(i in 1:1000){
  post.betas1[i,] <- betas[[i]][1,]
  post.betas2[i,] <- betas[[i]][2,]
  post.betas3[i,] <- betas[[i]][3,]
}
par(mfrow = c(3,3))
plot(density(post.betas1[,1]))
plot(density(post.betas1[,2]))
plot(density(post.betas1[,3]))
plot(density(post.betas2[,1]))
plot(density(post.betas2[,2]))
plot(density(post.betas2[,3]))
plot(density(post.betas3[,1]))
plot(density(post.betas3[,2]))
plot(density(post.betas3[,3]))



##############################################################################
#Applied lecture 3

library(MASS)
dim(Cars93)

keep = c('Price', 'MPG.city', 'MPG.highway', 'EngineSize',
         'Horsepower', 'RPM', 'Rev.per.mile', 'Fuel.tank.capacity',
         'Length', 'Wheelbase', 'Width', 'Turn.circle', 'Rear.seat.room')
cars93 = Cars93[,keep]
head(cars93,n=2)


#######Multiple Linear Regression

#Estimating the linear model

linear.model <- lm(Price ~ ., data = cars93)
summary(linear.model)

attach(cars93)

logPrice = log(Price)
par(mfrow = c(1,2))
hist(Price, probability = TRUE)
lines(density(Price))
hist(logPrice, probability = TRUE)
lines(density(logPrice))

cars93[, 'Price'] = logPrice
names(cars93)[1] = 'logPrice'
head(cars93, n = 2)

linear.model <- lm(logPrice ~ ., data = cars93)
summary(linear.model)


apply(cars93[,-1], 2, range)
scaled.cars93 = cars93
scaled.cars93[, -1] = scale(cars93[, -1])
linear.model = lm(logPrice ~ ., data = scaled.cars93)
summary(linear.model)

cov93 = cov(scaled.cars93[, -1])
cor93 = cov2cor(cov93)
library(pheatmap)
pheatmap(cor93, display_numbers = TRUE, cluster_rows = FALSE,
         cluster_cols = FALSE, fontsize_number = 10,)
         

#Dealing with multicollinerity

library(car)

vif(linear.model)

linear.model <- lm(logPrice ~. - EngineSize, data = scaled.cars93)
vif(linear.model)

linear.model = lm(logPrice ~ . - EngineSize - MPG.city - Length -
                    Width - Wheelbase, data = scaled.cars93)
vif(linear.model)

summary(linear.model)
Wald <- confint(linear.model)
Wald

#Bootstrap regression confidence intervals

#NPB

set.seed(1)
NPB.results <- Boot(linear.model, method  = 'case', R = 10000)
summary(NPB.results)

normal.NPB = confint(NPB.results, level = 0.95, type = "norm")
basic.NPB = confint(NPB.results, level = 0.95, type = "basic")
percent.NPB = confint(NPB.results, level = 0.95, type = "perc")
bc.NPB = confint(NPB.results, level = 0.95, type = "bca")
normal.NPB
basic.NPB
percent.NPB
bc.NPB

hist(NPB.results, legend = 'top', col = 'red', main = 'Nonparametric Bootstrap')

#SRB

set.seed(1)
SRB.results <- Boot(linear.model, method = 'residual', R = 10000)
summary(SRB.results)

normal.SRB = confint(SRB.results, level= 0.95, type="norm")
basic.SRB = confint(SRB.results, level= 0.95, type="basic")
percent.SRB = confint(SRB.results, level= 0.95, type="perc")
bc.SRB = confint(SRB.results, level= 0.95, type="bca")
hist(SRB.results, legend ="top", col = 'lightgreen',
     main = 'Semiparametric Residual Bootstrap')

#bayesian boot

library(bayesboot)
betas <- function(data){
  coef(lm(logPrice ~ . - EngineSize - MPG.city - Length - Width - Wheelbase,
          data = data))
}

betas(scaled.cars93)

set.seed(1)
BBR <- bayesboot(data = scaled.cars93, statistic = betas, 
                 R = 10000, use.weights = FALSE)

summary(BBR)


#Visualising multiples CI's

library(plotrix)

bayes.MPG <- c(-0.15734644, -0.005605435)

CI.MPG <- rbind(Wald[2,],normal.NPB[2,],basic.NPB[2,],percent.NPB[2,],bc.NPB[2,],
                normal.SRB[2,],basic.SRB[2,],percent.SRB[2,],bc.SRB[2,], bayes.MPG)
rownames(CI.MPG) = c('Wald.MPG', 'normal.NPB.MPG', 'basic.NPB.MPG', 'percent.NPB.MPG', 'bc.NPB.MPG',
                     'normal.SRB.MPG', 'basic.SRB.MPG', 'percent.SRB.MPG', 'bc.SRB.MPG','bayes.MPG')
CI.MPG

#MPG
bayes.mean.MPG <- -0.07975419
point.MPG <- c(rep(coef(linear.model)[2], 9), bayes.mean.MPG)
point.MPG

par(mfrow = c(1,1))

plotCI(x = 1:10, y = point.MPG,li = CI.MPG[,1], ui = CI.MPG[,2], lwd = 2, main = 'MPG.highway',
       col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',xlab = 'Methods', xaxt ='n')
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('topright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)

#Horsepower
bayes.horsepower = c(0.14448541, 0.392370433)
CI.horsepower = rbind(Wald[3,],normal.NPB[3,],basic.NPB[3,],percent.NPB[3,],bc.NPB[3,],
                      normal.SRB[3,],basic.SRB[3,],percent.SRB[3,],bc.SRB[3,],bayes.horsepower)
bayes.mean.horsepower= 0.26680123
point.horsepower = c(rep(coef(linear.model)[3],9), bayes.mean.horsepower)
plotCI(x = 1:10, y = point.horsepower,li = CI.horsepower[,1], ui = CI.horsepower[,2], lwd = 2,
       main = 'Horsepower', col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',
       xlab = 'Methods', xaxt ='n')
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('bottomright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)

#Fuel Tank
bayes.fuel = c(0.03098518, 0.235375479)
CI.fuel = rbind(Wald[6,],normal.NPB[6,],basic.NPB[6,],percent.NPB[6,],bc.NPB[6,],
                normal.SRB[6,],basic.SRB[6,],percent.SRB[6,],bc.SRB[6,],bayes.fuel)
bayes.mean.fuel= 0.13097722
point.fuel = c(rep(coef(linear.model)[6],9), bayes.mean.fuel)
plotCI(x = 1:10, y = point.fuel,li = CI.fuel[,1], ui = CI.fuel[,2], lwd = 2,
       main = 'Fuel Tank Capacity', col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',
       xlab = 'Methods', xaxt ='n')
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'),
     cex.axis = 0.8)
legend('bottomright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15),
       bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)


###A DPM regression model

library(BNPmix)
library(gridExtra)

grid_y = seq(min(logPrice), max(logPrice), length.out = 100)
X = cars93[, c('MPG.highway', 'Horsepower', 'Fuel.tank.capacity')]
grid_x1 = quantile(X[,1], probs = c(0.1, 0.9))
grid_x2 = c(50, 100, 150, 200, 300)
grid_x3 = quantile(X[,3], probs = c(0.50))

grid_X = expand.grid(grid_x1, grid_x2, grid_x3)

DPprior = PYcalibrate(Ek = 10, n = n, discount = 0)
prior = list(strength = DPprior$strength, discount = 0, hyper = FALSE,
             m0 = c(0,0,0,0), a0 = 3, b0 = 100)
output <- list(grid_x = grid_X, grid_y = grid_y, out_type = 'FULL',
               out_param = TRUE)
mcmc.total.n = 11000
mcmc.burnin = 1000
mcmc = list(niter = mcmc.total.n, nburn = mcmc.burnin, method = 'ICS')
set.seed(1)
fit.reg = PYregression(y = logPrice, x = X, prior = prior, mcmc = mcmc,
                       output = output)

summary(fit.reg)

#Multivariate unsupervised learning

Cars93[17,1] = 'Chrysler'
Cars93[,1] = factor(Cars93[,1])
table(Cars93[,1])

Y = cars93[, c('logPrice', 'MPG.highway', 'Horsepower', 'Fuel.tank.capacity')]
grid.length = 20
grid1 = seq(min(Y[,1]), max(Y[,1]), length.out = grid.length)
grid2 = seq(min(Y[,2]), max(Y[,2]), length.out = grid.length)
grid3 = seq(min(Y[,3]), max(Y[,3]), length.out = grid.length)
grid4 = seq(min(Y[,4]), max(Y[,4]), length.out = grid.length)
grid = expand.grid(grid1, grid2, grid3, grid4)
dim(grid)

DPprior = PYcalibrate(Ek = 10, n = n, discount = 0)
DPprior

mcmc.total.n = 5000
mcmc.burnin = 1000
mcmc = list(niter = mcmc.total.n, nburn = mcmc.burnin, method = 'ICS')
prior = list(strength = DPprior$strength, discount = 0, m0 = rep(0,4),
             Sigma0 = diag(1000,4), hyper = FALSE)
output = list(grid = grid, out_param = TRUE)
set.seed(1)
mult.fit = PYdensity(y = Y, mcmc = mcmc, prior = prior, output = output)

summary(mult.fit)

cluster.labels = mult.fit$clust
n.cluster = apply(cluster.labels, 1, unique)
n.cluster = lapply(n.cluster, length)
n.cluster = unlist(n.cluster)
counts = table(n.cluster)
counts

MCMC.allocations = mult.fit$clust
dim(MCMC.allocations)

allocations <- round(apply(MCMC.allocations, 2, mean))
allocations
table(allocations)
