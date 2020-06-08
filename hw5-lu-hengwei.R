# Hengwei Lu
# A99010013

# Q1
# Test Pearson test
ps_test <- function(B = 999){
  t <- c()
  for (i in 1: B){
    X <- runif(10000, min = -1, max = 1)
    Y <- X ^ 2
    t <- c(t, cor.test(X, Y)$p.value)
  }
  return(t)
}
t <- ps_test(999)

# Test Spearman test
sprm_test <- function(B = 999){
  S <- c()
  for (i in 1: B){
    X <- runif(10000, min = -1, max = 1)
    Y <- X ^ 2
    S <- c(S, cor.test(X, Y,  method="spearman", exact = FALSE)$p.value)
  }
  return(S)
}
S <- sprm_test(999)

# Test Kendall test
kd_test <- function(B = 999){
  z <- c()
  for (i in 1: B){
    X <- runif(10000, min = -1, max = 1)
    Y <- X ^ 2
    z <- c(z, cor.test(X, Y,  method="kendall", alternative = "two.sided")$p.value)
  }
  return(z)
}
z <- kd_test(10)

# make side-by-side boxlpot
boxplot(t,S,z)
# Median is around 0.4. But desending.
# From the plot, we can see that kendall is the most powerful and pearson is the least.

# Q2
hoeff.test <- function(z, B = 999){
  d0 <- max(hoeffd(z)$maxad)
  d <- c()
  # bootstrap
  for (i in 1:B){
    X <- z[,1]
    Y <- sample(z[,2], nrow(z))
    d <- c(d, max(hoeffd(X, Y)$maxad))
  }
  # return p-value
  return((sum(d > d0)+1)/(B+1))
}

# Q3
# A
attach(read.table("sea_ice_data.txt", header = TRUE))
# Linear regression
model1 <- lm(Ice ~ Year)
plot(Year, Ice, pch=16)
abline(model1, col='red', lwd=2)
model1_r2 <- summary(model1)$r.squared

# B
# degree 2 polynomial regression
model2 <- lm(Ice ~ Year + I(Year^2))
plot(Year, Ice, pch=16)
lines(Year, predict(model2, data.frame(Year = Year)), col='red', lwd=2)
model2_r2 <- summary(model2)$r.squared

# degree 3 polynomial regression
model3 <- lm(Ice ~ poly(Year, 3))
plot(Year, Ice, pch=16)
lines(Year, predict(model3, data.frame(Year = Year)), col='red', lwd=2)
model3_r2 <- summary(model3)$r.squared

# degree 4 polynomial regression
model4 <- lm(Ice ~ poly(Year, 4))
plot(Year, Ice, pch=16)
lines(Year, predict(model4, data.frame(Year = Year)), col='red', lwd=2)
model4_r2 <- summary(model4)$r.squared

# degree 3 and 4 almost have the same r^2 so we stop

# C
library(stats)
# monotone regression
mono_model <- loess(Ice ~ Year)
plot(Year, Ice, pch=16)
lines(Year, predict(mono_model, Year), col='red', lwd=2)

# D
par(new=TRUE)
plot(Year, Ice, pch=16)
# plot each line
lines(Year, predict(model1, data.frame(Year = Year)), col='red', lwd=2)
lines(Year, predict(model2, data.frame(Year = Year)), col='purple', lwd=2)
lines(Year, predict(model3, data.frame(Year = Year)), col='blue', lwd=2)
lines(Year, predict(model4, data.frame(Year = Year)), col='green', lwd=2)
lines(Year, predict(mono_model, Year), col='yellow', lwd=2)
# Signals
legend('topright', 
       legend = c("linear", "quardratic", "cubic", "quartic", "mono"),
       col = c("red", "purple", "blue", "green", "yellow"),
       lty = 1)

# Q4
boot.regression <- function(x, y, conf=0.95, residual=F, B=999) {
  # compute the linear model of the original data
  model <- lm(y~x)
  beta <- model$coefficients
  b0 <- beta[1]
  b1 <- beta[2]
  n <- length(x)
  # function to compute linear model from bootstrap samples 
  f0 <- function(ix) {
    l <- sample(ix, n, replace=T)
    xb <- x[l]
    yb <- y[l]
    model <- lm(yb ~ xb)
    beta <- model$coefficients
    res <- model$residuals
    b0h <- beta[1]
    b1h <- beta[2]
    # return beta0, beta1, residuals, bootstrap index
    return(c(b0h, b1h, res, l))
  }
  # function to compute standard error of beta
  f1 <- function(xb, res) {
    xbar <- mean(xb)
    ssx <- var(xb)
    sigma <- sqrt(Reduce(function(x, y){x + y^2}, res, 0)/(n-2))
    seb0h <- sigma*sqrt(1/n+xbar^2/ssx)
    seb1h <- sigma*sqrt(1/ssx)
    # return standard error of beta0 and beta1
    return(c(seb0h, seb1h))
  }
  # function to compute t statistics
  boot_t <- function(ix) {
    b <- f0(ix)
    b0h <- b[1]
    b1h <- b[2]
    res <- b[3:(n+2)]
    l <- b[-(1:(n+2))]
    if (residual) {
      bb <- sapply(1:B, function(x){f0(l)})
      seb0h <- sqrt(var(bb[1,]))
      seb1h <- sqrt(var(bb[2,]))
    }else {
      se <- f1(x[l], res)
      seb0h <- se[1]
      seb1h <- se[2]
    }
    return(c((b0h- b0)/seb0h, (b1h - b1)/seb1h, b0h, b1h))
  }
  # run function for B times
  t <- sapply(1:B, function(x){boot_t(1:n)})
  t0 <- t[1,]
  t1 <- t[2,]
  bb0 <- t[3,]
  bb1 <- t[4,]
  alpha <- 1 - conf
  # sort t statistics and get the corresponding quantile
  t0ha <- sort(t0)[as.integer(B*alpha/2)+1]
  t1ha <- sort(t1)[as.integer(B*alpha/2)+1]
  # return confidence interval of beta0 and beta1
  return(c(b0+t0ha*sqrt(var(bb0)), b0-t0ha*sqrt(var(bb0)),
           b1+t1ha*sqrt(var(bb1)), b1-t1ha*sqrt(var(bb1))))
}