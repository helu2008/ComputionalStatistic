# Hengwei Lu
# A99010013
# Q1:
# A: 
permutation_test <- function(m = 1000, n = 3000) {
  # generate two normal sample
  x <- rnorm(m, 0, 1)
  y <- rnorm(n, 0, 5)
  # D under null hypothesis
  D <- mean(x)-mean(y)
  D.sim <- numeric(10000)
  Z <- c(x,y)
  for (b in 1:10000) {
    Zperm <- sample(Z)
    D.sim[b] = mean(Zperm[1:m]) - mean(Zperm[(m+1):(m+n)])
  }
  hist(D.sim,breaks=50,freq=FALSE, main='Histogram of permutated differences')
  N <- dnorm(seq(-1, 1, 0.01), 0, 1/m+5/n)
  lines(seq(-1, 1, 0.01), N, col='red')
}
permutation_test(1000,3000)
# From the plots, we can see that they do not match well

# B:
permutation_test(2000,2000)
# From the plots, we can see that they do not match well

# Q2:
# A:
# H0: X's and Y's are from the same distribution
# B:
load("cloudseeding.rda")
test_dat <- sort(c(cloudseeding[,1], cloudseeding[,2]), index.return =TRUE)$ix <= 26
library(tseries)
runs.test(as.factor(test_dat))
# p-value = 0.7794 (alternative hypothesis: two.sided)
# So we do not reject Null hyphothesis

# C:
nb.runs.test <- function(x, y, B=999) {
  # Sort x and y
  test_dat <- sort(c(x, y), index.return =TRUE)$ix <= length(x)
  run_num <- Reduce(function(a,b){
    if (a[1] == b)
      a
    else 
      c(b, a[2] + 1)
  }, test_dat, c(test_dat[1],1))[2]
  x_length <- length(x)
  y_length <- length(y)
  # Run test formulas
  mu <- 2 * x_length * y_length / (x_length + y_length) + 1
  sigma <- sqrt((mu - 1) * (mu - 2)/(x_length + y_length - 1))
  num <- rnorm(B, mu , sigma)
  # Use Monte Carlo similation  
  if(run_num > 2 * mu - run_num){
    p1 <- (sum(num > run_num)+1)/(B+1)
    p2 <- (sum(num < (2 * mu - run_num))+1)/(B+1)
    p <- p1 + p2
  }
  else {
    p1 <- (sum(num < run_num)+1)/(B+1)
    p2 <- (sum(num > (2 * mu - run_num))+1)/(B+1)
    p <- p1 + p2
  }
  return(p)
}

nb.runs.test(cloudseeding[,1],cloudseeding[,2], B = 10000)
# p-value = 0.7754225, almost the same

# D:
nb.runs.test_onesided <- function(x, y, B=999) {
  # Sort x and y
  test_dat <- sort(c(x, y), index.return =TRUE)$ix <= length(x)
  run_num <- Reduce(function(a,b){
    if (a[1] == b)
      a
    else 
      c(b, a[2] + 1)
  }, test_dat, c(test_dat[1],1))[2]
  # run test formulas
  x_length <- length(x)
  y_length <- length(y)
  mu <- 2 * x_length * y_length / (x_length + y_length) + 1
  sigma <- sqrt((mu - 1) * (mu - 2)/(x_length + y_length - 1))
  num <- rnorm(B, mu , sigma)
  # MC similation
  p <- (sum(num > run_num)+1)/(B+1)
  return(p)
}

nb.runs.test_onesided(cloudseeding[,1],cloudseeding[,2], B = 10000)
# p = 0.6147385, almost the same as twice of (1- value from last question)

# Q3
nb.runs.sym.test <- function(x, B=999) {
  # Sort positive or negative numbers in x
  test_dat <- x > 0
  run_num <- Reduce(function(a,b){
    if (a[1] == b)
      a
    else 
      c(b, a[2] + 1)
  }, test_dat, c(test_dat[1],1))[2]
  # Run test formulas
  x_length <- sum(test_dat)
  y_length <- sum(!test_dat)
  mu <- 2 * x_length * y_length / (x_length + y_length) + 1
  sigma <- sqrt((mu - 1) * (mu - 2)/(x_length + y_length - 1))
  num <- rnorm(B, mu , sigma)
  # Use Monte Carlo similation  
  if(run_num > 2 * mu - run_num){
    p1 <- (sum(num > run_num)+1)/(B+1)
    p2 <- (sum(num < (2 * mu - run_num))+1)/(B+1)
    p <- p1 + p2
  }
  else {
    p1 <- (sum(num < run_num)+1)/(B+1)
    p2 <- (sum(num > (2 * mu - run_num))+1)/(B+1)
    p <- p1 + p2
  }
  return(p)
}

# Test
# 2 normal sample
# Null hypothesis: x and y are from the same distribution
x <- rnorm(1000, 0, 1)
y <- rnorm(1000, -100, 6)
nb.runs.sym.test((x-y), B = 10000)
# p-value = 0.00019998, so we reject null hypothesis

# 2 normal sample
# Null hypothesis: x and y are from the same distribution
x <- rnorm(1000, 5, 1)
y <- rnorm(1000, 5, 1)
nb.runs.sym.test((x-y), B = 10000)
# p-value = 0.805, so we accept null hypothesis