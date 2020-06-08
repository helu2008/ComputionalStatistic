# Hengwei Lu
# A99010013
# HW2:

# Q1a: Because the null hypothesis assumes X1,...,Xn are from a 
#      distribution which is symmetric, (-Xn) and Xn are from 
#      the same distribution. For every Xn, there are two choice: 
#      (-Xn) or Xn. So in total there are 2^n choices.

# Q1b:
flipSignTest1 <- function(x, B = 999) {
  count <- 0
  m <- mean(x)
  for (i in 1 : B){
    mean_test <- mean(x*(sample(c(1, -1), length(x), TRUE)))
    if (mean_test > m) {
      count <- count + 1
    }
  }
  return((count + 1)/(B + 1))
}

# Q1c:
flipSignTest2 <- function(x, B = 999) {
  count <- 0
  m <- abs(mean(x))
  for (i in 1 : B){
    mean_test <- mean(x*(sample(c(1, -1), length(x), TRUE)))
    if (mean_test >= m || mean_test <= -m) {
      count <- count + 1
    }
  }
  return((count + 1)/(B + 1))
}

# Q2a:
data <- father.son[,2] - father.son[,1]
hist(data, 100)
# It seems that the axis of symmetry is larger than 0, 
# so sons tend to be higher than fathers.

# Q2b:
flipSignTest1(data)
# p = 0.001. It's almost zero. So it is impossible that data is a 
# symmetric distribution.
# This test is in fact a permutation test because if we change x to 
# -x, where x is the height of son minus the height of father, we are
# permuting the data of father and son. 

# Q3:
one_side_bootstrap <- function(x, func = mean, alpha = 0.1, B = 1000){
  theta_hat <- c()
  for (i in 1:B) {
    theta_hat <- c(theta_hat, func(sample(x, length(x), TRUE)))
  }
  theta_hat <- sort(theta_hat)
  CI_index <- B * alpha
  return(theta_hat[floor(CI_index)] + 
           (theta_hat[ceil(CI_index)] - theta_hat[floor(CI_index)]) * 
           (CI_index - floor(CI_index)))
}
one_side_bootstrap(data)
# The result of the test is 0.8851017, so 90% of the differences
# between heights are larger than 0.885, which is larger than 0. 
# We use bootstrap so we draw samples from sample distribution.