# 1
# Calculate original FT
twowayTest <- function(dat){
  # Set Treatment and Block groups
  Trt <- as.factor(dat[, 2])
  Blk <- as.factor(dat[, 3])
  I <- length(levels(Trt))
  J <- length(levels(Blk))
  n <- length(dat[, 1])/I/J
  # Calculate grand mean
  grd <- mean(dat[, 1])
  # Calculate treat sum of squares
  SST <- J * n * sum((tapply(dat[, 1], Trt, mean) - grd)^2)
  # Calculate error sum of squares
  SSE <- sum((dat[, 1] - grd)^2)
  # FT under null hypothesis
  FT <- (SST/(I-1))/(SSE/((n-1)*J*I))
  return(FT)
}

# Permutation test
twowayPermTest <- function(dat, B=999){
  # Calculate original FT
  FT <- twowayTest(dat)
  n <- 0
  # Permutation
  for (i in 1:B){
    dat[, 2] <- sample(dat[, 2], length(dat[, 2]), FALSE)
    FT_new <- twowayTest(dat)
    if (FT_new > FT){
      n <- n + 1
    }
  }
  # Monte Carlo estimate
  p <- (n+1)/(B+1)
  return(p)
}

twowayPermTest(ToothGrowth, 2000)
# p-value = 0.0519, so it's unlikely that the null hypothesis is true

# 2
load("alon.RData")
hlth <- alon$x[alon$y == 'n',]
tumor <- alon$x[alon$y == 't',]
p_list <- c()
t <- c()
for (i in 1 : ncol(alon$x)){
  num <- t.test(hlth[,i], tumor[,i]) $ p.value
  t <- c(t, abs(t.test(hlth[,i], tumor[,i]) $ statistic))
  # Calculate two-sided p value of t test
  num <- min(num, 1-num) * 2
  p_list <- c(p_list, num)
}

# Apply all methods for controlling FWER
p_list_bon = p.adjust(p_list, "bon")
p_list_holm = p.adjust(p_list, "holm")
p_list_hoch = p.adjust(p_list, "hoch") 
# Try setting alpha = 0.05
sum(p_list_bon<0.05)
sum(p_list_holm<0.05)
sum(p_list_hoch<0.05)
# They all reject 8 instances
# Apply all methods for controlling FDR
p_list_bh = p.adjust(p_list, "BH")
sum(p_list_bh<0.05)
# It rejects 118 instances

# 3
B <- 2000
count <- rep(0, B)
for (i in 1 : B){
  y_new <- sample(alon$y, length(alon$y), FALSE)
  hlth <- alon$x[y_new == 'n',]
  tumor <- alon$x[y_new == 't',]
  t_perm <- c()
  # Calculate t test statistics
  for (j in 1 : ncol(alon$x)){
    t_perm <- c(t_perm, abs(t.test(hlth[,j], tumor[,j]) $ statistic))
  }
  count <- count + t_perm > t
}
p_perm <- (count+1)/(B+1)

# Apply all methods for controlling FWER
p_perm_bon = p.adjust(p_perm, "bon")
p_perm_holm = p.adjust(p_perm, "holm")
p_perm_hoch = p.adjust(p_perm, "hoch") 
# Try setting alpha = 0.05
sum(p_perm_bon<0.05)
sum(p_perm_holm<0.05)
sum(p_perm_hoch<0.05)
# Apply all methods for controlling FDR
p_perm_bh = p.adjust(p_perm, "BH")
sum(p_perm_bh<0.05)