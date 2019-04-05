<!--
```r
##########Project 5#############
set.seed(1)
X <- rnorm(m, 0, sigma); nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
while (m < (qt(1 - alpha/2, m - 1)/d)^2 * var(X)) {
  X <- c(X, rnorm(1, 0, sigma))
  m <- m + 1
}

library(doParallel)
library(foreach)


h <- function(k, x){
  if (k == 1) 1 else 
    sum((x - a[k])^(2:k - 1)/factorial(2:k - 1) * sapply(k:2 - 1, function(i) h(i, a[k])))
}
Ginf <- function(k) exp(-a[k]) * sum(sapply(1:k, function(i) h(i, a[k])))

### more precise ###
h <- function(k, x){
  if (k == 1) 1 else 
    sum(exp((2:k - 1) * log(x - a[k]) - sum(log(2:k - 1))) * sapply(k:2 - 1, function(i) h(i, a[k])))
}
Ginf <- function(k) exp(-a[k]) * sum(sapply(1:k, function(i) h(i, a[k])))
system.time(sapply(1:20, function(i) Ginf(i)))

### parallel computation ###
{
cl <- makeCluster(detectCores())      
registerDoParallel(cl)       
getDoParWorkers()

system.time(G <- foreach(i = 1:20, .combine = "c") %dopar% Ginf(i))

stopImplicitCluster()
stopCluster(cl)
}

m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.05

pureSeqEst <- function(m, d, sigma, alpha){
  A <- (qnorm(1 - alpha/2) * sigma/d)^2; nopt <- ceiling(A)
  
  kcandi <- ceiling((m - 1)/2):(nopt + 4)
  a <- (kcandi - 1) * (2 * kcandi - 1)/A
  
  dist <- -diff(sapply(kcandi, function(m) Ginf(m))) 
  
  Ncandi <- 2 * kcandi[1:25] + 1
  ch <- t(dist) %*% cbind(Ncandi, Ncandi^2, pnorm(d * sqrt(Ncandi)/sigma))
  
  list(DistN = dist, EN = ch[1], SigN = sqrt(ch[2] - ch[1]^2), CovProb = 2 * ch[3] - 1,
       hatNopt = N, Nopt = nopt, SampMean = mean(X))
}

exp(-a[2]) - exp(-a[3]) * (1 + a[3] - a[2])
p2 <- pgamma(a[3] - a[2], shape = 2) * factorial(1)
p1 <- pgamma(a[2] - a[1], shape = 1) * factorial(0)
(1 - p1) * p2


> dist
[1] 1.681932e-01 5.796593e-02 5.113208e-02 5.775334e-02
[5] 7.035698e-02 8.518424e-02 9.770710e-02 1.028833e-01
[9] 9.726518e-02 8.120275e-02 5.910416e-02 3.712438e-02
[13] 1.995566e-02 9.116229e-03 3.518392e-03 1.141433e-03
[17] 3.098972e-04 7.013897e-05 1.318799e-05


##########Project 5#############
set.seed(1)
X <- rnorm(m, 0, sigma); nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
while (m < (qt(1 - alpha/2, m - 1)/d)^2 * var(X)) {
  X <- c(X, rnorm(1, 0, sigma))
  m <- m + 1
}

library(doParallel)
library(foreach)

### maybe no precise ###
h <- function(k, x){
  if (k == 1) 1 else 
    sum((x - a[k])^(2:k - 1)/factorial(2:k - 1) * sapply(k:2 - 1, function(i) h(i, a[k])))
}

### more precise ###
h <- function(k, x){
  if (k == 1) 1 else 
    sum(exp((2:k - 1) * log(x - a[k]) - log(factorial(2:k - 1))) * sapply(k:2 - 1, function(i) h(i, a[k])))
}

### sapply mathod + more precise ###
h <- function(k, x){
  if (k == 1) 1 else
    sum(sapply(2:k - 1, function(i) {
      exp(i * log(x - a[k]) - log(factorial(i))) * h(k - i, a[k])
    }))
}


Ginf <- function(k) exp(-a[k]) * sum(sapply(1:k, function(i) h(i, a[k])))
system.time(print(-diff(sapply(1:20, function(i) Ginf(i)))))




### parallel computation ###
{
  cl <- makeCluster(detectCores())      
  registerDoParallel(cl)       
  getDoParWorkers()
  
  system.time(print(-diff(foreach(i = 1:20, .combine = "c") %dopar% Ginf(i))))
  
  stopImplicitCluster()
  stopCluster(cl)
}

m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.05

pureSeqEst <- function(m, d, sigma, alpha){
  A <- (qnorm(1 - alpha/2) * sigma/d)^2; nopt <- ceiling(A)
  
  kcandi <- ceiling((m - 1)/2):(nopt + 4)
  a <- (kcandi - 1) * (2 * kcandi - 1)/A
  
  dist <- -diff(sapply(kcandi, function(m) Ginf(m))) 
  dist <- -diff(foreach(i = 1:25, .combine = "c") %dopar% Ginf(i))
  
  Ncandi <- 2 * kcandi[1:19] + 1
  ch <- t(dist) %*% cbind(Ncandi, Ncandi^2, pnorm(d * sqrt(Ncandi)/sigma))
  
  list(DistN = dist, EN = ch[1], SigN = sqrt(ch[2] - ch[1]^2), CovProb = 2 * ch[3] - 1,
       hatNopt = N, Nopt = nopt, SampMean = mean(X))
}

Ix <- function(x, p) pgamma(x, shape = p + 1)
p1 <- Ix(a[2], 0)
p2 <- Ix(a[3] - a[2], 1)
B2 <- 1 - Ix(a[3] - a[2], 1)
p3 <- 1/(2 * B2 * exp(a[3] - a[2])) * (2 * Ix(a[4] - a[3], 2) + 2 * (a[3] - a[2]) * Ix(a[4] - a[3], 1))
p1
(1 - p1) * p2
(1 - p1) * (1 - p2) * p3






p3 <- 1/(2 * B2 * exp(a[2] - a[1])) * (2 * pgamma(a[3] - a[2], shape = 2) * factorial(1) + 2 * (a[2] - a[1]) * pgamma(a[3] - a[2], shape = 1))

(1 - p1) * (1 - p2) * p3


############################Revise at 04/02/2019###########################################333333
##########Project 5#############
set.seed(1)
X <- rnorm(m, 0, sigma); nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
while (m < (qt(1 - alpha/2, m - 1)/d)^2 * var(X)) {
  X <- c(X, rnorm(1, 0, sigma))
  m <- m + 1
}

library(doParallel)
library(foreach)


### maybe no precise ###
h <- function(k, x){
  if (k == 1) 1 else 
    sum((x - a[k])^(2:k - 1)/factorial(2:k - 1) * sapply(k:2 - 1, function(i) h(i, a[k])))
}

### more precise ###

k = 6;x = a[6]
h(6, a[6])
h <- function(k, x){
  if (k == 1) 1 else 
    sum(exp((2:k - 1) * log(x - a[k]) - log(factorial(2:k - 1))) * sapply(k:2 - 1, function(i) h(i, a[k])))
}
h(2, a[2])
a[7]
### sapply mathod + more precise ###
h <- function(k, x){
  if (k == 1) 1 else
    sum(sapply(2:k - 1, function(i) {
      exp(i * log(x - a[k]) - log(factorial(i))) * h(k - i, a[k])
    }))
}
Ginf <- function(k) exp(-a[k]) * sum(sapply(1:k, function(i) h(i, a[k])))

k = 5
h(3, a[5])

### parallel computation ###
{
  cl <- makeCluster(detectCores())      
  registerDoParallel(cl)       
  getDoParWorkers()
  
  system.time(print(-diff(foreach(i = 1:10, .combine = "c") %dopar% Ginf(i))))
  
  stopImplicitCluster()
  stopCluster(cl)
}
system.time(print(-diff(sapply(1:10, function(i) Ginf(i)))))
system.time(sapply(1:10, function(i) Ginf(i)))
-diff(sapply(1:10, function(i) Ginf(i)))


m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.05

pureSeqEst <- function(m, d, sigma, alpha){
  nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
  
  kcandi <- ceiling((m - 1)/2):(nopt + 10)
  a <- (kcandi - 1) * (2 * kcandi - 1) * (d/(qt(1 - alpha/2, 2 * kcandi) * sigma))^2
  {# A <- (qnorm(1 - alpha/2) * sigma/d)^2; nopt <- ceiling(A)
  # kcandi <- ceiling((m - 1)/2):(nopt + 4)
  # a <- (kcandi - 1) * (2 * kcandi - 1)/A
  }
  
  dist <- -diff(sapply(1:20, function(m) Ginf(m)))
  dist <- -diff(foreach(i = 1:20, .combine = "c") %dopar% Ginf(i))
  
  
  sum(dist)
  Ncandi <- 2 * kcandi[1:19] + 1
  ch <- t(dist) %*% cbind(Ncandi, Ncandi^2, 2 * pnorm(d * sqrt(Ncandi)/sigma) - 1)
  ch
  list(DistN = dist, EN = ch[1], SigN = sqrt(ch[2] - ch[1]^2), CovProb = 2 * ch[3] - 1,
       hatNopt = N, Nopt = nopt, SampMean = mean(X))
}

exp(-a[2]) - exp(-a[3]) * (1 + a[3] - a[2])
p2 <- pgamma(a[3] - a[2], shape = 2) * factorial(1)
p1 <- pgamma(a[2] - a[1], shape = 1) * factorial(0)
(1 - p1) * p2
h(5, a[6])

Ix <- function(x, p) pgamma(x, shape = p + 1)
p1 <- Ix(a[2], 0)
p2 <- Ix(a[3] - a[2], 1)
B2 <- 1 - Ix(a[3] - a[2], 1)
p3 <- 1/(2 * B2 * exp(a[3] - a[2])) * (2 * Ix(a[4] - a[3], 2) + 2 * (a[3] - a[2]) * Ix(a[4] - a[3], 1))
p1
(1 - p1) * p2
(1 - p1) * (1 - p2) * p3



```

############################## Revise at 04/03/2019 ###################################
set.seed(1)
X <- rnorm(m, 0, sigma); nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
while (m < (qt(1 - alpha/2, m - 1)/d)^2 * var(X)) {
  X <- c(X, rnorm(1, 0, sigma))
  m <- m + 1
}

library(doParallel)
library(foreach)


### maybe no precise ###
h <- function(k, x){
  if (k == 1) 1 else 
    sum((x - a[k])^(2:k - 1)/factorial(2:k - 1) * sapply(k:2 - 1, function(i) h(i, a[k])))
}

### more precise ###
log(factorial(1000))
sum(sapply(1:1000, function(i) log(i)))
h <- function(k, x){
  if (k == 1) 1 else 
    sum(exp((2:k - 1) * log(x - a[k]) - log(factorial(2:k - 1))) * sapply(k:2 - 1, function(i) h(i, a[k])))
}

### sapply mathod + more precise ###
h <- function(k, x){
  if (k == 1) 1 else
    sum(sapply(2:k - 1, function(i) {
      exp(i * log(x - a[k]) - log(factorial(i))) * h(k - i, a[k])
    }))
}
Ginf <- function(k) exp(-a[k]) * sum(sapply(1:k, function(i) h(i, a[k])))

### parallel computation ###
{
  cl <- makeCluster(detectCores())      
  registerDoParallel(cl)       
  getDoParWorkers()
  
  system.time(print(-diff(foreach(i = 1:10, .combine = "c") %dopar% Ginf(i))))
  
  stopImplicitCluster()
  stopCluster(cl)
}


m <- 2; d <- 0.8; sigma <- 1; alpha <- 0.05

pureSeqEst <- function(m, d, sigma, alpha){
  nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
  
  kcandi <- ceiling((m - 1)/2):(nopt + 25)
  a <- (kcandi - 1) * (2 * kcandi - 1) * c(0, (d/(qt(1 - alpha/2, 2 * kcandi[-1] - 2) * sigma))^2)
  dist <- -diff(foreach(i = 1:20, .combine = "c") %dopar% Ginf(i))
  
  
  sum(dist)
  Ncandi <- 2 * kcandi[1:19] + 1
  ch <- t(dist) %*% cbind(Ncandi, Ncandi^2, pnorm(d * sqrt(Ncandi)/sigma))
  c(EN = ch[1], SigN = sqrt(ch[2] - ch[1]^2), CovProb = 2 * ch[3] - 1)

}
########################## Revise at 4/4/2019##############################
######################## for loop ############################
h <- function(k, x){
  ifelse(k == 1, 1, ifelse(k == 2, x - a[k], {
    coef <- 1; mat <- matrix(1, k - 2, k - 1)
    hin <- function(j, y){
      sum(exp(((j - 1):1) * log(y - a[j]) - sapply((j - 1):1, function(l) sum(log(1:l)))) * coef)
    }
    for (i in 2:(k - 1)) {
      mat[(i - 1):(k - 2), i] <- sapply((i + 1):k, function(m) hin(i, a[m]))
      coef <- mat[i - 1, 1:i]
    }
    hin(k, x)
  }))
}
Ginf <- function(k) exp(-a[k]) * sum(sapply(1:k, function(i) h(i, a[k])))




pureSeqEst <- function(m, d, sigma, alpha){
  nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
  
  kcandi <- ceiling((m - 1)/2):(nopt + 20)
  a <- (kcandi - 1) * (2 * kcandi - 1) * c(0, (d/(qt(1 - alpha/2, 2 * kcandi[-1] - 2) * sigma))^2)
  
  dist <- -diff(foreach(i = kcandi, .combine = "c") %dopar% Ginf(i))
  dist <- -diff(sapply(kcandi, function(i) Ginf(i)))

  Ncandi <- 2 * kcandi[-length(kcandi)] + 1
  ch <- t(dist) %*% cbind(Ncandi, Ncandi^2, 2 * pnorm(d * sqrt(Ncandi)/sigma) - 1)
  c(EN = ch[1], SigN = sqrt(ch[2] - ch[1]^2), CovProb = ch[3], nOpt = nopt, cdf = sum(dist)) %>% round(3)
}
#### (a) ####
library(dplyr)
m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.05
pureSeqEst(2, 0.5, 1, 0.05) %>% round(3)

m <- 2; d <- 0.3; sigma <- 1; alpha <- 0.05
pureSeqEst(2, 0.3, 1, 0.05) %>% round(3)

m <- 2; d <- 0.5; sigma <- 2; alpha <- 0.05
pureSeqEst(2, 0.5, 2, 0.05) %>% round(3)

m <- 2; d <- 0.3; sigma <- 2; alpha <- 0.05
pureSeqEst(2, 0.3, 2, 0.05) %>% round(3)

m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.1
pureSeqEst(2, 0.5, 1, 0.1) %>% round(3)

m <- 2; d <- 0.3; sigma <- 1; alpha <- 0.1
pureSeqEst(2, 0.3, 1, 0.1) %>% round(3)


############### AR(2) process ##################

########2#########
## a ##
set.seed(3)
X <- arima.sim(list(ar = c(1.5, -0.75)), n = 50) + 400
plot(X[1:30] ~ c(1:30), pch = 19, xlab = "k", ylab = expression(X[k]), xlim = c(1,50), ylim = range(X) + c(-5, 5))
points(X[31:50] ~ c(31:50), pch = 1)


## b ##
matpow <- function(mat, n) if (n == 0) diag(dim(mat)[1]) else mat %*% matpow(mat, n - 1)
matpow(Fmat, 100)

Fmat <- matrix(c(1.5, 1, -0.75, 0), 2)
Xkl <- 400 + sapply(1:20, function(l) matpow(Fmat, l) %*% matrix(X[30:29] - 400))[1, ]
points(Xkl ~ c(31:50), pch = 0)


## c ##
sigm0 <- matrix(c(1, 0, 0, 0), 2)
mse <- function(l) if (l == 1) sigm0 else sigm0 + Fmat %*% mse(l - 1) %*% t(Fmat)
RMSE <- sqrt(unlist(lapply(seq(2, 20, 2), function(l) rev(diag(mse(l))))))



RMSE <- sqrt(cumsum(sapply(0:19, function(j) matpow(Fmat, j)[1, 1]^2)))
d <- RMSE %*% t(qnorm(c(0.95, 0.975, 0.995)))
sapply(c(-1, 1), function(i) {
  points(c(31:50), Xkl + i * d[, 1], type = "l", lty = 2, col = "tomato3")
  points(c(31:50), Xkl + i * d[, 2], type = "l", lty = 3, col = "cyan3")
  points(c(31:50), Xkl + i * d[, 3], type = "l", lty = 4, col = "green3")
})
legend("topleft", text.col = c("tomato3", "cyan3", "green3"), bty = "n", lty = 2:4, col 
       = c("tomato3", "cyan3", "green3"), legend = c("90% CI", "95% CI", "99% CI"))
sapply(1:3, function(i) sum(abs(X[31:50] - Xkl) > d[, i]))

###################################### Finally ##########################################################
library(doParallel)
library(foreach)
### parallel computation ###
{
  cl <- makeCluster(detectCores())      
  registerDoParallel(cl)       
  getDoParWorkers()
  
  system.time(print(-diff(foreach(i = 1:10, .combine = "c") %dopar% Ginf(i))))
  
  stopImplicitCluster()
  stopCluster(cl)
}

h <- function(k, x){
  ifelse(k == 1, 1, ifelse(k == 2, x - a[k], {
    coef <- 1; mat <- matrix(1, k - 2, k - 1)
    hin <- function(j, y){
      sum(exp(((j - 1):1) * log(y - a[j]) - sapply((j - 1):1, function(l) sum(log(1:l)))) * coef)
    }
    for (i in 2:(k - 1)) {
      mat[(i - 1):(k - 2), i] <- sapply((i + 1):k, function(m) hin(i, a[m]))
      coef <- mat[i - 1, 1:i]
    }
    hin(k, x)
  }))
}
Ginf <- function(k) exp(-a[k]) * sum(sapply(1:k, function(i) h(i, a[k])))



pureSeqEst <- function(m, d, sigma, alpha){
  h <- function(k, x){
    ifelse(k == 1, 1, ifelse(k == 2, x - a[k], {
      coef <- 1; mat <- matrix(1, k - 2, k - 1)
      hin <- function(j, y){
        sum(exp(((j - 1):1) * log(y - a[j]) - sapply((j - 1):1, function(l) sum(log(1:l)))) * coef)
      }
      for (i in 2:(k - 1)) {
        mat[(i - 1):(k - 2), i] <- sapply((i + 1):k, function(m) hin(i, a[m]))
        coef <- mat[i - 1, 1:i]
      }
      hin(k, x)
    }))
  }
  Ginf <- function(k) exp(-a[k]) * sum(sapply(1:k, function(i) h(i, a[k])))
  
  
  
  nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
  
  kcandi <- ceiling((m - 1)/2):(nopt + 20)
  a <- (kcandi - 1) * (2 * kcandi - 1) * c(0, (d/(qt(1 - alpha/2, 2 * kcandi[-1] - 2) * sigma))^2)
  
  dist <- -diff(foreach(i = kcandi, .combine = "c") %dopar% Ginf(i))
  #dist <- -diff(sapply(kcandi, function(i) Ginf(i)))
  
  Ncandi <- 2 * kcandi[-length(kcandi)] + 1
  ch <- t(dist) %*% cbind(Ncandi, Ncandi^2, 2 * pnorm(d * sqrt(Ncandi)/sigma) - 1)
  list(EN = ch[1], SigN = sqrt(ch[2] - ch[1]^2), CovProb = ch[3], nOpt = nopt, 
       cdf = sum(dist), DistN = dist, Ncan = Ncandi) 
}
#### (a) ####

library(dplyr)
m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.05
sec1 <- pureSeqEst(2, 0.5, 1, 0.05) %>% round(3)
#     EN    SigN CovProb    nOpt     cdf 
#16.726   6.291   0.931  16.000   1.000 


m <- 2; d <- 0.3; sigma <- 1; alpha <- 0.05
sec2 <- pureSeqEst(2, 0.3, 1, 0.05) %>% round(3)
#     EN    SigN CovProb    nOpt     cdf 
#43.502  11.445   0.935  43.000   1.000


m <- 2; d <- 0.5; sigma <- 2; alpha <- 0.05
sec3 <- pureSeqEst(2, 0.5, 2, 0.05) %>% round(3)
#     EN    SigN CovProb    nOpt     cdf 
#62.399  13.501   0.940  62.000   1.000

m <- 2; d <- 0.3; sigma <- 2; alpha <- 0.05
sec4 <- pureSeqEst(2, 0.3, 2, 0.05) %>% round(3)
#     EN    SigN CovProb    nOpt     cdf 
#171.930  21.523   0.947 171.000   1.000


m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.1
sec5 <- pureSeqEst(2, 0.5, 1, 0.1) %>% round(3)
#     EN    SigN CovProb    nOpt     cdf 
#11.827   5.028   0.879  11.000   1.000


m <- 2; d <- 0.3; sigma <- 1; alpha <- 0.1
sec6 <- pureSeqEst(2, 0.3, 1, 0.1) %>% round(3)
#     EN    SigN CovProb    nOpt     cdf 
#29.927  10.031   0.871  31.000   1.000


#### (b) ####
par(mai = c(0.9, 0.9, 0.3, 0.2))
plot(sec1$DistN ~ c(sec1$Ncan - 2), xlim = range(sec2$Ncan) - 2, type = "p", pch = 8, col = "tomato3", 
     xlab = "k", ylab = expression(P[mu][","][sigma](N == m+k)))
points(sec2$DistN ~ c(sec2$Ncan - 2), type = "p", pch = 18, col = "cyan3")
legend("topright", legend = c("d = 0.5", "d = 0.3"), pch = c(8, 18), 
       col = c("tomato3", "cyan3"), text.col = c("tomato3", "cyan3"))



plot(sec3$DistN ~ c(sec3$Ncan - 2), xlim = range(sec4$Ncan) - 2, type = "p", pch = 8, col = "tomato3", 
     xlab = "k", ylab = expression(P[mu][","][sigma](N == m+k)))
points(sec4$DistN ~ c(sec4$Ncan - 2), type = "p", pch = 18, col = "cyan3")
legend("topright", legend = c("d = 0.5", "d = 0.3"), pch = c(8, 18), 
       col = c("tomato3", "cyan3"), text.col = c("tomato3", "cyan3"))


plot(sec5$DistN ~ c(sec5$Ncan - 2), xlim = range(sec6$Ncan) - 2, type = "p", pch = 8, col = "tomato3", 
     xlab = "k", ylab = expression(P[mu][","][sigma](N == m+k)))
points(sec6$DistN ~ c(sec6$Ncan - 2), type = "p", pch = 18, col = "cyan3")
legend("topright", legend = c("d = 0.5", "d = 0.3"), pch = c(8, 18), 
       col = c("tomato3", "cyan3"), text.col = c("tomato3", "cyan3"))

#### (c) ####




#### (d) ####
-->

```r
exp(-a[2]) - exp(-a[3]) * (1 + a[3] - a[2])
p2 <- pgamma(a[3] - a[2], shape = 2) * factorial(1)
p1 <- pgamma(a[2] - a[1], shape = 1) * factorial(0)
(1 - p1) * p2
h(5, a[6])

Ix <- function(x, p) pgamma(x, shape = p + 1)
p1 <- Ix(a[2], 0)
p2 <- Ix(a[3] - a[2], 1)
B2 <- 1 - Ix(a[3] - a[2], 1)
p3 <- 1/(2 * B2 * exp(a[3] - a[2])) * (2 * Ix(a[4] - a[3], 2) + 2 * (a[3] - a[2]) * Ix(a[4] - a[3], 1))
p1
(1 - p1) * p2
(1 - p1) * (1 - p2) * p3



####### scenario4 #########
> sec4
$`EN`
[1] 171.93

$SigN
[1] 21.52301

$CovProb
[1] 0.9465981

$nOpt
[1] 171

$cdf
[1] 1

$DistN
  [1] 3.639483e-03 3.195251e-04 4.885214e-05 1.133864e-05 3.504875e-06 1.339307e-06 6.041281e-07
  [8] 3.118164e-07 1.800920e-07 1.144550e-07 7.900506e-08 5.861695e-08 4.634810e-08 3.877849e-08
 [15] 3.412592e-08 3.142459e-08 3.014403e-08 3.000380e-08 3.088090e-08 3.276452e-08 3.573705e-08
 [22] 3.997183e-08 4.574366e-08 5.345129e-08 6.365333e-08 7.712086e-08 9.491204e-08 1.184766e-07
 [29] 1.498011e-07 1.916106e-07 2.476477e-07 3.230581e-07 4.249234e-07 5.629935e-07 7.506948e-07
 [36] 1.006509e-06 1.355859e-06 1.833669e-06 2.487837e-06 3.383910e-06 4.611347e-06 6.291860e-06
 [43] 8.590444e-06 1.172986e-05 1.600946e-05 2.182953e-05 2.972229e-05 4.039114e-05 5.475958e-05
 [50] 7.403150e-05 9.976432e-05 1.339563e-04 1.791489e-04 2.385436e-04 3.161332e-04 4.168432e-04
 [57] 5.466792e-04 7.128718e-04 9.240089e-04 1.190139e-03 1.522833e-03 1.935169e-03 2.441638e-03
 [64] 3.057925e-03 3.800550e-03 4.686346e-03 5.731757e-03 6.951961e-03 8.359812e-03 9.964638e-03
 [71] 1.177095e-02 1.377709e-02 1.597402e-02 1.834415e-02 2.086058e-02 2.348666e-02 2.617613e-02
 [78] 2.887378e-02 3.151685e-02 3.403703e-02 3.636309e-02 3.842398e-02 4.015229e-02 4.148778e-02
 [85] 4.238079e-02 4.279523e-02 4.271094e-02 4.212517e-02 4.105308e-02 3.952715e-02 3.759550e-02
 [92] 3.531942e-02 3.276998e-02 3.002426e-02 2.716139e-02 2.425861e-02 2.138781e-02 1.861255e-02
 [99] 1.598596e-02 1.354938e-02 1.133194e-02 9.350807e-03 7.612229e-03 6.112950e-03 4.842003e-03
[106] 3.782635e-03 2.914202e-03 2.213921e-03 1.658379e-03 1.224751e-03 8.916970e-04 6.399655e-04
[113] 4.527208e-04 3.156489e-04 2.168923e-04 1.468644e-04 9.799159e-05 6.442117e-05 4.172577e-05
[120] 2.662470e-05 1.673552e-05 1.036186e-05 6.319030e-06 3.795320e-06 2.244930e-06 1.307632e-06
[127] 7.500126e-07 4.235697e-07 2.355198e-07 1.289288e-07 6.948096e-08 3.685954e-08 1.924760e-08
[134] 9.892845e-09 5.004481e-09 2.491530e-09 1.220726e-09 5.885623e-10 2.792321e-10 1.303509e-10
[141] 5.987096e-11 2.705513e-11 1.202798e-11 5.260458e-12 2.263191e-12 9.577775e-13 3.986880e-13
[148] 1.632324e-13 6.573013e-14 2.603079e-14 1.013807e-14 3.882841e-15 1.462348e-15 5.415517e-16
[155] 1.971963e-16 7.060069e-17 2.485148e-17 8.600250e-18 2.925955e-18 9.785983e-19 3.217396e-19
[162] 1.039802e-19 3.303136e-20 1.031371e-20 3.165194e-21 9.546977e-22 2.830062e-22 8.244695e-23
[169] 2.360404e-23 6.640731e-24 1.835895e-24 4.987320e-25 1.331249e-25 3.491477e-26 8.997099e-27
[176] 2.277848e-27 5.665809e-28 1.384522e-28 3.323722e-29 7.838316e-30 1.815852e-30 4.132229e-31
[183] 9.236792e-32 2.028046e-32 4.373622e-33 9.264004e-34 1.927247e-34 3.937720e-35 7.901482e-36
[190] 1.557099e-36

$Ncan
  [1]   3   5   7   9  11  13  15  17  19  21  23  25  27  29  31  33  35  37  39  41  43  45  47  49
 [25]  51  53  55  57  59  61  63  65  67  69  71  73  75  77  79  81  83  85  87  89  91  93  95  97
 [49]  99 101 103 105 107 109 111 113 115 117 119 121 123 125 127 129 131 133 135 137 139 141 143 145
 [73] 147 149 151 153 155 157 159 161 163 165 167 169 171 173 175 177 179 181 183 185 187 189 191 193
 [97] 195 197 199 201 203 205 207 209 211 213 215 217 219 221 223 225 227 229 231 233 235 237 239 241
[121] 243 245 247 249 251 253 255 257 259 261 263 265 267 269 271 273 275 277 279 281 283 285 287 289
[145] 291 293 295 297 299 301 303 305 307 309 311 313 315 317 319 321 323 325 327 329 331 333 335 337
[169] 339 341 343 345 347 349 351 353 355 357 359 361 363 365 367 369 371 373 375 377 379 381
```

