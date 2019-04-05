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


```

