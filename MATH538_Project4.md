```r
m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.05

options(digits = 6)

twoStageSamp <- function(m, d, sigma, alpha){
  X <- rnorm(m, 0, sigma); nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
  Ntilde <- ceiling((qt(1 - alpha/2, m - 1)/d)^2 * var(X))
  N <- max(m, Ntilde); X <- c(X, rnorm(N - m, 0, sigma)) 
  
  Ncandi <- m:(10^5)
  Qchi <- c(0, Ncandi) * (m - 1) * (d/(sigma * qt(1 - alpha/2, m - 1)))^2
  dist <- diff(pchisq(Qchi, m - 1))
  
  EN <- t(Ncandi) %*% dist
  VarN <- t(Ncandi^2) %*% dist - EN^2
  Pcov <- t(dist) %*% pnorm(d * sqrt(Ncandi)/sigma) 
  
  list(DistN = dist, EN = EN, SigN = sqrt(VarN), CovProb = Pcov, hatNopt = N, Nopt = nopt)
}

library(dplyr)
set.seed(1)
res <- twoStageSamp(2, 0.1, 1, 0.05)
sum(res$DistN)
max(res$DistN)

plot(res$DistN)
### table 1 ###

M <- c(2, 10, 20, 30, 500)
set.seed(1)
sapply(M, function(i) {
  res <- twoStageSamp(m = i, d = 0.5, sigma = 1, alpha = 0.05)
  c(E = res$EN, Sig = res$SigN, Pmu = res$CovProb, n_opt = res$Nopt)
}) 


### (c) ###
set.seed(1)
sce1 <- twoStageSamp(m = 10, d = 0.5, sigma = 1, alpha = 0.05)
set.seed(1)
sce2 <- twoStageSamp(m = 10, d = 0.3, sigma = 1, alpha = 0.05)
set.seed(1)
sce3 <- twoStageSamp(m = 10, d = 0.1, sigma = 1, alpha = 0.05)


set.seed(1)
sce4 <- twoStageSamp(m = 10, d = 0.5, sigma = 2, alpha = 0.05)
set.seed(1)
sce5 <- twoStageSamp(m = 10, d = 0.3, sigma = 2, alpha = 0.05)
set.seed(1)
sce6 <- twoStageSamp(m = 10, d = 0.1, sigma = 2, alpha = 0.05)


set.seed(1)
sce7 <- twoStageSamp(m = 10, d = 0.5, sigma = 1, alpha = 0.1)
set.seed(1)
sce8 <- twoStageSamp(m = 10, d = 0.3, sigma = 1, alpha = 0.1)
set.seed(1)
sce9 <- twoStageSamp(m = 10, d = 0.1, sigma = 1, alpha = 0.1)


par(mai = c(0.9, 0.9, 0.3, 0.2))
plot(sce1$DistN[1:1500] ~ c(0:1499), type = "l", col = "tomato3", 
     xlab = "k", pch = 20, lty = 1, ylab = expression(P[mu][","][sigma](N == m+k)))
points(sce2$DistN[1:1500], type = "l", lty = 2, col = "cyan3", pch = 18)
points(sce3$DistN[1:1500], type = "l", lty = 6, col = "gray", pch = 8)
legend("topright", legend = c("d = 0.5", "d = 0.3", "d = 0.1"), 
       lty = c(1, 2, 6), col = c("tomato3", "cyan3", "gray"), text.col = c("tomato3", "cyan3", "gray"))



plot(sce4$DistN[1:5000] ~ c(0:4999), type = "l", col = "tomato3", 
     xlab = "k", pch = 20, lty = 1, ylab = expression(P[mu][","][sigma](N == m+k)))
points(sce5$DistN[1:5000], type = "l", lty = 2, col = "cyan3", pch = 18)
points(sce6$DistN[1:5000], type = "l", lty = 6, col = "gray", pch = 8)
legend("topright", legend = c("d = 0.5", "d = 0.3", "d = 0.1"), 
       lty = c(1, 2, 6), col = c("tomato3", "cyan3", "gray"), text.col = c("tomato3", "cyan3", "gray"))




plot(sce7$DistN[1:800] ~ c(0:799), type = "l", col = "tomato3", 
     xlab = "k", pch = 20, lty = 1, ylab = expression(P[mu][","][sigma](N == m+k)))
points(sce8$DistN[1:800], type = "l", lty = 2, col = "cyan3", pch = 18)
points(sce9$DistN[1:800], type = "l", lty = 6, col = "gray", pch = 8)
legend("topright", legend = c("d = 0.5", "d = 0.3", "d = 0.1"), 
       lty = c(1, 2, 6), col = c("tomato3", "cyan3", "gray"), text.col = c("tomato3", "cyan3", "gray"))





### (d) ###

twoStageSamp <- function(m, d, sigma, alpha){
  X <- rnorm(m, 0, sigma); nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
  Ntilde <- ceiling((qt(1 - alpha/2, m - 1)/d)^2 * var(X))
  N <- max(m, Ntilde); X <- c(X, rnorm(N - m, 0, sigma)) 
  
  Ncandi <- m:(10^5)
  Qchi <- c(0, Ncandi) * (m - 1) * (d/(sigma * qt(1 - alpha/2, m - 1)))^2
  dist <- diff(pchisq(Qchi, m - 1))
  
  EN <- t(Ncandi) %*% dist
  # VarN <- t(Ncandi^2) %*% dist - EN^2
  # Pcov <- t(dist) %*% pnorm(d * sqrt(Ncandi)/sigma) 
  # 
  # list(DistN = dist, EN = EN, SigN = sqrt(VarN), CovProb = Pcov, hatNopt = N, Nopt = nopt)
  list(EN = EN)
}


set.seed(1)
dsce1 <- sapply(2:500, function(i) twoStageSamp(m = i, d = 0.5, sigma = 1, alpha = 0.05)$EN)
set.seed(1)
dsce2 <- sapply(2:500, function(i) twoStageSamp(m = i, d = 0.3, sigma = 1, alpha = 0.05)$EN)
set.seed(1)
dsce3 <- sapply(2:500, function(i) twoStageSamp(m = i, d = 0.1, sigma = 1, alpha = 0.05)$EN)


set.seed(1)
dsce4 <- sapply(2:2000, function(i) twoStageSamp(m = i, d = 0.5, sigma = 2, alpha = 0.05)$EN)
set.seed(1)
dsce5 <- sapply(2:2000, function(i) twoStageSamp(m = i, d = 0.3, sigma = 2, alpha = 0.05)$EN)
set.seed(1)
dsce6 <- sapply(2:2000, function(i) twoStageSamp(m = i, d = 0.1, sigma = 2, alpha = 0.05)$EN)

set.seed(1)
dsce7 <- sapply(2:400, function(i) twoStageSamp(m = i, d = 0.5, sigma = 1, alpha = 0.1)$EN)
set.seed(1)
dsce8 <- sapply(2:400, function(i) twoStageSamp(m = i, d = 0.3, sigma = 1, alpha = 0.1)$EN)
set.seed(1)
dsce9 <- sapply(2:400, function(i) twoStageSamp(m = i, d = 0.1, sigma = 1, alpha = 0.1)$EN)



plot(dsce1[1:500] ~ c(2:501), type = "l", col = "tomato3", 
     xlab = "m", pch = 20, lty = 1, ylab = expression(E[mu][","][sigma](N)))
points(dsce2[1:500], type = "l", lty = 2, col = "cyan3", pch = 18)
points(dsce3[1:500], type = "l", lty = 6, col = "gray", pch = 8)
legend("bottomright", legend = c("d = 0.5", "d = 0.3", "d = 0.1"), 
       lty = c(1, 2, 6), col = c("tomato3", "cyan3", "gray"), text.col = c("tomato3", "cyan3", "gray"))



plot(dsce4[1:3000] ~ c(2:3001), type = "l", col = "tomato3", 
     xlab = "m", pch = 20, lty = 1, ylab = expression(E[mu][","][sigma](N)))
points(dsce5[1:3000], type = "l", lty = 2, col = "cyan3", pch = 18)
points(dsce6[1:3000], type = "l", lty = 6, col = "gray", pch = 8)
legend("bottomright", legend = c("d = 0.5", "d = 0.3", "d = 0.1"), 
       lty = c(1, 2, 6), col = c("tomato3", "cyan3", "gray"), text.col = c("tomato3", "cyan3", "gray"))

plot(dsce9[1:400] ~ c(2:401), type = "l", lty = 6, col = "gray", pch = 8)



plot(dsce7[1:400] ~ c(2:401), type = "l", col = "tomato3", 
     xlab = "m", pch = 20, lty = 1, ylab = expression(E[mu][","][sigma](N)))
points(dsce8[1:400], type = "l", lty = 2, col = "cyan3", pch = 18)
points(dsce9[1:400], type = "l", lty = 6, col = "gray", pch = 8)
legend("bottomright", legend = c("d = 0.5", "d = 0.3", "d = 0.1"), 
       lty = c(1, 2, 6), col = c("tomato3", "cyan3", "gray"), text.col = c("tomato3", "cyan3", "gray"))


### (e) ###
plot(dsce7[2:50] ~ c(3:51), type = "p", col = "tomato3", 
     xlab = "m", pch = 20, lty = 1, ylab = expression(E[mu][","][sigma](N)))



```






https://www.jstor.org/stable/2236965?seq=4#metadata_info_tab_contents
