```r
m <- 2; d <- 0.5; sigma <- 1; alpha <- 0.05

options(digits = 6)

twoStageSamp <- function(m, d, sigma, alpha){
  X <- rnorm(m, 0, sigma); nopt <- ceiling((qnorm(1 - alpha/2) * sigma/d)^2)
  Ntilde <- ceiling((qt(1 - alpha/2, m - 1)/d)^2 * var(X))
  N <- max(m, Ntilde); X <- c(X, rnorm(N - m, 0, sigma)) 
  
  Ncandi <- m:(100^3)
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

### table 2 ###
set.seed(1)
sapply(M, function(i) {
  res <- twoStageSamp(m = i, d = 0.3, sigma = 1, alpha = 0.05)
  c(E = res$EN, Sig = res$SigN, Pmu = res$CovProb, n_opt = res$Nopt)
})

### table 3 ###
options(digits = 3)
set.seed(1)
sapply(M, function(i) {
  res <- twoStageSamp(m = i, d = 0.1, sigma = 1, alpha = 0.05)
  c(E = res$EN, Sig = res$SigN, Pmu = res$CovProb, n_opt = res$Nopt)
})
```






https://www.jstor.org/stable/2236965?seq=4#metadata_info_tab_contents
