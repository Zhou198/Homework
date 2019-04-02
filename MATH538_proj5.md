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
```
-->
