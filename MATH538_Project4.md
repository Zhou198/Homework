```r
m <- 30; d <- 0.5; sigma <- 1; alpha <- 0.05

set.seed(1); X <- rnorm(m, 0, sigma)
Ntilde <- ceiling(qt(1 - alpha/2, m - 1)^2 * var(X)/d^2)
nopt <- (qnorm(1 - alpha/2) * sigma/d)^2

Qchi <- m * (m - 1) * (d/(sigma * qt(1 - alpha/2, m - 1)))^2

Pm <- pchisq(Qchi, m - 1); dist <- c(Pm, 1 - Pm)
EN <- t(dist) %*% c(m, Ntilde)
VarN <- cbind(Pm, 1 - Pm) %*% c(m, Ntilde)^2 - EN^2
Pcov <- pnorm(d * sqrt(max(m, Ntilde))/sigma)

set.seed(1)
X <- c(X, rnorm(max(Ntilde - m, 0), 0, sigma)) 

list(DistN = dist, EN = EN, SigN = sqrt(VarN), CovProb = Pcov, Nopt = nopt)
```






https://www.jstor.org/stable/2236965?seq=4#metadata_info_tab_contents
