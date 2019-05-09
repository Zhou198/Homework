```{r}
##########Project 8#############
cusum <- function(truThe, theta, N, ARL){
  y <- seq(0, ARL, length = N + 1); x <- (y[-1] + y[-(N + 1)])/2 
  K <- matrix(0, N, N)
  
  unlist(lapply(truThe, function(l) {
    for (i in 1:N) {
      K[i, ] <- diff(pnorm(y - x[i], theta * (l - 0.5 * theta), theta))
      num <- solve(diag(N) - K, cbind(1, diff(pnorm(-x, theta * (l - 0.5 * theta), theta))))[1, ]
      }
    num[1]/(1 - num[2])
    }))
}

SPRT <- function(the, Hyp, sigma, N, a0, a1){
  y <- seq(-a0, a1, length = N + 1); x <- (y[-1] + y[-(N + 1)])/2 
  K <- matrix(0, N, N)
  
  unlist(lapply(the, function(l) {
    for (i in 1:N) {
      K[i, ] <- diff(pnorm(y - x[i], diff(Hyp)/sigma^2 * (l - mean(Hyp)), diff(Hyp)/sigma))
      solve(diag(N) - K, rep(1, N))[which.min(abs(x))]
      }}))
}
```
