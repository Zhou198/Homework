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




cusum <- function(truThe, theta, N, ARL){
  y <- seq(0, ARL, length = N + 1); x <- (y[-1] + y[-(N + 1)])/2 
  K <- matrix(0, N, N)
  
  unlist(lapply(truThe, function(l) {
    for (i in 1:N) {
      K[i, ] <- diff(pnorm(y - x[i], theta * (l - 0.5 * theta), theta))
    }
    num <- solve(diag(N) - K, cbind(rep(1, N), pnorm(-x, theta * (l - 0.5 * theta), theta)))[1, ]
    num[1]/(1 - num[2])
    }))
}

cusum(0.5, 0.5, 1000, 5)



#################### updated #############################
cusum <- function(theta, N, ARL){
  unlist(lapply(ARL, function(l) {
    y <- seq(0, l, length = N + 1); x <- (y[-1] + y[-(N + 1)])/2
    K <- matrix(0, N, N)
    
    for (i in 1:N) {
      K[i, ] <- diff(pnorm(y - x[i], 0.5 * theta^2 , theta))
    }
    
    num <- solve(diag(N) - K, cbind(rep(1, N), pnorm(-x, 0.5 * theta^2, theta)))[1, ]
    num[1]/(1 - num[2])
  }))
}

ARL <- seq(0.1, 100, 0.1)
length(ARL)
EN <- cusum(0.5, 1000, ARL)
plot(ARL, EN, type = "l", col = "cyan3")

```
