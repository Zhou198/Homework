```r
alpha0 <- 0.001; alpha1 <- 0.001; N <- 1000
sigma <- 1; theta0 <- 0; theta1 <- 0.5; theta <- 0.5
a <- log(alpha1/(1 - alpha0)) - 0.07
b <- log((1 - alpha1)/alpha0) + 0.07
a;b
y <- seq(a, b, length = N + 1)
x <- (y[-1] + y[-(N + 1)])/2
Ktheta <- matrix(0, N, N)

for (i in 1:N) {
  Ktheta[i, ] <- diff(pnorm(y - x[i], (theta1 - theta0)/sigma^2 * (theta - 0.5 * (theta0 + theta1)), (theta1 - theta0)/sigma)) 
}

ASN <- solve(diag(N) - Ktheta, rep(1, N))

OC <- solve(diag(N) - Ktheta, pnorm(a - x, (theta1 - theta0)/sigma^2 * (0.5 - 0.5 * (theta0 + theta1)), (theta1 - theta0)/sigma))

ASN[which.min(abs(x))]
OC[which.min(abs(x))]

1 - OC[which.min(abs(x))]
```
