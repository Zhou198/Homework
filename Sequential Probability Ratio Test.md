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

SPRT <- function(the, the0, the1, sigma, N, alpha0, alpha1){
  a0 <- log((1 - alpha0)/alpha1); a1 <- log((1 - alpha1)/alpha0)
  y <- seq(-a0, a1, length = N + 1); x <- (y[-1] + y[-(N + 1)])/2 
  K <- matrix(0, N, N)
  
  ASN <- lapply(the, function(l) {
    for (i in 1:N) { K[i, ] <- diff(pnorm(y - x[i], (the1 - the0)/sigma^2 * (l - (the0 + the1)/2), (the1 - the0)/sigma))}
    solve(diag(N) - K, rep(1, N))[which.min(abs(x))]
  })
  list(bound = c(a0 = a0, a1 = a1), ASN = unlist(ASN))
}

##### Scenario1 #####

theta <- seq(0, 0.5, length = 100)
sprt1 <- SPRT(theta, 0, 0.5, 1, 1000, 0.001, 0.001)
a <- sprt1$bound
a

par(mai = c(0.9, 0.9, 0.3, 0.2))
plot(sprt1$ASN ~ theta, ylim = c(0, 300), type = "l", xlab = expression(theta), ylab = expression(ASN(theta)), col = "tomato3")




##### Scenario2 #####
sprt2 <- SPRT(theta, 0, 0.5, 1, 1000, 0.001, 0.0001)
a <- sprt2$bound
a

plot(sprt2$ASN ~ theta, ylim = c(0, 300), type = "l", xlab = expression(theta), ylab = expression(ASN(theta)), col = "tomato3")



##### Scenario3 #####
sprt3 <- SPRT(theta, 0, 0.5, 1, 1000, 0.0001, 0.001)
a <- sprt3$bound
a

plot(sprt3$ASN ~ theta, ylim = c(0, 300), type = "l", xlab = expression(theta), ylab = expression(ASN(theta)), col = "tomato3")



log((1 - 0.0001)/0.001)

```