################  Homework1 01/23/2020  #################
#### T3 ####
T <- 10000; n <- 25; SigX <- 2; SigY <- 3; gam <- (SigY/SigX)^2
set.seed(1)
hatgam <- sapply(1:T, function(t) {
  X <- rnorm(n, 0, SigX); Y <- rnorm(n, 0, SigY)
  var(Y)/var(X)
})
CI <- hatgam %*% matrix(qf(c(0.025, 0.975), 24, 24), ncol = 2)

mean(CI[, 1] <= gam & gam <= CI[, 2])

############  ———————— ###############
################  Homework2 01/25/2020  #################
#### T3 ####
### Method 1 ###
########### Randomn Number Generator #############
rX <- function(n){
  u <- runif(n)
  ifelse(u < 1/3, 3 * u, 1.5 * u + 0.5)
}

## sample mean
N <- 5000; Y <- 1
set.seed(1)
#tilMu1 <- sum(X * (X <= Y))/sum(X <= Y)
tilMu1 <- sapply(10:N, function(n){
  X <- rX(n); sum(X * (X <= Y))/sum(X <= Y)
})

par(mai = c(0.83, 0.8, 0.2, 0.1)) ## 600 x 420
plot(10:N, tilMu1, type = "l", col = "cyan3", xlab = "Sample Size", ylab = expression(tilde(mu)))
abline(h = 1/2, col = "tomato3")

## modified edf
par(mai = c(0.83, 0.88, 0.2, 0.1))
x <- seq(-1, 3, length = 2000)
Fx <- x/3 * (0 < x & x <= 1) + (2 * x - 1)/3 * (1 < x & x <= 2) + (x > 2)
plot(x, Fx, type = "l", ylab = expression(tilde(F)(x)))


tilF1 <- function(x){
  sum(X <= Y & X <= x)/sum(X <= Y)
}

tilP1 <- sapply(x, function(t) tilF1(t))
points(x, tilP1, type = "l")

sam <- c(5000, 1000, 100, 50, 20)
col <- c("tomato3", "green3", "cyan3", "darkorchid3", "burlywood3")
set.seed(1)
for (i in 1:5) {
   X <- rX(sam[i]); Y <- 1
   tilP1 <- sapply(x, function(t) tilF1(t))
   points(x, tilP1, type = "l", col = col[i])
}

legend("topleft", legend = c("N = 5000", "N = 1000", "N = 100", "N = 50", "N = 20", "Exact CDF"), bty = "n", lty = rep(1, 6), col = c("tomato3", "green3", "cyan3", "darkorchid3", "burlywood3", "black"), text.col = c("tomato3", "green3", "cyan3", "darkorchid3", "burlywood3", "black"))

### Method 2 ###
## sample mean
#tilMu2 <- mean(X * (X <= Y) + Y * (X > Y))
set.seed(1)
tilMu2 <- sapply(10:N, function(n) {
  X <- rX(n); mean(pmin(X, Y))
})

plot(10:N, tilMu2, type = "l", col = "cyan3", xlab = "Sample Size", ylab = expression(tilde(mu)))
abline(h = 5/6, col = "tomato3")

## modified edf
tilF2 <- function(x) mean(pmin(X, Y) <= x) 

x <- seq(-1, 3, length = 2000)
Fx <- x/3 * (0 < x & x <= 1) + (2 * x - 1)/3 * (1 < x & x <= 2) + (x > 2)
plot(x, Fx, type = "l", ylab = expression(tilde(F)(x)))

set.seed(1)
for (i in 1:5) {
  X <- rX(sam[i]); Y <- 1
  tilP2 <- sapply(x, function(t) tilF2(t))
  points(x, tilP2, type = "l", col = col[i])
}

legend("topleft", legend = c("N = 5000", "N = 1000", "N = 100", "N = 50", "N = 20", "Exact CDF"), bty = "n", lty = rep(1, 6), col = c("tomato3", "green3", "cyan3", "darkorchid3", "burlywood3", "black"), text.col = c("tomato3", "green3", "cyan3", "darkorchid3", "burlywood3", "black"))


## change y into 1.5
x <- seq(-1, 3, length = 2000)
Fx <- x/3 * (0 < x & x <= 1) + (2 * x - 1)/3 * (1 < x & x <= 2) + (x > 2)
plot(x, Fx, type = "l", ylab = expression(tilde(F)(x)))


set.seed(2)
for (i in 1:5) {
  X <- rX(sam[i]); Y <- 1.5
  tilP2 <- sapply(x, function(t) tilF2(t))
  points(x, tilP2, type = "l", col = col[i])
}

legend("topleft", legend = c("N = 5000", "N = 1000", "N = 100", "N = 50", "N = 20", "Exact CDF"), bty = "n", lty = rep(1, 6), col = c("tomato3", "green3", "cyan3", "darkorchid3", "burlywood3", "black"), text.col = c("tomato3", "green3", "cyan3", "darkorchid3", "burlywood3", "black"))



############  ———————— ###############
################  Homework3 01/29/2020  #################
library(ggplot2)

#### T1 ####
rN <- function(n, p){
  Intvl <- cumsum(c(0, p))
  sapply(runif(n), function(u) sum(Intvl <= u))
} # p is a preset prob. vector for R.V. N


Time <- 100; p <- c(0.1, 0.2, 0.05, 0.2, 0.15, 0.15, 0.1, 0.05)

set.seed(1)
X <- rexp(Time, 1/3); N <- rN(Time, p)

XStar <- sapply(1:Time, function(i){
  Y <- runif(N[i], 0, 2) + 2 * (1:N[i])
  if (Y[1] < X[i] & X[i] <= Y[N[i]]) {k <- sum(Y < X[i]); 0.5 * (Y[k] + Y[k + 1])}
  else Y[1] * (X[i] <= Y[1]) + Y[N[i]] * (X[i] > Y[N[i]])
})

range(XStar)
x <- seq(min(XStar) - 5, max(XStar) + 5, length = 500)
cdf <- sapply(x, function(s) mean(XStar <= s))
excdf <- sapply(x, function(s) mean(X <= s))

par(mai = c(0.83, 0.8, 0.2, 0.1))
#plot(x, cdf, type = "l")

datCDF <- cbind(x = x, prob. = cdf, exCDF = excdf)
ggplot(data.frame(datCDF), aes(x = x)) + geom_line(aes(y = prob.), col = "tomato3") + geom_line(aes(y = exCDF), col = "cyan3") + theme(legend.position = c(0.1, 0.86))  + theme(legend.title = element_blank())

                

                

################ Revise #################                
Ftil <- function(p, Size, t, Time) {
  sapply(1:Time, function(Ti) {
    X <- rexp(Size, 1/3); N <- rN(Size, p)
    
    XStar <- sapply(1:Size, function(i){
      Y <- runif(N[i], 0, 2) + 2 * (1:N[i])
      if (Y[1] < X[i] & X[i] <= Y[N[i]]) {k <- sum(Y < X[i]); 0.5 * (Y[k] + Y[k + 1])}
      else Y[1] * (X[i] <= Y[1]) + Y[N[i]] * (X[i] > Y[N[i]])
    })
    
    mean(XStar <= t)
  })
}

1 - exp(-1/3 * 2.5)


set.seed(1)
F3 <- Ftil(p, 100, 2.5, 10)
plot(F3)


set.seed(1)
F3asy <- sapply(1:5000, function(i) Ftil(p, i, 3, 1))
plot(1:5000, F3asy, type = "l")
abline(h = 0.34)

set.seed(1)
F2.5asy <- sapply(1:5000, function(i) Ftil(p, i, 2.5, 1))
plot(1:5000, F2.5asy, type = "l", col = "tomato3")

par(mai = c(0.83, 0.8, 0.2, 0.1))
#plot(x, cdf, type = "l")

datCDF <- cbind(x = x, prob. = cdf, hatCDF = hatcdf, exCDF = excdf)
ggplot(data.frame(datCDF), aes(x = x)) + 
  geom_line(aes(y = prob.), col = "tomato3") + 
  geom_line(aes(y = hatCDF), col = "black") +
  geom_line(aes(y = exCDF), col = "cyan3") + 
  


                
                
                
#### T2 ####
set.seed(1)
U <- runif(400, 0, 100); V <- U + 40
Tf <- rexp(400, 1/56); ind <- which(Tf <= V)
TStar <- (U * (Tf <= U) + 0.5 * (U + V) * (U < Tf & Tf <= V))[ind]

x <- seq(min(Tf) - 40, max(Tf) + 20, length = 1000)
S <- (x <= 0) + exp(-1/56 * x) * (x > 0)
Stil <- sapply(x, function(t) mean(TStar > t))

datS <- cbind(x = x, exS = S, Stil = Stil)
ggplot(data.frame(datS), aes(x = x)) + geom_line(aes(y = exS, col = "tomato3")) + geom_line(aes(y = Stil, col = "cyan3")) + labs(x = "x", y = "Probability") + theme(legend.background = element_blank(), legend.title = element_blank(), legend.position = c(0.93, 0.9)) + scale_color_manual(values = c("tomato3" = "tomato3", "cyan3" = "cyan3"), breaks = c("tomato3", "cyan3"), labels = c("S(x)", expression(tilde(S)(x))))





























