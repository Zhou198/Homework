################  Homework5 02/14/2020  #################
library(splines)
library(survival)

#### T1 ####
set.seed(1)
Y <- rweibull(100, 2, 1); U <- runif(100, 0, 2)
d <- as.numeric(Y <= U)
M <- Y * (Y <= U) + U * (Y > U)

sum(d)

#### T2 ####
# set.seed(1)
# Y <- rweibull(100, 0.5, 1); U <- runif(100, 0, 2)
# d <- as.numeric(Y <= U)
# M <- Y * d + U * (1 - d)
# 
# sum(d)


model <- survreg(Surv(M, d) ~ 1)
val <- summary(model)
c(val$coef, sigma = val$scale)


der2 <- function(ka, ta, x, d) {
  xdt <- x/ta
  matrix(c(-sum(d)/(ka^2) - sum(xdt^ka * (log(xdt))^2),
           -sum(d)/ta + sum(xdt^ka/ta * (1 + ka * log(xdt))),
           -sum(d)/ta + sum(xdt^ka/ta * (1 + ka * log(xdt))),
           ka * sum(d)/(ta^2) - (ka^2 + ka)/(ta^2) * sum(xdt^ka)), nrow = 2)
}

der1 <- function(ka, ta, x, d) {
  xdt <- x/ta
  matrix(c((1/ka - log(ta)) * sum(d) + sum(d * log(x)) - sum(xdt^ka * log(xdt)),
           -ka/ta * sum(d) + ka/ta * sum(xdt^ka)), ncol = 1)
}

KT <- c(1, 1); S <- c(0, 0) 
while (max(abs(KT - S)) >= 0.00001) {
  S <- KT; inv <- solve(der2(KT[1], KT[2], M, d), diag(2))
  KT <- KT - inv %*% der1(KT[1], KT[2], M, d)
}

c(alpha = log(KT[2]), logsigma = -log(KT[1]), sigma = 1/KT[1])

KT <- c(1/val$scale, exp(val$coef))


####### T3 #######
val$var

A <- matrix(c(0, 1/KT[2], -1/KT[1], 0), nrow = 2)
t(A) %*% -solve(der2(KT[1], KT[2], M, d), diag(2)) %*% A



#### T4 ####
vrho <- t(c(0, -1/(KT[2]^2))) %*% -solve(der2(KT[1], KT[2], M, d), diag(2)) %*% c(0, -1/(KT[2]^2))
1/KT[2] + c(-1, 1) * qnorm(0.975) * as.vector(sqrt(vrho))



#### T5 ####
# set.seed(1)
# Y <- rweibull(100, 2, 1); U <- runif(100)
# d <- as.numeric(Y <= U) 
# M <- Y * d + U * (1 - d)

model5 <-  survreg(Surv(M, d) ~ 1)
est <- data.frame(ka1 = 1/model5$scale, ta1 = exp(model5$coef), ta0 = mean(M)/mean(d))


survreg(Surv(M, d) ~ 1, dist = "exponential")


2 * ((log(est[1] * est[3]) - est[1] * log(est[2])) * sum(d) + (est[1] - 1) * sum(d * log(M)) + sum(M)/est[3] - sum((M/est[2])^est[1]))

qchisq(0.95, df = 1)



#### T7 ####
set.seed(1)
Y <- rweibull(4, 2, 1); U <- runif(4)
d <- as.numeric(Y <= U) 
M <- Y * d + U * (1 - d)


set.seed(1)
Y <- rweibull(4, 2, 1); U <- runif(4, 0, 2)
d <- as.numeric(Y <= U)
M <- Y * (Y <= U) + U * (Y > U)

model5 <-  survreg(Surv(M, d) ~ 1)
est <- data.frame(ka1 = 1/model5$scale, ta1 = exp(model5$coef), ta0 = mean(M)/mean(d))

2 * ((log(est[1] * est[3]) - est[1] * log(est[2])) * sum(d) + (est[1] - 1) * sum(d * log(M)) + sum(M)/est[3] - est[2]^(-est[1]) * sum(M^est[1]))

