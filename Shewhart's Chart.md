<!--
```{r}
gam <- 10
theta <- 0.1

A <- exp(theta * qnorm(1 - 1/gam) - theta^2/2)
1/pnorm(1/theta * log(A) + theta/2, theta, 1, lower.tail = FALSE)



ADD <- function(ARL, theta){
  unlist(lapply(ARL, function(i) {quan <- qnorm(1 - 1/i); 1/(1 - pnorm(quan, theta, 1))}))
}


ADD1 <- ADD(1:10000, 0.1)
ADD2 <- ADD(1:10000, 0.5)
ADD3 <- ADD(1:10000, 1)

plot(1:10000, ADD1, type = "l", lty = 3, col = "tomato3", xlab = "ARL", ylab = "ADD")
points(ADD2, type = "l", lty = 5, col = "cyan3")
points(ADD3, type = "l", lty = 7, col = "gray")
legend("topleft", legend = c(expression(theta == 0.1), expression(theta == 0.5), expression(theta == 1)),
       lty = c(3, 5, 7), col = c("tomato3", "cyan3", "gray"), text.col = c("tomato3", "cyan3", "gray"), bty = "n")


```
-->
