###################### Golden section method #########################
##### Copyright (C) 2016.11.16 Jaeyoung Jung all rights reserved. ####

rm(list=ls(all=TRUE))
gc(reset=TRUE)


### Target function : Problem 4.3
z <- function(x) {
  y = 3 * x^2 - 4 * x + 1
  return(y)
}

# Draw function graph
curve(z, 0, 2)


### Initial value
a <- 0
b <- 2

# Reduction ratio - Golden ratio : 1/2 * (sqrt(5) - 1)
ratio <- 1/2 * (sqrt(5) - 1)

# Convergence ratio
epsilon <- 0.001


### Optimization
for (i in 1:100) {
  x_L <- (b - a) * (1 - ratio) + a
  x_R <- (b - a) * ratio + a
  
  if (z(x_L) <= z(x_R)) {
    b <- x_R
    if ((b - a)/2 <= epsilon) {
      x_hat <- (b + a) / 2
      print(paste("x_hat = ", x_hat, ", ", "Niter = ", i, ", ", "Convergence interval = ", (b - a) / 2, sep = ""))
      break()
    } else {
      x_R <- x_L
      x_L <- (b - a) * (1 - ratio) + a
      print(paste("x_hat = ", (b + a) / 2, ", ", "Niter = ", i, ", ", "Convergence interval = ", (b - a) / 2, sep = ""))
      next
    }
  } else {
    a <- x_L
    if ((b - a)/2 <= epsilon) {
      x_hat <- (b + a) / 2
      print(paste("x_hat = ", x_hat, ", ", "Niter = ", i, ", ", "Convergence interval = ", (b - a) / 2, sep = ""))
      break()
    } else {
      x_L <- x_R
      x_R <- (b - a) * ratio + a
      print(paste("x_hat = ", (b + a) / 2, ", ", "Niter = ", i, ", ", "Convergence interval = ", (b - a) / 2, sep = ""))
      next
    }
  }
} 
