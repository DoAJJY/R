######################### Bisection method ###########################
##### Copyright (C) 2016.11.16 Jaeyoung Jung all rights reserved. ####

rm(list=ls(all=TRUE))
gc(reset=TRUE)


### Target function : Problem 4.3
z <- function(x) {
  y <- 3 * x^2 - 4 * x + 1
  return(y)
}

# First order derivative
dz <- function(x) {
  y <- 6 * x - 4
  return(y)
}

# Draw function graph
curve(z, 0, 2)


### Initial value
a <- 0
b <- 2

# Reduction ratio - Bisection ratio : 1/2 ...
ratio <- 1/2

# Convergence ratio
epsilon <- 0.001


### Optimization
Niter <- 10000 # Sufficient number
x_B <- rep(0, Niter)
x_B[1] <- a

for (i in 2:Niter) {
  x_B[i] <- (b + a) * ratio
  
  if (dz(x_B[i]) > 0) {
    b <- x_B[i]
    
    if (abs(x_B[i] - x_B[i-1]) <= epsilon) {
      x_hat <- x_B[i]
      print(paste("x_hat = ", x_B[i], ", ", "Niter = ", i-1, ", ", "Convergence interval = ", abs(x_B[i] - x_B[i-1]), sep = ""))
      break()
    } else {
      print(paste("x_hat = ", x_B[i], ", ", "Niter = ", i-1, ", ", "Convergence interval = ", abs(x_B[i] - x_B[i-1]), sep = ""))
      next
    }
  } else {
    a <- x_B[i]
    
    if (abs(x_B[i] - x_B[i-1]) <= epsilon) {
      x_hat <- x_B[i]
      print(paste("x_hat = ", x_B[i], ", ", "Niter = ", i-1, ", ", "Convergence interval = ", abs(x_B[i] - x_B[i-1]), sep = ""))
      break()
    } else {
      print(paste("x_hat = ", x_B[i], ", ", "Niter = ", i-1, ", ", "Convergence interval = ", abs(x_B[i] - x_B[i-1]), sep = ""))
      next
    }
  }
} 