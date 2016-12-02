###################### Golden section method #########################
##### Copyright (C) 2016.11.21 Jaeyoung Jung all rights reserved. ####
library(mosaic)
library(limSolve)

rm(list=ls(all=TRUE))
gc(reset=TRUE)


### Define function : Problem 4.20
# Target function
f <- function(x) {
  y <- 4*(x[1]-10)^2 + (x[2]-4)^2
  return(y)
}

# First order derivative
df <- function(x) {
  y <- c(8*(x[1]-10), 2*(x[2]-4))
  return(y)
}

# Step size
alpha <- function(x, y) {
  numer <- 40*(y[1]-x[1]) + 4*(y[2]-x[2]) - 4*x[1]*(y[1]-x[1]) - x[2]*(y[2]-x[2])
  denom <- 4*(y[1]-x[1])^2 + (y[2]-x[2])^2
  h <- numer / denom
  return(h)
}


### Initial value
a <- c(14, 0) # Set random initial point x^1 in feasible region

# Convergence value
epsilon <- 0.0001


### Convex combination optimization algorithm
Ni <- 1 # Set initial number of iteration

# Start optimization
while (TRUE) {
  # Step 1 : Solve linear program
  G <- matrix(nrow = 4, data = c(1, 0, 0, 1, -1, 0, 0, -1), byrow = TRUE) # Inequality conditions
  H <- c(9, 0, -14, -6) # Range of x_1 and x_2
  (L <- linp(E = NULL, F = NULL, Cost = df(a), G = G, H = H)) # Solve linear program by using 'limSolve' package
  
  b <- L$X # Y value that minimize the direction, 'gradient*Y'
  print(paste("The point that minimizes the direction is ", "[", b[1], ", ", b[2], "]", sep = ""))
  
  # Step 2 : Calculate the step size
  alp <- alpha(a, b)
  
  # Step 3 : Calculate the next point
  x_n <- a + alp * (b - a)
  
  # Print the result of iteration
  print(paste("Niter = ", Ni, ", ", Ni+1, "th X is ", "[", x_n[1], ", ", x_n[2], "]", ", ", "Convergence is ",
              f(a)-f(x_n), sep = ""))
  
  # Step 4 : Convergence test
  if (f(a)-f(x_n) <= epsilon) break
  
  # Set x for the next iteration
  a <- x_n
  
  # Go next iteration
  Ni <- Ni + 1
} 
