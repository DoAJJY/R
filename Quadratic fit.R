####################### False position method ########################
##### Copyright (C) 2016.11.16 Jaeyoung Jung all rights reserved. ####

rm(list=ls(all=TRUE))
gc(reset=TRUE)

library(mosaic)


### Define function and derivatives
# Original function
a <- 0
b <- 4

z <- makeFun(2*x^3 + 3*x^2 - 12*x + 5 ~ x)

curve(z, a, b)

# Tolerence
epsilon <- 0.001


### Newton method
# Initial values(condition of problem 4.7)
x_0 <- 0
x_1 <- 0.5
x_2 <- 1

Niter <- 1000 # Sufficiently large number

# Define newton method as a function
quadratic_fit <- function(f, tol, x_0, x_1, x_2, N) {
    # First order derivative
  
    i <- 1
    p <- numeric(N)
    
    while (i <= N) {
        numer <- (x_1^2 - x_2^2)*f(x_0) + (x_2^2 - x_0^2)*f(x_1) + (x_0^2 - x_1^2)*f(x_2)
        denom <- (x_1 - x_2)*f(x_0) + (x_2 - x_0)*f(x_1) + (x_0 - x_1)*f(x_2)
        x_3 <- (numer / denom) / 2
        p[i] <- x_3
        i <- i+1
        if (abs(x_3 - x_2) < tol) break
        x_0 <- x_1
        x_1 <- x_2
        x_2 <- x_3
        print(paste("Niter = ", i, ", ", "x_(n+1) = ", x_2, sep = ""))
    }
    
    return(p[1:(i-1)])
}

# Optimization
quadratic_fit(z, epsilon, x_0, x_1, x_2, Niter)
